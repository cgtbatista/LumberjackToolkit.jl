"""
    densityprofile(pdb::String, trajectory::String; ...)

It computes the density profile for each one of the frames inside the trajectory.

# Arguments
- `pdb::String`: a String with the pdb file path.
- `trajectory::String`: a String with the trajectory file path.
- `profile::String`: a String with the profile type. It can be "number", "mass", "charge", or "electron".
- `selection::String`: a String with the selection of the atoms.
- `axis::String`: a String with the axis flag.
- `dimensions::Int`: an integer with the number of dimensions. It can be 1 or 2.
- `resolution::Float64`: a Float64 with the resolution of the bins.
- `center::Bool`: a Bool to center the profile.
- `symmetry::Bool`: a Bool to apply the symmetry protocol.
- `symmetry_selection::String`: a String with the selection of the atoms to apply the symmetry protocol.
- `first::Int`: an integer with the first frame to be considered.
- `step::Int`: an integer with the step to be considered.
- `last::Int`: an integer with the last frame to be considered.
- `charge_column::String`: a String with the charge column to be considered.
- `edp::String`: a String with the electron density profile.
- `charge_correction::Bool`: a Bool to apply the charge correction.
- `normalization::Bool`: a Bool to apply the normalization.
"""
function densityprofile(
        pdb::String,
        trajectory::String;
        profile="electron", selection="all",
        axis="z", dimensions=1, resolution=1, center=true, symmetry=true, symmetry_selection="not water", first=1, step=1, last=nothing,
        charge_column=nothing, edp=nothing, charge_correction=true, normalization=true
    )
    
    distances, densities = [], []

    println("$(dimensions)D -- $profile density distribution of the molecules inside the box:")
    println("   the selection is `$selection`")
    println("   resolution on $axis is equal to $resolution Å")
    println("")

    simulation = MolSimToolkit.Simulation(pdb, trajectory; first=first, last=last, step=step)
    selected_atoms = PDBTools.select(MolSimToolkit.atoms(simulation), selection)
    idx = PDBTools.index.(selected_atoms)

    bins, V_norm = binning(simulation, axis=axis, resolution=resolution)
    if symmetry || center;
        references = _get_reference(simulation, symmetry_selection)
        if lowercase(string(axis)) == "x"
            reference = [ ref[1] for ref in references ]
        elseif lowercase(string(axis)) == "y"
            reference = [ ref[2] for ref in references ]
        elseif lowercase(string(axis)) == "z"
            reference = [ ref[3] for ref in references ]
        else
            throw(ArgumentError("The axis flag should be coded as `x`, `y`, or `z` strings."))
        end
    end

    if profile == "number"
        println("   ... choosing the number of molecules density profile (NDP)")
        println("")
        property = ones(Float64, length(idx))
    elseif profile == "mass"
        property = PDBTools.atomic_mass.(selected_atoms)
    elseif profile == "charge"
        if isnothing(charge_column) || charge_column == "beta" || charge_column == "b-factor"
            property = PDBTools.beta.(selected_atoms)
        elseif charge_column == "occup" || charge_column == "occupancy"
            property = PDBTools.occup.(selected_atoms)
        end
    elseif profile == "electron"
        if charge_correction
            if isnothing(charge_column) || charge_column == "beta" || charge_column == "b-factor"
                atomic_charges = PDBTools.beta.(selected_atoms)
            elseif charge_column == "occup" || charge_column == "occupancy"
                atomic_charges = PDBTools.occup.(selected_atoms)
            end
            property = PDBTools.atomic_number.(selected_atoms) .- atomic_charges
        else
            property = PDBTools.atomic_number.(selected_atoms)
        end
    end

    for frame in simulation
        i = simulation.frame_index
        println("    - frame $i")
        coords = MolSimToolkit.positions(frame)[idx]
        new_coords = [ MolSimToolkit.Point3D(coords[n][1], coords[n][2], coords[n][3]) for n in eachindex(coords) ]

        d = ρ(bins, new_coords, V_norm, axis=axis, prop=property)

        if center && !symmetry
            new_bins = _correct_center(bins, reference[i])
        elseif symmetry
            new_bins = _correct_symmetry(bins, reference[i])
            new_bins, d = _ordering_symmetry(new_bins, d)
        else
            new_bins = average_bins(bins)
        end

        push!(distances, new_bins)
        push!(densities, d)
    end

    return convert(Vector{Vector{Float64}}, distances), convert(Vector{Vector{Float64}}, densities)
end


"""
    averaging_profile(distances::Vector{Vector{Float64}}, densities::Vector{Vector{Float64}})

This function takes the `distances` and `densities` from Vector{Vector{Float64}} to Vector{Float64} by means and standard deviation computations.

# Arguments
- `distances::Vector{Vector{Float64}}`: a Vector{Vector{Float64}} with the distances from the center of the profile.
- `densities::Vector{Vector{Float64}}`: a Vector{Vector{Float64}} with the densities of the profile.
"""
function averaging_profile(distances::Vector{Vector{Float64}}, densities::Vector{Vector{Float64}})

    avg_distances = Statistics.mean(hcat(distances...), dims=2)[:,1]
    avg_densities = Statistics.mean(hcat(densities...), dims=2)[:,1]
    std_densities = Statistics.std(hcat(densities...), dims=2)[:,1]

    return avg_distances, avg_densities, std_densities
end


"""
    _ordering_symmetry(bins::Vector{Float64}, ρ::Vector{Float64})

This function checks if the bins are ordered after the symmetry protocol and if not, it orders them.

# Arguments
- `bins::Vector{Float64}`: a Vector{Float64} with the bins.
- `ρ::Vector{Float64}`: a Vector{Float64} with the densities.
"""
function _ordering_symmetry(bins::Vector{Float64}, ρ::Vector{Float64})
    if issorted(bins)
        return bins, ρ
    else
        sorted_idx = sortperm(bins)
        new_bins, new_ρ = bins[sorted_idx], ρ[sorted_idx]
        return new_bins, new_ρ
    end
end


"""
    average_bins(bins::Vector{Float64})

It takes an average vector based on the binned box parameters. So it can be rightly compare with the ρ density with ``length(density) ≡ length(bins) - 1``.
"""
function average_bins(bins::Vector{Float64})
    return [ 0.5*(bins[i]+bins[i+1]) for i in 1:(length(bins)-1) ]
end


"""
    _get_reference(simulation::MolSimToolkit.Simulation, selection::String)

It aims to get the reference center for each frame along the trajectory.

# Arguments
- `simulation::MolSimToolkit.Simulation`: a MolSimToolkit.Simulation object.
- `selection::String`: a selection string to get the reference center. For exemple, "all" or "not water".
"""
function _get_reference(simulation::MolSimToolkit.Simulation, selection::String)

    xyz = MolSimToolkit.Point3D[]

    reference_atoms = PDBTools.select(MolSimToolkit.atoms(simulation), selection)
    idx = PDBTools.index.(reference_atoms)

    for frame in simulation
        append!(xyz, mean(MolSimToolkit.positions(frame)[idx], dims=1))
    end

    return xyz
end


"""
    _correct_center(bins::Vector{Float64}, reference::Float64)

It corrects the bin values based on the referenced center. In this way, the center will become zero and all the bins values will be reajusted to apply it.
Now, every distance is compared with the central value.

# Arguments
- `bins::Vector{Float64}`: a Vector{Float64} with the bins.
- `reference::Float64`: a Float64 with the reference center.
"""
function _correct_center(bins::Vector{Float64}, reference::Float64)
    
    idx = findfirst(x -> x >= reference, bins)
    center = 0.5*(bins[idx]+bins[idx-1])
    centered_bins = average_bins(bins) .- center

    return centered_bins
end

"""
    _correct_symmetry(bins::Vector{Float64}, reference::Float64)

Do the same thing that `_correct_center(...)`, but it extends the concept to adjust the data in such way that the center will become the center of the box in every frame.

# Arguments
- `bins::Vector{Float64}`: a Vector{Float64} with the bins.
- `reference::Float64`: a Float64 with the reference center.
"""
function _correct_symmetry(bins::Vector{Float64}, reference::Float64)
    
    symmetric_bins = Float64[]

    min_threshold, max_threshold = extrema(bins)
    resolution = (max_threshold - min_threshold) / (length(bins)-1)
    ΔL = max_threshold - min_threshold - resolution

    centered_bins = _correct_center(bins, reference)
   
    for i in eachindex(centered_bins)
        if (sign(centered_bins[i]) == -1.) && (abs(centered_bins[i]) / (0.5*ΔL) > 1.)
            push!(symmetric_bins, centered_bins[i] + ΔL)
        elseif (sign(centered_bins[i]) == 1.) && (abs(centered_bins[i]) / (0.5*ΔL) > 1.)
            push!(symmetric_bins, centered_bins[i] - ΔL)
        else
            push!(symmetric_bins, centered_bins[i])
        end
    end

    return symmetric_bins
end


"""
    ρ(bins::Vector{Float64}, positions::Vector{Point3D{Float64}}, N::Float64; axis="z", prop=nothing)

Computes the density for each bin interval given the frame positions. Here, we apply a normalization base on the bin dimensions, so all the frames will have the same `sum(ρ)`.
The `prop` argument is used to compute the density of a specific property, like the atomic number, mass, or charge.

# Arguments
- `bins::Vector{Float64}`: a Vector{Float64} with the bin intervals.
- `positions::Vector{Point3D{Float64}}`: a Vector{Point3D{Float64}} with the atomic positions.
- `N::Float64`: a Float64 with the normalization factor.
- `axis::String`: a String with the axis flag.
- `prop::Vector{Float64}`: a Vector{Float64} with the property to be computed.
"""
function ρ(bins::Vector{Float64}, positions::Vector{MolSimToolkit.Point3D{Float64}}, N::Float64; axis="z", prop=nothing)

    if lowercase(string(axis)) == "x" || axis == 1
        coords = [ positions[i][1] for i in eachindex(positions) ]
    elseif lowercase(string(axis)) == "y" || axis == 2
        coords = [ positions[i][2] for i in eachindex(positions) ]
    elseif lowercase(string(axis)) == "z" || axis == 3
        coords = [ positions[i][3] for i in eachindex(positions) ]
    else
        println("Error: The axis flag must be related to the x, y, and z cartesian axes codification. For example, on axis z, it can be `z`, `Z` or `3`.")
    end

    ρ = Float64[]
    
    for b in eachindex(bins)
        if b == lastindex(bins); break; end
        coded = δ.(coords, bins[b], bins[b+1])
        ith_ρ = sum(coded .* prop) / N
        append!(ρ, ith_ρ)
    end

    return ρ
end

"""
    δ(value::Float64, lower::Float64, upper::Float64)

It is an indicator function that maps the elements inside ``x_lower ≤ x_value ≤ x_upper`` interval.
"""
function δ(value::Float64, lower::Float64, upper::Float64)
    if lower <= value <= upper
        return 1
    else
        return 0
    end
end

"""
    binning(simulation::MolSimToolkit.Simulation; axis="z", resolution=1.)

It reports the bins and normalization factor for the loaded simulation box.

# Arguments
- `simulation::MolSimToolkit.Simulation`: a MolSimToolkit.Simulation object.
- `axis::String`: a String with the axis flag.
- `dimensions::Int`: an integer with the number of dimensions. It can be 1 or 2.
- `resolution::Float64`: a Float64 with the resolution of the bins.
"""
function binning(simulation::MolSimToolkit.Simulation; axis="z", dimensions=1, resolution=1.)
    
    xmin, xmax = [ 0., 0., 0. ] , [ 0., 0., 0. ]

    for frame in simulation
        atomic_coordinates = MolSimToolkit.positions(frame)
        x = [ atom[1] for atom in atomic_coordinates ]; y = [ atom[2] for atom in atomic_coordinates ]; z = [ atom[3] for atom in atomic_coordinates ];
        xmin[1], xmax[1] = ifelse(minimum(x) < xmin[1], minimum(x), xmin[1]), ifelse(maximum(x) > xmax[1], maximum(x), xmax[1])
        xmin[2], xmax[2] = ifelse(minimum(y) < xmin[2], minimum(y), xmin[2]), ifelse(maximum(y) > xmax[2], maximum(y), xmax[2])
        xmin[3], xmax[3] = ifelse(minimum(z) < xmin[3], minimum(z), xmin[3]), ifelse(maximum(z) > xmax[3], maximum(z), xmax[3])
    end

    if lowercase(string(axis)) == "x" || axis == 1
        i_min = xmin[1]; i_max = xmax[1];
        L1 = xmax[2]-xmin[2]; L2 = xmax[3]-xmin[3];
    elseif lowercase(string(axis)) == "y" || axis == 2
        i_min = xmin[2]; i_max = xmax[2];
        L1 = xmax[1]-xmin[1]; L2 = xmax[3]-xmin[3];
    elseif lowercase(string(axis)) == "z" || axis == 3
        i_min = xmin[3]; i_max = xmax[3];
        L1 = xmax[1]-xmin[1]; L2 = xmax[2]-xmin[2];
    else
        println("Error: The axis flag must be related to the x, y, and z cartesian axes codification. For example, on axis z, it can be `z`, `Z` or `3`.")
    end

    return collect(i_min:resolution:i_max), Float64(resolution*L1*L2)
end

