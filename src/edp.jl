using PDBTools, MolSimToolkit, StatsPlots, Plots

function densityprofile(
        pdb::String,
        trajectory::String;
        profile="electron", selection="all", axis="z", resolution=1, center=true, symmetry=true, symmetry_selection="not water", first=1, step=1, last=nothing,
        charge_column=nothing, edp=nothing, charge_correction=true
    )
    
    distances, densities = [], []

    println("DP -- $profile density distribution of the molecules inside the box:")
    println("   the selection is `$selection`")
    println("   resolution on $axis is equal to $resolution Å")
    println("")

    simulation = MolSimToolkit.Simulation(pdb, trajectory; first=first, last=last, step=step)
    selected_atoms = PDBTools.select(MolSimToolkit.atoms(simulation), selection)
    idx = PDBTools.index.(selected_atoms)

    bins, V_norm = _unidimensional_binning(simulation, axis=axis, resolution=resolution)
    if symmetry || center;
        references = _get_reference(simulation, symmetry_selection)
        if lowercase(string(axis)) == "x" || axis == 1
            reference = [ ref[1] for ref in references ]
        elseif lowercase(string(axis)) == "y" || axis == 2
            reference = [ ref[2] for ref in references ]
        elseif lowercase(string(axis)) == "z" || axis == 3
            reference = [ ref[3] for ref in references ]
        else
            println("Error: The axis flag must be related to the x, y, and z cartesian axes codification. For example, on axis z, it can be `z`, `Z` or `3`.")
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
        coords = positions(frame)[idx]
        new_coords = [ Point3D(coords[n][1], coords[n][2], coords[n][3]) for n in eachindex(coords) ]

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

function averaging_profile(distances::Vector{Vector{Float64}}, densities::Vector{Vector{Float64}})

    avg_distances = Statistics.mean(hcat(distances...), dims=2)[:,1]
    avg_densities = Statistics.mean(hcat(densities...), dims=2)[:,1]
    std_densities = Statistics.std(hcat(densities...), dims=2)[:,1]

    return avg_distances, avg_densities, std_densities
end

function averaging_plotting(distances::Vector{Vector{Float64}}, densities::Vector{Vector{Float64}}; profile=nothing, type=nothing)
    
    avg_distances, avg_densities, std_densities = averaging_profile(distances, densities)
    
    if profile == "number"
        y = "number density (molecules/Å³)"
    elseif profile == "mass"
        y = "mass density (Da/Å³)"
    elseif profile == "charge"
        y = "charge density (e/Å³)"
    elseif isnothing(profile) || profile == "electron"
        y = "electron density (e/Å³)"
    end

    plotting = plot(avg_distances, avg_densities, label=:none, xlabel="distance (Å)", ylabel=y, linewidth=2)
    if type == "error"
        plotting = plot!(avg_distances, avg_densities, yerr=std_densities, marker=:circle)
    elseif type == "box"
        merging_densities = Float64[]
        for i in eachindex(densities)
            append!(merging_densities, densities[i])
        end
        merging_distances = Float64[]
        for i in eachindex(distances)
            append!(merging_distances, distances[i])
        end
        plotting = boxplot!(merging_distances, merging_densities, label=:none, outliers=false, bar_width=0.25)
    end

    return plotting
end

function _ordering_symmetry(bins::Vector{Float64}, ρ::Vector{Float64})
    if issorted(bins)
        return bins, ρ
    else
        sorted_idx = sortperm(bins)
        new_bins, new_ρ = bins[sorted_idx], ρ[sorted_idx]
        return new_bins, new_ρ
    end
end

function average_bins(bins::Vector{Float64})
    return [ 0.5*(bins[i]+bins[i+1]) for i in 1:(length(bins)-1) ]
end

function _get_reference(simulation::Simulation, selection::String)

    xyz = Point3D[]

    reference_atoms = PDBTools.select(MolSimToolkit.atoms(simulation), selection)
    idx = PDBTools.index.(reference_atoms)

    for frame in simulation
        append!(xyz, mean(positions(frame)[idx], dims=1))
    end

    return xyz
end

function _correct_center(bins::Vector{Float64}, reference::Float64)
    
    idx = findfirst(x -> x >= reference, bins)
    center = 0.5*(bins[idx]+bins[idx-1])
    centered_bins = average_bins(bins) .- center

    return centered_bins
end

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

function ρ(bins::Vector{Float64}, positions::Vector{Point3D{Float64}}, Nfactor::Float64; axis="z", prop=nothing)

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
        δ = δ_bin.(coords, bins[b], bins[b+1])
        ith_ρ = sum(δ .* prop) / Nfactor
        append!(ρ, ith_ρ)
    end

    return ρ
end

function δ_bin(value::Float64, lower::Float64, upper::Float64)
    if lower <= value <= upper
        return 1
    else
        return 0
    end
end

function _unidimensional_binning(simulation::Simulation; axis="z", resolution=1)
    
    xmin, xmax = [ 0., 0., 0. ] , [ 0., 0., 0. ]

    for frame in simulation
        atomic_coordinates = positions(frame)
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


function vmd_get_charges(rawfile::String; newpdb="new.pdb", pdb_column="beta", vmd="vmd")
    
    psf = rawfile * ".psf"
    pdb = rawfile * ".pdb"

    vmdinput_file = tempname() * ".tcl"
    
    vmdinput = open(vmdinput_file, "w")
    Base.write(vmdinput, "mol new \"$psf\" \n")
    Base.write(vmdinput, "mol addfile \"$pdb\" \n")
    Base.write(vmdinput, "set sel [ atomselect top \"all\" ] \n")
    Base.write(vmdinput, "\$sel set $pdb_column [\$sel get charge] \n")
    Base.write(vmdinput, "\$sel writepdb $newpdb \n")
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "exit")
    Base.close(vmdinput)
    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")

    Base.rm(vmdinput_file)
    
    return vmdoutput
end