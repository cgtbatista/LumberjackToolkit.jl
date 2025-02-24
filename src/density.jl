function edp(
    pdbname::String, trjname::String; first=1, last=nothing, step=1,
    psfname=nothing, selection="not water", isdimensionless=true,
    axis="z", cutoff=4.0, resolution=0.6, e=nothing, σ=1.0,
    hascenter=true, hassymmetry=true
)
    if !in(lowercase(axis), Set(["x", "y", "z", "xy", "xz", "yz", "yx", "zx", "zy"]))
        throw(ArgumentError("""
        The axis should be associated with the labels of cartesian axes or plane.
        Please insert it such as: `x`, `y`, `z`, `xy`, `xz`, or `yz`.
        """))
    end
    atoms = PDBTools.select(PDBTools.readPDB(pdbname), selection)
    sim = MolSimToolkit.Simulation(atoms, trjname, first=first, last=last, step=step)
    psfname = isnothing(psfname) ? replace(pdbname, ".pdb" => ".psf") : psfname
    #e = isnothing(e) ? chargesPSF(psfname)[PDBTools.index.(atoms)] : e
    e = isnothing(e) ? abs.(chargesPSF(psfname)) : e
    hascenter = hassymmetry ? true : hascenter
    if length(axis) == 1
        if isdimensionless
            return edp_1d_dimensionless(
                sim, e,
                axis=axis, cutoff=cutoff, hascenter=hascenter, hassymmetry=hassymmetry
            )
        else
            return edp_1d_gaussian(
                sim, e,
                axis=axis, cutoff=cutoff, σ=σ, hascenter=hascenter, hassymmetry=hassymmetry
            )
        end
    elseif length(axis) == 2
        if isdimensionless
            return edp_2d_dimensionless(
                sim, e,
                axis=axis, cutoff=cutoff, hascenter=hascenter, hassymmetry=hassymmetry
            )
        else
            return edp_2d_gaussian(
                sim, e,
                axis=axis, cutoff=cutoff, σ=σ, hascenter=hascenter, hassymmetry=hassymmetry
            )
        end
    else
        throw(ArgumentError("The axis should be a string with one or two characters."))
    end
end

function edp_1d_dimensionless(
    sim::MolSimToolkit.Simulation, electrons::Vector{Float64};
    axis="z", cutoff=1.0, hascenter=true, hassymmetry=true
)
    axescode = Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3)
    if haskey(axescode, lowercase(axis))
        dim = axescode[axis]
        V, bin = binning(sim, axis=axis, cutoff=cutoff)
        bins = Matrix{Float64}(undef, length(bin)-1, length(sim))
        densities = deepcopy(bins)
    else
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    iatoms = PDBTools.index.(sim.atoms)
    for (iframe, frame) in enumerate(sim)
        coords = [ coord[dim] for coord in MolSimToolkit.positions(frame)[iatoms] ]
        for i in 1:length(bin)-1
            idx = findall(coord -> bin[i] <= coord <= bin[i+1], coords)
            densities[i, iframe] = !isempty(idx) ? V \ sum(electrons[iatoms][idx]) : 0.0
        end
        if !hascenter
            bins[:,iframe] = avgbins(bin)
        else
            avg = mean(coords)
            if hassymmetry
                bins[:,iframe] = _fix_binning_symmetry(bin, avg)
                if !issorted(bins[:,iframe])
                    isorted = sortperm(bins[:,iframe])
                    bins[:,iframe], densities[:,iframe] = bins[isorted,iframe], densities[isorted,iframe]
                end
            else
                bins[:,iframe] = _fix_binning_center(bin, avg)
            end
        end
    end
    return bins, densities
end

function edp_2d_dimensionless(
    sim::MolSimToolkit.Simulation, electrons::Vector{Float64};
    axis="xz", cutoff=1.0, hascenter=true, hassymmetry=true
)
    axescode = Dict{String, Function}([
        ["xy", "yx"] .=> () -> [1,2]; ["xz", "zx"] .=> () -> [1,3]; ["yz", "zy"] .=> () -> [2,3]
    ])
    if haskey(axescode, lowercase(axis))
        dims = axescode[axis]()
        V, bin1, bin2 = binning(sim, axis=axis, cutoff=cutoff)
        bins1 = Matrix{Float64}(undef, length(bin1)-1, length(sim))
        bins2 = Matrix{Float64}(undef, length(bin2)-1, length(sim))
        densities = Array{Float64, 3}(undef, length(bin1)-1, length(bin2)-1, length(sim))
        #density = Matrix{Float64}(undef, length(bin1)-1, length(bin2)-1)
    else
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    iatoms = PDBTools.index.(sim.atoms)
    #densities = Vector{Matrix{Float64}}(undef, length(sim))
    for (iframe, frame) in enumerate(sim)
        coords = [ coord[dims] for coord in MolSimToolkit.positions(frame)[iatoms] ]
        for i in 1:length(bin1)-1, j in 1:length(bin2)-1
            idx = findall(
                coord -> (bin1[i] <= coord[1] <= bin1[i+1]) && (bin2[j] <= coord[2] <= bin2[j+1]),
                coords
            )
            densities[i,j,iframe] = !isempty(idx) ? V \ sum(electrons[iatoms][idx]) : 0.0
        end
        #densities[iframe] = deepcopy(density)
        if !hascenter
            bins1[:,iframe], bins2[:,iframe] = avgbins(bin1), avgbins(bin2)
        else
            avg = mean(coords)
            if hassymmetry
                bins1[:,iframe] = _fix_binning_symmetry(bin1, avg[1])
                bins2[:,iframe] = _fix_binning_symmetry(bin2, avg[2])
                if !issorted(bins1[:,iframe]) || !issorted(bins2[:,iframe])
                    isorted, jsorted = sortperm(bins1[:,iframe]), sortperm(bins2[:,iframe])
                    bins1[:,iframe], bins2[:,iframe] = bins1[isorted,iframe], bins2[jsorted,iframe]
                    densities[:,:,iframe] = densities[isorted,jsorted,iframe]
                end
            else
                bins1[:,iframe] = _fix_binning_center(bin1, avg[1])
                bins2[:,iframe] = _fix_binning_center(bin2, avg[2])
            end
        end
    end
    return bins1, bins2, densities
end

## sqrt((3V/4π)/m)
function density_grid(coords::Vector{SVector{3, Float64}}, resolution::Float64)
    xlower, xupper = extrema([ coord[1] for coord in coords ])
    ylower, yupper = extrema([ coord[2] for coord in coords ])
    zlower, zupper = extrema([ coord[3] for coord in coords ])
    points = SVector{3, Float64}[]
    for x in xlower:resolution:xupper, y in ylower:resolution:yupper, z in zlower:resolution:zupper
        push!(points, SVector(x, y, z))
    end
    return points
end

function density_grid(coords::MolSimToolkit.FramePositions, resolution::Float64)
    coords = [ SVector(coord) for coord in coords ]
    return density_grid(coords, resolution) 
end

function edp_1d_gaussian(
    sim::MolSimToolkit.Simulation, electrons::Vector{Float64};
    axis="z", cutoff=1.0, σ=1.0, hascenter=true, hassymmetry=true
)
    axescode = Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3)
    if haskey(axescode, lowercase(axis))
        dim = axescode[axis]
        V, bin = binning(sim, axis=axis, cutoff=cutoff)
        bins = Matrix{Float64}(undef, length(bin)-1, length(sim))
        densities = deepcopy(bins)
    else
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    iatoms = PDBTools.index.(sim.atoms)
    for (iframe, frame) in enumerate(sim)
        coords = MolSimToolkit.positions(frame)[iatoms]
        σ = σ .* ones(Float64, length(iatoms))
        A = electrons[iatoms] ./ (sqrt(2π) .* σ).^3
        α = inv.(-2σ.^2)
        ρ(r) = sum(
            A .* exp.(α .* [ sum((r .- coord).^2) for coord in coords ])
        )
        points = density_grid(coords, cutoff)
        density = ρ.(points)
        for i in 1:length(bin)-1
            idx = findall(
                coord -> bin[i] <= coord[dim] <= bin[i+1],
                points
            )
            densities[i,iframe] = !isempty(idx) ? V \ sum(density[idx]) : 0.0
        end
        if !hascenter
            bins[:,iframe] = avgbins(bin)
        else
            avg = mean([ coord[dim] for coord in coords ])
            if hassymmetry
                bins[:,iframe] = _fix_binning_symmetry(bin, avg)
                if !issorted(bins[:,iframe])
                    isorted = sortperm(bins[:,iframe])
                    bins[:,iframe], densities[:,iframe] = bins[isorted,iframe], densities[isorted,iframe]
                end
            else
                bins[:,iframe] = _fix_binning_center(bin, avg)
            end
        end
    end
    return bins, densities
end

function edp_2d_gaussian(
    sim::MolSimToolkit.Simulation, electrons::Vector{Float64};
    axis="xz", cutoff=1.0, σ=1.0, hascenter=true, hassymmetry=true
)
    axescode = Dict{String, Function}([
        ["xy", "yx"] .=> () -> [1,2]; ["xz", "zx"] .=> () -> [1,3]; ["yz", "zy"] .=> () -> [2,3]
    ])
    if haskey(axescode, lowercase(axis))
        dims = axescode[axis]()
        V, bin1, bin2 = binning(sim, axis=axis, cutoff=cutoff)
        bins1 = Matrix{Float64}(undef, length(bin1)-1, length(sim))
        bins2 = Matrix{Float64}(undef, length(bin2)-1, length(sim))
        densities = zeros(Float64, length(bin1)-1, length(bin2)-1, length(sim))
    else
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    iatoms = PDBTools.index.(sim.atoms)
    for (iframe, frame) in enumerate(sim)
        coords = MolSimToolkit.positions(frame)[iatoms]
        σ = σ .* ones(Float64, length(iatoms))
        A = electrons[iatoms] ./ (sqrt(2π) .* σ).^3
        α = inv.(-2σ.^2)
        ρ(r) = sum(
            A .* exp.(α .* [ sum((r .- coord).^2) for coord in coords ])
        )
        points = density_grid(coords, cutoff)
        density = ρ.(points)
        #fill!(densities[:,:,iframe], 0.0)
        ibin = [ searchsortedfirst(bin1, point[dims[1]]) for point in points ]
        jbin = [ searchsortedfirst(bin2, point[dims[2]]) for point in points ]
        @inbounds for p in 1:length(points)
            i, j = ibin[p], jbin[p]
            if 1 <= i <= size(densities, 1) && 1 <= j <= size(densities, 2)
                densities[i, j, iframe] += density[p]
            end
        end
        densities[:,:,iframe] ./= V
        # for i in 1:length(bin1)-1, j in 1:length(bin2)-1
        #     idx = findall(
        #         coord -> (bin1[i] <= coord[1] <= bin1[i+1]) && (bin2[j] <= coord[2] <= bin2[j+1]),
        #         points
        #     )
        #     densities[i,j,iframe] = !isempty(idx) ? V \ sum(density[idx]) : 0.0
        # end
        if !hascenter
            bins1[:,iframe], bins2[:,iframe] = avgbins(bin1), avgbins(bin2)
        else
            avg = mean([ coord[dims] for coord in coords ])
            if hassymmetry
                bins1[:,iframe] = _fix_binning_symmetry(bin1, avg[1])
                bins2[:,iframe] = _fix_binning_symmetry(bin2, avg[2])
                if !issorted(bins1[:,iframe]) || !issorted(bins2[:,iframe])
                    isorted, jsorted = sortperm(bins1[:,iframe]), sortperm(bins2[:,iframe])
                    bins1[:,iframe], bins2[:,iframe] = bins1[isorted,iframe], bins2[jsorted,iframe]
                    densities[:,:,iframe] = densities[isorted,jsorted,iframe]
                end
            else
                bins1[:,iframe] = _fix_binning_center(bin1, avg[1])
                bins2[:,iframe] = _fix_binning_center(bin2, avg[2])
            end
        end
    end
    return bins1, bins2, densities
end

function avgbins(bins::Vector{Float64})
    return 0.5 .* (bins[1:end-1] + bins[2:end])
end

function _fix_binning_center(bins::Vector{Float64}, center::Float64)
    return avgbins(bins .- center)
end

function _fix_binning_symmetry(bins::Vector{Float64}, center::Float64)
    lower, upper = extrema(bins)
    ΔL = upper - lower
    avgbins = _fix_binning_center(bins, center)
    return mod.(avgbins .+ 0.5*ΔL, ΔL) .- 0.5*ΔL
end

function binning(avgbins::Vector{Float64})
    bins = Vector{Float64}(undef, length(avgbins)+1)
    Δbin = avgbins[2] - avgbins[1]
    for (i, bin) in enumerate(avgbins)
        bins[i] = bin - 0.5*Δbin
    end
    bins[end] = avgbins[end] + 0.5*Δbin
    return bins
end

function binning(sim::MolSimToolkit.Simulation; axis="z", cutoff=1.0)
    lower, upper = zeros(3), zeros(3)
    for frame in sim
        xyz = MolSimToolkit.positions(frame)
        if sim.frame_index == 1
            lower[1], upper[1] = extrema([ ijk[1] for ijk in xyz ])
            lower[2], upper[2] = extrema([ ijk[2] for ijk in xyz ])
            lower[3], upper[3] = extrema([ ijk[3] for ijk in xyz ])
        else
            lower[1], upper[1] = min(lower[1], minimum([ ijk[1] for ijk in xyz ])), max(upper[1], maximum([ ijk[1] for ijk in xyz ]))
            lower[2], upper[2] = min(lower[2], minimum([ ijk[2] for ijk in xyz ])), max(upper[2], maximum([ ijk[2] for ijk in xyz ]))
            lower[3], upper[3] = min(lower[3], minimum([ ijk[3] for ijk in xyz ])), max(upper[3], maximum([ ijk[3] for ijk in xyz ]))
        end
    end
    return binspecs(axis, lower=lower, upper=upper, cutoff=cutoff)
end

function binspecs(axis::String; lower=zeros(3), upper=ones(3), cutoff=1.0)
    axis = split(lowercase(axis), "")
    if length(axis) == 1
        dict = Dict{String, Function}(
            "x" => () -> (cutoff*(upper[2]-lower[2])*(upper[3]-lower[3]), collect(lower[1]:cutoff:upper[1])),
            "y" => () -> (cutoff*(upper[1]-lower[1])*(upper[3]-lower[3]), collect(lower[2]:cutoff:upper[2])),
            "z" => () -> (cutoff*(upper[1]-lower[1])*(upper[2]-lower[2]), collect(lower[3]:cutoff:upper[3]))
        )
        return dict[axis[1]]()
    elseif length(axis) == 2
        dict = Dict{Tuple{String, String}, Function}([
            [("x", "y"), ("y", "x")] .=> () -> (cutoff*(upper[3]-lower[3]), collect(lower[1]:cutoff:upper[1]), collect(lower[2]:cutoff:upper[2]));
            [("x", "z"), ("z", "x")] .=> () -> (cutoff*(upper[2]-lower[2]), collect(lower[1]:cutoff:upper[1]), collect(lower[3]:cutoff:upper[3]));
            [("y", "z"), ("z", "y")] .=> () -> (cutoff*(upper[1]-lower[1]), collect(lower[2]:cutoff:upper[2]), collect(lower[3]:cutoff:upper[3]))
        ])
        return dict[(axis[1], axis[2])]()
    end
end


function fibrilwidth(
    bin1::Matrix{Float64},
    bin2::Matrix{Float64},
    densities::Array{Float64, 3},
    step=1.0, tol=0.005
)   
    ϕ = collect(0:step:180-step)

    d = Matrix{Float64}(undef, length(ϕ), size(densities, 3))
    for iframe in 1:size(densities, 3)
        x_coords = bin1[:, iframe]
        y_coords = bin2[:, iframe]
        densities_frame = densities[:, :, iframe]
        
        x_min, x_max = extrema(x_coords)
        y_min, y_max = extrema(y_coords)
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        
        sorted_x_perm = sortperm(x_coords)
        sorted_x = x_coords[sorted_x_perm]
        sorted_y_perm = sortperm(y_coords)
        sorted_y = y_coords[sorted_y_perm]
        
        dx = length(sorted_x) > 1 ? minimum(diff(sorted_x)) : Inf
        dy = length(sorted_y) > 1 ? minimum(diff(sorted_y)) : Inf
        dr = min(dx, dy) / 10.0
        dr = isfinite(dr) ? dr : 0.1

        function radii(theta_rad)
            max_r = 0.0
            t = 0.0
            while true
                x = x_center + t * cos(theta_rad)
                y = y_center + t * sin(theta_rad)
                
                if x < x_min || x > x_max || y < y_min || y > y_max
                    break
                end
                
                idx_x_sorted = searchsortednearest(sorted_x, x)
                idx_x = sorted_x_perm[idx_x_sorted]
                idx_y_sorted = searchsortednearest(sorted_y, y)
                idx_y = sorted_y_perm[idx_y_sorted]
                
                if densities_frame[idx_x, idx_y] < tol
                    break
                end
                
                max_r = t
                t += dr
            end
            return max_r
        end

        for (i, ϕi) in enumerate(ϕ)
            ϕi = deg2rad(ϕi)
            d[i, iframe] = sum(
                radii.([ϕi, ϕi + π])
            )
        end
    end
    return d
end