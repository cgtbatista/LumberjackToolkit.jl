function electrons(pdbname::String; charges=true, psfname=nothing)
    @assert isfile(pdbname) "The PDB file does not exist."
    psfname = isnothing(psfname) ? replace(pdbname, ".pdb" => ".psf") : psfname
    @assert isfile(psfname) || (isnothing(psfname) && !charges) "The PSF file does not exist."
    el = PDBTools.atomic_number.(PDBTools.readPDB(pdbname))
    if charges
        q = chargesPSF(psfname)
        @assert length(el) == length(q) "The number of charges should be the same as the number of atoms."
        return el .- q
    else
        return el
    end
end

## sqrt((3V/4π)/m)
function density_grid(coords::Vector{SVector{3, Float64}}, resolution::Float64)
    xlower, xupper = extrema([ coord[1] for coord in coords ])
    ylower, yupper = extrema([ coord[2] for coord in coords ])
    zlower, zupper = extrema([ coord[3] for coord in coords ])
    points = SVector{3, Float64}[]
    for x in xlower-4.0:resolution:xupper+4.0, y in ylower-4.0:resolution:yupper+4.0, z in zlower:resolution:zupper
        push!(points, SVector(x, y, z))
    end
    return points
end

function density_grid(coords::MolSimToolkit.FramePositions, resolution::Float64)
    coords = [ SVector(coord) for coord in coords ]
    return density_grid(coords, resolution) 
end

function density_profile(
    pdbname::String, trjname::String, property::Vector{Float64};
    selection="not water", isdimensionless=true,  axis="z", cutoff=1.0, resolution=0.6, σ=1.0,
    cs=4.0, hascenter=true, hassymmetry=true, first=1, last=nothing, step=1
)
    @assert isnothing(last) || last > first "The last frame should be greater than the first frame."
    @assert isfile(pdbname) && isfile(trjname) "The PDB or trajectory file does not exist."
    @assert in(lowercase(axis), Set(["x", "y", "z", "xy", "xz", "yz", "yx", "zx", "zy"])) "This axis is not available."
    selection = isnothing(selection) || selection == "" ? "all" : selection
    atoms = PDBTools.select(PDBTools.readPDB(pdbname), selection)
    sim = MolSimToolkit.Simulation(atoms, trjname, first=first, last=last, step=step)
    return density_profile(
        sim, property, axis=axis, cutoff=cutoff, isdimensionless=isdimensionless, σ=σ, resolution=resolution, cs=cs, hascenter=hascenter, hassymmetry=hassymmetry
    )
end

function density_profile(
    sim::MolSimToolkit.Simulation, property::Vector{Float64}; isdimensionless=true,
    axis="z", cutoff=1.0,
    resolution=0.6, σ=1.0,          ## only for gaussian profile
    cs=4.0,
    hascenter=true, hassymmetry=true
)
    @assert in(lowercase(axis), Set(["x", "y", "z", "xy", "xz", "yz", "yx", "zx", "zy"])) "This axis should be cartesian (e.g. `x` or `xz`)."
    #@assert length(sim.atoms) == length(property) "The number of atoms should be the same as the number of properties."
    hascenter = hassymmetry ? true : hascenter
    if isdimensionless
        return density_profile_dimensionless(
            sim, property, axis=axis, cutoff=cutoff, cs=cs, hascenter=hascenter, hassymmetry=hassymmetry
        )
    else
        return density_profile_gaussian(
            sim, property, axis=axis, cutoff=cutoff, σ=σ, resolution=resolution, cs=cs, hascenter=hascenter, hassymmetry=hassymmetry
        )
    end
end

function density_profile_dimensionless(sim::MolSimToolkit.Simulation, property::Vector{Float64}; axis="z", cutoff=1.0, cs=1.0, hascenter=true, hassymmetry=true)
    if length(axis) == 1
        return dimensionless_profile_1D(sim, property, axis=axis, cutoff=cutoff, hascenter=hascenter, hassymmetry=hassymmetry)
    elseif length(axis) == 2
        return dimensionless_profile_2D(sim, property, axis=axis, cutoff=cutoff, cs=cs, hascenter=hascenter, hassymmetry=hassymmetry)
    else
        throw(ArgumentError("The axis should be a string with one or two characters (e.g. `z`, `x`, `xy`)."))
    end
end

function density_profile_gaussian(
    sim::MolSimToolkit.Simulation, property::Vector{Float64}; axis="z", resolution=0.5, cutoff=1.0, σ=1.0, cs=4.0, hascenter=true, hassymmetry=true)
    if length(axis) == 1
        return gaussian_profile_1D(sim, property, axis=axis, cutoff=cutoff, resolution=resolution, σ=σ, hascenter=hascenter, hassymmetry=hassymmetry)
    elseif length(axis) == 2
        return gaussian_profile_2D(sim, property, axis=axis, cutoff=cutoff, resolution=resolution, σ=σ, cs=cs, hascenter=hascenter, hassymmetry=hassymmetry)
    else
        throw(ArgumentError("The axis should be a string with one or two characters (e.g. `z`, `x`, `xy`)."))
    end
end

function dimensionless_profile_1D(sim::MolSimToolkit.Simulation, Q::Vector{Float64}; axis="z", cutoff=1.0, hascenter=true, hassymmetry=true)
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
            densities[i, iframe] = !isempty(idx) ? V \ sum(Q[iatoms][idx]) : 0.0
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

function dimensionless_profile_2D(
    sim::MolSimToolkit.Simulation, Q::Vector{Float64};
    axis="xz", cutoff=1.0, cs=4.0, hascenter=true, hassymmetry=true
)
    BIN1, BIN2, DENS = [], [], []
    axescode = Dict{String, Function}([
        ["xy", "yx"] .=> () -> [1,2]; ["xz", "zx"] .=> () -> [1,3]; ["yz", "zy"] .=> () -> [2,3]
    ])
    if !haskey(axescode, lowercase(axis))
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    dims = axis2dims(axis)
    dim = axis2dims(notaxis(axis))
    V, bin1, bin2 = binning(sim, axis=axis, cutoff=cutoff)
    cs = binning(sim, axis=notaxis(axis), cutoff=cs)[2]
    bins1 = Matrix{Float64}(undef, length(bin1)-1, length(sim))
    bins2 = Matrix{Float64}(undef, length(bin2)-1, length(sim))
    densities = Array{Float64, 3}(undef, length(bin1)-1, length(bin2)-1, length(sim))
    iatoms = PDBTools.index.(sim.atoms)
    for (iframe, frame) in enumerate(sim)
        xyz = [ SVector(ijk) for ijk in MolSimToolkit.positions(frame)[iatoms] ]
        for slice in 1:length(cs)-1
            coords = filter(coord -> cs[slice] <= coord[dim] <= cs[slice+1], xyz)
            for i in 1:length(bin1)-1, j in 1:length(bin2)-1
                idx = findall(
                    coord -> (bin1[i] <= coord[dims[1]] <= bin1[i+1]) && (bin2[j] <= coord[dims[2]] <= bin2[j+1]),
                    coords
                )
                densities[i,j,iframe] = !isempty(idx) ? V \ sum(Q[iatoms][idx]) : 0.0
            end
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
            push!(BIN1, bins1)
            push!(BIN2, bins2)
            push!(DENS, densities)
        end
    end
    return BIN1, BIN2, DENS
end

function gaussian_profile_1D(
    sim::MolSimToolkit.Simulation, Q::Vector{Float64};
    axis="z", cutoff=1.0, resolution=0.5, σ=1.0, hascenter=true, hassymmetry=true
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
        A = Q[iatoms] ./ (sqrt(2π) .* σ).^3
        α = inv.(-2σ.^2)
        ρ(r) = sum(
            A .* exp.(α .* [ sum((r .- coord).^2) for coord in coords ])
        )
        points = density_grid(coords, resolution)
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

function gaussian_profile_2D(
    sim::MolSimToolkit.Simulation, Q::Vector{Float64};
    axis="xz", cutoff=1.0, cs=4.0, resolution=0.5, σ=1.0, hascenter=true, hassymmetry=true
)
    BIN1, BIN2, DENS = [], [], []
    axescode = Dict{String, Function}([
        ["xy", "yx"] .=> () -> [1,2]; ["xz", "zx"] .=> () -> [1,3]; ["yz", "zy"] .=> () -> [2,3]
    ])
    if !haskey(axescode, lowercase(axis))
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    dims = axis2dims(axis)
    dim = axis2dims(notaxis(axis))
    V, bin1, bin2 = binning(sim, axis=axis, cutoff=cutoff)
    cs = binning(sim, axis=notaxis(axis), cutoff=cs)[2]
    bins1 = Matrix{Float64}(undef, length(bin1)-1, length(sim))
    bins2 = Matrix{Float64}(undef, length(bin2)-1, length(sim))
    densities = zeros(Float64, length(bin1)-1, length(bin2)-1, length(sim))
    iatoms = PDBTools.index.(sim.atoms)
    for (iframe, frame) in enumerate(sim)
        xyz = MolSimToolkit.positions(frame)[iatoms]
        σ = σ .* ones(Float64, length(iatoms))
        A = Q[iatoms] ./ (sqrt(2π) .* σ).^3
        α = inv.(-2 .* σ.^ 2)
        ρ(r) = sum(
            A .* exp.(α .* [ sum((r .- ijk).^2) for ijk in xyz ])
        )
        for slice in 1:length(cs)-1
            coords = filter(coord -> cs[slice] <= coord[dim] <= cs[slice+1], xyz)
            if isempty(coords)
                continue
            end
            points = density_grid(coords, resolution)
            density = ρ.(points)
            ibin = [ searchsortedfirst(bin1, point[dims[1]]) for point in points ]
            jbin = [ searchsortedfirst(bin2, point[dims[2]]) for point in points ]
            @inbounds for p in 1:length(points)
                i, j = ibin[p], jbin[p]
                if 1 <= i <= size(densities, 1) && 1 <= j <= size(densities, 2)
                    densities[i, j, iframe] += density[p]
                end
            end
            densities[:,:,iframe] ./= V
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
            push!(BIN1, deepcopy(bins1))
            push!(BIN2, deepcopy(bins2))
            push!(DENS, deepcopy(densities))
            densities .= 0.0
        end
    end
    return BIN1, BIN2, DENS
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

function notaxis(axis::String)
    axis = split(lowercase(axis), "")
    notaxis = filter(i -> !in(i, axis), ["x", "y", "z"])
    return join(notaxis)
end

function axis2dims(axis::String)
    axis = split(lowercase(axis), "")
    if length(axis) == 1
        return Dict{String, Int}(
            "x" => 1,
            "y" => 2,
            "z" => 3
        )[axis[1]]
    elseif length(axis) == 2
        return Dict{String, Function}([
            ["xy", "yx"] .=> () -> [1, 2];
            ["xz", "zx"] .=> () -> [1, 3];
            ["yz", "zy"] .=> () -> [2, 3]
        ])[join(axis)]()
    end
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

function mapping_density(density::Array{Float64, 3}; tol=0.005)
    checklist = zeros(Bool, size(density, 1), size(density, 2), size(density, 3))
    for iframe in 1:size(density, 3)
        for i in 1:size(density, 1), j in 1:size(density, 2)
            checklist[i,j,iframe] = density[i,j,iframe] > tol
        end
    end
    return BitArray(checklist)
end

function fibrilwidth(
    bin1::Matrix{Float64},
    bin2::Matrix{Float64},
    densities::Array{Float64, 3};
    step=1.0,
    tol=0.005
)
    valid_regions = mapping_density(densities, tol=tol)   
    n_frames = size(densities, 3)
    angles = collect(0:step:180-step)
    results = Matrix{Float64}(undef, n_frames, length(angles))
    for frame in 1:n_frames
        x = bin1[:, frame]
        y = bin2[:, frame]
        mask = valid_regions[:, :, frame]
        xc = (maximum(x) + minimum(x)) / 2
        yc = (maximum(y) + minimum(y)) / 2
        sorted_x, idx_x = sortperm_with_values(x)
        sorted_y, idx_y = sortperm_with_values(y)
        dr = minimum([minimum(diff(sorted_x)), minimum(diff(sorted_y))]) / 10
        for (i, θ) in enumerate(angles)
            dir1 = deg2rad(θ)
            dir2 = dir1 + π
            r1 = radial_scan(xc, yc, dir1, sorted_x, sorted_y, idx_x, idx_y, mask, dr)
            r2 = radial_scan(xc, yc, dir2, sorted_x, sorted_y, idx_x, idx_y, mask, dr)
            results[frame, i] = r1 + r2
        end
    end
    return results
end

function sortperm_with_values(v)
    perm = sortperm(v)
    (v[perm], perm)
end

function radial_scan(xc, yc, θ, sorted_x, sorted_y, idx_x, idx_y, mask, dr)
    t = 0.0
    max_r = 0.0
    while true
        x = xc + t * cos(θ)
        y = yc + t * sin(θ)
        (x < sorted_x[1] || x > sorted_x[end] || 
         y < sorted_y[1] || y > sorted_y[end]) && break
        i = find_nearest_index(sorted_x, x)
        j = find_nearest_index(sorted_y, y)
        mask[idx_x[i], idx_y[j]] || break
        max_r = t
        t += dr
    end
    return max_r
end

function find_nearest_index(sorted_v, value)
    idx = searchsortedfirst(sorted_v, value)
    clamp(idx - (abs(sorted_v[idx] - value) > abs(sorted_v[idx-1] - value)), 1, length(sorted_v))
end

function interpol(x::Vector{Float64}, y::Vector{Float64}, z::Matrix{Float64}; bins=200)
    itp = Interpolations.linear_interpolation((x, y), z, extrapolation_bc=NaN)
    xnew = range(minimum(x), maximum(x), length=bins)
    ynew = range(minimum(y), maximum(y), length=bins)
    znew = [ itp(i,j) for i in xnew, j in ynew ]
    return xnew, ynew, znew
end
