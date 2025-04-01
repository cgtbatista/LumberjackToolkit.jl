const _bglc_charmm_partialcharges = Dict{String, Float64}([
    # carbons                                ##            O6-HO6                 
    [ "C1" ]                    .=>  0.340;  ##            |                
    [ "C2", "C3", "C4" ]        .=>  0.140;  ##        H61-C6-H62
    [ "C5" ]                    .=>  0.110;  ##            |
    [ "C6" ]                    .=>  0.050;  ##         H5-C5---O5
    # oxygens                                ##      H4   /       \    O1-HO1
    [ "O1", "O2", "O3", "O4",                ##        \ / HO3     \  /
      "O6" ]                    .=> -0.650;  ##         C4 |        C1
    [ "O5" ]                    .=> -0.400;  ##        / \ O3   H2 /  \
    # hydrogens                              ##  HO4-O4   \|    | /    H1
    [ "H1", "H2", "H3", "H4",                ##            C3---C2
      "H5", "H61", "H62" ]      .=>  0.090;  ##            |    |
    [ "HO1", "HO2", "HO3",                   ##            H3   O2-HO2
      "HO4", "HO6" ]            .=>  0.420;  ##                             β-D-Glucose – CHARMM36 Partial Charges
])

function charmm_partialcharges(atom::String, residue::String)
    if in(residue, ["BGLC", "BGC", "GLC"])
        if haskey(_bglc_charmm_partialcharges, atom)
            return round(_bglc_charmm_partialcharges[atom], digits=3)
        end
    end
    return nothing
end

function charges(pdbname::String; selection="all")
    if !isfile(pdbname) || filesize(pdbname) == 0
        throw(ArgumentError("The $pdbname file does not exist or is empty."))
    end
    atoms = PDBTools.readPDB(pdbname, selection)
    return charges(atoms)
end

function charges(atoms::Vector{<:PDBTools.Atom})
    q = zeros(Float64, length(atoms))
    for (i, at) in enumerate(atoms)
        q[i] = charmm_partialcharges(
            PDBTools.name(at),
            PDBTools.resname(at)
        )
    end
    return q
end

function electrons(pdbname::String; selection="all", qcorrection=true)
    if !isfile(pdbname) || filesize(pdbname) == 0
        throw(ArgumentError("The $pdbname file does not exist or is empty."))
    end
    e = PDBTools.atomic_number.(
        PDBTools.readPDB(pdbname, selection)
    )
    if qcorrection
        q = charges(pdbname, selection=selection)
        return round.(e .- q, digits=3)
    end
    return round.(e, digits=3)
end

function electrons(atoms::Vector{<:PDBTools.Atom}; qcorrection=true)
    e = PDBTools.atomic_number.(atoms)
    if qcorrection
        q = charges(atoms)
        return round.(e .- q, digits=3)
    end
    return round.(e, digits=3)
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

function binning(
    simulation::MolSimToolkit.Simulation;
    axis="z",
    binwidth::Float64=1.0,
    margin::Float64=0.0
)
    isempty(axis) && throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    axis = lowercase(axis)
    if !in(axis, Set(["x", "y", "z", "xy", "xz", "yz", "yx", "zx", "zy"]))
        throw(ArgumentError("The axis $axis is not valid. Please use x, y, z, xy, xz, yz, yx, zx or zy."))
    end
    lower, upper = extremacoordinates(simulation)
    lower .-= margin
    upper .+= margin
    if length(axis) == 1
        dim = findfirst(isequal(axis), ["x", "y", "z"])
        return binspecs1D(dim, lower, upper, binwidth)
    end
    if length(axis) == 2
        dims = [ findfirst(isequal(string(code)), ["x", "y", "z"]) for code in axis ]
        dims = sort(dims)
        return binspecs2D(dims[1], dims[2], lower, upper, binwidth)
    end
end

function binspecs1D(dim::Int, lower::AbstractVector, upper::AbstractVector, binwidth::Float64)
    edges = range(lower[dim], upper[dim], step=binwidth)
    @fastmath i, j = dim%3 + 1, (dim+1)%3 + 1
    V = vcat(upper[j] - lower[j], upper[i] - lower[i], binwidth)
    V = prod(V)
    return V, collect(edges)
end

function binspecs2D(dim1::Int, dim2::Int, lower::AbstractVector, upper::AbstractVector, binwidth::Float64)
    edges1 = range(lower[dim1], upper[dim1], step=binwidth)
    edges2 = range(lower[dim2], upper[dim2], step=binwidth)
    @fastmath k = 6 - dim1 - dim2
    V = vcat(upper[k] - lower[k], binwidth^2)
    V = prod(V)
    return V, collect(edges1), collect(edges2)
end

@inline function extremacoordinates(simulation::MolSimToolkit.Simulation)
    lower, upper = fill(Inf, 3), fill(-Inf, 3)
    @inbounds for frame in simulation
        coords = MolSimToolkit.positions(frame)
        mincoord, maxcoord = extremacoordinates(coords)
        lower .= min.(lower, mincoord)
        upper .= max.(upper, maxcoord)
    end
    return lower, upper
end

@inline function extremacoordinates(coords::MolSimToolkit.FramePositions)
    lower, upper = zeros(3), zeros(3)
    @inbounds @simd for dim in 1:3
        lower[dim], upper[dim] = extrema(coord[dim] for coord in coords)
    end
    return lower, upper
end

function interpol(x::AbstractVector{T}, y::AbstractVector{T}, z::Matrix{T}; bins=200) where T
    itp = Interpolations.linear_interpolation((x, y), z, extrapolation_bc=NaN)
    xnew = range(minimum(x), maximum(x), length=bins)
    ynew = range(minimum(y), maximum(y), length=bins)
    znew = [ itp(i,j) for i in xnew, j in ynew ]
    return xnew, ynew, znew
end

function gridpoints(
    coords::MolSimToolkit.FramePositions,
    resolution::Float64;
    margin::Float64=0.0
)
    xlower, xupper = extrema(coord[1] for coord in coords) .+ (-1, 1) .* margin
    ylower, yupper = extrema(coord[2] for coord in coords) .+ (-1, 1) .* margin
    zlower, zupper = extrema(coord[3] for coord in coords) .+ (-1, 1) .* margin
    points = [
        SVector(x, y, z) for
            x in xlower:resolution:xupper,
            y in ylower:resolution:yupper,
            z in zlower:resolution:zupper
    ]
    return points
end

function gridpoints(
    sim::MolSimToolkit.Simulation,
    resolution::Float64;
    margin::Float64=0.0
)
    lower, upper = extremacoordinates(sim)
    lower .-= margin
    upper .+= margin
    points = [
        SVector(x, y, z) for
            x in lower[1]:resolution:upper[1],
            y in lower[2]:resolution:upper[2],
            z in lower[3]:resolution:upper[3]
    ]
    return points
end

function density_profile(
    sim::MolSimToolkit.Simulation;
    axis="z",
    binwidth=0.5,
    cutoff=0.05,
    resolution=0.25,
    σ=1.0,
    cs=4.0,
    margin=2.0,
    hascenter=true,
    hassymmetry=true
)
    if !in(lowercase(axis), Set(["x", "y", "z", "xy", "xz", "yz", "yx", "zx", "zy"]))
        throw(ArgumentError("The axis $axis is not valid."))
    end
    hascenter = hassymmetry ? true : hascenter
    property = electrons(sim.atoms)
    if length(axis) == 1
        return gaussian_profile_1D(
            sim, property,
            axis=axis, binwidth=binwidth, cutoff=cutoff, resolution=resolution, σ=σ, margin=margin,
            hascenter=hascenter, hassymmetry=hassymmetry
        )
    end
    if length(axis) == 2
        return gaussian_profile_2D(
            sim, property,
            axis=axis, cutoff=cutoff, resolution=resolution, σ=σ, cs=cs,
            hascenter=hascenter, hassymmetry=hassymmetry
        )
    end
end

function gaussmap(coords::MolSimToolkit.FramePositions, Q::Vector{Float64}, σ::Vector{Float64})
    A = Q ./ (sqrt(2π) .* σ).^3
    α = 2σ.^2 .\ -1
    ρ(r) = sum(
        A .* exp.(α .* [ sum((r .- coord).^2) for coord in coords ])
    )
    return ρ
end

function gaussian_profile_1D(
    sim::MolSimToolkit.Simulation, Q::Vector{Float64};
    axis="z",
    cutoff=0.05,
    binwidth=0.25,
    resolution=0.5,
    σ=1.0,
    margin=0.0,
    hascenter=true, hassymmetry=true
)
    dim = get(
        Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3), lowercase(axis), nothing
    )
    isnothing(dim) && throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    Vbin, bin_edges = binning(sim; axis=axis, binwidth=binwidth, margin=margin)
    nbins, nframes = length(bin_edges)-1, length(sim)
    Mbins = zeros(Float64, nbins, nframes)
    densities = deepcopy(Mbins)
    ## gaussian data
    #points = gridpoints(sim, resolution; margin=margin)
    iatoms = PDBTools.index.(sim.atoms)
    σ = σ .* ones(Float64, length(iatoms))
    Q = Q[iatoms]
    for (iframe, frame) in enumerate(sim)
        coords = MolSimToolkit.positions(frame)[iatoms]
        surface = gaussmap(coords, Q, σ)
        points = gridpoints(coords, resolution, margin=margin)
        axis_coords = [ point[dim] for point in points ]
        density = surface.(points)
        # @threads :static for i in eachindex(bins[1:end-1])
        #     idx = findall(
        #         coord -> bins[i] <= coord[dim] < bins[i+1],
        #         points
        #     )
        #     isempty(idx) && continue
        #     @fastmath ρbin = Vbin \ sum(density[idx])
        #     densities[i,iframe] = isnothing(cutoff) || ρbin >= cutoff ? ρbin : 0.0
        # end
        bin_indexes = clamp.(
            searchsortedfirst.(Ref(bin_edges), axis_coords) .- 1, 1, nbins
        )
        hist = zeros(nbins)
        @inbounds for (idx, d) in zip(bin_indexes, density)
            hist[idx] += d
        end
        hist ./= Vbin
        hist = ifelse.(hist .>= cutoff, hist, 0.0)
        densities[:,iframe] = hist
        if !hascenter
            Mbins[:,iframe] = avgbins(bin_edges)
        else
            avg = mean(coord[dim] for coord in coords)
            if hassymmetry
                Mbins[:,iframe] = _fix_binning_symmetry(bin_edges, avg)
                if !issorted(Mbins[:,iframe])
                    isorted = sortperm(Mbins[:,iframe])
                    Mbins[:,iframe], densities[:,iframe] = Mbins[isorted,iframe], densities[isorted,iframe]
                end
            else
                Mbins[:,iframe] = _fix_binning_center(bin_edges, avg)
            end
        end
    end
    return Mbins, densities
end


# function gaussian_profile_1D(
#     sim::MolSimToolkit.Simulation, Q::Vector{Float64};
#     axis="z",
#     cutoff=0.05,
#     binwidth=0.25,
#     resolution=0.5,
#     σ=1.0,
#     margin=0.0,
#     hascenter=true, hassymmetry=true
# )
#     axescode = Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3)
#     if !haskey(axescode, lowercase(axis))
#         throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
#     end
#     dim = axescode[axis]
#     Vbin, bins = binning(sim, axis=axis, binwidth=binwidth, margin=margin)
#     Mbins = zeros(Float64, length(bins)-1, length(sim))
#     densities = deepcopy(Mbins)
#     iatoms = PDBTools.index.(sim.atoms)
#     for (iframe, frame) in enumerate(sim)
#         coords = MolSimToolkit.positions(frame)[iatoms]
#         Q = Q[iatoms]
#         σ = σ .* ones(Float64, length(iatoms))
#         surface = gaussmap(coords, Q, σ)
#         points = gridpoints(coords, resolution, margin=margin)
#         density = surface.(points)
#         @threads :static for i in eachindex(bins[1:end-1])
#             idx = findall(
#                 coord -> bins[i] <= coord[dim] < bins[i+1],
#                 points
#             )
#             isempty(idx) && continue
#             @fastmath ρbin = Vbin \ sum(density[idx])
#             densities[i,iframe] = isnothing(cutoff) || ρbin >= cutoff ? ρbin : 0.0
#         end
#         if !hascenter
#             Mbins[:,iframe] = avgbins(bins)
#         else
#             avg = mean([ coord[dim] for coord in coords ])
#             if hassymmetry
#                 Mbins[:,iframe] = _fix_binning_symmetry(bins, avg)
#                 if !issorted(Mbins[:,iframe])
#                     isorted = sortperm(Mbins[:,iframe])
#                     Mbins[:,iframe], densities[:,iframe] = Mbins[isorted,iframe], densities[isorted,iframe]
#                 end
#             else
#                 Mbins[:,iframe] = _fix_binning_center(bins, avg)
#             end
#         end
#     end
#     return Mbins, densities
# end

################################################################ Fibril width

function mapping_density(density::Array{T, 3}; tol=0.005) where T
    checklist = zeros(Bool, size(density, 1), size(density, 2), size(density, 3))
    for iframe in 1:size(density, 3)
        for i in 1:size(density, 1), j in 1:size(density, 2)
            checklist[i,j,iframe] = density[i,j,iframe] > tol
        end
    end
    return BitArray(checklist)
end

function sortperm_indexes_and_values(vector::AbstractVector{T}) where T
    perm = sortperm(vector)
    return vector[perm], perm
end

function searchsortednearest(sorted_values, value)
    idx = searchsortedfirst(sorted_values, value)
    isnearest = abs(sorted_values[idx] - value) > abs(sorted_values[idx-1] - value)
    return clamp(
        idx - Int(isnearest),
        1,
        length(sorted_values)
    )
end

function radial_scan(xc, yc, θ, sorted_x, sorted_y, idx_x, idx_y, mask, dr)
    t = 0.0
    max_r = 0.0
    while true
        x = xc + t * cos(θ)
        y = yc + t * sin(θ)
        (x < sorted_x[1] || x > sorted_x[end] || 
         y < sorted_y[1] || y > sorted_y[end]) && break
        i = searchsortednearest(sorted_x, x)
        j = searchsortednearest(sorted_y, y)
        mask[idx_x[i], idx_y[j]] || break
        max_r = t
        t += dr
    end
    return max_r
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
        sorted_x, idx_x = sortperm_indexes_and_values(x)
        sorted_y, idx_y = sortperm_indexes_and_values(y)
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