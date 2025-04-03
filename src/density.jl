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

function charmm_partialcharges(atom::AbstractString, residue::AbstractString)
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

function interpol(x::AbstractVector{T}, y::AbstractVector{T}, z::Matrix{T}; bins=200) where T
    itp = Interpolations.linear_interpolation((x, y), z, extrapolation_bc=NaN)
    xnew = range(minimum(x), maximum(x), length=bins)
    ynew = range(minimum(y), maximum(y), length=bins)
    znew = [ itp(i,j) for i in xnew, j in ynew ]
    return xnew, ynew, znew
end

function notaxis(axis::String)
    notaxis = filter(
        i -> !in(i, split(lowercase(axis), "")),
        ["x", "y", "z"]
    )
    return join(notaxis)
end

# function axis2dims(axis::String)
#     axis = split(lowercase(axis), "")
#     if length(axis) == 1
#         return Dict{String, Int}(
#             "x" => 1,
#             "y" => 2,
#             "z" => 3
#         )[axis[1]]
#     elseif length(axis) == 2
#         return Dict{String, Function}([
#             ["xy", "yx"] .=> () -> [1, 2];
#             ["xz", "zx"] .=> () -> [1, 3];
#             ["yz", "zy"] .=> () -> [2, 3]
#         ])[join(axis)]()
#     end
# end

@inline function avgbins(bin_edges::Vector{Float64})
    return 0.5 .* (bin_edges[1:end-1] + bin_edges[2:end])
end

@inline function avgbins(bin_edges::Vector{Float64}, average::Float64; hascenter::Bool=true, hassymmetry::Bool=true)
    if !hascenter
        return avgbins(bin_edges)
    end
    return hassymmetry ? _fix_binning_symmetry(bin_edges, average) : _fix_binning_center(bin_edges, average)
end

@inline function _fix_binning_center(bin_edges::Vector{Float64}, center::Float64)
    return avgbins(bin_edges .- center)
end

@inline function _fix_binning_symmetry(bin_edges::Vector{Float64}, center::Float64)
    lower, upper = extrema(bin_edges)
    ΔL = upper - lower
    avg = _fix_binning_center(bin_edges, center)
    return mod.(avg .+ 0.5*ΔL, ΔL) .- 0.5*ΔL
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

@inline function extremacoordinates(simulation::MolSimToolkit.Simulation, indices::AbstractVector)
    lower, upper = fill(Inf, 3), fill(-Inf, 3)
    @inbounds for frame in simulation
        coords = MolSimToolkit.positions(frame)[indices]
        mincoord, maxcoord = extremacoordinates(coords)
        lower .= min.(lower, mincoord)
        upper .= max.(upper, maxcoord)
    end
    return lower, upper
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
    simulation::MolSimToolkit.Simulation,
    resolution::Float64;
    margin::Float64=0.0
)
    lower, upper = extremacoordinates(simulation)
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

function gridpoints(
    simulation::MolSimToolkit.Simulation,
    indices::AbstractVector,
    resolution::Float64;
    margin::Float64=0.0
)
    lower, upper = extremacoordinates(simulation, indices)
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
            axis=axis, binwidth=binwidth, cutoff=cutoff, resolution=resolution, σ=σ, cswidth=cs, margin=margin,
            hascenter=hascenter, hassymmetry=hassymmetry     
        )
    end
end

# function gaussmap(coords::MolSimToolkit.FramePositions, Q::Vector{Float64}, σ::Vector{Float64})
#     A = Q ./ (sqrt(2π) .* σ).^3
#     α = 2σ.^2 .\ -1
#     ρ(r) = sum(
#         A .* exp.(α .* [ sum((r .- coord).^2) for coord in coords ])
#     )
#     return ρ
# end

function ρ(coords::MolSimToolkit.FramePositions, A::Vector{Float64}, β::Vector{Float64})
    density(r) = sum(
        A .* exp.(β .* [ sum((r .- coord).^2) for coord in coords ])
    )
    return density
end

function gaussian_density_function!(
    density::Vector{T},
    points::AbstractArray{<:SVector{3, T}},
    coords::MolSimToolkit.FramePositions,
    A::Vector{T},
    β::Vector{T}
) where T
    npoints = length(points)
    natoms = length(coords)
    Threads.@threads for i in 1:npoints
        r = points[i]
        sum_d = 0.0
        @inbounds @fastmath for j in 1:natoms
            Δx = r[1] - coords[j][1]
            Δy = r[2] - coords[j][2]
            Δz = r[3] - coords[j][3]
            sum_d += A[j] * exp(β[j] * (Δx^2 + Δy^2 + Δz^2))
        end
        density[i] = sum_d
    end
end

function gaussian_density_function(
    gridpoints::AbstractArray{<:SVector{3, T}},
    coords::MolSimToolkit.FramePositions,
    A::Vector{T},
    β::Vector{T}
) where T
    density = zeros(T, length(gridpoints))
    @inbounds @fastmath Threads.@threads for i in eachindex(gridpoints)
        r = gridpoints[i]
        sum_val = 0.0
        for j in eachindex(coords)
            dx = r[1] - coords[j][1]
            dy = r[2] - coords[j][2]
            dz = r[3] - coords[j][3]
            sum_val += A[j] * exp(β[j] * (dx^2 + dy^2 + dz^2))
        end
        density[i] = sum_val
    end
    return density
end

function searching_slabs(coords::MolSimToolkit.FramePositions, dim::Int, bin_edges::AbstractVector)
    indices = searchsortedfirst.(Ref(bin_edges), coord[dim] for coord in coords) .- 1
    return clamp.(
        indices, 1, length(bin_edges)-1
    )
end

function get_slab_indices(coords, z_dim, edges)
    z_coords = [c[z_dim] for c in coords]
    return clamp.(searchsortedfirst.(Ref(edges), z_coords) .- 1, 1, length(edges)-1)
end

function gaussian_profile_1D(
    simulation::MolSimToolkit.Simulation, Q::Vector{Float64};
    axis="z",
    cutoff=0.05,
    binwidth=0.25,
    resolution=0.5,
    σ=1.0,
    margin=0.0,
    hascenter=true, hassymmetry=true
)
    dim = get(Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3), lowercase(axis), nothing)
    isnothing(dim) && throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    Vbin, bin_edges = binning(
        simulation;
        axis=axis, binwidth=binwidth, margin=margin
    )
    nbins, nframes = length(bin_edges)-1, length(simulation)
    bincenters = zeros(Float64, nframes)
    densities = zeros(Float64, nbins, nframes)
    iatoms = PDBTools.index.(simulation.atoms)
    ## gaussian data
    A = Q[iatoms] ./ (sqrt(2π) .* σ).^3
    β = inv.(-2σ.^2) .* ones(Float64, length(iatoms))
    points = gridpoints(simulation, iatoms, resolution; margin=margin)
    dimcoords = [ point[dim] for point in points ]
    density = Vector{Float64}(undef, length(points))
    for (iframe, frame) in enumerate(simulation)
        coords = MolSimToolkit.positions(frame)[iatoms]
        bincenters[iframe] = mean(coords[dim] for coord in coords)[1]
        gaussian_density_function!(density, points, coords, A, β)
        bin_indices = clamp.(
            searchsortedfirst.(Ref(bin_edges), dimcoords) .- 1, 1, nbins
        )
        hist = zeros(nbins)
        @fastmath @inbounds for (idx, d) in zip(bin_indices, density)
            hist[idx] += d
        end
        densities[:,iframe] = ifelse.(hist .>= cutoff, hist / Vbin, 0.0)
    end
    Mbins = Matrix{Float64}(undef, nbins, nframes)
    @threads for iframe in 1:nframes
        bins = avgbins(
            bin_edges,
            bincenters[iframe];
            hascenter=hascenter,
            hassymmetry=hassymmetry
        )
        Mbins[:, iframe] = sort(bins)
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
#     dim = get(
#         Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3), lowercase(axis), nothing
#     )
#     isnothing(dim) && throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
#     Vbin, bin_edges = binning(sim; axis=axis, binwidth=binwidth, margin=margin)
#     nbins, nframes = length(bin_edges)-1, length(sim)
#     Mbins = zeros(Float64, nbins, nframes)
#     densities = deepcopy(Mbins)
#     iatoms = PDBTools.index.(sim.atoms)
#     ## gaussian data
#     A = Q[iatoms] ./ (sqrt(2π) .* σ).^3
#     β = inv.(-2σ.^2) .* ones(Float64, length(iatoms))
#     points = gridpoints(sim, iatoms, resolution; margin=margin)
#     axis_coords = [ point[dim] for point in points ]
#     density = Vector{Float64}(undef, length(points))
#     for (iframe, frame) in enumerate(sim)
#         coords = MolSimToolkit.positions(frame)[iatoms]
#         gaussian_density_function!(density, points, coords, A, β)
#         bin_indices = clamp.(
#             searchsortedfirst.(Ref(bin_edges), axis_coords) .- 1, 1, nbins
#         )
#         hist = zeros(nbins)
#         @fastmath @inbounds for (idx, d) in zip(bin_indices, density)
#             hist[idx] += d
#         end
#         densities[:,iframe] = ifelse.(hist .>= cutoff, hist / Vbin, 0.0)
#         Mbins[:, iframe] = avgbins(
#             bin_edges,
#             mean(coord[dim] for coord in coords);
#             hascenter=hascenter,
#             hassymmetry=hassymmetry
#         )
#     end
#     return Mbins, densities
# end


function gaussian_profile_2D(
    simulation::MolSimToolkit.Simulation, Q::Vector{Float64};
    axis="xy",
    cutoff=1.0,
    cswidth=4.0,
    binwidth=0.5,
    resolution=0.5,
    σ=1.0,
    margin=0.0,
    hascenter=true, hassymmetry=true
)
    dims = get(
        Dict{String, Vector{Int}}(
            "xy" => [1, 2], "yx" => [1, 2], "xz" => [1, 3],
            "zx" => [1, 3], "yz" => [2, 3], "zy" => [2, 3]),
        lowercase("xy"), nothing
    )
    isnothing(dims) && throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    dim = filter(t -> !in(t, dims), 1:3)[1]
    iatoms = PDBTools.index.(simulation.atoms)
    ## binning and pre-allocation
    Vbin, bin_edges_1, bin_edges_2 = binning(simulation, axis=axis, binwidth=binwidth, margin=margin)
    cs_bin_edges = binning(simulation, axis=notaxis(axis), binwidth=cswidth)[2]
    nbins1, nbins2, nframes = length(bin_edges_1)-1, length(bin_edges_2)-1, length(simulation)
    nslabs = length(cs_bin_edges)-1
    densities = [ zeros(Float64, nbins1, nbins2, nslabs) for _ in 1:nframes ]
    bincenters1 = zeros(Float64, nframes)
    bincenters2 = zeros(Float64, nframes)
    ## gaussian data
    A = Q[iatoms] ./ (sqrt(2π) .* σ).^3
    β = inv.(-2σ.^2) .* ones(Float64, length(iatoms))
    ## gridpoints to calculate density
    points = gridpoints(simulation, resolution; margin=margin)
    dimcoords1, dimcoords2 = [ point[dims[1]] for point in points ], [ point[dims[2]] for point in points ]
    for (iframe, frame) in enumerate(simulation)
        coords = MolSimToolkit.positions(frame)[iatoms]
        slabindices =get_slab_indices(coords, dim, cs_bin_edges)
        bincenters1 = mean(coord[dims[1]] for coord in coords)
        bincenters2 = mean(coord[dims[2]] for coord in coords)
        for islab in 1:nslabs
            mask = slabindices .== islab
            sum(mask) == 0 && continue
            density = zeros(Float64, length(points))
            gaussian_density_function!(density, points, coords[mask], A[mask], β[mask])
            ibin = clamp.(
                searchsortedfirst.(Ref(bin_edges_1), dimcoords1) .- 1, 1, nbins1
            )
            jbin = clamp.(
                searchsortedfirst.(Ref(bin_edges_2), dimcoords2) .- 1, 1, nbins2
            )
            @inbounds @fastmath for p in 1:length(points)
                i, j = ibin[p], jbin[p]
                densities[iframe][i, j, islab] += density[p]
            end
            densities[iframe][:, :, islab] ./= Vbin
            densities[:,iframe] = ifelse.(hist .>= cutoff, hist / Vbin, 0.0)
        end
    end
    bins1 = [ zeros(Float64, nbins1, nslabs) for _ in 1:nframes ]
    bins2 = [ zeros(Float64, nbins2, nslabs) for _ in 1:nframes ]
    @threads for iframe in 1:nframes
        tmp_bins1 = avgbins(bin_edges_1, bincenters1[iframe]; hascenter=hascenter, hassymmetry=hassymmetry)
        tmp_bins2 = avgbins(bin_edges_2, bincenters2[iframe]; hascenter=hascenter, hassymmetry=hassymmetry)
        tmp_bins1, tmp_bins2 = sort(tmp_bins1), sort(tmp_bins2)
        for islab in 1:nslabs
            bins1[iframe][:, islab] = sort(tmp_bins1)
            bins2[iframe][:, islab] = sort(tmp_bins2)
        end
    end
    return bins1, bins2, densities
end

# function gaussian_profile_2D(
#     sim::MolSimToolkit.Simulation, Q::Vector{Float64};
#     axis="xz", cutoff=1.0, cs=4.0, resolution=0.5, σ=1.0, hascenter=true, hassymmetry=true
# )
#     BIN1, BIN2, DENS = [], [], []
#     axescode = Dict{String, Function}([
#         ["xy", "yx"] .=> () -> [1,2]; ["xz", "zx"] .=> () -> [1,3]; ["yz", "zy"] .=> () -> [2,3]
#     ])
#     if !haskey(axescode, lowercase(axis))
#         throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
#     end
#     dims = axis2dims(axis)
#     dim = axis2dims(notaxis(axis))
#     V, bin1, bin2 = binning(sim, axis=axis, cutoff=cutoff)
#     cs = binning(sim, axis=notaxis(axis), cutoff=cs)[2]
#     bins1 = Matrix{Float64}(undef, length(bin1)-1, length(sim))
#     bins2 = Matrix{Float64}(undef, length(bin2)-1, length(sim))
#     densities = zeros(Float64, length(bin1)-1, length(bin2)-1, length(sim))
#     iatoms = PDBTools.index.(sim.atoms)
#     for (iframe, frame) in enumerate(sim)
#         xyz = MolSimToolkit.positions(frame)[iatoms]
#         σ = σ .* ones(Float64, length(iatoms))
#         A = Q[iatoms] ./ (sqrt(2π) .* σ).^3
#         α = inv.(-2 .* σ.^ 2)
#         ρ(r) = sum(
#             A .* exp.(α .* [ sum((r .- ijk).^2) for ijk in xyz ])
#         )
#         for slice in 1:length(cs)-1
#             coords = filter(coord -> cs[slice] <= coord[dim] <= cs[slice+1], xyz)
#             if isempty(coords)
#                 continue
#             end
#             points = density_grid(coords, resolution)
#             density = ρ.(points)
#             ibin = [ searchsortedfirst(bin1, point[dims[1]]) for point in points ]
#             jbin = [ searchsortedfirst(bin2, point[dims[2]]) for point in points ]
#             @inbounds for p in 1:length(points)
#                 i, j = ibin[p], jbin[p]
#                 if 1 <= i <= size(densities, 1) && 1 <= j <= size(densities, 2)
#                     densities[i, j, iframe] += density[p]
#                 end
#             end
#             densities[:,:,iframe] ./= V
#             if !hascenter
#                 bins1[:,iframe], bins2[:,iframe] = avgbins(bin1), avgbins(bin2)
#             else
#                 avg = mean([ coord[dims] for coord in coords ])
#                 if hassymmetry
#                     bins1[:,iframe] = _fix_binning_symmetry(bin1, avg[1])
#                     bins2[:,iframe] = _fix_binning_symmetry(bin2, avg[2])
#                     if !issorted(bins1[:,iframe]) || !issorted(bins2[:,iframe])
#                         isorted, jsorted = sortperm(bins1[:,iframe]), sortperm(bins2[:,iframe])
#                         bins1[:,iframe], bins2[:,iframe] = bins1[isorted,iframe], bins2[jsorted,iframe]
#                         densities[:,:,iframe] = densities[isorted,jsorted,iframe]
#                     end
#                 else
#                     bins1[:,iframe] = _fix_binning_center(bin1, avg[1])
#                     bins2[:,iframe] = _fix_binning_center(bin2, avg[2])
#                 end
#             end
#             push!(BIN1, deepcopy(bins1))
#             push!(BIN2, deepcopy(bins2))
#             push!(DENS, deepcopy(densities))
#             densities .= 0.0
#         end
#     end
#     return BIN1, BIN2, DENS
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