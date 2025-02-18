const cellulose_atoms_vdw_radii = Dict{String, Float64}(
    ## carbons
    "C1" => 2.000, "C2" => 2.000, "C3" => 2.000, "C4" => 2.000, "C5" => 2.000, "C6" => 2.010,
    ## oxygens
    "O1" => 1.765, "O2" => 1.765, "O3" => 1.765, "O4" => 1.765, "O5" => 1.650, "O6" => 1.765,
    ## hydrogens
    "H1" => 1.340, "H2" => 1.340, "H3" => 1.340, "H4" => 1.340, "H5" => 1.340, "H61" => 1.340,
    "H62" => 1.340, "HO1" => 1.800, "HO2" => 1.800, "HO3" => 1.800, "HO4" => 1.800, "HO6" => 1.800
)   ## they are the CHARMM36 force field standard radii! Please check top_all36_carb.rtf file

"""

    surface(rawfile::String; surface_method="distance", selection="all", parameters="-res 0.6 -cutoff 4.0", vmd="vmd")

Generate a surface from a VMD input file using the `volmap` plugin. The surface can be generated using different methods, such as `distance`, `density`, among others.
The default method is `distance`. The selection can be defined by the user, and the default is `all`.

The parameters can be set by the user, and the default is `-res 0.6 -cutoff 4.0`. The VMD executable can be provided as an argument whether is not on default PATH.

## Arguments

- `rawfile::String`: The name of the input file without extension. The script will look for `rawfile.pdb` and `rawfile.psf`.

- `surface_method::String="distance"`: The method to be used to generate the surface using the volmap tools. The default is `distance`.
    `density`
        creates a map of the weighted atomic density at each gridpoint. This is done by replacing each atom in the selection with a normalized gaussian distribution
        of width (standard deviation) equal to its atomic radius. The gaussian distribution for each atom is then weighted using an optional -weight argument, and
        defaults to a weight of one (i.e, the number density). The various gaussians are then additively distributed on a grid.
    `interp`
        creates a map with the atomic weights interpolated onto a grid. For each atom, its weight is distributed to the 8 nearest voxels via a trilinear interpolation.
        The optional -weight argument defaults to a weight of one.
    `distance`
        creates a map for which each gridpoint contains the distance between that point and the edge of the nearest atom. In other words, each gridpoint specifies the
        maximum radius of a sphere cnetered at that point which does not intersect with the spheres of any other atoms. All atoms are treated as spheres using the
        atoms' VMD radii.
    `coulomb`, `coulombmsm`
        Creates a map of the electrostatic field of the atom selection, made by computing the non-bonded Coulomb potential from each atom in the selection (in units
        of **k_B.T/e**). The coulomb map generation is optimized to take advantage of multi-core CPUs and programmable GPUs if they are available.
    `ils`
        Creates a free energy map of the distribution of a weakly-interacting monoatomic or diatomic gas ligand throughout the system using the Implicit Ligand
        Sampling (ILS) technique. See additional information about ILS below.
    `mask`
        Creates a map which is set to 0 or 1 depending on whether they are within a specified -cutoff distance argument of any atoms in the selection. The mask map is
        typically used in combination with other maps in order to hide/mask data that is far from a region of interest.
    `occupancy`
        Each grid point is set to either 0 or 1, depending on whether it contains onbe or more atoms or not. When averaged over many frames, this will provide the
        fractional occupancy of that grid point. By default, atoms are treated as spheres using the atomic radii and a gridpoint is considered to be "occupied" if
        it lies inside that sphere. Use the -points argument to treat atoms as points (a grid point is "occupied" if its grid cube contains an atom's center). 

- `selection::String="all"`: The selection to be used to generate the surface. The default is `all`.

- `parameters::String="-res 0.6 -cutoff 4.0"`: The parameters to be used in the `volmap` plugin. The default is `-res 0.6 -cutoff 4.0`. To understand the parameters,
   see the VMD documentation https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node157.html.

- `vmd::String="vmd"`: The VMD executable. The default is `vmd`.

### Examples

```julia-repl

julia > surface("new_crystal")

```

"""
function cellulose_surface(
    pdbname::String;
    psfname=nothing,
    vdw_radii=nothing,
    surface_method="distance", selection="all", parameters="-res 0.6 -cutoff 4.0",
    vmd="vmd", vmdDebug=false
)
    pdbname = abspath(pdbname)
    psfname = isnothing(psfname) ? replace(pdbname, ".pdb" => ".psf") : psfname
    if !isfile(pdbname) || !isfile(psfname)
        error("File not found: Be sure that exist $psfname and $pdbname files.")
    end
    if !isnothing(vdw_radii)
        if (typeof(vdw_radii) != Vector{Float64}) || (length(vdw_radii) != length(cellulose_atoms_vdw_radii))
            throw(ArgumentError("The vdw_radii should be a `Vector{Float64}` and it must have 12 elements."))
        end
    end
    tcl = tempname() * ".tcl"
    dxname, stlname = replace(pdbname, ".pdb" => ".dx"), replace(pdbname, ".pdb" => ".stl")         ## the surface files
    open(tcl, "w") do file
        println(file, """
        mol new     "$psfname"
        mol addfile "$pdbname"
        """)
        for (i, key) in enumerate(keys(cellulose_atoms_vdw_radii))
            r = isnothing(vdw_radii) ? cellulose_atoms_vdw_radii[key] : vdw_radii[i]
            println(file, "[atomselect top \"name $key\"] set radius $r")
        end
        println(file, """
        volmap $surface_method [atomselect top "$selection"] $parameters -o $dxname
        
        mol delete top
        
        mol new "$dxname"
        
        axes location off
        mol modstyle 0 top isosurface 0.5 0 0 0 1 1
        set unitary_space {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}
        
        molinfo top set center_matrix [list \$unitary_space]
        molinfo top set rotate_matrix [list \$unitary_space]
        molinfo top set scale_matrix  [list \$unitary_space]
        
        render STL $stlname true
        
        mol delete all

        exit
        """)
    end
    vmdoutput = execVMD(tcl, vmd=vmd)
    if vmdDebug
        return vmdoutput
    else
        return dxname, stlname
    end
end

function filterSTL(stlname::String)
    STL = Vector{SVector{3, Float64}}()
    open(stlname) do file
        for line in eachline(file)
            if occursin("vertex", line)
                vertex = string.(split(line))[2:end]
                push!(STL, SVector{3, Float64}(parse.(Float64, vertex)))
            end
        end
    end
    return unique(STL)
end

function chain_centers(pdbname::String; step::Float64=1.0, cutoff::Float64=0.5, isreference::Bool=true)
    center = Vector{SVector{3, Float64}}()
    atoms = PDBTools.readPDB(pdbname)
    monomers = PDBTools.resnum.(atoms)
    p_first = mean(
            PDBTools.coor(PDBTools.select(atoms, atom -> atom.resnum == minimum(monomers)))
        )
    p_last = mean(
            PDBTools.coor(PDBTools.select(atoms, atom -> atom.resnum == maximum(monomers)))
        )
    displace = p_last - p_first
    norm_displace = norm(displace)
    for i in 0:step:norm_displace
        icenter = p_first .+ i * (displace/norm_displace)
        if isreference
            push!(center, icenter)
        else
            iatoms = PDBTools.select(atoms, atom -> norm(PDBTools.coor(atom) - icenter) <= cutoff)
            if isempty(iatoms)
                push!(center, icenter)
            else
                push!(center, mean(PDBTools.coor(iatoms)))
            end
        end
    end
    return center
end

function fibrilwidth(
    pdbname::String;
    selection="not water",
    step=1.0, tol=0.5, isreference=false,
    psfname=nothing, vdw_radii=nothing,
    surface_parameters="-res 0.6 -cutoff 4.0", vmd="vmd"
)
    stl = cellulose_surface(
        pdbname, selection=selection, psfname=psfname, vdw_radii=vdw_radii, parameters=surface_parameters, vmd=vmd
    )[2]
    fibril_centers = chain_centers(
        pdbname, step=step, cutoff=0.5*step, isreference=isreference
    )
    if 2*tol > step
        throw(ArgumentError("The tolerance should be less than half of the step, at least, equal to it."))
    end
    xyz = filterSTL(stl)
    return scan_diameter(xyz, fibril_centers, step=step, tol=tol)
end

function scan_diameter(xyz::Vector{SVector{3, Float64}}, centers::Vector{SVector{3, Float64}}; step=1.0, tol=0.05)
    d, Φ, Θ = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
    for center in centers
        isurface = findall(k -> norm(k - center[3]) <= tol, [ ijk[3] for ijk in xyz ])
        surface = [ surf - center for surf in xyz[isurface] ]
        r = norm.(surface)
        x, y, z = [ s[1] for s in surface ], [ s[2] for s in surface ], [ s[3] for s in surface ]
        norm_xy = norm.([ s[1:2] for s in surface ])
        ϕ = mod.(
            sign.(y) .* acosd.(x ./ norm_xy),
            360.0
        )
        θ = mod.(
            acosd.(z ./ r),
            360.0
        )
        r, azimuth, inclination = getradii(r, ϕ, θ, step=step, tol=tol)
        push!(d, r); push!(Φ, azimuth); push!(Θ, inclination)
    end
    #return d, Φ, Θ  ## I will remove the inlcination later, cause I only need the azimuthal angle to calculate the diameter
    return radii2diameter(d, Φ)
end

function getradii(
    radii::Vector{Float64},
    azimuthal::Vector{Float64},
    polar::Vector{Float64};
    step=1.0, tol=0.05
)
    ρ, ϕ, θ = Float64[], Float64[], Float64[]
    for angle in 0.0:step:360.0
        a1, a2 = mod(angle - .5*step, 360.0), mod(angle + .5*step, 360.0)
        idxa = a1 < a2 ? findall(i -> (i >= a1) && (i <= a2), azimuthal) : findall(i -> (i >= a1) || (i <= a2), azimuthal)
        r, inclination = Vector{Float64}(undef, length(idxa)), Vector{Float64}(undef, length(idxa))
        token = 1
        for idx in idxa
            p1, p2 = mod(polar[idx] - tol, 360.0), mod(polar[idx] + tol, 360.0)
            idxp = p1 < p2 ? findall(i -> (i >= p1) && (i <= p2), polar) : findall(i -> (i >= p1) || (i <= p2), polar)
            if isempty(idxp)
                continue
            end
            imax = argmax(radii[idxp])
            r[token], inclination[token] = radii[imax], polar[imax]
            token += 1
        end
        if isempty(r)
            continue
        else
            push!(ρ, median(r)); push!(ϕ, angle); push!(θ, median(inclination))
        end
        println("")
    end
    return ρ, ϕ, θ
end

function radii2diameter(radii::Vector{Vector{Float64}}, angles::Vector{Vector{Float64}})
    diameters = Float64[]
    for istep in eachindex(radii, angles)
        r_slices, ϕ_slices = radii[istep], angles[istep]
        scanning = Dict{Tuple{Float64, Float64}, Float64}()
        for (ri, ϕi) in zip(r_slices, ϕ_slices)
            ϕi_norm = mod(ϕi, 360.0)    ## to avoid double counting the same angle (i.e. 0 and 360)
            ϕf = mod(ϕi_norm + 180.0, 360.0)
            if haskey(scanning, (ϕi_norm, ϕf)) || haskey(scanning, (ϕf, ϕi_norm))
                continue
            end
            idx = findfirst(ϕ -> ϕ == ϕf, ϕ_slices)
            if !isnothing(idx)
                scanning[(ϕi_norm, ϕf)] = ri + r_slices[idx]
            else
                idx1, idx2 = if ϕf ≈ 0.0 || ϕf ≈ 360.0
                        findfirst(ϕ -> ϕ > 0.0, ϕ_slices), findlast(ϕ -> ϕ < 360.0, ϕ_slices)
                    else
                        findfirst(ϕ -> ϕ > ϕf, ϕ_slices), findlast(ϕ -> ϕ < ϕf, ϕ_slices)
                end
                if !isnothing(idx1) && !isnothing(idx2)
                    ϕ1, ϕ2 = ϕ_slices[idx1], ϕ_slices[idx2]
                    r1, r2 = r_slices[idx1], r_slices[idx2]
                    scanning[(ϕi_norm, ϕf)] = ri + r1 + (r2 - r1) * (ϕf - ϕ1) / (ϕ2 - ϕ1)
                end
            end    
        end
        if !isempty(scanning)
            append!(diameters, values(scanning))
        end
    end    
    return diameters
end

##################################################################################################################################################
## pdb = "/home/user/Documents/FibrilWidth/setup/system.pdb"
## dcd = "/home/user/Documents/FibrilWidth/output/prod_NpT.dcd"

### Varrer apenas as cargas Σqi
### Com a gausiana seria r_i - r_bin
function edp(
    pdbname::String, trjname::String;
    psfname=nothing, profile_selection="not water", align_selection="not water",
    axis="z", cutoff=4.0, resolution=0.6, q=nothing, first=1, last=nothing, step=1
)
    if !in(lowercase(z), Set(["x", "y", "z", "xy", "xz", "yz", "yx", "zx", "zy"]))
        throw(ArgumentError("""
        The axis should be associated with the labels of cartesian axes or plane.
        Please insert it such as: `x`, `y`, `z`, `xy`, `xz`, or `yz`.
        """))
    end
    new_trjname = align_frames(psfname, trjname, selection=align_selection)
    atoms = PDBTools.select(PDBTools.readPDB(pdbname), profile_selection)
    sim = MolSimToolkit.Simulation(
        atoms, new_trjname,
        first=first, last=last, step=step
    )
    psfname = isnothing(psfname) ? replace(pdbname, ".pdb" => ".psf") : psfname
    q = isnothing(q) ? chargesPSF(psfname)[PDBTools.index.(atoms)] : q
    ndims, bins = length(axis), binning(sim, axis=axis, cutoff=cutoff)
    if ndims == 1
        return edp_1d(
            sim, bins, q,
            axis=axis, cutoff=cutoff
        )
    else
        throw(ArgumentError("The axis should be a string with one or two characters."))
    end
end

function edp_1d(
    sim::MolSimToolkit.Simulation,
    bins::Vector{Float64}, q::Vector{Float64};
    axis="z", cutoff=4.0, resolution=0.6, isdimensionless=true,
    hascenter=true, hassymmetry=true
)
    densities = Vector{Vector{Float64}}(undef, length(sim))
    byframes = Vector{Vector{Float64}}(undef, length(sim))
    axesdecode = Dict{String, Vector{Int64}}("x" => [1], "y" => [2], "z" => [3])
    dim = if haskey(axesdecode, lowercase(axis))
            axescode[axis]
        else
            throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    hascenter = hassymmetry ? true : hascenter
    iatoms = hascenter ? PDBTools.index.(sim.atoms) : nothing
    V, bins = binning(sim, axis=axis, cutoff=cutoff)
    ρtemp = zeros(Float64, length(bins)-1)
    for (iframe, frame) in enumerate(sim)
        coords = MolSimToolkit.positions(frame)
        ibin = 1
        while ibin < length(bins)
            idx = findall(coord -> bins[ibin] <= coord[dim] <= bins[ibin+1], coords)
            ρtemp[ibin] = isdimensionless ? ρ(q[idx], V)
            ibin += 1
        end
        densities[iframe] = deepcopy(ρtemp)
        if !hascenter
            byframes[iframe] = avgbins(bins)
        else
            avg = mean(coords[iatoms])
            if hassymmetry
                byframes[iframe] = _fix_binning_symmetry(bins, avg)
                if !issorted(byframes[iframe])
                    isorted = sortperm(byframes[iframe])
                    byframes[iframe], densities[iframe] = byframes[iframe][isorted], densities[iframe][isorted]
                end
            else
                byframes[iframe] = _fix_binning_center(bins, avg)
            end
        end
    end
    return byframes, densities
end

# for (iframe, frame) in enumerate(sim)
#     coords = MolSimToolkit.positions(frame)
#     xyz = [ xyz[dim] for xyz in coords ]
#     densities[iframe] = ρbins(bins, xyz, q, N=Vbin)
#     if isnothing(idx)
#         byframes[iframe] = avgbins(bins)
#     else
#         avg = mean(xyz[idx])
#         if hassymmetry
#             byframes[iframe] = _fix_binning_symmetry(bins, avg)
#             if !issorted(byframes[iframe])
#                 isorted = sortperm(byframes[iframe])
#                 byframes[iframe], densities[iframe] = byframes[iframe][isorted], densities[iframe][isorted]
#             end
#         else
#             byframes[iframe] = _fix_binning_center(bins, avg)
#         end
#     end
# end

function ρbins(bins::Vector{Float64}, coords::Vector{Float64}, q::Vector{Float64}; N=nothing, σ=nothing)
    N = isnothing(N) ? length(coords) : N
    σ = isnothing(σ) ? ones(length(coords)) : σ
    density = zeros(Float64, length(bins))
    for (i, bin) in enumerate(bins)
        if lastindex(bins) == i
            break
        end
        δ = δbin.(coords, bin, bins[i+1])
        electron = eprofile(coords, q, σ)
        density[i] = N \ sum(δ .* electron)
    end
    return density
end

function ρ(q::Vector{Float64}, N::Float64)
    return N \ sum(q)
end

function ρ(coords, q, N; σ=nothing)
    σ = isnothing(σ) ? ones(length(coords)) : σ
    A = inv.(
        sqrt.(2π) .* σ
    ).^3
    return N \ sum(δbin(coords, bin, bins[i+1]) .* eprofile(coords, q, σ) for (i, bin) in enumerate(bins))
end

function σ_edp(m::)
end

function eprofile(coords::Vector{Float64}, q::Vector{Float64}, σ::Vector{Float64}; isdimensionless=true)
    if isdimensionless
        return q
    end
    Ng = inv.(
        sqrt.(2π) .* σ
    ) ## gaussian normalization
    score = (coords .- mean(coords)) ./ σ
    return q .* Ng .* exp.(-0.5 .* score.^2)
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

# function binning_reference(
#     sim::MolSimToolkit.Simulation;
#     selection="not water", axis="z", Δ=1.0,
#     hascenter=true, hassymmetry=true
# )
#     atoms = PDBTools.select(
#         MolSimToolkit.atoms(sim), selection
#     )
#     idx, axis = PDBTools.index.(atoms), split(lowercase(axis), "")
#     reference = length(axis) == 1 ? Vector{SVector(1, Float64)}(undef, length(sim)) : Vector{SVector(2, Float64)}(undef, length(sim))
#     for (iframe, frame) in enumerate(sim)
#         coords = MolSimToolkit.positions(frame)[idx]
#         if length(axis) == 1
#             dict = Dict{String, Vector{Int64}}(
#                 "x" => [1], "y" => [2], "z" => [3]
#             )
#             dim = dict[axis[1]]
#         elseif length(axis) == 2
#             dict = Dict{Tuple{String, String}, Vector{Int64}}([
#                 [("x", "y"), ("y", "x")] .=> [1, 2]; [("x", "z"), ("z", "x")] .=> [1, 3]; [("y", "z"), ("z", "y")] .=> [2, 3]
#             ])
#             dim = dict[(axis[1], axis[2])]
#         end
#         reference[iframe] = mean(coords[dim])
#     end
#     return reference
# end

function chargesPSF(psfname::String)
    q = Vector{Float64}()
    natoms = 0.0
    open(psfname) do file
        for line in eachline(file)
            if occursin("!NATOM", line)
                natoms = parse(Float64, string(split(line)[begin]))
                continue
            end
            if iszero(natoms)
                continue
            else
                qi = string(split(line)[7])
                push!(q, parse(Float64, qi))
                natoms -= 1
            end
        end
    end
    return q
end

function align_frames(
    psfname::String,
    trjname::String; new_trajectory=nothing,
    selection="not water", reference=0,
    vmd="vmd", DebugVMD=false
)    
    new_trajectory = isnothing(new_trajectory) ? tempname() * ".dcd" : new_trajectory
    tcl = tempname() * ".tcl"
    Base.open(tcl, "w") do file
        println(file, """
        package require pbctools

        mol new     $psfname
        mol addfile $trjname waitfor all
        
        animate goto 0
        pbc wrap -centersel \"$selection\" -center com -compound residue -all        
        set ref_frame $reference
        set ref_selection [atomselect top \"$selection\" frame \$ref_frame]

        set num_frames [molinfo top get numframes]

        for {set i 0} {\$i < \$num_frames} {incr i} {
            set cur_selection [atomselect top \"$selection\" frame \$i]
            set mat_trans [measure fit \$cur_selection \$ref_selection]
            set sel_all [atomselect top "all" frame \$i]
            \$sel_all move \$mat_trans
        }

        animate write dcd $new_trajectory
        """)
    end
    vmdoutput = Base.read(`$vmd -dispdev text -e $tcl`, String)
    return DebugVMD ? vmdoutput : new_trajectory
end