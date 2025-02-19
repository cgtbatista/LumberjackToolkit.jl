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

### align_frames(psf, dcd, selection="not water")
### edp()
function plot_average_edp(nbins::Vector{Vector{Float64}}, densities::Vector{Vector{Float64}})
    avg = deepcopy(nbins[1])
    for i in 2:length(nbins)
        avg .+= nbins[i]
    end
    return avg ./ length(nbins)
end
#at -> PDBTools.element(at) != "H"
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
    e = isnothing(e) ? chargesPSF(psfname) : e
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
    else
        throw(ArgumentError("The axis should be a string with one or two characters."))
    end
end

function edp_1d_dimensionless(
    sim::MolSimToolkit.Simulation, e::Vector{Float64};
    axis="z", cutoff=1.0, hascenter=true, hassymmetry=true
)
    axescode = Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3)
    if haskey(axescode, lowercase(axis))
        dim = axescode[axis]
        V, bins = binning(sim, axis=axis, cutoff=cutoff)
        density = Vector{Float64}(undef, length(bins)-1)
    else
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    iatoms = PDBTools.index.(sim.atoms)    
    electrons = e[iatoms]
    densities = Vector{Vector{Float64}}(undef, length(sim))
    byframes = deepcopy(densities)
    for (iframe, frame) in enumerate(sim)
        coords = [ coord[dim] for coord in MolSimToolkit.positions(frame)[iatoms] ]
        ibin = 1
        while ibin < length(bins)
            idx = findall(
                coord -> bins[ibin] <= coord <= bins[ibin+1],
                coords
            )
            density[ibin] = !isempty(idx) ? V \ sum(abs.(electrons[idx])) : 0.0
            ibin += 1
        end
        densities[iframe] = deepcopy(density)
        if !hascenter
            byframes[iframe] = avgbins(bins)
        else
            avg = mean(coords)
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

function egrid(coords, resolution::Float64)
    xlower, xupper = extrema([ coord[1] for coord in coords ])
    ylower, yupper = extrema([ coord[2] for coord in coords ])
    zlower, zupper = extrema([ coord[3] for coord in coords ])
    points = SVector{3, Float64}[]
    for x in xlower:resolution:xupper, y in ylower:resolution:yupper, z in zlower:resolution:zupper
        push!(points, SVector(x, y, z))
    end
    return points
end

function ρ(coords, q, N; σ=nothing, resolution=0.6)
    σ = isnothing(σ) ? σ_edp(length(coords), N) : σ
    A = q ./ (sqrt(2π) .* σ).^3
    α = inv(-2σ^2)
    ρ(r) = sum(A .* exp.(α .* sum((r .- coord).^2 for coord in coords)))
    return ρ.(grid(coords, resolution))
end

function σ_edp(m::Float64, V::Float64)
    return sqrt(
        (3 * V / 4π) / m
    )
end

function edp_1d_gaussian(
    sim::MolSimToolkit.Simulation, e::Vector{Float64};
    axis="z", cutoff=1.0, σ=1.0, hascenter=true, hassymmetry=true
)
    axescode = Dict{String, Int64}("x" => 1, "y" => 2, "z" => 3)
    if haskey(axescode, lowercase(axis))
        dim = axescode[axis]
        V, bins = binning(sim, axis=axis, cutoff=cutoff)
        density = Vector{Float64}(undef, length(bins)-1)
    else
        throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
    end
    iatoms = PDBTools.index.(sim.atoms)    
    electrons = e[iatoms]
    densities = Vector{Vector{Float64}}(undef, length(sim))
    byframes = deepcopy(densities)
    for (iframe, frame) in enumerate(sim)
        coords = MolSimToolkit.positions(frame)[iatoms]
        A = electrons ./ (sqrt(2π) * σ)^3
        α = inv(-2σ^2)
        ρ(r) = sum(
            A .* exp.(α .* [ sum((r .- coord).^2) for coord in coords ])
        )
        points = egrid(coords, cutoff)
        ρbins = ρ.(points)
        ibin = 1
        while ibin < length(bins)
            idx = findall(
                coord -> bins[ibin] <= coord[dim] <= bins[ibin+1],
                points
            )
            density[ibin] = !isempty(idx) ? V \ sum(abs.(ρbins[idx])) : 0.0
            ibin += 1
        end
        densities[iframe] = deepcopy(density)
        if !hascenter
            byframes[iframe] = avgbins(bins)
        else
            avg = mean([ coord[dim] for coord in coords ])
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
#     coords = MolSimToolkit.positions(frame)[iatoms]
#     if isdimensionless 
#         ibin = 1
#         while ibin < length(bins)
#             idx = findall(coord -> bins[ibin] <= coord[dim] <= bins[ibin+1], coords)
#             ρtemp[ibin] = ρ(electrons[idx], V) : ρ(coords[idx], electrons[idx], V, resolution=resolution)
#             ibin += 1
#         end
#     end
#     densities[iframe] = deepcopy(ρtemp)
#     if !hascenter
#         byframes[iframe] = avgbins(bins)
#     else
#         avg = mean(coords[iatoms])
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
# return byframes, densities

function ρ(q::Vector{Float64}, N::Float64)
    return N \ sum(abs.(q))
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


"""
    chargesPSF(psfname::String)

Extract the partial charges from a PSF file.
"""
function chargesPSF(psfname::String; strfield=7)
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
                qi = string(split(line)[strfield])
                push!(q, parse(Float64, qi))
                natoms -= 1
            end
        end
    end
    return q
end

"""
    align_frames(psfname::String, trjname::String; new_trajectory=nothing, selection="not water", reference=0, vmd="vmd", DebugVMD=false)

Align the frames of a trajectory using the VMD `measure fit` command. The selection can be defined by the user,
and the default is `not water`. The reference frame can be defined by the user, and the default is `0`.

### Arguments
- `psfname::String`: The name of the PSF file.
- `trjname::String`: The name of the trajectory file. It can be in any format that VMD can read.
- `new_trajectory::String=nothing`: The name of the new DCD trajectory file. The default creates a temporary file.
- `selection::String="not water"`: The selection to be used to align the frames.
- `reference=0`: The reference frame to be used to align the frames.
- `vmd="vmd"`: The VMD executable. The default is `vmd`.
- `DebugVMD=false`: If `true`, the output of VMD will be printed.
"""
function align_frames(
    psfname::String, trjname::String; new_trajectory=nothing,
    selection="not water", reference::Int64=0,
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

# function edp_1d(
#     sim::MolSimToolkit.Simulation,
#     q::Vector{Float64};
#     axis="z", cutoff=1.0, resolution=0.6,
#     isdimensionless=true, hascenter=true, hassymmetry=true
# )
#     densities = Vector{Vector{Float64}}(undef, length(sim))
#     byframes = Vector{Vector{Float64}}(undef, length(sim))
#     axesdecode = Dict{String, Vector{Int64}}("x" => [1], "y" => [2], "z" => [3])
#     dim = if haskey(axesdecode, lowercase(axis))
#             axescode[axis]
#         else
#             throw(ArgumentError("The axis should be associated with the labels of cartesian axes."))
#     end
#     hascenter = hassymmetry ? true : hascenter
#     iatoms = hascenter ? PDBTools.index.(sim.atoms) : nothing
#     V, bins = binning(sim, axis=axis, cutoff=cutoff)
#     ρtemp = zeros(Float64, length(bins)-1)
#     for (iframe, frame) in enumerate(sim)
#         coords = MolSimToolkit.positions(frame)
#         ibin = 1
#         while ibin < length(bins)
#             idx = findall(coord -> bins[ibin] <= coord[dim] <= bins[ibin+1], coords)
#             ρtemp[ibin] = isdimensionless ? ρ(q[idx], V) : ρ(coords[idx], q[idx], V, resolution=resolution)
#             ibin += 1
#         end
#         densities[iframe] = deepcopy(ρtemp)
#         if !hascenter
#             byframes[iframe] = avgbins(bins)
#         else
#             avg = mean(coords[iatoms])
#             if hassymmetry
#                 byframes[iframe] = _fix_binning_symmetry(bins, avg)
#                 if !issorted(byframes[iframe])
#                     isorted = sortperm(byframes[iframe])
#                     byframes[iframe], densities[iframe] = byframes[iframe][isorted], densities[iframe][isorted]
#                 end
#             else
#                 byframes[iframe] = _fix_binning_center(bins, avg)
#             end
#         end
#     end
#     return byframes, densities
# end

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

# function ρbins(bins::Vector{Float64}, coords::Vector{Float64}, q::Vector{Float64}; N=nothing, σ=nothing)
#     N = isnothing(N) ? length(coords) : N
#     σ = isnothing(σ) ? ones(length(coords)) : σ
#     density = zeros(Float64, length(bins))
#     for (i, bin) in enumerate(bins)
#         if lastindex(bins) == i
#             break
#         end
#         δ = δbin.(coords, bin, bins[i+1])
#         electron = eprofile(coords, q, σ)
#         density[i] = N \ sum(δ .* electron)
#     end
#     return density
# end

# function eprofile(coords::Vector{Float64}, q::Vector{Float64}, σ::Vector{Float64}; isdimensionless=true)
#     if isdimensionless
#         return q
#     end
#     Ng = inv.(
#         sqrt.(2π) .* σ
#     ) ## gaussian normalization
#     score = (coords .- mean(coords)) ./ σ
#     return q .* Ng .* exp.(-0.5 .* score.^2)
# end

