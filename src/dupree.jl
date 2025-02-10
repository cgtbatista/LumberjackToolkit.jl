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

function teste2(;data="/home/carlos/Documents/Sandbox/paul/cellulose.dat", nbins=100)
    xyz = readdlm(data)
    x, y, z = xyz[:,1], xyz[:,2], xyz[:,3]
    zmin, zmax = extrema(z)
    bin_edges = range(
        zmin,
        zmax,
        length=nbins+1
    )
    diameters = Float64[]
    for i in 1:nbins
        lower = bin_edges[i]
        upper = i < nbins ? bin_edges[i+1] : zmax + eps()
        mask = (z .>= lower) .& (z .< upper)
        slice_x = x[mask]
        slice_y = y[mask]
        length(slice_x) < 2 && continue ## skip if there are less than 2 points
        d = 0.0
        n_pontos = length(slice_x)
        for j in 1:n_pontos
            for k in (j+1):n_pontos
                dist = hypot(slice_x[j] - slice_x[k], slice_y[j] - slice_y[k])
                d = max(d, dist)
            end
        end
        
        push!(diameters, d)
    end
    
    # Gerar histograma
    Plots.histogram(diameters, bins=20, xlabel="Diâmetro", ylabel="Frequência",
             title="Distribuição de Diâmetros", legend=false)
    
    return (d_mean = mean(diameters),
            d_media = median(diameters),
            d_sd = std(diameters),
            hist = Plots.current())
end

function diameter_analysis(;fibrildata="/home/carlos/Documents/Sandbox/paul/m23432_data.dat")
    data = readdlm(fibrildata, ',')
    Θ, diameter, dresidues = [], [], convert(Vector{Int64}, data[:,1])
    dΘ = missing
    for i in unique(dresidues)
        idx = findall(x -> x == i, dresidues)
        radii_info = data[idx,2]
        angle_info = data[idx,3]
        for idx_angle in eachindex(angle_info)
            if angle_info[idx_angle] < 0.0
                angle_info[idx_angle] = angle_info[idx_angle] + 360.0
            end
        end
        dΘ = std(diff(angle_info))
        for j in eachindex(idx)
            tmp_radius = radii_info[j]; tmp_angle = angle_info[j];
            if tmp_angle < 180.0
                reciprocal_angle = tmp_angle + 180.0
            else
                reciprocal_angle = tmp_angle - 180.0
            end
            lower_bound = reciprocal_angle-dΘ
            upper_bound = reciprocal_angle+dΘ
            nearest_angles = findall(x ->
                (x >= lower_bound) && (x <= upper_bound), angle_info)
            median_idx = Int64(round(median(nearest_angles)))
            tmp_radius2 = radii_info[median_idx]
            push!(Θ, round(tmp_angle, digits=3))
            push!(diameter, round(tmp_radius+tmp_radius2, digits=3))
        end
    end

    Θ = convert(Vector{Float64}, Θ)
    diameter = convert(Vector{Float64}, diameter)

    return dresidues, diameter, Θ

end

function chain_centers(pdbname::String; Δz::Float64=0.1, cutoff::Float64=0.5, isreference::Bool=true)
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
    for i in 0:Δz:norm_displace
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
    Δz=1.0, Δθ=1.0, tol=0.5, isreference=false,
    psfname=nothing, vdw_radii=nothing,
    surface_parameters="-res 0.6 -cutoff 4.0", vmd="vmd"
)
    stl = cellulose_surface(
        pdbname, psfname=psfname, vdw_radii=vdw_radii, parameters=surface_parameters, vmd=vmd
    )[2]
    fibril_centers = chain_centers(
        pdbname, Δz=Δz, cutoff=0.5*Δz, isreference=isreference
    )
    if 2*tol > Δz
        throw(ArgumentError("The tolerance should be less than half of the Δz."))
    end
    xyz = filterSTL(stl)
    return scan_diameter(xyz, fibril_centers, Δθ=Δθ, tol=tol)
end

function scan_diameter(xyz::Vector{SVector{3, Float64}}, centers::Vector{SVector{3, Float64}}; Δθ=1.0, tol=0.05)
    d = Float64[]
    for center in centers
        isurface = findall(k -> norm(k - center[3]) <= tol, [ ijk[3] for ijk in xyz ])
        surface = [ surf - center for surf in xyz[isurface] ]
        r = norm.(surface)
        x, y, z = [ s[1] for s in surface ], [ s[2] for s in surface ], [ s[3] for s in surface ]
        norm_xy = norm.([ s[1:2] for s in surface ])
        θ = sign.(y) .* acosd.(x ./ norm_xy)
        ϕ = acosd.(z ./ r)
        θ, ϕ = mod.(θ, 360.), mod.(ϕ, 360.)
        append!(d, getwidth(r, θ, ϕ, Δθ=Δθ, tol=tol))
    end
    return d
end

function getwidth(r::Vector{Float64}, θ::Vector{Float64}, ϕ::Vector{Float64}; Δθ=1.0, tol=0.05)
    radii = Float64[]
    for angle in 0.0:Δθ:360.0
        θi, θf = mod(angle - 0.5*Δθ, 360.0), mod(angle + 0.5*Δθ, 360.0)
        idxθ = θi < θf ? findall(i -> (i >= θi) && (i <= θf), θ) : findall(i -> (i >= θi) || (i <= θf), θ)
        radius, theta, phi = Vector{Float64}(undef, length(idxθ)), Vector{Float64}(undef, length(idxθ)), Vector{Float64}(undef, length(idxθ))
        for idx in idxθ
            ϕ1, ϕ2 = mod(ϕ[idx] - tol, 360.0), mod(ϕ[idx] + tol, 360.0)
            idxϕ = ϕ1 < ϕ2 ? findall(i -> (i >= ϕ1) && (i <= ϕ2), ϕ) : findall(i -> (i >= ϕ1) || (i <= ϕ2), ϕ)
            if isempty(idxϕ)
                continue
            end
            iradii = argmax(r[idxϕ])
            push!(radii, r[iradii])
            push!(theta, θ[i])
            push!(phi, ϕ[i])
        end
    end
    return d
end

# dΘ = std(diff(angle_info))
# for j in eachindex(idx)
#     tmp_radius = radii_info[j]; tmp_angle = angle_info[j];
#     if tmp_angle < 180.0
#         reciprocal_angle = tmp_angle + 180.0
#     else
#         reciprocal_angle = tmp_angle - 180.0
#     end
#     lower_bound = reciprocal_angle-dΘ
#     upper_bound = reciprocal_angle+dΘ
#     nearest_angles = findall(x ->
#         (x >= lower_bound) && (x <= upper_bound), angle_info)
#     median_idx = Int64(round(median(nearest_angles)))
#     tmp_radius2 = radii_info[median_idx]
#     push!(Θ, round(tmp_angle, digits=3))
#     push!(diameter, round(tmp_radius+tmp_radius2, digits=3))
# end

function fibrilradii(
    pdb::String,
    file::String;
    Δz=0.01,
    steps=false, Δstep=0.1
)
    if (steps == true) && (Δstep == 0.0)
        throw(ArgumentError("The Δstep must be greater than zero."))
    end
    if Δstep <= Δz
        throw(ArgumentError("The Δstep ($Δstep Å) should be greater than Δz ($Δz Å) to avoid duplicate radius quantification."))
    end
    atoms = PDBTools.readPDB(pdb)
    lower_cutoff = mean([
            atom.z for atom in PDBTools.select(atoms, atom -> atom.resnum == minimum(PDBTools.resnum.(atoms)))
        ])
    upper_cutoff = mean([
            atom.z for atom in PDBTools.select(atoms, atom -> atom.resnum == maximum(PDBTools.resnum.(atoms)))
        ])
    data = readdlm(file)
    xyz = unique(
            [ [irow[1], irow[2], irow[3]] for irow in eachrow(data) ]
        )
    zaxis = [ ijk[3] for ijk in xyz ]
    xyz, zaxis, ztoken = xyz[sortperm(zaxis)], zaxis[sortperm(zaxis)], 0.
    r = Float64[]; ϕ = Float64[]; θ = Float64[]; central_vectors = Vector{Float64}[]; labels = Int64[];
    for (i, z) in enumerate(zaxis)
        if (z < lower_cutoff) || (z > upper_cutoff)
            continue
        end
        if any(test -> test == true, steps)
            if (i != 1)
                if abs(z-ztoken) < Δstep
                    continue
                else
                    ztoken = z
                end
            else
                ztoken = z
            end
        end
        Δatoms = PDBTools.select(atoms, atom -> (atom.z >= z - Δz) && (atom.z <= z + Δz))
        xyz0 = [ mean([ xyz[1] for xyz in PDBTools.coor(Δatoms) ]), mean([ xyz[2] for xyz in PDBTools.coor(Δatoms) ]), mean([ xyz[3] for xyz in PDBTools.coor(Δatoms) ]) ]
        if isnan(mean(xyz0)); continue; end
        push!(central_vectors, xyz0)
        for idx in findall(z -> (z >= z - Δz) && (z <= z + Δz), zaxis)
            ## picking {r, θ, ϕ} -- following the physical spherical coordinates classification as { radius, azimuthal, polar }
            vector_xyz = xyz[idx] - xyz0
            idx_r = norm(vector_xyz)
            idx_θ = sign(vector_xyz[2])*acosd(vector_xyz[1]/norm(vector_xyz[1:2]))
            idx_ϕ = acosd(vector_xyz[3]/norm(vector_xyz))
            if idx_θ < 0; idx_θ += 360; end
            if idx_ϕ < 0; idx_ϕ += 360; end
            ## pushing the data
            push!(r, idx_r); push!(θ, idx_θ); push!(ϕ, idx_ϕ); push!(labels, i)
        end
    end
    return xyz, central_vectors, labels, [ r, θ, ϕ ]
end

# if centermethods == "real"
#     selstructure = select(structure, by = atoms -> (atoms.z >= lower_bound) && (atoms.z <= upper_bound))
#     xyz0 = [ mean([ x[1] for x in coor(selstructure) ]), mean([ y[2] for y in coor(selstructure) ]), mean([ z[3] for z in coor(selstructure) ]) ];        
# elseif centermethods == "imaginary"
#     pdbresids = resnum.(structure); coordinates = coor(structure);
#     x = [ i[1] for i in coordinates ]; y = [ j[2] for j in coordinates ]; z = [ k[3] for k in coordinates ];
#     initial_idx = findall(resids -> resids == minimum(pdbresids), pdbresids); x0 = mean(x[initial_idx]); y0 = mean(y[initial_idx]); z0 = mean(z[initial_idx]);
#     final_idx   = findall(resids -> resids == maximum(pdbresids), pdbresids); x1 = mean(x[final_idx]); y1 = mean(y[final_idx]); z1 = mean(z[final_idx]);
#     vectorx = [ x1-x0, y1-y0, z1-z0 ]; norm_vector = norm(vectorx);
#     dteste = collect(0.:dz:norm_vector); xteste = randn(length(dteste)); yteste = randn(length(dteste)); zteste = randn(length(dteste));
#     for i in eachindex(dteste)
#         xteste[i] = x0 + dteste[i]*vectorx[1]/norm_vector
#         yteste[i] = y0 + dteste[i]*vectorx[2]/norm_vector
#         zteste[i] = z0 + dteste[i]*vectorx[3]/norm_vector
#     end
#     selstructure = findall(ithz -> (ithz >= lower_bound) && (ithz <= upper_bound), zteste)
#     xyz0 = [ mean(xteste[selstructure]), mean(yteste[selstructure]), mean(zteste[selstructure]) ]
# else
#     throw(ArgumentError("The method $method is not implemented."))
# end


function get_diameters(radiidata::Vector{Vector{Float64}}; threshold=1.5, tolerance=0.05)
    ## The tolerance is related to the angle θ and the marginal error to be outsisde the plane of the fibril is related to the angle ϕ
    dθ = threshold; dϕ = tolerance;
    diameters = Float64[]
    for ith in eachindex(radiidata[1])
        r = radiidata[1][ith]; θ = radiidata[2][ith]; ϕ = radiidata[3][ith];
        new_θ = θ + 180.; lower_bound = mod(new_θ - dθ, 360.); upper_bound = mod(new_θ + dθ, 360.);
        if lower_bound < upper_bound
            raw_idx = findall(θi -> θi >= lower_bound && θi <= upper_bound, radiidata[2])
        else
            raw_idx = findall(θi -> θi >= lower_bound || θi <= upper_bound, radiidata[2])
        end
        idx = findall(ϕi -> ϕi >= mod(ϕ - dϕ, 360) && ϕi <= mod(ϕ + dϕ, 360), radiidata[3][raw_idx])
        if length(idx) == 0
            continue
        else
            inv_r = radiidata[1][idx]; inv_θ = radiidata[2][idx]; inv_ϕ = radiidata[3][idx]
            distances = @. sqrt((r)^2 + (inv_r)^2 - 2*(r)*(inv_r)*(sind(inv_ϕ)*sind(ϕ)*cosd(inv_θ-θ)+cosd(inv_ϕ)*cosd(ϕ)))
            push!(diameters, median(distances))
        end
    end
    return diameters
end
