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

function fibrilradii(
    pdb::String,
    file::String;
    rounding_xyz=3,
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
    lower_cutoff = mean([ atom.z for atom in PDBTools.select(atoms, atom -> atom.resnum == minimum(PDBTools.resnum.(atoms))) ])
    upper_cutoff = mean([ atom.z for atom in PDBTools.select(atoms, atom -> atom.resnum == maximum(PDBTools.resnum.(atoms))) ])
    data = readdlm(file)
    rawxyz = [ 
        [ round(
            i[1], digits=rounding_xyz), round(i[2], digits=rounding_xyz), round(i[3], digits=rounding_xyz
        ) ] for i in eachrow(data) ]
    xyz = unique(rawxyz)
    zaxis = [ i[3] for i in xyz ]
    xyz, zaxis = xyz[sortperm(zaxis)], zaxis[sortperm(zaxis)]
    ztoken = 0.
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
        lower_bound = z-Δz; upper_bound = z+Δz;
        Δatoms = PDBTools.select(atoms, atom -> (atom.z >= lower_bound) && (atom.z <= upper_bound))
        xyz0 = [ mean([ xyz[1] for xyz in PDBTools.coor(Δatoms) ]), mean([ xyz[2] for xyz in PDBTools.coor(Δatoms) ]), mean([ xyz[3] for xyz in PDBTools.coor(Δatoms) ]) ]
        if isnan(mean(xyz0)); continue; end
        push!(central_vectors, xyz0)
        for idx in findall(z -> (z >= lower_bound) && (z <= upper_bound), zaxis)
            ## picking {r, θ, ϕ} -- following the physical spherical coordinates classification as { radius, azimuthal, polar }
            vector_xyz = xyz[idx] - xyz0
            idx_r = norm(vector_xyz)
            idx_θ = sign(vector_xyz[2])*acosd(vector_xyz[1]/norm(vector_xyz[1:2]))
            idx_ϕ = acosd(vector_xyz[3]/norm(vector_xyz))
            if idx_θ < 0; idx_θ += 360; end
            if idx_ϕ < 0; idx_ϕ += 360; end
            ## pushing the data
            push!(r, round(idx_r, digits=rounding_xyz)); push!(θ, round(idx_θ, digits=1)); push!(ϕ, round(idx_ϕ, digits=1)); push!(labels, i)
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
