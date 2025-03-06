## Loading the PDB and TRAJECTORY information
function loading_files()
    pdbfile = "/home/user/Documents/phd/sandbox/pcw.pdb"
    trajectory = "/home/user/Documents/phd/sandbox/equilibration.dcd"
    return pdbfile, trajectory
end

## Setting the scanned residues and segments
function pdb_dummy_selection(pdbfile::String, monitored_segment::String; nonhydrogen=:false)

    reference = PDBTools.readPDB(pdbfile, only = atom -> (
            atom.segname == "A24" || atom.segname == "A23" || atom.segname == "A30" || atom.segname == "A29" || atom.segname == "A28" || atom.segname == "A19" ||
            atom.segname == "A11" || atom.segname == "A12" || atom.segname == "A5"  || atom.segname == "A6"  || atom.segname == "A7"  || atom.segname == "A16" ||
            atom.segname == "B24" || atom.segname == "B23" || atom.segname == "B30" || atom.segname == "B29" || atom.segname == "B28" || atom.segname == "B19" ||
            atom.segname == "B11" || atom.segname == "B12" || atom.segname == "B5"  || atom.segname == "B6"  || atom.segname == "B7"  || atom.segname == "B16" ||
            atom.segname == "C24" || atom.segname == "C23" || atom.segname == "C30" || atom.segname == "C29" || atom.segname == "C28" || atom.segname == "C19" ||
            atom.segname == "C11" || atom.segname == "C12" || atom.segname == "C5"  || atom.segname == "C6"  || atom.segname == "C7"  || atom.segname == "C16" ||
            atom.segname == "D24" || atom.segname == "D23" || atom.segname == "D30" || atom.segname == "D29" || atom.segname == "D28" || atom.segname == "D19" ||
            atom.segname == "D11" || atom.segname == "D12" || atom.segname == "D5"  || atom.segname == "D6"  || atom.segname == "D7"  || atom.segname == "D16" ||
            atom.segname == "E24" || atom.segname == "E23" || atom.segname == "E30" || atom.segname == "E29" || atom.segname == "E28" || atom.segname == "E19" ||
            atom.segname == "E11" || atom.segname == "E12" || atom.segname == "E5"  || atom.segname == "E6"  || atom.segname == "E7"  || atom.segname == "E16" ||
            atom.segname == "F24" || atom.segname == "F23" || atom.segname == "F30" || atom.segname == "F29" || atom.segname == "F28" || atom.segname == "F19" ||
            atom.segname == "F11" || atom.segname == "F12" || atom.segname == "F5"  || atom.segname == "F6"  || atom.segname == "F7"  || atom.segname == "F16"
            )
        )
        
    monitored = PDBTools.readPDB(pdbfile, only = atom -> (
            atom.segname == monitored_segment)
        )

    if nonhydrogen
        reference = PDBTools.select(reference, by = atom -> (
            atom.name != "H1" && atom.name != "H2" && atom.name != "H3" && atom.name != "H4" && atom.name != "H5" && atom.name != "H6" && atom.name != "HO1" &&
            atom.name != "HO2" && atom.name != "HO3" && atom.name != "HO4" && atom.name != "HO5" &&  atom.name != "HO6" && atom.name != "H52" &&
            atom.name != "H51" && atom.name != "H61" && atom.name != "H62" && atom.name != "HB1" && atom.name != "HB2" && atom.name != "HB3")
        )
        monitored = PDBTools.select(monitored, by = atom -> (
            atom.name != "H1" && atom.name != "H2" && atom.name != "H3" && atom.name != "H4" && atom.name != "H5" && atom.name != "H6" && atom.name != "HO1" &&
            atom.name != "HO2" && atom.name != "HO3" && atom.name != "HO4" && atom.name != "HO5" &&  atom.name != "HO6" && atom.name != "H52" &&
            atom.name != "H51" && atom.name != "H61" && atom.name != "H62" && atom.name != "HB1" && atom.name != "HB2" && atom.name != "HB3")
        )            
    end

    ireference = index.(reference); namereference = name.(reference); segnamereference = segname.(reference);
    imonitored = index.(monitored); namemonitored = name.(monitored); segnamemonitored = segname.(monitored); monitored_resnums = resnum.(monitored);
    return ireference, namereference, segnamereference, imonitored, namemonitored, segnamemonitored, monitored_resnums
end

## Evaluating the minimum distances between the monitored atoms and the reference atoms, and writing the output dat file
function mindist(pdbfile::String, trajectory::String, segment::String; ffirst=1, flast=nothing, fstep=1, dist_cutoff=20.0, nonhydrogen=:false)
    
    println(" ~~ Loading the PDB and TRAJECTORY information..."); println("")
    ireference, namereference, segnamereference, imonitored, namemonitored, segnamemonitored, monitored_resnums = pdb_dummy_selection(pdbfile, segment, nonhydrogen=nonhydrogen)
    unique_resnums = unique(monitored_resnums)

    println(" ~~ Setting the monitored residues..."); println("")
    function mol_indices(ith_atom)
        value = monitored_resnums[ith_atom] ## You need to set the value of the monitored_resnums before the function
        selected_residues = unique(monitored_resnums)
        
        findfirst(x -> x == value, selected_residues)
    end
    
    println(" ~~ Calculating the shortest distances between the monitored atoms and the reference atoms:")
    simulation = MolSimToolkit.Simulation(pdbfile, trajectory; first=ffirst, last=flast, step=fstep)
    mindist = [] ## empty array to store the minimum distances
    open("./closest_$(segment).dat", "w") do dat
        ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
        for frame in simulation
            ith_frame = simulation.frame_index; println("    - frame $ith_frame")
            coor = MolSimToolkit.positions(frame)
            uc = diag(unitcell(frame))
            crosspairs_system = CrossPairs(
                xpositions = [ SVector(coor[i]) for i in imonitored ], # solvent - the monitored atoms around the reference (e.g. xylan and mannans)
                ypositions = [ SVector(coor[i]) for i in ireference ], # solute  - the reference atoms (e.g. cellulsoe microfibril)
                xmol_indices = mol_indices,
                cutoff = dist_cutoff,
                unitcell = uc
            ); mindist = minimum_distances!(crosspairs_system)
            ## wrinting the output dat file
            j = 0; for resid in mindist
                if (resid.i == 0) && (resid.j == 0)
                    continue
                else
                    wfile_3, wfile_4 = namemonitored[resid.i], segnamemonitored[resid.i]
                    wfile_5, wfile_6 = namereference[resid.j], segnamereference[resid.j]
                    wfile_7 = resid.d
                end; j += 1; write(dat, "$ith_frame $(getindex(unique_resnums, j)) $wfile_3 $wfile_4 $wfile_5 $wfile_6 $wfile_7\n")
            end
        end
    end; println(""); println(" ~~ Done!"); return nothing

end

function closest2fibril(
    pdbname::String, trjname::String, segcode::String; first=1, last=nothing, step=1,
    selection="resname BGLC", without_hydrogens=true, cutoff=5.0
)   
    simulation = MolSimToolkit.Simulation(pdbname, trjname, first=first, step=step, last=last)
    return closest2fibril(
        simulation,
        segcode,
        selection = selection,
        without_hydrogens = without_hydrogens,
        cutoff = cutoff
    )
    
end

function closest2fibril(simulation::MolSimToolkit.Simulation, segcode::String; selection="resname BGLC", without_hydrogens=true, cutoff=5.0)
    reference = PDBTools.select(simulation.atoms, selection) 
    monitored = PDBTools.select(simulation.atoms, at -> at.segname == segcode)
    if without_hydrogens
        reference = PDBTools.select(reference, at -> PDBTools.element(at) != "H")
        monitored = PDBTools.select(monitored, at -> PDBTools.element(at) != "H")
    end
    imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)
    monitored_resnums = PDBTools.resnum.(monitored)
    molindexes(iatom) = findfirst(x -> x == monitored_resnums[iatom], unique(monitored_resnums))
    iresid, distances = Vector{Vector{Int64}}(undef, length(simulation)), Vector{Vector{Float64}}(undef, length(simulation))
    for (iframe, frame) in enumerate(simulation)
        xyz, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))
        crosspairs = MolSimToolkit.CrossPairs(
            xpositions = [ SVector(xyz[i]) for i in imonitored ], # solvent - the monitored atoms around the reference (e.g. xylan and mannans)
            ypositions = [ SVector(xyz[j]) for j in jreference ], # solute  - the reference atoms (e.g. cellulsoe microfibril)
            xmol_indices = molindexes,
            cutoff = cutoff,
            unitcell = uc
        )
        mindistances = MolSimToolkit.minimum_distances!(crosspairs)
        iresid[iframe] = [ monitored_resnums[md.i] for md in mindistances ]
        distances[iframe] = [ md.d for md in mindistances ]
    end
    return iresid, distances
end

function closest2fibril(files::Vector{String}, first::Int64, last::Int64)
    distances = Vector{Matrix{Float64}}(undef, length(files))
    for (iseg, file) in enumerate(files)
        if filesize(file) == 0
            continue
        end
        distances[iseg] = closest2fibril(file, first, last)
    end
    return distances
end

function closest2fibril(file::String, first::Int64, last::Int64)
    data = readdlm(file)
    nframes = minimum(data[:,1]):maximum(data[:,1])
    distance = zeros(Float64, length(first:last), length(nframes))
    for iframe in nframes
        idx = findall((data[:,1] .== iframe) .& (first .<= data[:,2] .<= last))
        distance[:, iframe] .= data[idx,7]
    end
    return distance
end

function closest2fibril(data1::Vector{Matrix{Float64}}, data2::Vector{Matrix{Float64}}; colors=["coral2", "skyblue2"])
    # Plots.gr(size=(1000,800), dpi=900, fmt=:png)
    # Plots.plot(
    #     title="", legend=:topleft, xlabel="# backbone residue", ylabel="Minimum distances (Å)", fontfamily=:arial,
    #     ## Axis configs
    #     framestyle=:box, grid=true, minorgrid=true, minorticks=5, thick_direction=:out,
    #     ## Font configs
    #     titlefontsize=18, guidefontsize=16, tickfontsize=16, labelfontsize=18, legendfontsize=16, guidefonthalign=:center,
    #     ## Margins
    #     left_margin=5Plots.Measures.mm, right_margin=10Plots.Measures.mm, top_margin=10Plots.Measures.mm, bottom_margin=1Plots.Measures.mm
    # )
    data1 = hcat(data1...)
    data2 = hcat(data2...)
    average1, average2 = mean(data1, dims=2), mean(data2, dims=2)
    #extrema1, extrema2 = extrema(data1, dims=2), extrema(data2, dims=2)
    println("""
    $(round(mean(average1), sigdigits=2)) ± $(round(std(average1), sigdigits=2)) Å
    $(round(mean(average2), sigdigits=2)) ± $(round(std(average2), sigdigits=2)) Å
    """)
    # Plots.plot!(
    #     average1,
    #     #ribbon=(average1 - extrema1[1], extrema1[2] - average1),
    #     #fillalpha=0.3,
    #     color=colors[1],
    #     linewidth=2,
    #     label=label1
    # )
    # Plots.plot!(
    #     average2,
    #     #ribbon=(average2 - extrema2[1], extrema2[2] - average2),
    #     #fillalpha=0.3,
    #     color=colors[2],
    #     linewidth=2,
    #     label=label2
    # )
    # return Plots.current()
end