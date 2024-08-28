## FUNCTION 01 - DISTANCE BETWEEN WATER AND LIGNIN  
function water_lignin(pdbfile::String, trajectory::String, segment::String; nonhydrogen=:false, ffirst=1, flast=nothing, fstep=1, dist_cutoff=20.0)
    
    println(" ~~ Loading the PDB and TRAJECTORY information..."); println("")
    simulation = MolSimToolkit.Simulation(pdbfile, trajectory; first=ffirst, last=flast, step=fstep)

    println(" ~~ Setting the monitored and reference residues..."); println("")
    water_selection = PDBTools.readPDB(pdbfile, only = atom -> (atom.resname == "TIP3")); imonitored = index.(water_selection)
    lignin_selection = PDBTools.readPDB(pdbfile, only = atom -> (atom.segname == segment)); ireference = index.(lignin_selection)
    if nonhydrogen
        water_selection = PDBTools.select(water_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H2")
        )
        lignin_selection = PDBTools.select(lignin_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H2" && atom.name != "H3" && atom.name != "H4" && atom.name != "H5" && atom.name != "H6" && atom.name != "H7" &&
            atom.name != "H8" && atom.name != "H91" && atom.name != "H92" && atom.name != "H93" && atom.name != "H101" && atom.name != "H102" && atom.name != "H103" &&
            atom.name != "H111" && atom.name != "H112" && atom.name != "H113" && atom.name != "HO4" && atom.name != "HO7" && atom.name != "HO9")
        )
    end
    namemonitored = name.(water_selection); resnummonitored = resnum.(water_selection); segnamemonitored = segname.(water_selection)
    namereference = name.(lignin_selection); segnamereference = segname.(lignin_selection)

    println(" ~~ Calculating the shortest distances between the monitored atoms and the reference atoms:")
    
    water_lignin = [] ## empty array to store the minimum distances
    open("./water_$(segment).dat", "w") do dat
        ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
        for frame in simulation
            ith_frame = simulation.frame_index; println("    - frame $ith_frame")
            coor = MolSimToolkit.positions(frame)
            uc = diag(unitcell(frame))
            water_lignin = minimum_distances(
                xpositions = [ SVector(coor[i]) for i in imonitored ], # solvent - the monitored atoms around the reference (e.g. xylan and mannans)
                ypositions = [ SVector(coor[i]) for i in ireference ], # solute  - the reference atoms (e.g. cellulsoe microfibril)
                xn_atoms_per_molecule = 3,
                cutoff = dist_cutoff,
                unitcell = uc
            )
            ## wrinting the output dat file from water_lignin variable
            j = 0; for resid in water_lignin
                if (resid.i == 0) && (resid.j == 0)
                    wfile_3, wfile_4, wfile_5 = NaN, NaN, NaN
                    wfile_6, wfile_7 = NaN, NaN, NaN
                    wfile_8 = NaN
                else
                    wfile_3, wfile_4, wfile_5 = namemonitored[resid.i], resnummonitored[resid.i], segnamemonitored[resid.i]
                    wfile_6, wfile_7 = namereference[resid.j], segnamereference[resid.j]
                    wfile_8 = resid.d
                end; j += 1; write(dat, "$ith_frame $wfile_3 $wfile_4 $wfile_5 $wfile_6 $wfile_7 $wfile_8\n")
            end
        end
    end; println(""); println(" ~~ Done!"); return nothing

end

## FUNCTION 02 - DISTANCE BETWEEN WATER AND CENTER OF 2-LIGNINS
function water_lignin(pdbfile::String, trajectory::String, segment1::String, segment2::String; nonhydrogen=:false, ffirst=1, flast=nothing, fstep=1, dist_cutoff=20.0)
    
    println(" ~~ Loading the PDB and TRAJECTORY information..."); println("")
    simulation = MolSimToolkit.Simulation(pdbfile, trajectory; first=ffirst, last=flast, step=fstep)

    println(" ~~ Setting the monitored and reference residues..."); println("")
    water_selection = PDBTools.readPDB(pdbfile, only = atom -> (atom.resname == "TIP3")); imonitored = index.(water_selection)
    lignin_selection = PDBTools.readPDB(pdbfile, only = atom -> (atom.segname == segment1 || atom.segname == segment2)); ireference = index.(lignin_selection)
    if nonhydrogen
        water_selection = PDBTools.select(water_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H2")
        )
        lignin_selection = PDBTools.select(lignin_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H2" && atom.name != "H3" && atom.name != "H4" && atom.name != "H5" && atom.name != "H6" && atom.name != "H7" &&
            atom.name != "H8" && atom.name != "H91" && atom.name != "H92" && atom.name != "H93" && atom.name != "H101" && atom.name != "H102" && atom.name != "H103" &&
            atom.name != "H111" && atom.name != "H112" && atom.name != "H113" && atom.name != "HO4" && atom.name != "HO7" && atom.name != "HO9")
        )
    end
    namemonitored = name.(water_selection); resnummonitored = resnum.(water_selection); segnamemonitored = segname.(water_selection)
    namereference = name.(lignin_selection); segnamereference = segname.(lignin_selection)

    println(" ~~ Calculating the shortest distances between the monitored atoms and the reference atoms:")
    water_lignin = [] ## empty array to store the minimum distances
    open("./water_$(segment1)and$(segment2).dat", "w") do dat
        ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
        for frame in simulation
            ith_frame = simulation.frame_index; println("    - frame $ith_frame")
            coor = MolSimToolkit.positions(frame)
            uc = diag(unitcell(frame))
            water_lignin = minimum_distances(
                xpositions = [ SVector(coor[i]) for i in imonitored ], # solvent - the monitored atoms around the reference (e.g. xylan and mannans)
                ypositions = [ SVector(coor[i]) for i in ireference ], # solute  - the reference atoms (e.g. cellulsoe microfibril)
                xn_atoms_per_molecule = 3,
                cutoff = dist_cutoff,
                unitcell = uc
            )
            ## wrinting the output dat file from water_lignin variable
            j = 0; for resid in water_lignin
                if (resid.i == 0) && (resid.j == 0)
                    continue
                else
                    wfile_3, wfile_4, wfile_5 = namemonitored[resid.i], resnummonitored[resid.i], segnamemonitored[resid.i]
                    wfile_6, wfile_7 = namereference[resid.j], segnamereference[resid.j]
                    wfile_8 = resid.d; write(dat, "$ith_frame $wfile_3 $wfile_4 $wfile_5 $wfile_6 $wfile_7 $wfile_8\n")
                end; j += 1;
            end
        end
    end; println(""); println(" ~~ Done!"); return nothing

end

## FUNCTION 03 - DISTANCE BETWEEN THE CENTER OF LIGNIN-HEMICELLULOSE BONDED RESIDUES
function water_lcc(pdbfile::String, trajectory::String, segment1::String, segment2::String; nonhydrogen=:false, ffirst=1, flast=nothing, fstep=1, dist_cutoff=20.0)
    
    println(" ~~ Loading the PDB and TRAJECTORY information..."); println("")
    simulation = MolSimToolkit.Simulation(pdbfile, trajectory; first=ffirst, last=flast, step=fstep)

    println(" ~~ Setting the monitored and reference residues..."); println("")
    water_selection = PDBTools.readPDB(pdbfile, only = atom -> (atom.resname == "TIP3")); imonitored = index.(water_selection)
    segment1_options = split(segment1, " "); segment1_segname = segment1_options[1]; segment1_resnum = parse(Int64, segment1_options[2])
    segment2_options = split(segment2, " "); segment2_segname = segment2_options[1]; segment2_resnum = parse(Int64, segment2_options[2])
    lcc_selection = PDBTools.readPDB(pdbfile, only = atom -> (
            (atom.segname == segment1_segname && atom.resnum == segment1_resnum) || (atom.segname == segment2_segname && atom.resnum == segment2_resnum)
        )); ireference = index.(lcc_selection)
    if nonhydrogen
        water_selection = PDBTools.select(water_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H2")
        )
        lcc_selection = PDBTools.select(lignin_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H2" && atom.name != "H3" && atom.name != "H4" && atom.name != "H5" && atom.name != "H6" && atom.name != "H7" &&
            atom.name != "H8" && atom.name != "H91" && atom.name != "H92" && atom.name != "H93" && atom.name != "H101" && atom.name != "H102" &&
            atom.name != "H103" && atom.name != "H111" && atom.name != "H112" && atom.name != "H113" && atom.name != "HO4" && atom.name != "HO7" && atom.name != "HO9" &&
            atom.name != "H1" && atom.name != "H2" && atom.name != "H3" && atom.name != "H4" && atom.name != "H5" && atom.name != "H6" && atom.name != "HO1" &&
            atom.name != "HO2" && atom.name != "HO3" && atom.name != "HO4" && atom.name != "HO5" &&  atom.name != "HO6" && atom.name != "H52" &&
            atom.name != "H51" && atom.name != "H61" && atom.name != "H62" && atom.name != "HB1" && atom.name != "HB2" && atom.name != "HB3")
        )
    end
    namemonitored = name.(water_selection); resnummonitored = resnum.(water_selection); segnamemonitored = segname.(water_selection)
    namereference = name.(lcc_selection); segnamereference = segname.(lcc_selection)

    println(" ~~ Calculating the shortest distances between the monitored atoms and the reference atoms:")
    water_lcc = [] ## empty array to store the minimum distances
    open("./water_$(segment1_segname)and$(segment2_segname).dat", "w") do dat
        ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
        for frame in simulation
            ith_frame = simulation.frame_index; println("    - frame $ith_frame")
            coor = MolSimToolkit.positions(frame)
            uc = diag(unitcell(frame))
            water_lcc = minimum_distances(
                xpositions = [ SVector(coor[i]) for i in imonitored ], # solvent - the monitored atoms around the reference (e.g. xylan and mannans)
                ypositions = [ SVector(coor[i]) for i in ireference ], # solute  - the reference atoms (e.g. cellulsoe microfibril)
                xn_atoms_per_molecule = 3,
                cutoff = dist_cutoff,
                unitcell = uc
            )
            ## wrinting the output dat file from water_lignin variable
            j = 0; for resid in water_lcc
                if (resid.i == 0) && (resid.j == 0)
                    continue
                else
                    wfile_3, wfile_4, wfile_5 = namemonitored[resid.i], resnummonitored[resid.i], segnamemonitored[resid.i]
                    wfile_6, wfile_7 = namereference[resid.j], segnamereference[resid.j]
                    wfile_8 = resid.d; write(dat, "$ith_frame $wfile_3 $wfile_4 $wfile_5 $wfile_6 $wfile_7 $wfile_8\n")
                end; j += 1;
            end
        end
    end; println(""); println(" ~~ Done!"); return nothing

end

## ----------------------------------------

## FUNCTION 04 - DISTANCE BETWEEN THE WATER AND A DUMMY POLYMER
function water_segment(pdbfile::String, trajectory::String, segment::String; nonhydrogen=:false, ffirst=1, flast=nothing, fstep=1, dist_cutoff=20.0)
    
    println(" ~~ Loading the PDB and TRAJECTORY information..."); println("")
    simulation = MolSimToolkit.Simulation(pdbfile, trajectory; first=ffirst, last=flast, step=fstep)

    println(" ~~ Setting the monitored and reference residues..."); println("")
    water_selection = PDBTools.readPDB(pdbfile, only = atom -> (atom.resname == "TIP3"))
    polymer_selection = PDBTools.readPDB(pdbfile, only = atom -> (atom.segname == segment))

    if nonhydrogen
        water_selection = PDBTools.select(water_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H2")
        )
        polymer_selection = PDBTools.select(polymer_selection, by = atom -> (
            atom.name != "H1" && atom.name != "H101" && atom.name != "H102" && atom.name != "H103" && atom.name != "H111" && atom.name != "H112" &&
            atom.name != "H113" && atom.name != "H2" && atom.name != "H3" && atom.name != "H4" && atom.name != "H5" && atom.name != "H51" && atom.name != "H52" &&
            atom.name != "H6" && atom.name != "H61" && atom.name != "H62" && atom.name != "H7" && atom.name != "H8" && atom.name != "H91" && atom.name != "HB1" &&
            atom.name != "HB2" && atom.name != "HB3" && atom.name != "HO1" && atom.name != "HO2" && atom.name != "HO3" && atom.name != "HO4" && atom.name != "HO5" &&
            atom.name != "HO6" && atom.name != "HO7" && atom.name != "HO9")
        )
    end
    namemonitored = name.(water_selection); resnummonitored = resnum.(water_selection); segnamemonitored = segname.(water_selection)
    namereference = name.(polymer_selection); resnumreference = resnum.(polymer_selection); segnamereference = segname.(polymer_selection)
    imonitored = index.(water_selection); ireference = index.(polymer_selection);

    println(" ~~ Calculating the shortest distances between the monitored atoms and the reference atoms:")
    
    water_polymer = [] ## empty array to store the minimum distances
    open("./water_$(segment).dat", "w") do dat
        ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
        for frame in simulation
            ith_frame = simulation.frame_index; println("    - frame $ith_frame")
            coor = MolSimToolkit.positions(frame)
            uc = diag(unitcell(frame))
            water_polymer = minimum_distances(
                xpositions = [ SVector(coor[i]) for i in imonitored ], # solvent - the monitored atoms around the reference (e.g. water)
                ypositions = [ SVector(coor[i]) for i in ireference ], # solute  - the reference atoms (e.g. polymeric matrix)
                xn_atoms_per_molecule = 3,
                cutoff = dist_cutoff,
                unitcell = uc
            )
            ## wrinting the output dat file from water_lignin variable
            j = 0; for resid in water_polymer
                if (resid.i == 0) && (resid.j == 0)
                    continue
                else
                    wfile_3, wfile_4, wfile_5 = namemonitored[resid.i], resnummonitored[resid.i], segnamemonitored[resid.i]
                    wfile_6, wfile_7, wfile_8 = namereference[resid.j], resnumreference[resid.j], segnamereference[resid.j]
                    wfile_9 = resid.d
                end; j += 1; write(dat, "$ith_frame $wfile_3 $wfile_4 $wfile_5 $wfile_6 $wfile_7 $wfile_8 $wfile_9\n")
            end
        end
    end; println(""); println(" ~~ Done!"); return nothing

end