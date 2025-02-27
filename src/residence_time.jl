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

function mapwater(
    pdbname::String, trjname::String;
    first=1, step=1, last=nothing,
    selection="not water", water_selection="water",
    without_hydrogens=true, cutoff=5.0
)
    simulation = MolSimToolkit.Simulation(pdbname, trjname, first=first, step=step, last=last)
    return mapwater(
        simulation,
        selection = selection, water_selection = water_selection,
        without_hydrogens = without_hydrogens,
        cutoff = cutoff
    )
end

function mapwater(
    simulation::MolSimToolkit.Simulation;
    selection="not water",                      # selection of the reference atoms -- it's the solute! It can be the carbohydrates, the ions, etc.
    water_selection="water",                    # selection of the water molecules, ie. the solvent of the PCW matrix
    without_hydrogens=true,                     # if the hydrogens should be removed from the analysis. Very useful to approach some experimental evidences.
    cutoff=5.0                                  # the cutoff distance to consider the water molecule around the reference atoms
)
    reference = PDBTools.select(simulation.atoms, selection) 
    monitored = PDBTools.select(simulation.atoms, water_selection)
    if without_hydrogens
        reference = PDBTools.select(reference, at -> PDBTools.element(at) != "H")
        monitored = PDBTools.select(monitored, at -> PDBTools.element(at) != "H")
        natoms = 1
    else
        natoms = 3
    end
    imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)     
    M = Matrix{Bool}(undef, length(imonitored), length(simulation.frame_range))
    for (iframe, frame) in enumerate(simulation)
        xyz, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))
        mindistances = MolSimToolkit.minimum_distances(
            xpositions = [ SVector(xyz[i]) for i in imonitored ],     # solvent - the monitored atoms around the reference (e.g. water)
            ypositions = [ SVector(xyz[j]) for j in jreference ],     # solute  - the reference atoms (e.g. polymeric matrix)
            xn_atoms_per_molecule = natoms,
            cutoff = cutoff,
            unitcell = uc
        )
        for (iwater, md) in enumerate(mindistances)
           M[iwater, iframe] = md.within_cutoff
        end
    end
    return BitMatrix(M)
end

function mapwater(M1::BitMatrix, M2::BitMatrix, operation::Int64)
    M = M1 + M2
    if operation in Set([0,1,2])
        return ifelse.(M .== operation, true, false)
    else
        throw(ArgumentError("The operation must be 0 (complement), 1 (disjoint), or 2 (intersection)."))
    end
end

function t_residence(M::BitMatrix; timestep=0.1)
    t = Float64[]
    sizehint!(t, size(M,1) * 10)
    for row in eachrow(M)
        trow = checking_residence(
            BitVector(row)
        )
        append!(t, timestep .* trow)
    end
    return filter(!iszero, t)
end

function checking_residence(checklist::BitVector)
    if iszero(sum(checklist))
        return 0
    end
    counting = Int64[]
    counter = 0
    for present in checklist
        if present
            counter += 1
        else
            if counter > 0
                push!(counting, counter)
                counter = 0
            end
        end
    end
    if counter > 0
        push!(counting, counter)
    end
    return counting
end

# function t_residence(M::BitMatrix; timestep=0.1)
#     t = Float64[]
#     for i in 1:size(M,1)
#         append!(
#             t,
#             timestep * checking_residence(M[i,:])
#         )
#     end
#     return t[findall(it -> !iszero(it), t)]
# end

# function checking_residence(checklist::BitVector)
#     if iszero(sum(checklist))
#         return 0
#     end
#     counting, counter = Int64[], 0
#     itsnot = findall(frame -> !frame, checklist)
#     for iframe in 1:length(checklist)
#         if !in(iframe, Set(itsnot))
#             counter += 1
#         else
#             push!(counting, counter)
#             counter = 0
#         end
#         if (iframe == length(checklist)) && (counter != 0)
#             push!(counting, counter)
#         end
#     end
#     return counting
# end

# function mapwater(
#     simulation::MolSimToolkit.Simulation;
#     selection="not water", water_selection="water",
#     without_hydrogens=true,
#     cutoff=10.0,
#     filename=nothing
# )
#     waters, reference = PDBTools.select(simulation.atoms, water_selection), PDBTools.select(simulation.atoms, selection)
#     if without_hydrogens
#         waters = PDBTools.select(waters, at -> PDBTools.element(at) != "H")
#         reference = PDBTools.select(reference, at -> PDBTools.element(at) != "H")
#         water_natoms = 1
#     else
#         water_natoms = 3
#     end
#     iwater, jreference = PDBTools.index.(waters), PDBTools.index.(reference)
#     filename = isnothing(filename) ? tempname() * ".dat" : filename
#     println("""
#     ~~ Calculating the shortest distances for all water molecules and the reference atoms based on $(length(simulation.frame_range)) frames:
#     """)
#     Base.open(filename, "w") do file
#         ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
#         distances = []
#         for frame in simulation
#             iframe = simulation.frame_index
#             println(" frame $iframe")
#             coor, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))
#             distances = MolSimToolkit.minimum_distances(
#                 xpositions = [ SVector(coor[i]) for i in iwater ],     # solvent - the monitored atoms around the reference (e.g. water)
#                 ypositions = [ SVector(coor[i]) for i in jreference ], # solute  - the reference atoms (e.g. polymeric matrix)
#                 xn_atoms_per_molecule = water_natoms,
#                 cutoff = cutoff,
#                 unitcell = uc
#             )
#             for distance in distances
#                 if !iszero([distance.i , distance.j])
#                     resnum_water, segname_water = PDBTools.resnum(waters[distance.i]), PDBTools.segname(waters[distance.i])
#                     println(file, "$iframe $resnum_water $segname_water $(distance.d)")
#                 end
#             end
#         end
#     end
#     println("""

#     ~~ It's done! The information was saved in $(abspath(filename)).
#     """)
#     return filename
# end


# function retention_time(filename::String; nframes = 10, timestep=1, md_time="ns", dist_cutoff=10.0)
    
#     closest_water_time = Dict{Tuple{Int64, String}, Vector{Int64}}()
#     closest_water_dist = Dict{Tuple{Int64, String}, Vector{Float64}}()

#     water_data = readdlm(filename)

#     for irow in eachindex(water_data[:,1])
#         if isnan(water_data[irow, end]) || water_data[irow, end] > dist_cutoff
#             continue
#         end
#         resid, segid = Int64(water_data[irow, 3]), String(water_data[irow, 4]) 
#         key = (resid, segid)
#         if haskey(closest_water_time, key) && haskey(closest_water_dist, key)
#             push!(closest_water_time[key], water_data[irow,1])
#             push!(closest_water_dist[key], water_data[irow,end])
#         else
#             closest_water_time[key] = [water_data[irow,1]]
#             closest_water_dist[key] = [water_data[irow,end]]
#         end
#     end
    
#     ## Residence time
#     water_molecules = collect(1:length(closest_water_time)); unique_water_molecules = [ key for key in keys(closest_water_time) ]
#     filename_label = split(split(filename, ".")[1], "_")[2]
#     #a
#     println("~ $filename_label")
#     println("  $(size(water_molecules)[1]) water molecules in range")
#     filename_array = []; ith_water_array = [];
#     min_rtime_array = []; max_rtime_array = []; contact_array = []; fluctuation_array = []; mean_dist_array = []; std_dist_array = []; closest_dist_array = [];
#     for iwater in eachindex(unique_water_molecules)
#         ikey = unique_water_molecules[iwater]; iwater_label = string( ikey[2], ikey[1])
#         water_time = closest_water_time[ikey]; water_dist = closest_water_dist[ikey]
#         if iwater != 0; ## == 1 just to dummy the simulation
#             println("  $(iwater) - the water molecule $iwater_label appeared $(size(water_time)[1])/$nframes:")
#             min_rtime=0; max_rtime=0; contact_time=0; fluctuation=0;
#             mean_dist=mean(water_dist); std_dist=std(water_dist); closest_dist=minimum(water_dist);
#             el=1; dt=0; lost_contact=0;
#             while el < length(water_time)
#                 if water_time[el+1]-water_time[el] == 1
#                     dt += 1; contact_time += 1; max_rtime += 1;
#                 else
#                     lost_contact += 1
#                     if min_rtime > dt; min_rtime = dt; end
#                     if max_rtime < dt; min_rtime = max_rtime; max_rtime = dt; end
#                     dt = 0
#                 end
#                 if length(water_time)-max_rtime == 1; min_rtime = max_rtime; end
#                 el += 1
#             end
#             max_rtime = max_rtime * timestep; min_rtime = min_rtime * timestep; contact_time = contact_time * timestep; fluctuation = lost_contact;
#             push!(filename_array, filename_label); push!(ith_water_array, iwater_label);
#             push!(min_rtime_array, min_rtime); push!(max_rtime_array, max_rtime); push!(contact_array, contact_time); push!(fluctuation_array, fluctuation);
#             push!(mean_dist_array, mean_dist); push!(std_dist_array, std_dist); push!(closest_dist_array, closest_dist);
#             println("  + the water retention time reached a minimum of $min_rtime $md_time and a maximum of $max_rtime $md_time.")
#             println("  + the water spent $contact_time $md_time on surroundings, but it lost the contact $lost_contact time(s).")
#             println("  + the water molecule mean distance was ($(round(mean_dist, digits=2)) ± $(round(std_dist, digits=2))) Å, with $(round(closest_dist, digits=2)) Å as the closest distance.")
#             println("")
#         end
#         ## ikey = iésima água do vetor unique_water_molecules, onde ikey[2] é o segmento e ikey[1] é o número do resíduo
#     end
#     return filename_array, ith_water_array, min_rtime_array, max_rtime_array, contact_array, fluctuation_array, mean_dist_array, std_dist_array, closest_dist_array
#     #println("- water molecules mean distance (histogram like 2D to see better effect of time in the distribution)")
#     #println("- relantionship between retention time and distance (linear regression??)")       

# end