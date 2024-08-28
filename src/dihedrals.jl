## dihedral

function dihedral(atom_1::SVector, atom_2::SVector, atom_3::SVector, atom_4::SVector)
    # getting the vectors
    v1 = atom_2 - atom_1; v2 = atom_3 - atom_2; v3 = atom_4 - atom_3
    # getting the plane data - cross product of the vectors
    n1 = cross(v1,v2); n2 = cross(v2,v3)
    u1 = cross(n1,n2); unitary_v2 = v2/norm(v2)
    # getting the dihedral angle using the -- atan function
    m1 = dot(u1, unitary_v2)
    m2 = dot(n1,n2)
    dihedral_angle = atan(m1,m2) * (180 / π)
    # standartization of the dihedral angle to be positive
    if dihedral_angle < 0
        dihedral_angle = 360 + dihedral_angle
    end; return dihedral_angle
end

# ERROR!!!
#function dihedral2(atom_1::SVector, atom_2::SVector, atom_3::SVector, atom_4::SVector)
#    # getting the vectors
#    v1 = atom_2 - atom_1; v2 = atom_3 - atom_2; v3 = atom_4 - atom_3
#    # getting the plane data - cross product of the vectors
#    n1 = cross(v1,v2)/norm(cross(v1,v2))
#    n2 = cross(v2,v3)/norm(cross(v2,v3))
#    # getting the dihedral angle using the -- atan function
#    m1 = cross(n1,v2/norm(v2))
#    x = dot(n1,n2)
#    y = dot(m1,n2)
#    dihedral_angle = atan(y,x) * (180 / π)
#    # standartization of the dihedral angle to be positive
#    if dihedral_angle < 0
#        dihedral_angle = 360 + dihedral_angle
#    end; return dihedral_angle
#end

## select_pdb_atom
## aim to get the atom selection by the PDBTools using Julia syntax
function select_pdb_atom(atom, segname, resnum, name)
    if atom.segname == segname && atom.resnum == resnum && atom.name == name
        return true
    else
        return false
    end
end

function hemicellulose_dihedrals(pdbfile::String, trajectory::String, segment::String, first_xylan_residue::Int64, last_xylan_residue::Int64;
    ffirst=1, flast=nothing, fstep=1)
    ##
    println(" ~~ Loading the PDB and TRAJECTORY information..."); println("")
    simulation = MolSimToolkit.Simulation(pdbfile, trajectory; first=ffirst, last=flast, step=fstep)

    println(" ~~ Calculating the shortest distances between the monitored atoms and the reference atoms:")
    segment_atoms = PDBTools.readPDB(pdbfile, only = atom -> (
        atom.segname == segment)
    ); segment_resnums = resnum.(segment_atoms)
    
    println(" ~~ Calculating the dihedrals associated to each backbone residue pair:")

    open("./dihedral_$(segment).dat", "w") do dat
        ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
        for frame in simulation
            ith_frame = simulation.frame_index; println("    - frame $ith_frame")
            coor = MolSimToolkit.positions(frame)
            ## wrinting the output dat file
            resid_token = first_xylan_residue
            while resid_token < last_xylan_residue
                # setting the dummies index variables
                i = resid_token; j = resid_token + 1; resid_token = resid_token + 1
                println("      ... for the residue pair $i and $j")
                # getting the phi dihedral angle -- O5'-C1'-O4-C4
                da1_atom1 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == i && atom.name == "O5"))[1]
                da1_atom2 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == i && atom.name == "C1"))[1]
                da1_atom3 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "O4"))[1]
                da1_atom4 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "C4"))[1]
                ϕ = dihedral(@SVector([ coor[da1_atom1][1], coor[da1_atom1][2], coor[da1_atom1][3] ]),
                             @SVector([ coor[da1_atom2][1], coor[da1_atom2][2], coor[da1_atom2][3] ]),
                             @SVector([ coor[da1_atom3][1], coor[da1_atom3][2], coor[da1_atom3][3] ]),
                             @SVector([ coor[da1_atom4][1], coor[da1_atom4][2], coor[da1_atom4][3] ]))
                # getting the psi dihedral angle -- C1'-O4-C4-C3
                da2_atom1 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == i && atom.name == "C1"))[1]
                da2_atom2 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "O4"))[1]
                da2_atom3 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "C4"))[1]
                da2_atom4 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "C3"))[1] 
                ψ = dihedral(@SVector([ coor[da2_atom1][1], coor[da2_atom1][2], coor[da2_atom1][3] ]),
                             @SVector([ coor[da2_atom2][1], coor[da2_atom2][2], coor[da2_atom2][3] ]),
                             @SVector([ coor[da2_atom3][1], coor[da2_atom3][2], coor[da2_atom3][3] ]),
                             @SVector([ coor[da2_atom4][1], coor[da2_atom4][2], coor[da2_atom4][3] ]))
                # sum of the angles phi + psi
                sum_angles = ϕ + ψ
                if sum_angles > 360 && sum_angles < 720
                    sum_angles = sum_angles - 360
                elseif sum_angles > 720
                    sum_angles = sum_angles - 720
                end
                write(dat, "$ith_frame $i $j $ϕ $ψ $sum_angles\n")
            end
            
        end
    end; println(""); println(" ~~ Done!"); return nothing

end

function hemicellulose_dihedrals2(pdbfile::String, trajectory::String, segment::String, first_xylan_residue::Int64, last_xylan_residue::Int64;
    ffirst=1, flast=nothing, fstep=1)
    ##
    println(" ~~ Loading the PDB and TRAJECTORY information..."); println("")
    simulation = MolSimToolkit.Simulation(pdbfile, trajectory; first=ffirst, last=flast, step=fstep)

    println(" ~~ Calculating the shortest distances between the monitored atoms and the reference atoms:")
    segment_atoms = PDBTools.readPDB(pdbfile, only = atom -> (
        atom.segname == segment)
    ); segment_resnums = resnum.(segment_atoms)
    
    println(" ~~ Calculating the dihedrals associated to each backbone residue pair:")

    open("./dihedral2_$(segment).dat", "w") do dat
        ## ith frame, monitored residue, monitored atom name, monitored segment name, reference atom name, reference segment name, minimum distance
        for frame in simulation
            ith_frame = simulation.frame_index; println("    - frame $ith_frame")
            coor = MolSimToolkit.positions(frame)
            ## wrinting the output dat file
            resid_token = first_xylan_residue
            while resid_token < last_xylan_residue
                # setting the dummies index variables
                i = resid_token; j = resid_token + 1; resid_token = resid_token + 1
                println("      ... for the residue pair $i and $j")
                # getting the phi dihedral angle -- O5'-C1'-O4-C4
                da1_atom1 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == i && atom.name == "O5"))[1]
                da1_atom2 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == i && atom.name == "C1"))[1]
                da1_atom3 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "O4"))[1]
                da1_atom4 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "C4"))[1]
                ϕ = dihedral(@SVector([ coor[da1_atom1][1], coor[da1_atom1][2], coor[da1_atom1][3] ]),
                             @SVector([ coor[da1_atom2][1], coor[da1_atom2][2], coor[da1_atom2][3] ]),
                             @SVector([ coor[da1_atom3][1], coor[da1_atom3][2], coor[da1_atom3][3] ]),
                             @SVector([ coor[da1_atom4][1], coor[da1_atom4][2], coor[da1_atom4][3] ]))
                # getting the psi dihedral angle -- C1'-O4-C4-C5
                da2_atom1 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == i && atom.name == "C1"))[1]
                da2_atom2 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "O4"))[1]
                da2_atom3 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "C4"))[1]
                da2_atom4 = PDBTools.index.(select(segment_atoms, by = atom -> atom.segname == segment && atom.resnum == j && atom.name == "C5"))[1] 
                ψ = dihedral(@SVector([ coor[da2_atom1][1], coor[da2_atom1][2], coor[da2_atom1][3] ]),
                             @SVector([ coor[da2_atom2][1], coor[da2_atom2][2], coor[da2_atom2][3] ]),
                             @SVector([ coor[da2_atom3][1], coor[da2_atom3][2], coor[da2_atom3][3] ]),
                             @SVector([ coor[da2_atom4][1], coor[da2_atom4][2], coor[da2_atom4][3] ]))
                # sum of the angles phi + psi
                sum_angles = ϕ + ψ
                if sum_angles > 360 && sum_angles < 720
                    sum_angles = sum_angles - 360
                elseif sum_angles > 720
                    sum_angles = sum_angles - 720
                end
                write(dat, "$ith_frame $i $j $ϕ $ψ $sum_angles\n")
            end
            
        end
    end; println(""); println(" ~~ Done!"); return nothing

end