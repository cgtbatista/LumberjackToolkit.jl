struct CarbohydrateDihedrals
    frame::Vector{Int64}
    i::Vector{Int64}
    j::Vector{Int64}
    ϕ::Vector{Float64}
    ψ::Vector{Float64}
    sum_dihedrals::Vector{Float64}
end


function dihedral_atoms(dihedral::String, resnum::Int64; sep=" ")
    
    names, resnums = String[], Int64[]
    atoms = String.(split(dihedral, sep))

    for atom in atoms
        if occursin("'", atom)
            push!(names, replace(atom, "'" => ""))
            push!(resnums, resnum)
        else
            push!(names, atom)
            push!(resnums, resnum+1)
        end
    end

    return names, resnums
end


function dihedral_indexes(atoms::Vector{PDBTools.Atom}, segname::String, resnum::Int64; dihedral="O5'-C1'-O4-C4", sep="-")
    
    indexes = Int64[]
    names, resnums = dihedral_atoms(dihedral, resnum, sep=sep)

    for (name, resnum) in zip(names, resnums)
        idx = PDBTools.index.(
                            PDBTools.select(atoms, atom -> (atom.segname == segname) && (atom.resnum == resnum) && (atom.name == name))
                        )[1]
        push!(indexes, idx)
    end

    return indexes    
end


"""
    dihedral(atom1::StaticArrays.SVector, atom2::StaticArrays.SVector, atom3::StaticArrays.SVector, atom4::StaticArrays.SVector)

Calculate the dihedral angle between four atoms. The dihedral angle is the angle between the planes defined by the atoms (atom1, atom2, atom3)
and (atom2, atom3, atom4). The dihedral angle is calculated using the atan method.
"""
function dihedral(atom1::SVector, atom2::SVector, atom3::SVector, atom4::SVector)
    
    v1, v2, v3 = atom2 - atom1, atom3 - atom2, atom4 - atom3    
    w1, w2 = cross(v1, v2), cross(v2, v3)
    
    dihedral = atand(
                dot(cross(w1, w2), v2/norm(v2)),
                dot(w1, w2)
            )
    
    return mod(dihedral, 360.0)
end


"""
   dihedrals(pdbname::String, trajectory::String, segname::String; ...)

Calculate the dihedral angles and their sum based on the carbohydrate chain. The `pdbname` carries the atom information, and the `trajectory` holds the coordinates
for each atoms in the simulation. The carbohydrate chain is defined by the `segname` parameter, that is the same segname present on PDB.

Wohlert uses a lot the ϕ(O5'-C1'-O4-C4) and ψ(C1'-O4-C4-C3) definition, while Carol uses ϕ(O5'-C1'-O4-C4) and ψ(C1'-O4-C4-C5)
"""
function dihedrals(
    pdbname::String, trajectory::String, segname::String;
                first_resid=3, last_resid=98, dihedral1="O5'-C1'-O4-C4", dihedral2="C1'-O4-C4-C5", sep="-",
                first=1, last=nothing, step=1, 
                outfile=nothing
            )

    println("")
    println(" ~~ Loading the simulation data...")
    simulation = MolSimToolkit.Simulation(
                            pdbname,
                            trajectory;
                            first=first, last=last, step=step
                        )
    println("")

    println(" ~~ Picking all of the $segname atoms:")
    atoms = PDBTools.readPDB(pdbname, only = (atom -> atom.segname == segname))
    println("")
    
    println(" ~~ Calculating the dihedrals associated to each residue pair:")
    println("")
    
    if !isnothing(outfile)
        out = Base.open(outfile, "w")
    else
        f = Int64[]
        index1 = Int64[]; index2 = Int64[]
        angle1 = Float64[]; angle2 = Float64[]; sum_dihedrals = Float64[]
    end

    for frame in simulation
        
        i = simulation.frame_index
        
        coor = MolSimToolkit.positions(frame)

        println("""

        ~~~~ Frame $i ~~~~

        frame           i      j          ϕ          ψ      ϕ + ψ
        ---------------------------------------------------------
        """)        

        resid = first_resid
        
        while resid < last_resid
            
            # getting the ϕ dihedral angle
            ϕ_atoms = dihedral_indexes(atoms, segname, resid, dihedral=dihedral1, sep=sep)
            ϕ = dihedral(
                    @SVector([ coor[ϕ_atoms[1]][1], coor[ϕ_atoms[1]][2], coor[ϕ_atoms[1]][3] ]),
                    @SVector([ coor[ϕ_atoms[2]][1], coor[ϕ_atoms[2]][2], coor[ϕ_atoms[2]][3] ]),
                    @SVector([ coor[ϕ_atoms[3]][1], coor[ϕ_atoms[3]][2], coor[ϕ_atoms[3]][3] ]),
                    @SVector([ coor[ϕ_atoms[4]][1], coor[ϕ_atoms[4]][2], coor[ϕ_atoms[4]][3] ])
                )
            
            # getting the ψ dihedral angle
            ψ_atoms = dihedral_indexes(atoms, segname, resid, dihedral=dihedral2, sep=sep)
            ψ = dihedral(
                    @SVector([ coor[ψ_atoms[1]][1], coor[ψ_atoms[1]][2], coor[ψ_atoms[1]][3] ]),
                    @SVector([ coor[ψ_atoms[2]][1], coor[ψ_atoms[2]][2], coor[ψ_atoms[2]][3] ]),
                    @SVector([ coor[ψ_atoms[3]][1], coor[ψ_atoms[3]][2], coor[ψ_atoms[3]][3] ]),
                    @SVector([ coor[ψ_atoms[4]][1], coor[ψ_atoms[4]][2], coor[ψ_atoms[4]][3] ])
                )
            
            ϕ = mod(ϕ, 360.0)
            ψ = mod(ψ, 360.0)

            info = (i, resid, resid+1, ϕ, ψ, ϕ+ψ)
            @printf("%-10d %6d %6d %10.2f %10.2f %10.2f\n", info...)

            if !isnothing(outfile)
                Base.write(out, "$i $resid $(resid+1) $ϕ $ψ $(ϕ+ψ)\n")
            else
                push!(f, i)
                push!(index1, resid); push!(index2, resid+1)
                push!(angle1, ϕ); push!(angle2, ψ); push!(sum_dihedrals, ϕ+ψ)
            end

            resid += 1

        end

        println("")
        
    end

    if !isnothing(outfile)
        Base.close(out)
    end

    println("")
    println(" ~~ Done!")

    Base.GC.gc()

    if !isnothing(outfile)
        return outfile
    else
        return CarbohydrateDihedrals(f, index1, index2, angle1, angle2, sum_dihedrals)
    end

end