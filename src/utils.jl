"""
    select_atoms(atom::Vector{PDBTools.Atom}, name::String, resname::String, resnum::Int64, segname::String)

    Select the atoms that match the given name, residue name, residue number, and segment name.
"""
function select_atoms(atom::Vector{PDBTools.Atom}, name::String, resname::String, resnum::Int64, segname::String)
    if (atom.segname == segname) && (atom.resnum == resnum) && (atom.resname == resname) && (atom.name == name)
        return true
    else
        return false
    end
end

"""
    simulation_steps(t_simulation::Int64; unit="ns", timestep=2.0, output=1000)

    This function calculates the number of steps you need to set and the number of frames that you will get in the MD simulation.

# Arguments
- `t_simulation`: The total simulation time in the specified unit.
- `unit`: The unit of the simulation time. It can be "ns" (nan
- `timestep`: The time step of the simulation in the specified unit.
- `output`: The frequency of frame storage during the MD simulation.
"""
function simulation_steps(t_simulation::Int64; unit="ns", timestep=2.0, output=1000)
    
    timestep = timestep * 1e-3                                                                  # converting to ps (will be easier to work)

    if unit == "ps"
        total_time = t_simulation
    elseif unit == "ns"
        total_time = t_simulation * 1e3
    elseif unit == "μs"
        total_time = t_simulation * 1e6
    elseif unit == "ms"
        total_time = t_simulation * 1e9
    else
        throw(ArgumentError("Unit $unit not found (try to use `ps`, `ns`, `μs` or `ms`)."))
    end

    steps  = round(Int64, (total_time / timestep))
    frames = round(Int64, steps / output)

    println(
        raw"""Simulation:
        You should run""" * " $steps steps to complete $t_simulation $unit" * raw""".
        In the end, you will have""" * " $frames frames in the output, with each frame covering" * " $(t_simulation / Float64(frames)) $unit " * raw"""
        of your simulation.""")

    return steps, frames

end

function center_of_mass(pdbname::String, trajectory::String; selection="all", first=1, last=nothing, step=1)
    
    return center_of_mass(MolSimToolkit.Simulation(
                                    pdbname,
                                    trajectory;
                                    first=first, last=last, step=step
                                ), selection=selection)
end

function center_of_mass(simulation::MolSimToolkit.Simulation; selection="all")
    
    com = Vector{MolSimToolkit.Point3D}()

    idxs = PDBTools.selindex(MolSimToolkit.atoms(simulation), selection)

    for frame in simulation
        coor = MolSimToolkit.positions(frame)
        push!(com, MolSimToolkit.center_of_mass(idxs, simulation, coor))
    end

    return com
end

function rmsd(pdbname::String, trajectory::String; selection="all", mass=nothing, reference_frame=nothing, show_progress=false, first=1, last=nothing, step=1)
    
    return rmsd(MolSimToolkit.Simulation(
                        pdbname,
                        trajectory;
                        first=first, last=last, step=step
                    ), selection=selection, mass=mass, reference_frame=reference_frame, show_progress=show_progress)
end

function rmsd(simulation::MolSimToolkit.Simulation; selection="all", mass=nothing, reference_frame=nothing, show_progress=false)
    
    return MolSimToolkit.rmsd(
                    simulation,
                    PDBTools.selindex(MolSimToolkit.atoms(simulation), selection);
                    mass=mass, reference_frame=reference_frame, show_progress=show_progress
                )
end

"""
    writepdb_trajectory(pdbname::String, trajectory::String; selection="all", first=1, last=nothing, step=1)

Write a trajectory in a temporary PDB file. This is very useful when you want to visualize the edited trajectory in VMD.
The MDLovoFit program uses this function to write the trajectory in a PDB file before the analysis.
"""
function writepdb_trajectory(
                        pdbname::String,
                        trajectory::String;
                        selection="all",
                        first=1, last=nothing, step=1
                    )

    atoms = PDBTools.readPDB(pdbname, selection)
    simulation = MolSimToolkit.Simulation(atoms, trajectory, first=first, last=last, step=step)
    
    tempPDB = tempname() * ".pdb"
    
    io = Base.open(tempPDB, "w")
    for frame in simulation
        write_frame!(io, atoms, frame)
    end
    Base.close(io)

    return tempPDB
end

"""
    readcoords(pdbname::String, trjname::String; selection="protein and name CA")

Read the coordinates of a trajectory in a PDB file. Returns a vector of vectors of coordinates of selected atoms.
"""
function readcoords(pdbname::String, trjname::String; selection="protein and name CA")
    return readcoords(
                    PDBTools.readPDB(pdbname, selection),
                    trjname
                )
end

"""
    readcoords(atoms::AbstractVector{<:PDBTools.Atom}, pdbname::String)
"""
function readcoords(atoms::AbstractVector{<:PDBTools.Atom}, pdbname::String)

    frames = Vector{Vector{SVector{3, Float64}}}()
    frame = Vector{SVector{3, Float64}}(undef, length(atoms))

    file = Base.open(pdbname, "r")

    idx = 1
    for line in eachline(file)
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            atom = PDBTools.read_atom_pdb(line)
            frame[idx] = @SVector(Float64[atom.x, atom.y, atom.z])
            idx += 1
        end        
        if startswith(line, "END")
            push!(frames, deepcopy(frame))
            idx = 1
        end
    end
    Base.close(file)

    return frames
end

"""
    write_frame!(trajectory_pdb_file, atoms, frame

Write a frame of a trajectory in the temporary PDB trajectory file.
"""
function write_frame!(pdbfile::IO, atoms::AbstractVector{<:PDBTools.Atom}, frame::MolSimToolkit.Chemfiles.Frame) 
    
    idx_map = Dict(atom.index => idx for (idx, atom) in enumerate(atoms))
    coor = MolSimToolkit.positions(frame)

    for (i,xyz) in enumerate(coor)
        
        idx = get(idx_map, i, nothing)

        if !isnothing(idx)
           atoms[idx].x = xyz[1]
           atoms[idx].y = xyz[2]
           atoms[idx].z = xyz[3]
           Base.write(pdbfile, PDBTools.write_pdb_atom(atoms[idx]) * "\n")
        end

    end

    Base.write(pdbfile, "ENDMDL\n")

    return nothing
end

function pdb2trajectory(pdbname::String; trajectory=nothing, vmd="vmd", format="dcd", DebugVMD=false)
    
    trajectory = isnothing(trajectory) ? tempname() * "." * format : trajectory

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")
    
    Base.write(vmdinput,"""
        mol new     $pdbname
        
        animate goto 0
        animate write $format $trajectory
        """
        )

    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $tcl`, String), "\n")

    if DebugVMD
        return vmdoutput, tcl
    else
        return trajectory
    end
end

function frame_coordinates(
                    pdbname::String,
                    trajectory::String;
                    selection="protein and name CA",
                )                
    
    atoms = PDBTools.readPDB(pdbname, selection)

    simulation = MolSimToolkit.Simulation(atoms, trajectory)

    ##coor = zeros(SVector{length(atoms), MolSimToolkit.Point3D{Float64}}, length(simulation.frame_range))
    coor = []

    for (i,frame) in enumerate(simulation)
        push!(coor, MolSimToolkit.positions(frame))
    end

    return coor
    
end