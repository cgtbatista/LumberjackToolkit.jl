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

function center_of_mass(
    pdbname::String, trajectory::String;
    selection="all",
    first=1, last=nothing, step=1
)
    simulation = MolSimToolkit.Simulation(
        pdbname, trajectory; first=first, last=last, step=step
    )
    return center_of_mass(simulation, selection=selection)
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
    selection="all", filename=nothing,
    first=1, last=nothing, step=1
)
    atoms = PDBTools.readPDB(pdbname, selection)
    simulation = MolSimToolkit.Simulation(
            atoms, trajectory,
            first=first, last=last, step=step
        )
    filename = isnothing(filename) ? tempname() * ".pdb" : filename
    
    if filename[end-3:end] != ".pdb"
        throw(ArgumentError("The filename must have the .pdb extension."))
    end
    
    Base.open(filename, "w") do io
        for frame in simulation
            write_frame!(io, atoms, frame)
        end
    end
    return filename
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
           Base.write(pdbfile, PDBTools.write_atom(atoms[idx]) * "\n")
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

function molindexes(pdbname::String; selection="water")    
    return molindexes(
                    PDBTools.readPDB(pdbname, selection)
                )
end

function molindexes(atoms::Vector{<:PDBTools.Atom})
    
    indexes = Dict{Tuple{Int64, String}, Vector{Int64}}()

    for (index, atom) in enumerate(atoms)
        key = (atom.resnum, atom.segname)
        
        if !haskey(indexes, key)
            indexes[key] = [index]
        else
            push!(indexes[key], index)
        end
    end

    return collect(values(indexes))
end

"""
    readcoords(atoms::AbstractVector{<:PDBTools.Atom}, pdbname::String)

Read the coordinates of a trajectory in a PDB file. Returns a vector of vectors of coordinates of selected atoms.
"""
function readcoords(atoms::AbstractVector{<:PDBTools.Atom}, pdbname::String)
    trajectory, iframe = Vector{Vector{SVector{3, Float64}}}(), Vector{SVector{3, Float64}}(undef, length(atoms))
    idx = 1
    open(pdbname, "r") do file
        for line in eachline(file)
            if startswith(line, "ATOM") || startswith(line, "HETATM")
                atom = PDBTools.read_atom_pdb(line)
                iframe[idx] = @SVector(Float64[atom.x, atom.y, atom.z])
                idx += 1
            end        
            if startswith(line, "END")
                push!(trajectory, deepcopy(iframe))
                idx = 1
            end
        end
    end
    return trajectory
end

"""
    readcoords(pdbname::String, trjname::String; selection="protein and name CA")

"""
# function readcoords2(
#     pdbname::String,
#     trjname::String;
#     selection="water and name OH2", averaging=false, bymol=false,
#     first=1, last=nothing, step=1
# )
#     atoms = PDBTools.readPDB(pdbname, selection)
#     simulation = MolSimToolkit.Simulation(atoms, trjname, first=first, last=last, step=step)

#     if averaging
#         iavg = molindexes(atoms)
#         N = length(iavg)
#     else
#         N = length(atoms)
#     end
    
#     coordinates = if bymol
#             [Vector{SVector{3, Float64}}(undef, length(simulation.frame_range)) for _ in 1:N]
#         else
#             Vector{Vector{SVector{3, Float64}}}(undef, length(simulation.frame_range))
#     end
#     for (i, frame) in enumerate(simulation)        
#         xyz = if averaging
#                 avgpositions(frame, iavg)
#             else
#                 MolSimToolkit.positions(frame)
#         end
#         if bymol
#             for mol in 1:N
#                 coordinates[mol][i] = SVector{3, Float64}(xyz[mol]...)
#             end
#         else
#             coordinates[i] = [SVector{3, Float64}(c...) for c in xyz]
#         end
#     end
#     return coordinates
# end


# function readcoords(
#                 pdbname::String,
#                 trjname::String;
#                 selection="water and name OH2", averaging=false, bymol=false,
#                 first=1, last=nothing, step=1
#             )

#     atoms = PDBTools.readPDB(pdbname, selection)
#     simulation = MolSimToolkit.Simulation(atoms, trjname, first=first, last=last, step=step)

#     iavg = if averaging
#         molindexes(atoms)
#     end

#     coords = Vector{Vector{MolSimToolkit.Point3D{Float64}}}(undef, length(simulation.frame_range))

#     for (i, frame) in enumerate(simulation)
        
#         coords[i] = if averaging
#             avgpositions(frame, iavg)
#         else
#             MolSimToolkit.positions(frame)
#         end

#     end

#     return coords
# end


function writecoords(
    pdbname::String,
    trjname::String;
    selection="water and name OH2", filename=nothing, isaverage=false,
    first=1, last=nothing, step=1
)
    atoms = PDBTools.readPDB(pdbname, selection)
    if isaverage
        idx = molindexes(atoms)
        N = length(idx)
    else
        N = length(atoms)
    end

    simulation = MolSimToolkit.Simulation(
        atoms, trjname,
        first=first, last=last, step=step
    )
    F = length(simulation.frame_range)

    filename = isnothing(filename) ? tempname() * ".bin" : filename

    Base.open(filename, "w") do io
        Base.write(io, Int32(N), Int32(F))
        bincoords = Vector{Float32}(undef, 3 * N)
        for frame in simulation
            coords = if isaverage
                    avgpositions(frame, idx)
                else
                    MolSimToolkit.positions(frame)
            end
            @inbounds for i in 1:N
                bincoords[3 * (i - 1) + 1] = Float32(coords[i][1])
                bincoords[3 * (i - 1) + 2] = Float32(coords[i][2])
                bincoords[3 * (i - 1) + 3] = Float32(coords[i][3])
            end
            # bincoords = reduce(vcat, [Float32.(coord) for coord in coords])   ## não é muito eficiente, pois reatribui a variável a cada iteração
            Base.write(io, bincoords)
        end
    end
    return filename
end

function avgpositions(frame::MolSimToolkit.Chemfiles.Frame, iavg::Vector{Vector{Int64}})
    coords = MolSimToolkit.positions(frame)
    return MolSimToolkit.Point3D{Float64}[ mean(coords[ijk]) for ijk in iavg ]
end

function unwrapping(
    psfname::String,
    pdbname::String,
    trajectory::String;
    new_trajectory=nothing, vmd="vmd", DebugVMD=false
)
    
    new_trajectory = isnothing(new_trajectory) ? tempname() * ".dcd" : new_trajectory

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")    
    Base.write(vmdinput,"""
        package require pbctools

        mol new     $psfname
        mol addfile $pdbname
        mol addfile $trajectory waitfor all
        
        animate goto 0
        pbc unwrap -all

        animate write dcd $new_trajectory

        exit
        """
        )
    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $tcl`, String), "\n")

    if DebugVMD
        return vmdoutput, tcl
    else
        return new_trajectory
    end
end

function molframes(
    pdbname::String,
    trjname::String;
    selection="water and name OH2",
    first=1, last=nothing, step=1
)
    atoms = PDBTools.readPDB(pdbname, selection)
    sim = MolSimToolkit.Simulation(
        atoms,
        trjname,
        first=first, last=last, step=step
    )    
    return molframes(sim)
end

function molframes(sim::MolSimToolkit.Simulation)
    natoms = length(sim.atoms)
    nframes = length(sim.frame_range)
    trajectory = [Vector{SVector{3, Float64}}(undef, nframes) for _ in 1:natoms]
    for (iframe, frame) in enumerate(sim)        
        coords =  MolSimToolkit.positions(frame)
        for atom in 1:natoms
            trajectory[atom][iframe] = SVector{3, Float64}(coords[atom]...)
        end
    end
    return trajectory
end

# function molframes(sim::MolSimToolkit.Simulation)
#     natoms, nframes = length(sim.atoms), length(sim.frame_range)
#     trajectory = [Vector{SVector{3, Float64}}(undef, nframes) for _ in 1:natoms]
#     MolSimToolkit.first_frame!(sim)
#     plast = deepcopy(
#             MolSimToolkit.positions(MolSimToolkit.current_frame(sim))
#         )
#     for iframe in 2:length(sim)
#         frame = MolSimToolkit.next_frame!(sim)
#         p =  MolSimToolkit.positions(frame)
#         uc = MolSimToolkit.unitcell(frame)
#         for iatom in 1:natoms
#             pat = MolSimToolkit.wrap(p[iatom], plast[iatom], uc)
#             trajectory[iatom][iframe] = SVector{3, Float64}(pat...)
#         end
#         plast .= p
#     end
#     return trajectory
# end

### VMD

function execVMD(input::String; vmd::String="vmd")
    
    if !hasVMD(vmd)
        error("VMD executable not found. Please, check if VMD is installed and in your PATH.")
    end

    if !isfile(input)
        error("File not found: $input. Be sure to provide a valid file.")
    end

    return Base.read(`$vmd -dispdev text -e $input`, String)
end

function hasVMD(vmd::String)
    return Sys.which(vmd) !== nothing
end