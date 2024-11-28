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