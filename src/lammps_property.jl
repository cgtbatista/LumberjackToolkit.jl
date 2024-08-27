## (TODO) colocar uma função para ler seleções e só pegar as propriedades de interesse

#using MolSimToolkit
#using Chemfiles
#using StaticArrays

function get_property(pdbfile::String, lammps_trajectory::String; lammps_property="v_vonmises", first=1, last=nothing, step=1)

    simulation = MolSimToolkit.Simulation(pdbfile, lammps_trajectory; first=first, last=last, step=step)
    natoms = size(simulation.atoms)[1]

    lammpstrj_data = []

    for frame in simulation
        ith_frame = simulation.frame_index
        println(" ~~ frame $ith_frame")
        for j_atom in 1:natoms
            println("    $j_atom")
            atomdata = Chemfiles.Atom(frame, j_atom-1)
            if count(x -> x == lammps_property, Chemfiles.list_properties(atomdata)) == 1
                push!(lammpstrj_data, property(atomdata, lammps_property))
            end; atomdata = nothing
        end
        Base.GC.gc()
    end

    return convert(Vector{Float64}, lammpstrj_data)

end