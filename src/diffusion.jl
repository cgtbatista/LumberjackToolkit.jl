function testfiles()
    return "/home/user/Documents/Softwoods/full/sys.pdb", "/home/user/Documents/Softwoods/full/full.dcd"
end

# function MSD(trajectories)
#     n_particles = length(trajectories)
#     n_frames = length(trajectories[1])
    
#     msd = zeros(n_frames)  # Array para armazenar o MSD para cada passo de tempo
    
#     for t in 1:n_frames
#         total_displacement = 0.0
#         for particle in trajectories
#             for t0 in 1:(n_frames - t)
#                 displacement = particle[t0 + t] - particle[t0]
#                 total_displacement += norm(displacement)^2
#             end
#         end
#         msd[t] = total_displacement / (n_particles * (n_frames - t))
#     end
    
#     return msd
# end

function coef_diffusion(pdbname::String, trajname::String; selection="water", d=3, first=1, last=nothing, step=1)
    
    atoms = PDBTools.index.(PDBTools.readPDB(pdbname, selection))
    simulation = MolSimToolkit.Simulation(pdbname, trajname; first=first, last=last, step=step)

    msd = MSD(simulation, atoms)

    fit = EasyFit.fitlinear(
                        collect(simulation.frame_range),
                        msd
                    )
    D = fit.a * inv(2d)
    println("""
    Diffusion coefficient: $D Å² / ns
    
    Based on the fit:
           y = $(fit.a) * x + $(fit.b);     R² = $(fit.R)
    """)
    return D
end

function MSD(simulation::MolSimToolkit.Simulation, indexes::Vector{Int64})
    
    natoms = length(indexes)
    
    msd = zeros(
                length(simulation.frame_range)
            )

    ref = MolSimToolkit.positions(simulation.frame)[indexes]

    for (i, frame) in enumerate(simulation)
        
        pos = MolSimToolkit.positions(frame)[indexes]
        
        Δr = 0.
        for (r0, r) in zip(ref, pos)
            Δr += norm(r - r0)^2
        end
        msd[i] = inv(natoms) * Δr
    end

   return msd
end
