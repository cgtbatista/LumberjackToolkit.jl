function testfiles()
    return "/home/user/Documents/Softwoods/dried/sys.pdb", "/home/user/Documents/Softwoods/dried/1000.dcd"
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

# function coef_diffusion(pdbname::String, trajname::String; selection="water", d=3, first=1, last=nothing, step=1)
    
#     atoms = PDBTools.index.(PDBTools.readPDB(pdbname, selection))
#     simulation = MolSimToolkit.Simulation(pdbname, trajname; first=first, last=last, step=step)

#     msd = MSD(simulation, atoms)

#     fit = EasyFit.fitlinear(
#                         collect(simulation.frame_range),
#                         msd
#                     )
#     D = fit.a * inv(2d)
#     println("""
#     Diffusion coefficient: $D Å² / ns
    
#     Based on the fit:
#            y = $(fit.a) * x + $(fit.b);     R² = $(fit.R)
#     """)
#     return D
# end

function coef_diffusion(t, msd; d=3)

    fit = EasyFit.fitlinear(
                        t,
                        msd
                    )
    D = fit.a * inv(2d)
    println("""
    Diffusion coefficient: $D Å² / ns
    
    Based on the fit:
           y = $(fit.a) * x + $(fit.b);     R² = $(fit.R)
    """)
    return D * 1.0E+9 * (1.0E-8)^2 # to cm²/s
end


function diffusion(pdbname::String, trjname::String)
    pdbname, trjname = testfiles()
    coef_diffusion(pdbname, trjname)
end

"""
    msd(trajectory::Vector{Vector{SVector{3, Float64}}}; step=1, dim=collect(1:3))

Calculate the mean square displacement (MSD) of a trajectory. The MSD is calculated using Green-Kubo relation: ```<r²(t)> = MSD(t) = 1/N * Σ <[r_i(t) - r_i(0)]²>```.
It is important to notice that trajectories must have been adjusted to approach:

        traj = [
                    [coord(t1), coord(t2), coord(t3), ..., coord(tmax)],    ## molecule 1
                    [coord(t1), coord(t2), coord(t3), ..., coord(tmax)],    ## molecule 2
                    [coord(t1), coord(t2), coord(t3), ..., coord(tmax)],    ## molecule 3
                    ...
                    [coord(t1), coord(t2), coord(t3), ..., coord(tmax)],    ## molecule N
            ]

### Arguments
- `trajectory::Vector{Vector{SVector{3, Float64}}}`: A vector of vectors of 3D coordinates.
- `step::Int=1`: The time step to calculate the MSD.
- `dim::Vector{Int}=collect(1:3)`: The dimensions to calculate the MSD.
"""
function msd(trajectory::Vector{Vector{SVector{3, Float64}}}; step=1, dims=collect(1:3))
    if step <= 0
        throw(ArgumentError("step must be greater than 0"))
    end

    N, tmax = length(trajectory), length(trajectory[1])
        
    τ = collect(
            0:step:(tmax-step)
        )
    MSD = zeros(Float64, length(τ)) 
    for (i, Δt) in enumerate(τ)
        Σ = 0.0
        for molecule in trajectory
            for t0 in 1:(tmax - Δt)
                Δr = molecule[t0 + Δt] - molecule[t0]
                Σ += sum(Δr[dim]^2 for dim in dims)
            end
        end
        MSD[i] = inv(N * (tmax - Δt)) * Σ
    end    
    return τ, MSD
end

# function msd(trajectory::Vector{Vector{SVector{3, Float64}}}; timestep=1, dims=collect(1:3))
    
#     if timestep <= 0
#         throw(ArgumentError("step must be greater than 0"))
#     end

#     N, tmax = length(trajectory[1]), length(trajectory)
        
#     τ = collect(
#             1:timestep:(tmax-timestep)
#         )

#     MSD = zeros(Float64, length(τ))
    
#     for (i, Δt) in enumerate(τ)
#         Σ = 0.0
#         println("$i / $(length(τ))")
#         for t0 in 1:(tmax - Δt)
#             r0, rt = trajectory[t0], trajectory[t0 + Δt]
#             for molecule in 1:N
#                 Δr = rt[molecule] - r0[molecule]
#                 Σ += norm(Δr[dims])^2
#             end
#         end

#         MSD[i] = inv(N * (tmax-Δt)) * Σ
#     end
    
#     return τ, MSD
# end

# function msd(trajectory::Vector{Vector{SVector{3, Float64}}}; timestep=1, dims=1:3)
#     T = length(trajectory)        
#     N = length(trajectory[1])     
#     τ = 1:timestep:(T-1)          

#     coords = zeros(T, N, 3)
#     for t in 1:T
#         for n in 1:N
#             coords[t, n, :] = trajectory[t][n][dims]
#         end
#     end

#     MSD = zeros(length(τ))
    
#     # Calcular MSD via FFT para cada molécula
#     @views for n in 1:N
#         # Extrair coordenadas da molécula n
#         x = coords[:, n, 1]
#         y = coords[:, n, 2]
#         z = coords[:, n, 3]
        
#         # Calcular autocorrelação via FFT
#         Fx = rfft(x)       # FFT real (mais eficiente)
#         Fy = rfft(y)
#         Fz = rfft(z)
        
#         S = abs2.(Fx) .+ abs2.(Fy) .+ abs2.(Fz)  # |FFT|² para cada componente
#         autocorr = irfft(S, length(x))           # FFT inversa
        
#         # Acumular MSD para todos os Δt
#         msd_n = 2 .* (autocorr[1] .- autocorr[1 .+ τ])
#         MSD .+= msd_n
#     end

#     # Média sobre as moléculas e ajuste de normalização
#     MSD ./= N

#     return τ, MSD
# end