function diffusion(t, msd; dim=3, timestep=0.2)
    fit = EasyFit.fitlinear(timestep .* t, msd)
    D = fit.a * inv(2*dim)
    println("""
    Diffusion coefficient: $D Å² / ns
    
    Based on the fit:
           y = $(fit.a) * x + $(fit.b);     R² = $(fit.R)
    """)
    return D * 1.0E+9 * (1.0E-8)^2 # to cm²/s
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
    τ = timelag(tmax, step)
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

function timelag(nframes::Int64, Δf::Int64)
    return collect(
        0:Δf:(nframes-Δf)
    )
end

