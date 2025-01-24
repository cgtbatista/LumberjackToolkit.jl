"""
    
"""
function simsteps(t::Float64, Δt::Float64; unit="ns", timestep=2.0)
    
    conversion = Dict("ps" => 1, "ns" => 1e3, "μs" => 1e6, "ms" => 1e9)
    
    if !(unit in keys(conversion))
        throw(ArgumentError("Unit $unit not found (try to use `ps`, `ns`, `μs` or `ms`)."))
    end
    
    timestep = timestep * 1e-3

    t_ps, Δt_ps = t * conversion[unit], Δt * conversion[unit]
    
    freq = round(Int64, Δt_ps / timestep)
    
    nsteps, nframes = round(Int64, t_ps / timestep), round(Int64, t_ps / (freq * timestep))
    
    println("""
    Simulation: 
    You should run $nsteps steps to complete $t $unit.
    In the end, you will have $nframes frames that will be saved after $freq simulation steps, covering the simulation time by $Δt $unit.
    """)

    return nsteps, freq
end

function realtime(nsteps::Int64, freq::Int64; n=1, unit="ns", timestep=2.0)

    conversion = Dict("ps" => 1, "ns" => 1e-3, "μs" => 1e-6, "ms" => 1e-9)
    
    if !(unit in keys(conversion))
        throw(ArgumentError("Unit $unit not found (try to use `ps`, `ns`, `μs` or `ms`)."))
    end

    timestep = timestep * 1e-3
    
    t , Δt = nsteps * timestep * conversion[unit], n * freq * timestep * conversion[unit]
    nframes = round(Int64, nsteps / freq)

    println("""
    Real time:
    The simulation will take $t $unit to complete $nsteps steps.
    If you pick the $nframes frames by $n frame(s), you will cover $Δt $unit every step with only $(round(Int64, nframes / n)) frames.
    """)
end