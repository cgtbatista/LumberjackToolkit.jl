function closest2fibril(
    pdbname::String, trjname::String, segcode::String; first=1, last=nothing, step=1,
    selection="resname BGLC", without_hydrogens=true, cutoff=5.0
)   
    simulation = MolSimToolkit.Simulation(pdbname, trjname, first=first, step=step, last=last)
    return closest2fibril(
        simulation,
        segcode,
        selection = selection,
        without_hydrogens = without_hydrogens,
        cutoff = cutoff
    )
    
end

function closest2fibril(simulation::MolSimToolkit.Simulation, segcode::String; selection="resname BGLC", without_hydrogens=true, cutoff=5.0)
    reference = PDBTools.select(simulation.atoms, selection) 
    monitored = PDBTools.select(simulation.atoms, at -> at.segname == segcode)
    if without_hydrogens
        reference = PDBTools.select(reference, at -> PDBTools.element(at) != "H")
        monitored = PDBTools.select(monitored, at -> PDBTools.element(at) != "H")
    end
    imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)
    monitored_resnums = PDBTools.resnum.(monitored)
    molindexes(iatom) = findfirst(x -> x == monitored_resnums[iatom], unique(monitored_resnums))
    iresid, distances = Vector{Vector{Int64}}(undef, length(simulation)), Vector{Vector{Float64}}(undef, length(simulation))
    for (iframe, frame) in enumerate(simulation)
        xyz, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))
        crosspairs = MolSimToolkit.CrossPairs(
            xpositions = [ SVector(xyz[i]) for i in imonitored ], # solvent - the monitored atoms around the reference (e.g. xylan and mannans)
            ypositions = [ SVector(xyz[j]) for j in jreference ], # solute  - the reference atoms (e.g. cellulsoe microfibril)
            xmol_indices = molindexes,
            cutoff = cutoff,
            unitcell = uc
        )
        mindistances = MolSimToolkit.minimum_distances!(crosspairs)
        iresid[iframe] = [ monitored_resnums[md.i] for md in mindistances ]
        distances[iframe] = [ md.d for md in mindistances ]
    end
    return iresid, distances
end

function closest2fibril(files::Vector{String}, first::Int64, last::Int64)
    distances = Vector{Matrix{Float64}}(undef, length(files))
    for (iseg, file) in enumerate(files)
        if filesize(file) == 0
            continue
        end
        distances[iseg] = closest2fibril(file, first, last)
    end
    return distances
end

function closest2fibril(file::String, first::Int64, last::Int64)
    data = readdlm(file)
    nframes = minimum(data[:,1]):maximum(data[:,1])
    distance = zeros(Float64, length(first:last), length(nframes))
    for iframe in nframes
        idx = findall((data[:,1] .== iframe) .& (first .<= data[:,2] .<= last))
        distance[:, iframe] .= data[idx,7]
    end
    return distance
end

function closest2fibril(
    data1::Vector{Matrix{Float64}}, data2::Vector{Matrix{Float64}};
    colors=["coral2", "skyblue2"], xlim=(1,96), ylim=(0,40), labels=["close", "far"]
)
    Plots.gr(size=(1000,800), dpi=900, fmt=:png)
    Plots.plot(
        title="", legend=:topleft, xlabel="# backbone residue", ylabel="Minimum distances (Å)",                                                     ## Basics
        framestyle=:box, grid=true, minorgrid=true, minorticks=5, thick_direction=:out, xlim=xlim, ylim=ylim,                                       ## Axis
        fontfamily=:arial, titlefontsize=18, guidefontsize=16, tickfontsize=16, labelfontsize=18, legendfontsize=16, guidefonthalign=:center,       ## Font
        left_margin=5Plots.Measures.mm, right_margin=10Plots.Measures.mm, top_margin=10Plots.Measures.mm, bottom_margin=1Plots.Measures.mm          ## Margins
    )
    data1 = hcat(data1...)
    data2 = hcat(data2...)
    average1, average2 = mean(data1, dims=2), mean(data2, dims=2)
    extrema1, extrema2 = extrema(data1, dims=2), extrema(data2, dims=2)
    lower1, upper1 = [ ext[1] for ext in extrema1 ], [ ext[2] for ext in extrema1 ]
    lower2, upper2 = [ ext[1] for ext in extrema2 ], [ ext[2] for ext in extrema2 ]
    Plots.plot!(
        average1,
        ribbon=(average1 .- lower1, upper1 .- average1),
        fillalpha=0.2,
        color=colors[1],
        linewidth=5,
        label=labels[1]
    )
    Plots.plot!(lower1, linestyle=:dashdot, color=colors[1], labels = "")
    Plots.plot!(upper1, linestyle=:dashdot, color=colors[1], labels = "")
    Plots.plot!(
        average2,
        ribbon=(average2 .- lower2, upper2 .- average2),
        fillalpha=0.2,
        color=colors[2],
        linewidth=5,
        label=labels[2]
    )
    Plots.plot!(lower2, linestyle=:dashdot, color=colors[2], labels = "")
    Plots.plot!(upper2, linestyle=:dashdot, color=colors[2], labels = "")
    println("""
    Closest distances to the fibril:

    $(labels[1])
    The average distance between the monitored atoms and the fibril is $(round(mean(average1), sigdigits=2)) ± $(round(std(average1), sigdigits=2)) Å
    The average minimum distance is $(round(mean(lower1[1]), sigdigits=2)) Å, with the lowest distance being $(round(minimum(lower1[1]), sigdigits=2)) Å
    The average maximum distance is $(round(upper1[1], sigdigits=2)) Å, with the highest distance being $(round(maximum(upper1[1]), sigdigits=2)) Å
    
    $(labels[2])
    The average distance between the monitored atoms and the fibril is $(round(mean(average2), sigdigits=2)) ± $(round(std(average2), sigdigits=2)) Å.
    The average minimum distance is $(round(mean(lower2[1]), sigdigits=2)) Å, with the lowest distance being $(round(minimum(lower2[1]), sigdigits=2)) Å
    The average maximum distance is $(round(upper2[1], sigdigits=2)) Å, with the highest distance being $(round(maximum(upper2[1]), sigdigits=2)) Å
    """)
    return Plots.current()
end

function catalytic_distances(
    pdbname::String, trjname::String; selection0="not water", selection1="water",
    first=1, last=nothing, step=1
)
    simulation = MolSimToolkit.Simulation(
        pdbname, trjname; first=first, last=last, step=step
    )
    idx0 = PDBTools.index(PDBTools.select(simulation.atoms, selection0))
    idx1 = PDBTools.index(PDBTools.select(simulation.atoms, selection1))
    return MolSimToolkit.distances(simulation, idx0, idx1)
end

function catalytic_distances(
    pdbname::String, trjname::String, idx0::AbstractVector{<:Integer}, idx1::AbstractVector{<:Integer};
    first=1, last=nothing, step=1
)
    simulation = MolSimToolkit.Simulation(
        pdbname, trjname; first=first, last=last, step=step
    )
    return MolSimToolkit.distances(simulation, idx0, idx1)
end

function catalytic_distances(dist1::Vector{Float64}, dist2::Vector{Float64})
    Plots.gr(size=(1000,800), dpi=900, fmt=:png)
    Plots.scatter(
        dist1, dist2,
        xlabel="distância O-Asp140 (Å)", ylabel="distância Asp140-Asp94 (Å)", labels="por frame",
        title="", fontfamily=:arial, alfa=0.5, color = "#0173B2",
        ## Axis configs
        framestyle=:box,
        grid=true,
        minorgrid=true,
        minorticks=5,
        thick_direction=:out,
        ## Font configs
        titlefontsize=18,
        guidefontsize=16,
        tickfontsize=16,
        labelfontsize=18,
        legendfontsize=16,
        guidefonthalign=:center,
        ## Margins
        left_margin=5Plots.Measures.mm,
        right_margin=10Plots.Measures.mm,
        top_margin=10Plots.Measures.mm,
        bottom_margin=1Plots.Measures.mm
    )
    avg_d1, avg_d2 = mean(dist1), mean(dist2)
    std_d1, std_d2 = std(dist1), std(dist2)
    println("Average distance O-Asp140: $avg_d1 ± $std_d1 Å")
    println("Average distance Asp140-Asp94: $avg_d2 ± $std_d2 Å")
    Plots.scatter!(
        [avg_d1], [avg_d2],
        label="média", color=:black, markersize=30, marker=:star
    )
    return Plots.current()
end ### colocar 