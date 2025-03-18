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

function min_distances(
    pdbname::String, trjname::String; selection1="not water", selection2="water",
    first=1, last=nothing, step=1
)
    simulation = MolSimToolkit.Simulation(
        pdbname, trjname; first=first, last=last, step=step
    )
    idx1 = PDBTools.index.(PDBTools.select(simulation.atoms, selection1))
    idx2 = PDBTools.index.(PDBTools.select(simulation.atoms, selection2))
    mindistances = zeros(Float64, length(simulation))
    for (iframe, frame) in enumerate(simulation)
        p1, p2 = MolSimToolkit.positions(frame)[idx1], MolSimToolkit.positions(frame)[idx2]
        uc = MolSimToolkit.unitcell(frame)
        d12 = +Inf
        for x1 in p1, x2 in p2
            x2_wrapped = MolSimToolkit.wrap(x2, x1, uc)
            d12 = min(
                d12,
                norm(x2_wrapped .- x1)
            )
        end
        mindistances[iframe] = d12
    end
    return mindistances
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

function water_hbonding(
    pdbname::String, trjname::String; first=1, step=1, last=nothing,
    selection="not water", water="TIP3", OO=3.5, HO=2.5, HOO=30
)
    simulation = MolSimToolkit.Simulation(pdbname, trjname, first=first, step=step, last=last)
    return water_hbonding(
        simulation,
        selection = selection, water = water, OO = OO, HO = HO, HOO = HOO
    )
end

function water_hbonding(
    simulation::MolSimToolkit.Simulation;
    selection="not water",                          # selection of the reference atoms -- it's the solute! It can be the carbohydrates, the ions, etc.
    water="TIP3", OO=3.5, HO=2.5, HOO=30            # the cutoff distance to consider the water molecule around the reference atoms
)
    monitored = PDBTools.select(simulation.atoms, at -> at.resname == water && PDBTools.element(at) != "H")
    reference = PDBTools.select(
        PDBTools.select(simulation.atoms, selection),
        at -> PDBTools.element(at) != "H"
    )
    imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)
    M = Matrix{Bool}(undef, length(imonitored), length(simulation.frame_range))
    for (iframe, frame) in enumerate(simulation)
        xyz, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))
        mindist = MolSimToolkit.minimum_distances(
            xpositions = [ SVector(xyz[i]) for i in imonitored ],     # solvent - the monitored atoms around the reference (e.g. water)
            ypositions = [ SVector(xyz[j]) for j in jreference ],     # solute  - the reference atoms (e.g. polymeric matrix)
            xn_atoms_per_molecule = 1,
            cutoff = OO,
            unitcell = uc
        )
        for (iwater, md) in enumerate(mindist)
            M[iwater, iframe] = md.within_cutoff ? hbond_extended_geomcriteria(
                simulation.atoms, xyz, uc=uc, i=md.i, j=md.j, HO=HO, HOO=HOO
            ) : false
        end
    end
    return BitMatrix(M)
end

function water_hbonding_parallel(
    simulation::MolSimToolkit.Simulation;
    selection="not water",                          
    water="TIP3", OO=3.5, HO=2.5, HOO=30
)
    monitored = PDBTools.select(simulation.atoms, at -> at.resname == water && PDBTools.element(at) != "H")
    reference = PDBTools.select(
        PDBTools.select(simulation.atoms, selection),
        at -> PDBTools.element(at) != "H"
    )
    imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)
    M = Matrix{Bool}(undef, length(imonitored), length(simulation.frame_range))
    frames = [ frame for frame in simulation ]
    BLAS.set_num_threads(1)
    @threads for iframe in eachindex(frames)
        frame = frames[iframe]
        xyz, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))        
        mindist = MolSimToolkit.minimum_distances(
            xpositions = [ SVector(xyz[i]) for i in imonitored ],
            ypositions = [ SVector(xyz[j]) for j in jreference ],
            xn_atoms_per_molecule = 1,
            cutoff = OO,
            unitcell = uc
        )
        for (iwater, md) in enumerate(mindist)
            M[iwater, iframe] = md.within_cutoff ? hbond_extended_geomcriteria(
                simulation.atoms, xyz, uc=uc, i=md.i, j=md.j, HO=HO, HOO=HOO
            ) : false
        end
    end
    
    return BitMatrix(M)
end

function hbond_extended_geomcriteria(
    atoms::AbstractVector{<:PDBTools.Atom}, xyz::MolSimToolkit.FramePositions;
    uc=zeros(Float64, 3, 3), i=nothing, j=nothing, HO=2.5, HOO=30
)
    if isnothing(HO) || isnothing(HOO)
        return true
    end
    ## checking water as donor...
    OA = SVector(xyz[j])
    OD = SVector(xyz[i])
    hydrogens = PDBTools.index.(PDBTools.select(atoms, at -> at.resnum == atoms[i].resnum && PDBTools.element(at) == "H"))
    isdonor = false
    for idx in hydrogens
        HD = SVector(xyz[idx])
        if hbond_angle_checking(HD, OD, OA, cutoff=HOO, uc=uc) && hbond_bond2_checking(HD, OA, cutoff=HO, uc=uc)
            isdonor = true
            break
        end
    end
    ## checking water as acceptor...
    OA = SVector(xyz[i])
    OD = SVector(xyz[j])
    hydrogens = PDBTools.index.(PDBTools.select(atoms, at -> at.resnum == atoms[j].resnum && at.name == string("H", atoms[j].name)))
    if length(hydrogens) == 0
        isacepptor = false
    else
        idx = hydrogens[1]
        HD = SVector(xyz[idx])
        isacepptor = hbond_angle_checking(HD, OD, OA, cutoff=HOO, uc=uc) && hbond_bond2_checking(HD, OA, cutoff=HO, uc=uc)
    end
    if isdonor || isacepptor
        return true
    end
    return false
end

function hbond_angle_checking(
    HD::SVector{3, Float64}, OD::SVector{3, Float64}, OA::SVector{3, Float64};
    cutoff=30.0, uc=zeros(Float64, 3, 3)
)
    HD, OD = if uc == zeros(Float64, 3, 3)
            HD, OD
        else
            MolSimToolkit.wrap(HD, OA, uc), MolSimToolkit.wrap(OD, OA, uc)
    end
    v1 = HD .- OD
    v2 = OA .- OD
    ratio = clamp(
        dot(v1, v2) / (norm(v1) * norm(v2)),
        -1.0, 1.0    
    )
    return acosd(ratio) < cutoff
end

function hbond_bond2_checking(
    HD::SVector{3, Float64}, OA::SVector{3, Float64};
    cutoff=2.5, uc=zeros(Float64, 3, 3)
)
    HD = if uc == zeros(Float64, 3, 3)
            HD
        else
            MolSimToolkit.wrap(HD, OA, uc)
    end
    return norm(HD .- OA) < cutoff
end

function mapwater(
    pdbname::String, trjname::String;
    first=1, step=1, last=nothing,
    selection="not water", water_selection="water",
    without_hydrogens=true, cutoff=5.0
)
    simulation = MolSimToolkit.Simulation(pdbname, trjname, first=first, step=step, last=last)
    return mapwater(
        simulation,
        selection = selection, water_selection = water_selection,
        without_hydrogens = without_hydrogens,
        cutoff = cutoff
    )
end

function mapwater(
    simulation::MolSimToolkit.Simulation;
    selection="not water",                      # selection of the reference atoms -- it's the solute! It can be the carbohydrates, the ions, etc.
    water_selection="water",                    # selection of the water molecules, ie. the solvent of the PCW matrix
    without_hydrogens=true,                     # if the hydrogens should be removed from the analysis. Very useful to approach some experimental evidences.
    cutoff=5.0                                  # the cutoff distance to consider the water molecule around the reference atoms
)
    reference = PDBTools.select(simulation.atoms, selection) 
    monitored = PDBTools.select(simulation.atoms, water_selection)
    if without_hydrogens
        reference = PDBTools.select(reference, at -> PDBTools.element(at) != "H")
        monitored = PDBTools.select(monitored, at -> PDBTools.element(at) != "H")
        natoms = 1
    else
        natoms = 3
    end
    imonitored, jreference = PDBTools.index.(monitored), PDBTools.index.(reference)     
    M = Matrix{Bool}(undef, length(imonitored), length(simulation.frame_range))
    for (iframe, frame) in enumerate(simulation)
        xyz, uc = MolSimToolkit.positions(frame), diag(MolSimToolkit.unitcell(frame))
        mindistances = MolSimToolkit.minimum_distances(
            xpositions = [ SVector(xyz[i]) for i in imonitored ],     # solvent - the monitored atoms around the reference (e.g. water)
            ypositions = [ SVector(xyz[j]) for j in jreference ],     # solute  - the reference atoms (e.g. polymeric matrix)
            xn_atoms_per_molecule = natoms,
            cutoff = cutoff,
            unitcell = uc
        )
        for (iwater, md) in enumerate(mindistances)
           M[iwater, iframe] = md.within_cutoff
        end
    end
    return BitMatrix(M)
end

function mapwater(M1::BitMatrix, M2::BitMatrix, operation::Int64)
    M = M1 + M2
    if operation in Set([0,1,2])
        return ifelse.(M .== operation, true, false)
    else
        throw(ArgumentError("The operation must be 0 (complement), 1 (disjoint), or 2 (intersection)."))
    end
end

function mapwater(files::Vector{String}; threshold=nothing, nbins=40)
    t = Float64[]
    for file in files
        if filesize(file) == 0
            continue
        end
        data = isnothing(threshold) ? readdlm(file)[:,1] : filter(>(threshold), readdlm(file)[:,1])
        append!(t, data)
    end
    Plots.gr(size=(1000,800), dpi=900, fmt=:png)
    fitting = normalize(fit(Histogram, t, nbins=nbins), mode=:probability)
    Plots.plot(fitting, bins=nbins, label="", normalize=true)
    Plots.plot!(
        xlabel="residence time (ns)", ylabel="Normalized frequency", title="", fontfamily=:arial, alfa=1.0, color=:skyblue2,
        ## Axis configs
        xlims=extrema(t), ylims=(0.0, 1.0),
        framestyle=:box, grid=true, minorgrid=true, minorticks=5, thick_direction=:out,
        ## Font configs
        titlefontsize=20, guidefontsize=20, tickfontsize=18, labelfontsize=22, legendfontsize=16,
        guidefonthalign=:center,
        ## Margins
        left_margin=5Plots.Measures.mm,
        right_margin=10Plots.Measures.mm,
        top_margin=10Plots.Measures.mm,
        bottom_margin=1Plots.Measures.mm
    )
    println("""
    -------------------
    Time histogram data
    -------------------

    The average time is $(round(mean(t), sigdigits=2)) ⨦ $(round(std(t), sigdigits=2)) ns.
    The minimum time is $(round(minimum(t), sigdigits=2)) ns, and the maximum time is $(round(maximum(t), sigdigits=2)) ns.
    
    Filtering the values above 0.2 ns with average time of $(round(median(filter(>(0.2), t)), sigdigits=2)) ns.
    t > 1.0 ns: $(round(100 * length(filter(>(1.0), t)) / length(filter(>(0.2), t)), sigdigits=2)) of the total.
    t > 2.0 ns: $(round(100 * length(filter(>(2.0), t)) / length(filter(>(0.2), t)), sigdigits=2)) of the total.
    t > 3.0 ns: $(round(100 * length(filter(>(3.0), t)) / length(filter(>(0.2), t)), sigdigits=2)) of the total.
    """)
    return Plots.current()
end

function t_residence(t::Vector{Float64}; n=1)
    times = unique(t)
    prob = zeros(Float64, length(times))
    for (i, time) in enumerate(times)
        prob[i] = sum(t .== time) / length(t)
    end
    fit = EasyFit.fitexp(times, prob, n=n)
    return fit.b, fit.R
end

function t_residence(files::Vector{String}, nbins::Int64; lower=nothing, upper=nothing, n=1)
    t = Float64[]
    for file in files
        if filesize(file) == 0
            continue
        end
        data = readdlm(file)[:,1]
        if !isnothing(lower)
            data = filter(>(lower), data)
        end
        if !isnothing(upper)
            data = filter(<(upper), data)
        end
        append!(t, data)
    end
    #println(length(unique(t)))
    nbins = iszero(nbins) ? length(unique(t)) : nbins
    hist = fit(Histogram, t, nbins=nbins)
    idx = findall(h -> !iszero(h), hist.weights)
    prob = hist.weights[idx] / sum(hist.weights[idx])
    tbin = (hist.edges[1][1:end-1] .+ hist.edges[1][2:end]) ./ 2
    tbin = tbin[idx]
    fitting = EasyFit.fitexp(tbin, prob, n=n)
    return fitting.b, fitting.R
end

function t_residence(t, prob; n=1)
    fit = EasyFit.fitexp(t, prob, n=n)
    return fit.b, fit.R
end

function t_residence(M::BitMatrix; timestep=0.1)
    t = Float64[]
    sizehint!(t, size(M,1) * 10)
    for row in eachrow(M)
        trow = checking_residence(
            BitVector(row)
        )
        append!(t, timestep .* trow)
    end
    return filter(!iszero, t)
end

function residence(M::BitMatrix; step=1)
    if step <= 0
        throw(ArgumentError("step must be greater than 0"))
    end
    tmax = size(M, 2)
    τ = timelag(tmax, step)
    C = zeros(Float64, length(τ))
    for (i, Δt) in enumerate(τ)
        Σ, N = 0, 0
        for waters in eachrow(M)
            for t0 in 1:(tmax - Δt)
                Σ += waters[t0] && waters[t0 + Δt]
                N += waters[t0]
            end
        end
        C[i] = iszero(N) ? 0.0 : inv(N) * Σ
    end    
    return τ, C
end

function checking_residence(checklist::BitVector)
    if iszero(sum(checklist))
        return 0
    end
    counting = Int64[]
    counter = 0
    for present in checklist
        if present
            counter += 1
        else
            if counter > 0
                push!(counting, counter)
                counter = 0
            end
        end
    end
    if counter > 0
        push!(counting, counter)
    end
    return counting
end