function avg_densityprofile(distances::Vector{Vector{Float64}}, densities::Vector{Vector{Float64}}; profile=nothing, type=nothing,
    xlabel="distance (Å)", ylabel="electron density (e/Å³)", label=:none, linewidth=2, bar_width=0.25, marker=:circle, ms=4, lcol=:blue, pcol=:red)
    
    avg_distances, avg_densities, std_densities = averaging_profile(distances, densities)
    
    plotting = Plots.plot(avg_distances, avg_densities, label=label, xlabel=xlabel, ylabel=ylabel, linewidth=linewidth, color=lcol)
    if type == "error"
        plotting = Plots.scatter!(avg_distances, avg_densities, yerr=std_densities, marker=marker, color=pcol, ms=ms)
    elseif type == "box"
        merging_densities = Float64[]
        for i in eachindex(densities)
            append!(merging_densities, densities[i])
        end
        merging_distances = Float64[]
        for i in eachindex(distances)
            append!(merging_distances, distances[i])
        end
        plotting = StatsPlots.boxplot!(merging_distances, merging_densities, label=label, outliers=false, bar_width=bar_width)
    end

    return plotting
end

## RMSD

function rmsd_plot(inputfile::String; timestep=1, t_max=400, xlabel="Time (ns)", ylim=(0, 10), label="cellulose", color="black", outputfile=nothing)
    
    outputfile = isnothing(outputfile) ? tempname() * ".png" : outputfile
    
    t = collect(0:timestep:t_max)

    rmsd = readdlm(inputfile)[Base.OneTo(t_max+1),1]

    if length(t) != length(rmsd)
        throw(ArgumentError("The time vector must have the same length of the RMSD vector."))
    end
    
    Plots.plot(t, rmsd, label=label, color=color, linewidth=2)
    Plots.plot!(xlabel=xlabel, ylabel="\n RMSD (Å)")
    Plots.plot!(frame=:box, grid=false)
    Plots.plot!(xlim=(t[begin], t[end]), ylim=ylim)
    Plots.plot!(size=(1000, 800))
    Plots.savefig(outputfile)

    return outputfile
end
## preto, vermelho, azul e verde

function rmsd_plot(inputfiles::Vector{String}; timestep=1, t_max=400, xlabel="Time (ns)", ylim=(0, 10), label=nothing, color=nothing, outputfile=nothing)
    
    outputfile = isnothing(outputfile) ? tempname() * ".png" : outputfile
    
    t = collect(0:timestep:t_max)

    Plots.gr()

    Plots.plot(xlabel=xlabel, ylabel="RMSD (Å)", legend=:topright)
    Plots.plot!(frame=:box, framestyle=:box, grid=false, tickfontsize=12, legendfontsize=12, guidefontsize=12, thickness_scaling=2)
    Plots.plot!(topmargin=0.5Plots.Measures.cm, rightmargin=0.5Plots.Measures.cm, bottommargin=-0.5Plots.Measures.cm, leftmargin=-0.5Plots.Measures.cm)
    Plots.plot!(xlim=(t[begin], t[end]), ylim=ylim)
    
    for input in enumerate(inputfiles)
        
        i, file = input[1], input[2]

        rmsd = readdlm(file)[Base.OneTo(t_max+1),1]
        
        tmp_label = isnothing(label) ? "$i" : label[i]
        tmp_color = isnothing(color) ? "black" : color[i]

        Plots.plot!(t, rmsd, label=tmp_label, color=tmp_color, linewidth=2)
    end

    Plots.plot!(size=(1200, 800), dpi=900)
    Plots.savefig(outputfile)

    return outputfile
end

## PLOT BY RESIDS
function byresid_plot(resids::Vector{Int64}, property::Vector{Float64}; title="", color=:coral2, output_filename="byresid.png")
       
    xresids = collect(1:length(property))
    Plots.plot(xresids, property, title=title, legend=:none, label=:none, lw=5, color=color)
    Plots.plot!(xlabel="# biose units", ylabel="Distance (Å)")
    Plots.plot!(size=(800,600))

    Plots.savefig(output_filename)
    
end


## DIHEDRAL

function hist_dihedral(filename::String; bins=50, normalize=:true, label="Dihedral distribution", output_filename="histogram.png")
    
    data = readdlm(filename)
    
    sums_dict = Dict{Tuple{Int64, Int64}, Vector{Float64}}()

    for row in eachindex(data[:,1])
        res_i, res_j = Int64(data[row, 2]), Int64(data[row, 3])
        key = (res_i, res_j)
        if haskey(sums_dict, key)
            push!(sums_dict[key], data[row,end])
        else
            sums_dict[key] = [data[row,end]]
        end
    end
    

    mean_sums = [ mean(sums_dict[key]) for key in keys(sums_dict) ]
    
    GR.histogram(mean_sums, bins=bins, normalize=normalize, label=label, legend=:topright)
    GR.histogram!(xlabel="ϕ + ψ (degrees)", ylabel="Frequencies")
    GR.histogram!(size=(800,600))

    GR.savefig(output_filename)
end

function hist_dihedral(files::Vector{String}; bins=50, label="xylan", normalize=:true, color=:blue, title="Distribution", output_filename="histogram.png")
 
    all_mean_sums = []
    
    for file in files

        data = readdlm(file)
        sums_dict = Dict{Tuple{Int64, Int64}, Vector{Float64}}()

        for row in eachindex(data[:,1])
            res_i, res_j = Int64(data[row, 2]), Int64(data[row, 3])
            key = (res_i, res_j)
            if haskey(sums_dict, key)
                push!(sums_dict[key], data[row,end])
            else
                sums_dict[key] = [data[row,end]]
            end
        end

        mean_sums = [ mean(sums_dict[key]) for key in keys(sums_dict) ]
        append!(all_mean_sums, mean_sums)

    end

    Plots.histogram(all_mean_sums, bins=bins, normalize=normalize, label=label, color=color, title=title, legend=:topright)
    Plots.histogram!(xlabel="ϕ + ψ (degrees)", ylabel="Frequencies")
    Plots.histogram!(size=(800,600))

    Plots.savefig(output_filename)  
end

function hist_gradient(files::Vector{String}; bins=50, label="xylan", normalize=:true, color=:blue, title="Distribution", output_filename="histogram.png")
 
    all_mean_sums = []
    
    for file in files

        data = readdlm(file)
        sums_dict = Dict{Tuple{Int64, Int64}, Vector{Float64}}()

        for row in eachindex(data[:,1])
            res_i, res_j = Int64(data[row, 2]), Int64(data[row, 3])
            key = (res_i, res_j)
            if haskey(sums_dict, key)
                push!(sums_dict[key], data[row,end])
            else
                sums_dict[key] = [data[row,end]]
            end
        end

        mean_sums = [ mean(sums_dict[key]) for key in keys(sums_dict) ]
        append!(all_mean_sums, mean_sums)

    end

    Plots.histogram(all_mean_sums, bins=bins, normalize=normalize, label=label, color=color, title=title, legend=:topright)
    Plots.histogram!(xlabel="ϕ + ψ (degrees)", ylabel="Frequencies")
    Plots.histogram!(size=(800,600))

    Plots.savefig(output_filename)  
end

## DIHEDRAL 2D

function hist_dihedral2D(filename::String; bins=30, normalize=:true, color=:plasma, title="Dihedral distribution 2D", output_filename="histogram2D.png")
    
    data = readdlm(filename)
    
    sums_dict1 = Dict{Tuple{Int64, Int64}, Vector{Float64}}()
    sums_dict2 = Dict{Tuple{Int64, Int64}, Vector{Float64}}()

    for row in eachindex(data[:,1])
        res_i, res_j = Int64(data[row, 2]), Int64(data[row, 3])
        key = (res_i, res_j)
        if haskey(sums_dict1, key) && haskey(sums_dict2, key)
            push!(sums_dict1[key], data[row,4])
            push!(sums_dict2[key], data[row,5])
        else
            sums_dict1[key] = [data[row,4]]
            sums_dict2[key] = [data[row,5]]
        end
    end
    

    mean_sums1 = [ mean(sums_dict1[key]) for key in keys(sums_dict1) ]
    mean_sums2 = [ mean(sums_dict2[key]) for key in keys(sums_dict2) ]
    
    GR.histogram2d(mean_sums1, mean_sums2, bins=(bins,bins), normalize=normalize, color=color, title=title, legend=:topright, show_empty_bins=:true)
    GR.histogram2d!(xlabel="ϕ (degrees)", ylabel="ψ (degrees)")
    GR.histogram2d!(size=(800,600))
    
    GR.savefig(output_filename)
end

function hist_dihedral2D(files::Vector{String}; bins=80, normalize=:true, color=:plasma, title="Dihedral distribution 2D", output_filename="histogram2D.png")
    
    all_mean_sums1 = []
    all_mean_sums2 = []

    for file in files

        data = readdlm(file)
        
        sums_dict1 = Dict{Tuple{Int64, Int64}, Vector{Float64}}()
        sums_dict2 = Dict{Tuple{Int64, Int64}, Vector{Float64}}()

        for row in eachindex(data[:,1])
            res_i, res_j = Int64(data[row, 2]), Int64(data[row, 3])
            key = (res_i, res_j)
            if haskey(sums_dict1, key) && haskey(sums_dict2, key)
                push!(sums_dict1[key], data[row,4])
                push!(sums_dict2[key], data[row,5])
            else
                sums_dict1[key] = [data[row,4]]
                sums_dict2[key] = [data[row,5]]
            end
        end  

        mean_sums1 = [ mean(sums_dict1[key]) for key in keys(sums_dict1) ]; append!(all_mean_sums1, mean_sums1)
        mean_sums2 = [ mean(sums_dict2[key]) for key in keys(sums_dict2) ]; append!(all_mean_sums2, mean_sums2)
    
    end
    
    Plots.histogram2d(all_mean_sums1, all_mean_sums2, bins=(bins,bins), normalize=normalize, color=color, title=title, legend=:topright, show_empty_bins=:true)
    Plots.histogram2d!(xlabel="ϕ (degrees)", ylabel="ψ (degrees)")
    Plots.histogram2d!(size=(800,600))
    Plots.savefig(output_filename)
end

function dihedrals_plot(phi::Vector{Float64}, psi::Vector{Float64}; title="Dihedrals Scatter", output_filename="dihedrals.png")
       
    Plots.plot(phi, psi, title=title, legend=:none, seriestype=:scatter, marker=2)
    Plots.plot!(xlabel="ϕ (degrees)", ylabel="ψ (degrees)")
    Plots.plot!(size=(800,600))

    Plots.savefig(output_filename)
    
end

function dist_by_resids(property::Vector{Float64}; title="", color=:coral2, output_filename="dist_resids.png")
       
    xresids = collect(1:length(property))
    Plots.plot(xresids, property, title=title, legend=:none, label=:none, lw=5, color=color)
    Plots.plot!(xlabel="# biose units", ylabel="Distance (Å)")
    Plots.plot!(size=(800,600))

    Plots.savefig(output_filename)
    
end

function distance_plot(phi::Vector{Float64}, psi::Vector{Float64}, property::Vector{Float64}; bins=100, normalize=:true, color=:turbo, output_filename="dihedrals_dist.png")
       
    sum_angles_vec = [];
    for i in eachindex(phi)
        sum_angles = phi[i] + psi[i]
        if (sum_angles > 360) && (sum_angles < 720)
            sum_angles = sum_angles - 360
        elseif sum_angles > 720
            sum_angles = sum_angles - 720
        end
        push!(sum_angles_vec, sum_angles)
    end
    Plots.histogram2d(sum_angles_vec, property, bins=(bins,bins), normalize=normalize, color=color, legend=:none, show_empty_bins=:true, colorbar=:true)
    Plots.histogram2d!(xlabel="ϕ + ψ (degrees)", ylabel="Closest distance (Å)")
    #Plots.histogram2d!(ylims=(0, 20))
    Plots.histogram2d!(size=(800,600))

    Plots.savefig(output_filename)
    
end

function dihedrals_contour(phi::Vector{Float64}, psi::Vector{Float64}, property::Vector{Float64}; gridsize=50, color=:viridis, title="Contourn", output_filename="contourn.png")
    
    newphi, newpsi, newproperty = GR.gridit(phi, psi, property, gridsize, gridsize)
    countorn_levels = gridsize/10; countorn_levels = convert(Int64, countorn_levels);
    Plots.contourf(newphi, newpsi, newproperty, levels=countorn_levels, color=color, title=title, legend=:none, lw=0.5)
    Plots.contourf!(xlabel="ϕ (degrees)", ylabel="ψ (degrees)")
    Plots.contourf!(size=(800,600))

    Plots.savefig(output_filename)
    
end


function rt_plot(filename::String; timelength = 10, md_time="ns", output_filename="rtime_plot.png")
    
    closest_water_time = Dict{Tuple{Int64, String}, Vector{Int64}}()
    closest_water_dist = Dict{Tuple{Int64, String}, Vector{Float64}}()

    water_data = readdlm(filename)

    for irow in eachindex(water_data[:,1])
        if isnan(water_data[irow, end])
            continue
        end
        resid, segid = Int64(water_data[irow, 3]), String(water_data[irow, 4]) 
        key = (resid, segid)
        if haskey(closest_water_time, key) && haskey(closest_water_dist, key)
            push!(closest_water_time[key], water_data[irow,1])
            push!(closest_water_dist[key], water_data[irow,end])
        else
            closest_water_time[key] = [water_data[irow,1]]
            closest_water_dist[key] = [water_data[irow,end]]
        end
    end
    
    ## Residence time
    water_molecules = collect(1:length(closest_water_time)); unique_water_molecules = [ key for key in keys(closest_water_time) ]
    lignin_label = split(split(filename, ".")[1], "_")[2]
    #a
    println("~ $lignin_label")
    println("  $(size(water_molecules)[1]) water molecules in range")
    for iwater in eachindex(unique_water_molecules)
        ikey = unique_water_molecules[iwater]; water_residence = closest_water_time[ikey]
        if iwater == 1;
            println("  $(size(water_residence)[1]) water molecules in $(ikey[2])$(ikey[1])")
        end
    end
    println("  max contact time: on n molecules"); println("  min contact time: on n molecules")
    println("  retention time")

    #b plotting
    Plots.plot(xlabel="MD simulation time ($md_time)", ylabel="# water molecule", xlims=(1, timelength), ylims=(0, maximum(water_molecules)+1),
               legend=:none, size=(800, 600))
    for iwater in eachindex(unique_water_molecules)
        ikey = unique_water_molecules[iwater]; water_residence = closest_water_time[ikey]
        imolecule = [ water_molecules[iwater] for i in 1:length(water_residence) ]
        Plots.plot!(water_residence, imolecule, colour=:black, linewidth=4)
    end
    savefig("rtime_$output_filename")
    ## Plotting the distance
    
    println("- water molecules mean distance (histogram like 2D to see better effect of time in the distribution)")
    println("- relantionship between retention time and distance (linear regression??)")    
    
    #plot(
    #        timestep, rmsd, xlabel="MD simulation time ($md_time)", ylabel="RMSD", label=data_label[1], linewidth=2,
    #        xlim=(minimum(timestep), maximum(timestep)+60), ylim=(minimum(rmsd), maximum(rmsd)+5),
    #        size=(800,600)
    #    )

end


function rtime_plot(files::Vector{String}; timelength = 10, md_time="ns", output_filename="full_plot.png")
    
    
    water_code = []; closest_time = []; delta_water = 0; last_water = 0;
    Plots.plot(xlabel="MD simulation time ($md_time)", ylabel="# water molecule", legend=:none)

    for file in files
        closest_water_time = Dict{Tuple{Int64, String}, Vector{Int64}}()
        closest_water_dist = Dict{Tuple{Int64, String}, Vector{Float64}}()

        water_data = readdlm(file)

        for irow in eachindex(water_data[:,1])
            if isnan(water_data[irow, end])
                continue
            end
            resid, segid = Int64(water_data[irow, 3]), String(water_data[irow, 4]) 
            key = (resid, segid)
            if haskey(closest_water_time, key) && haskey(closest_water_dist, key)
                push!(closest_water_time[key], water_data[irow,1])
                push!(closest_water_dist[key], water_data[irow,end])
            else
                closest_water_time[key] = [water_data[irow,1]]
                closest_water_dist[key] = [water_data[irow,end]]
            end
        end

        dummy_water_vector = collect(1:length(closest_water_time))
        water_molecules = copy(dummy_water_vector) .+ delta_water; last_water=maximum(water_molecules)
        delta_water += length(closest_water_time);
        unique_water_molecules = [ key for key in keys(closest_water_time) ]
    
        #b plotting

        for iwater in eachindex(unique_water_molecules)
            ikey = unique_water_molecules[iwater]; water_residence = closest_water_time[ikey]
            imolecule = [ water_molecules[iwater] for i in 1:length(water_residence) ]
            Plots.plot!(water_residence, imolecule, colour=:black, linewidth=1)
        end
    end

    Plots.plot!(xlims=(1, timelength), ylims=(0, last_water+1))
    Plots.plot!(size=(800, 600))

    Plots.savefig("rtime_$output_filename")

    
    #plot(
    #        timestep, rmsd, xlabel="MD simulation time ($md_time)", ylabel="RMSD", label=data_label[1], linewidth=2,
    #        xlim=(minimum(timestep), maximum(timestep)+60), ylim=(minimum(rmsd), maximum(rmsd)+5),
    #        size=(800,600)
    #    )

end


function pmf1D(
                pmf::Vector{Float64}; coord=0:20:360,
                xaxis=0:60:360, title="PMF at kcal/mol", xlabel="ϕ (degrees)", ylabel="Potential of Mean Force",
                color=:blue, legend=:none, size=(800,600),
                output_filename=nothing
            )

    output_filename = isnothing(output_filename) ? tempname() * ".png" : output_filename

    new_coord = coord[1:end-1] .+ diff(coord) ./ 2

    Plots.scatter(new_coord, pmf, color=color, title=title, legend=legend)
    Plots.scatter!(
            xlabel=xlabel, xlims=(xaxis[begin], xaxis[end]), xticks=xaxis,
            ylabel=ylabel
        )
    Plots.scatter!(
                frame=:box, grid=false,
                tickfontsize=12, legendfontsize=12, guidefontsize=12,
                thickness_scaling=1.5
            )
    Plots.scatter!(size=size)

    Plots.savefig(output_filename)

    return output_filename
end


function pmf2D(
                pmf::Matrix{Float64}; type="heatmap", coord1=0:20:360, coord2=0:20:360,
                xaxis=0:60:360, yaxis=0:60:360, title="PMF at kcal/mol", xlabel="ϕ (degrees)", ylabel="ψ (degrees)",
                levels=5, lw=2.0, color=:plasma, color2=[:black], cbar=true, clabels=true, legend=:none, size=(800,600),
                output_filename=nothing
            )

    output_filename = isnothing(output_filename) ? tempname() * ".png" : output_filename

    if type == "heatmap"
        
        maxpmf = maximum(pmf[.!isnan.(pmf)])
        minpmf = minimum(pmf[.!isnan.(pmf)])

        Plots.heatmap(coord1, coord2, pmf', color=color, title=title, legend=legend, clims=(minpmf, maxpmf))
        Plots.heatmap!(
                xlabel=xlabel, xlims=(xaxis[begin], xaxis[end]), xticks=xaxis,
                ylabel=ylabel, ylims=(yaxis[begin], yaxis[end]), yticks=yaxis
            )
        Plots.heatmap!(
                frame=:box, grid=false,
                tickfontsize=12, legendfontsize=12, guidefontsize=12,
                thickness_scaling=1.5
            )
        Plots.heatmap!(size=size) 

    elseif type == "contour"

        new_coord1 = coord1[1:end-1] .+ diff(coord1) ./ 2
        new_coord2 = coord2[1:end-1] .+ diff(coord2) ./ 2

        Plots.contourf(
                new_coord1, new_coord2, pmf', levels=levels, color=color, title=title, legend=legend, cbar=cbar, clabels=clabels,
                lw=lw, lc=color2
            )
        Plots.contourf!(
                xlabel=xlabel, xlims=(xaxis[begin], xaxis[end]), xticks=xaxis,
                ylabel=ylabel, ylims=(yaxis[begin], yaxis[end]), yticks=yaxis
            )
        Plots.contourf!(
                frame=:box, grid=false,
                tickfontsize=12, legendfontsize=12, guidefontsize=12,
                thickness_scaling=1.5
            )
        Plots.contour!(size=size)

    else
        throw(ArgumentError("Invalid type of plot: $type. Try to use 'heatmap' or 'contour'"))
    end

    Plots.savefig(output_filename)

    return output_filename
end

# vmd_makesystems("crystal", cellulose_habits=12345432, pkgload=packages_list)
# vmd_getsurface("new_crystal")
# rdata = radii_distribution("new_crystal.pdb", "new_crystal.dat", dstep=0.2)[4]
# idx_cutoff = findall(ϕ -> ϕ >= 90.0-dϕ && ϕ <= 90.0+dϕ, rdata[3])
# idx = findall(x -> x > 0., rdata[1][idx_cutoff])
# d4 = 2*rdata[1][idx_cutoff][idx]

# vmd_makesystems("crystal", cellulose_habits=33333333, pkgload=packages_list)
# vmd_getsurface("new_crystal")
# rdata = radii_distribution("new_crystal.pdb", "new_crystal.dat", dstep=0.2)[4]
# idx_cutoff = findall(ϕ -> ϕ >= 90.0-dϕ && ϕ <= 90.0+dϕ, rdata[3])
# idx = findall(x -> x > 0., rdata[1][idx_cutoff])
# d5 = 2*rdata[1][idx_cutoff][idx]


# xyz = radii_distribution("new_crystal.pdb", "new_crystal.dat", dstep=0.2)[1]
# ## 333333 234432 12333321 33333333 12345432
# println("The mean radius is $(mean(d3)) Å with a standard deviation of $(std(d3)) Å and a total of $(length(d3)) points.")
# println("The mean radius is $(mean(d2)) Å with a standard deviation of $(std(d2)) Å and a total of $(length(d2)) points.")
# println("The mean radius is $(mean(d1)) Å with a standard deviation of $(std(d1)) Å and a total of $(length(d1)) points.")
# println("The mean radius is $(mean(d5)) Å with a standard deviation of $(std(d5)) Å and a total of $(length(d5)) points.")
# println("The mean radius is $(mean(d4)) Å with a standard deviation of $(std(d4)) Å and a total of $(length(d4)) points.")

# density(d3, bandwidth=0.8, label="333333", xlabel="Fibril diameter (Å)", ylabel="Frequency", title="18-chain", color="black", lw=3)
# density!(d2, bandwidth=0.8, label="234432", color="lightblue", lw=3)
# density!(d1, bandwidth=0.8, label="12333321", color="red", lw=3)
# plot!(size=(800, 800))
# savefig("18chain.png")
# density(d5, bandwidth=0.8, label="33333333", xlabel="Fibril diameter (Å)", ylabel="Frequency", title="24-chain", color="black", lw=3)
# density!(d4, bandwidth=0.8, label="12345432", color="lightblue", lw=3)
# plot!(size=(800, 800))
# savefig("24chain.png")


# ## violin plot for the all systems
# #function custom_kde(x; bw)
# #    kde(x, bandwidth=bw)
# #end

# d = [d1; d2; d3; d4; d5]
# smoothed_d = kde(d, bandwidth=0.8)
# r = [fill("12333321", length(d1)); fill("333333", length(d2)); fill("234432", length(d3)); fill("12345432", length(d4)); fill("33333333", length(d5))]
# violin(r, d, bw=0.5, xlabel="Cellulose habits", ylabel="Fibril diameter (Å)", title="Fibril diameter distribution", color=:steelblue, lw=2, label=:none)
# boxplot!(r, d, color=:coral2, lw=2, bar_width=0.25, label=:none, outliers=:false)
# plot!(size=(800,800), tickfontsize=11, guidefontsize=14)
# savefig("violin.png")


# csvdata = open("microfibrils.csv", "w")
# Base.write(csvdata, "habit,diameter\n")
# for i in d3
#     Base.write(csvdata, "333333,$i\n")
# end
# for i in d2
#     Base.write(csvdata, "234432,$i\n")
# end
# for i in d1
#     Base.write(csvdata, "12333321,$i\n")
# end
# for i in d5
#     Base.write(csvdata, "33333333,$i\n")
# end
# for i in d4
#     Base.write(csvdata, "12345432,$i\n")
# end

function fibril_surface(
    xyz::Vector{SVector{3, Float64}};
    centers=nothing, filename=nothing,
    legend=:none, size=(800, 800), color1=:blue, color2=:red, marker=1
)
    Plots.plot(xlabel="y (Å)", ylabel="z (Å)", zlabel="x (Å)", legend=legend, size=size) 
    Plots.scatter!(
        [ vec[2] for vec in xyz ], [ vec[3] for vec in xyz ], [ vec[1] for vec in xyz ],
        color=color1, marker=marker
    )
    if !isnothing(centers)
        Plots.scatter!(
            [ vec[2] for vec in centers ], [ vec[3] for vec in centers ], [ vec[1] for vec in centers ],
            color=color2, marker=marker
        )
    end
    #Plots.savefig(output_filename)
    return Plots.current()
end

function fibril_slice(
    xyz::Vector{SVector{3, Float64}}, step::Int64;
    Δz=1.0, legend=:none, size=(800, 800), color1=:blue, marker=1,
    filename=nothing
)
    Plots.plot(xlabel="y (Å)", ylabel="x (Å)", legend=legend, size=size)
    isorted = sortperm([ vec[3] for vec in xyz ])
    x, y, z = [ vec[2] for vec in xyz[isorted] ], [ vec[1] for vec in xyz[isorted] ], [ vec[3] for vec in xyz[isorted] ]
    idx = findall(k -> norm(k - unique(z)[step]) <= Δz, z)
    Plots.scatter!(x[idx], y[idx], color=color1, marker=marker)
    #Plots.savefig(output_filename)
    return Plots.current()
end


function density_profile(
    bin1::Matrix{Float64}, bin2::Matrix{Float64}, densities::Array{Float64, 3};
    frame=nothing, smoothing=false, bins=200, colors=:Greens,
    xlabel="x-axis (Å)", ylabel="y-axis (Å)", xlims=nothing, ylims=nothing
)
    @assert size(bin1, 2) == size(bin2, 2) == size(densities, 3) "The number of frames should be the same."
    xlims = !isnothing(xlims) ? (-xlims, xlims) : (minimum(bin1), maximum(bin1))
    ylims = !isnothing(ylims) ? (-ylims, ylims) : (minimum(bin2), maximum(bin2))
    if isnothing(frame)
        edp_plotting_averages(
            bin1, bin2, densities; smoothing=smoothing, bins=bins, colors=colors,
            xlabel=xlabel, ylabel=ylabel, xlims=xlims, ylims=ylims
        )
    end
    return Plots.current()
end

function edp_plotting_averages(
    bin1, bin2, densities; smoothing=true, bins=200, colors=:Greens,
    xlabel="x-axis (Å)", ylabel="y-axis (Å)", xlims=nothing, ylims=nothing
)
    nframes = size(densities, 3)
    projx = smoothing ? zeros(Float64, bins) : zeros(Float64, size(bin1, 1))
    projy = smoothing ? zeros(Float64, bins) : zeros(Float64, size(bin2, 1))
    ρ_new = smoothing ? zeros(Float64, bins, bins) : zeros(Float64, size(densities, 1), size(densities, 2))
    for iframe in 1:nframes
        x, y, z = if smoothing
                interpol(bin1[:,iframe], bin2[:,iframe], densities[:,:,iframe], bins=bins)
            else
                bin1[:,iframe], bin2[:,iframe], densities[:,:,iframe]
        end
        projx .+= x
        projy .+= y
        ρ_new .+= z
    end
    projx ./= nframes
    projy ./= nframes
    ρ_new ./= nframes
    Plots.gr(size=(1000,800), dpi=900, fmt=:png)
    Plots.heatmap(
        projx, projy, ρ_new',
        color=colors,
        clims=(minimum(ρ_new), maximum(ρ_new)),
        xlabel=xlabel, ylabel=ylabel, fontfamily=:arial,
        ## Axis configs
        xlims=xlims, ylims=ylims,
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
end