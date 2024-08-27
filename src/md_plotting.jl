## RMSD

function rmsd_plot(rmsd_data::String; timestep=1, md_time="ns", output_filename="rmsdplot.png")
    # Load data
    all_data = readdlm(rmsd_data)
    simulation_time = collect(0:size(all_data)[1]) * timestep; rmsd = [ 0.0; all_data[:,1] ]
    # Plot
    data_label=split(rmsd_data, ".")
    Plots.plot(simulation_time, rmsd, label=data_label[1], linewidth=2)
    Plots.plot!(xlabel="MD simulation time ($md_time)", ylabel="RMSD")
    Plots.plot!(xlim=(minimum(simulation_time), maximum(simulation_time)+60), ylim=(minimum(rmsd), maximum(rmsd)+5))
    Plots.plot!(size=(800, 600))
    
    Plots.savefig(output_filename)
end

function rmsd_plot(files::Vector{String}; timestep=1, md_time="ns", output_filename="rmsdplot.png")
    
    Plots.plot(xlabel="MD simulation time ($md_time)", ylabel="RMSD", legend=:bottomright)
    Plots.plot!(xlims=(Inf, -Inf), ylims=(Inf, -Inf))
    Plots.plot!(size=(800, 600))
    
    for file in files
        all_data = readdlm(file)
        simulation_time = collect(0:size(all_data)[1]) * timestep
        rmsd = [0.0; all_data[:, 1]]
        
        data_label = split(file, ".")[1]
        
        Plots.plot!(simulation_time, rmsd, label=data_label, linewidth=2)
        Plots.plot!(xlims=(minimum(simulation_time), maximum(simulation_time) + 60), ylims=(minimum(rmsd), maximum(rmsd) + 5))
    end
    
    Plots.savefig("full_$output_filename")
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