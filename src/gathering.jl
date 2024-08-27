## GET SOME PROPERTY/DATA PAIRWISE TO THE DIHEDRAL ANGLES
function property_by_resids(filename::String; framebyframe=:true, lower_limit=3, greater_limit=98, first_frame=0)

    property_data = readdlm(filename)

    if framebyframe

        extended_resids = []; extended_property = []
        for irow in eachindex(property_data[:,end])
            if property_data[irow, 1] >= first_frame
                push!(extended_resids, property_data[irow, 2])
                push!(extended_property, property_data[irow, end])
            end
        end

        extended_resids = convert(Vector{Int64}, extended_resids)
        extended_property = convert(Vector{Float64}, extended_property)
        return extended_resids, extended_property

    else

        property = Dict{Int64, Vector{Float64}}()

        for irow in eachindex(property_data[:,end])
            if property_data[irow, 1] >= first_frame
                ith_resid = Int64(property_data[irow, 2])
                if (ith_resid >= lower_limit) && (ith_resid < greater_limit)
                    catched_property = property_data[irow, end]
                    if haskey(property, ith_resid)
                        push!(property[ith_resid], catched_property)
                    else
                        property[ith_resid] = [catched_property]
                    end
                end
            end
        end

        mean_resids = [ i for i in keys(property)]
        mean_property = [ mean(property[key]) for key in keys(property) ]
        mean_resids = convert(Vector{Int64}, mean_resids); mean_property = convert(Vector{Float64}, mean_property);
        return mean_resids, mean_property

    end  

end


function property_in_dihedrals(filename1::String, filename2::String; framebyframe=:true, property_option="distance", lower_limit=3, greater_limit=98)
    
    # filename1 = "dihedral data" && filename2 = "property data"
    dihedral_data = readdlm(filename1); length_data=size(dihedral_data)[1]
    property_data = readdlm(filename2)
    
    local pairwise_property, extended_pairwise_property, mean_property
    ## getting the property dictionary
    if property_option == "distance"
        
        pairwise_property = Dict{Int64, Vector{Float64}}()
        extended_pairwise_property = []

        for irow in eachindex(property_data[:,1])
            iresnum = Int64(property_data[irow, 2])
            if (iresnum >= lower_limit) && (iresnum < greater_limit)
                jrow = irow + 1; mean = 0.5 * (property_data[irow, end] + property_data[jrow, end])
                push!(extended_pairwise_property, mean)
                if haskey(pairwise_property, iresnum)
                    push!(pairwise_property[iresnum], mean)
                else
                    pairwise_property[iresnum] = [mean]
                end
            end
        end
        mean_property = [ mean(pairwise_property[key]) for key in keys(pairwise_property) ]
        extended_pairwise_property = [ extended_pairwise_property[i] for i in eachindex(extended_pairwise_property) ]

    end

    local extendeded_phi, extendeded_psi
    ## getting the dihedrals dictionary for ϕ and ψ angles
    sums_dict1 = Dict{Tuple{Int64, Int64}, Vector{Float64}}(); extendeded_phi = []
    sums_dict2 = Dict{Tuple{Int64, Int64}, Vector{Float64}}(); extendeded_psi = []

    for row in eachindex(dihedral_data[:,1])
        res_i, res_j = Int64(dihedral_data[row, 2]), Int64(dihedral_data[row, 3])
        key = (res_i, res_j)
        push!(extendeded_phi, dihedral_data[row,4]); push!(extendeded_psi, dihedral_data[row,5])
        if haskey(sums_dict1, key) && haskey(sums_dict2, key)
            push!(sums_dict1[key], dihedral_data[row,4])
            push!(sums_dict2[key], dihedral_data[row,5])
        else
            sums_dict1[key] = [dihedral_data[row,4]]
            sums_dict2[key] = [dihedral_data[row,5]]
        end
    end  
    mean_sums1 = [ mean(sums_dict1[key]) for key in keys(sums_dict1) ]; extendeded_phi = [ extendeded_phi[i] for i in eachindex(extendeded_phi) ]
    mean_sums2 = [ mean(sums_dict2[key]) for key in keys(sums_dict2) ]; extendeded_psi = [ extendeded_psi[i] for i in eachindex(extendeded_psi) ]
    
    if framebyframe
        extendeded_phi = convert(Vector{Float64}, extendeded_phi)
        extendeded_psi = convert(Vector{Float64}, extendeded_psi)
        extended_pairwise_property = convert(Vector{Float64}, extended_pairwise_property)
        return extendeded_phi, extendeded_psi, extended_pairwise_property
    else
        mean_sums1 = convert(Vector{Float64}, mean_sums1)
        mean_sums2 = convert(Vector{Float64}, mean_sums2)
        mean_property = convert(Vector{Float64}, mean_property)
        return mean_sums1, mean_sums2, mean_property
    end

end

#################################################################################################################################

function property_by_resids(segname::Vector{String}; framebyframe=:true, lower_limit=3, greater_limit=98, first_frame=0)
    
    resids = []; property = [];
    
    for isegname in segname
        filename = "closest_$(isegname).dat"
        a, b = property_by_resids(filename; framebyframe=framebyframe, lower_limit=lower_limit, greater_limit=greater_limit, first_frame=first_frame)
        append!(resids, a); append!(property, b)
    end
       
    return convert(Vector{Int64}, resids), convert(Vector{Float64}, property)
end

function property_in_dihedrals(segname::Vector{String}; framebyframe=:true, property_option="distance", lower_limit=3, greater_limit=98)
    
    phi_angles = []; psi_angles = []; pairwise_property = [];
    
    for isegname in segname
        filename1 = "dihedral2_$(isegname).dat"; filename2 = "closest_$(isegname).dat"
        a, b, c = property_in_dihedrals(filename1, filename2; framebyframe=framebyframe, property_option=property_option, lower_limit=lower_limit, greater_limit=greater_limit)
        append!(phi_angles, a); append!(psi_angles, b); append!(pairwise_property, c)
    end
    phi_angles = convert(Vector{Float64}, phi_angles); psi_angles = convert(Vector{Float64}, psi_angles);
    pairwise_property = convert(Vector{Float64}, pairwise_property)
    return phi_angles, psi_angles, pairwise_property    
end

################################################################################################################################

function retention_time(filename::String; nframes = 10, timestep=1, md_time="ns", dist_cutoff=10.0)
    
    closest_water_time = Dict{Tuple{Int64, String}, Vector{Int64}}()
    closest_water_dist = Dict{Tuple{Int64, String}, Vector{Float64}}()

    water_data = readdlm(filename)

    for irow in eachindex(water_data[:,1])
        if isnan(water_data[irow, end]) || water_data[irow, end] > dist_cutoff
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
    filename_label = split(split(filename, ".")[1], "_")[2]
    #a
    println("~ $filename_label")
    println("  $(size(water_molecules)[1]) water molecules in range")
    filename_array = []; ith_water_array = [];
    min_rtime_array = []; max_rtime_array = []; contact_array = []; fluctuation_array = []; mean_dist_array = []; std_dist_array = []; closest_dist_array = [];
    for iwater in eachindex(unique_water_molecules)
        ikey = unique_water_molecules[iwater]; iwater_label = string( ikey[2], ikey[1])
        water_time = closest_water_time[ikey]; water_dist = closest_water_dist[ikey]
        if iwater != 0; ## == 1 just to dummy the simulation
            println("  $(iwater) - the water molecule $iwater_label appeared $(size(water_time)[1])/$nframes:")
            min_rtime=0; max_rtime=0; contact_time=0; fluctuation=0;
            mean_dist=mean(water_dist); std_dist=std(water_dist); closest_dist=minimum(water_dist);
            el=1; dt=0; lost_contact=0;
            while el < length(water_time)
                if water_time[el+1]-water_time[el] == 1
                    dt += 1; contact_time += 1; max_rtime += 1;
                else
                    lost_contact += 1
                    if min_rtime > dt; min_rtime = dt; end
                    if max_rtime < dt; min_rtime = max_rtime; max_rtime = dt; end
                    dt = 0
                end
                if length(water_time)-max_rtime == 1; min_rtime = max_rtime; end
                el += 1
            end
            max_rtime = max_rtime * timestep; min_rtime = min_rtime * timestep; contact_time = contact_time * timestep; fluctuation = lost_contact;
            push!(filename_array, filename_label); push!(ith_water_array, iwater_label);
            push!(min_rtime_array, min_rtime); push!(max_rtime_array, max_rtime); push!(contact_array, contact_time); push!(fluctuation_array, fluctuation);
            push!(mean_dist_array, mean_dist); push!(std_dist_array, std_dist); push!(closest_dist_array, closest_dist);
            println("  + the water retention time reached a minimum of $min_rtime $md_time and a maximum of $max_rtime $md_time.")
            println("  + the water spent $contact_time $md_time on surroundings, but it lost the contact $lost_contact time(s).")
            println("  + the water molecule mean distance was ($(round(mean_dist, digits=2)) ± $(round(std_dist, digits=2))) Å, with $(round(closest_dist, digits=2)) Å as the closest distance.")
            println("")
        end
        ## ikey = iésima água do vetor unique_water_molecules, onde ikey[2] é o segmento e ikey[1] é o número do resíduo
    end
    return filename_array, ith_water_array, min_rtime_array, max_rtime_array, contact_array, fluctuation_array, mean_dist_array, std_dist_array, closest_dist_array
    #println("- water molecules mean distance (histogram like 2D to see better effect of time in the distribution)")
    #println("- relantionship between retention time and distance (linear regression??)")       

end


function retention_time(files::Vector{String}; nframes = 10, timestep=1, md_time="ns", dist_cutoff=10.0)

    filename_array = []; ith_water_array = [];
    min_rtime_array = []; max_rtime_array = []; contact_array = []; fluctuation_array = []; mean_dist_array = []; std_dist_array = []; closest_dist_array = []

    for file in files
        a, b, c, d, e, f, g, h, i = retention_time(file; nframes=nframes, timestep=timestep, md_time=md_time, dist_cutoff=dist_cutoff)
        append!(filename_array, a); append!(ith_water_array, b);
        append!(min_rtime_array, c); append!(max_rtime_array, d); append!(contact_array, e); append!(fluctuation_array, f);
        append!(mean_dist_array, g); append!(std_dist_array, h); append!(closest_dist_array, i);

    end
    filename_array = convert(Vector{String}, filename_array); ith_water_array = convert(Vector{String}, ith_water_array);
    min_rtime_array = convert(Vector{Float64}, min_rtime_array); max_rtime_array = convert(Vector{Float64}, max_rtime_array);
    contact_array = convert(Vector{Float64}, contact_array); fluctuation_array = convert(Vector{Float64}, fluctuation_array);
    mean_dist_array = convert(Vector{Float64}, mean_dist_array); std_dist_array = convert(Vector{Float64}, std_dist_array);
    closest_dist_array = convert(Vector{Float64}, closest_dist_array);
    return filename_array, ith_water_array, min_rtime_array, max_rtime_array, contact_array, fluctuation_array, mean_dist_array, std_dist_array, closest_dist_array

end
