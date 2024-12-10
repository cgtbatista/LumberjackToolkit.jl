function WHAM(coord1::Vector{Float64}, coord2::Vector{Float64}; bins1=0.0:20.0:360.0, bins2=0.0:20.0:360.0, kB=1.9872041e-3, T=298.15)
    
    n_bins1, Δbins1 = length(bins1) - 1, step(bins1)
    n_bins2, Δbins2 = length(bins2) - 1, step(bins2)
    
    bins_coord1 = floor.(Int64, (coord1 .- minimum(bins1)) ./ Δbins1) .+ 1
    bins_coord2 = floor.(Int64, (coord2 .- minimum(bins2)) ./ Δbins2) .+ 1

    counting = zeros(Int64, n_bins1, n_bins2)

    for (i,j) in zip(bins_coord1, bins_coord2)
        counting[i,j] += 1
    end

    prob = counting ./ sum(counting)
    β = inv(kB * T)

    pmf = -log.(prob) / β

    for i in eachindex(pmf)
        if pmf[i] == Inf
            pmf[i] = NaN
        end
    end ## dealing with de Inf values from log(0), since it is indefinite

    pmf .-= minimum(pmf[.!isnan.(pmf)]) # normalizing the PMF by the minimum value

    return pmf
end

function WHAM(coord::Vector{Float64}; bins=0.0:20.0:360.0, kB=1.9872041e-3, T=298.15)
    
    n_bins, Δbins = length(bins) - 1, step(bins)
    
    bins_coord = floor.(Int64, (coord .- minimum(bins)) ./ Δbins) .+ 1

    counting = zeros(Int64, n_bins)

    for i in bins_coord
        counting[i] += 1
    end

    prob = counting ./ sum(counting)
    β = inv(kB * T)

    pmf = -log.(prob) / β

    for i in eachindex(pmf)
        if pmf[i] == Inf
            pmf[i] = NaN
        end
    end ## dealing with de Inf values from log(0), since it is indefinite

    pmf .-= minimum(pmf[.!isnan.(pmf)]) # normalizing the PMF by the minimum value

    return pmf
end