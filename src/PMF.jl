function WHAM(ϕ::Vector{Float64}, ψ::Vector{Float64}; bins=0.0:10.0:360.0, R=1.9872041e-3, T=298.15)
    
    n_bins = length(bins) - 1
    
    Δbins = step(bins)
    ϕ_bins = floor.(Int64, (ϕ .- minimum(bins)) ./ Δbins) .+ 1
    ψ_bins = floor.(Int64, (ψ .- minimum(bins)) ./ Δbins) .+ 1

    histogram = zeros(Float64, n_bins, n_bins)

    for (i,j) in zip(ϕ_bins, ψ_bins)
        histogram[i,j] += 1
    end

    histogram .= histogram ./ sum(histogram)

    β = 1.0 / (R * T)

    pmf = -β * log.(histogram)
    pmf .-= minimum(pmf)

    for i in eachindex(pmf)
        if pmf[i] == Inf
            pmf[i] = NaN
        end
    end ## adding a pseudo-count to avoid Inf values

    return pmf
end