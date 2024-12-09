function WHAM(ϕ::Vector{Float64}, ψ::Vector{Float64}; bins=(0.0,360.0), Δbins=10.0, R=1.9872041e-3, T=298.15)
    
    angle_range = collect(bins[1]:Δbins:bins[2])
    n_bins = length(angle_range) - 1

    ϕ_bins = floor.(Int64, ϕ ./ Δbins) .+ 1
    ψ_bins = floor.(Int64, ψ ./ Δbins) .+ 1

    histogram = zeros(Int64, n_bins, n_bins)
    for (i,j) in zip(eachindex(ϕ_bins), eachindex(ψ_bins))
        if ϕ_bins[i] > 0 && ϕ_bins[i] <= n_bins && ψ_bins[j] > 0 && ψ_bins[j] <= n_bins
            histogram[i,j] += 1
        end
    end
    histogram .= histogram ./ sum(histogram)

    β = 1.0 / (R * T)

    pmf = -β * log.(histogram)
    pmf .-= minimum(pmf)

    return pmf
end