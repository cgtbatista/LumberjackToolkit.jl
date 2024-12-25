@kwdef struct PCA
    λ::Vector{Float64}
    η::Matrix{Float64}
    variance::Vector{Float64}
    projection::Matrix{Float64}
end

function r_avg(frames::Vector{<:Vector{SVector{3, Float64}}})
    
    N, F = length(frames[1]), length(frames)

    r_avg = Vector{SVector{3, Float64}}(undef, N)

    for i in Base.OneTo(N)
        avg = zero(SVector{3, Float64})
        for j in Base.OneTo(F)
            avg += frames[j][i]
        end
        avg /= F
        r_avg[i] = avg
    end
    
    return r_avg
end


function X_disp(frames::Vector{<:Vector{SVector{3, Float64}}})
    
    N, F = length(frames[1]), length(frames)
    avg = r_avg(frames)

    X = zeros(Float64, 3*N, F)
    for j in Base.OneTo(F), i in Base.OneTo(N)
        for k in 1:3
            X[(i-1)*3+k, j] = frames[j][i][k] - avg[i][k]
        end
    end
    
    return X
end

function covariance(frames::Vector{<:Vector{SVector{3, Float64}}})
    
    X = X_disp(frames)
    bessel_correction = length(frames) - 1

    return (X * X') / bessel_correction
end

function pca(frames::Vector{<:Vector{SVector{3, Float64}}})
    C = covariance(frames)
    λ, η = pca(C)
    return PCA(
            λ,
            η,
            λ ./ sum(λ),
            η' * X_disp(frames)
        )
end

function pca(C::Matrix{Float64})    
    λ, η = eigen(C)
    idx = sortperm(λ, rev=true)
    return λ[idx], η[:, idx]
end