function _pca_data()
    pdbname = "/home/carlos/Documents/Sandbox/Cel6A/petncel/setup/system.pdb"
    trjname = "/tmp/jl_azqdEpyxch.pdb"
    return pdbname, trjname
end

@kwdef struct PCA
    λ::Vector{Float64}                      ## Eigenvalues
    η::Matrix{Float64}                      ## Eigenvectors
    variance::Vector{Float64}               ## Explained Variance
    projection::Matrix{Float64}             ## Projection of the trajectory in the PCA space
    ω::Vector{Float64}                      ## Quasiharmonic Frequencies -- on cm^-1 or rad/s
    Δx::Matrix{Float64}                     ## Quasiharmonic modes of vibration -- on Å
end

"""
    r_avg(frames::Vector{<:Vector{SVector{3, Float64}}})

Calculate the average position of each atom inside the trajectory.
"""
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


"""
    dispm(frames::Vector{<:Vector{SVector{3, Float64}}}, weights::Vector{Float64}; weight_method=nothing)

Calculate the displacement matrix of a trajectory with weights.
"""
function dispm(frames::Vector{<:Vector{SVector{3, Float64}}}; method=nothing, weights=nothing)
    avg = r_avg(frames)
    
    N, F = length(frames[1]), length(frames)
    if !isnothing(method) && (length(weights) != N)
        throw(ArgumentError("The length of the weights vector must be equal to the number of atoms in the frames."))
    end

    X = zeros(Float64, 3*N, F)
    for j in Base.OneTo(F), i in Base.OneTo(N)
        for k in 1:3
            if isnothing(method)
                weight = 1.
            elseif lowercase(method) == "sqrtmass"
                weight = sqrt(weights[i])
            elseif lowercase(method) == "mass"
                weight = weights[i]
            else
                throw(ArgumentError("The weight method $method is not implemented."))
            end
            X[(i-1)*3+k, j] = (frames[j][i][k] - avg[i][k]) * inv(weight)
        end
    end
    
    return X
end

function covarm(X::Matrix{Float64}; bessel=true)
    if bessel
        return (X * X') / (size(X, 2) - 1)
    else
        return (X * X') / size(X, 2)
    end
end

function covarm(frames::Vector{<:Vector{SVector{3, Float64}}}; bessel=true)    
    return covarm(dispm(frames), bessel=bessel)
end

function pca(
            pdbname::String,
            trjname::String;
            selection="protein and name CA",
            mass_weighting=true, quasiharmonic=true, negative=false,
            kB=1.9872041e-3, T=298.15
        )
    
    frames = readcoords(
                    PDBTools.readPDB(pdbname, selection),
                    trjname
                )

    sqrt_M, m = sqrt(massm(pdbname, selection=selection)), mass(pdbname, selection=selection)

    X = if mass_weighting
            dispm(frames, weights=m, method="sqrtmass")
        else
            dispm(frames)
    end
    
    C = if quasiharmonic
            sqrt_M * covarm(X) * sqrt_M
        else
            covarm(X)
    end

    λ, η = pca(C, negative=negative)

    β = inv(kB * T)
    return PCA(
            λ,
            η,
            explained_variance(λ),
            projections(η, X),
            quasiharmonic_frequency(λ, β),
            quasiharmonic_modes(η, sqrt_M)
        ) 
end


function pca(C::Matrix{Float64}; negative=false)
    eigenvalue, eigenvector = eigen(C)

    eigenvalue = if negative
            eigenvalue
        else
            clamp.(eigenvalue, 0, Inf)
    end

    idx = sortperm(eigenvalue, rev=true)
    return eigenvalue[idx], eigenvector[:, idx]
end

function massm(atoms::AbstractVector{<:PDBTools.Atom}; n=3)
    if n == 1
        return diagm(PDBTools.mass.(atoms))
    elseif n > 1
        return diagm(repeat(
                        PDBTools.mass.(atoms),
                        inner=n)
                    )
    else
        throw(ArgumentError("The number of dimensions must be greater than 0."))
    end
end

function massm(pdbname::String; selection="protein and name CA", n=3)
    return massm(PDBTools.readPDB(pdbname, selection), n=n)
end

function mass(pdbname::String; selection="protein and name CA")
    return PDBTools.mass.(PDBTools.readPDB(pdbname, selection))
end

function explained_variance(eigenvalues::Vector{Float64}; n=nothing)
    if isnothing(n)
        return eigenvalues ./ sum(eigenvalues)
    else
        return round.(eigenvalues ./ sum(eigenvalues), digits=n)
    end
end

function quasiharmonic_frequency(λ::Vector{Float64}, β::Float64, spectroscopic_unit=true)
    if spectroscopic_unit
        c = 2.99792458e10
        return sqrt.(inv(β) ./ λ) * inv(2 * π * c)
    else
        return sqrt.(inv(β) ./ λ)
    end
end

function quasiharmonic_modes(η::Matrix{Float64}, sqrt_M::Matrix{Float64})
    return inv(sqrt_M) * η
end

function projections(η::Matrix{Float64}, X::Matrix{Float64})
    return η' * X
end

function P_quasiharmonic(
                        X::Matrix{Float64},
                        C::Matrix{Float64};
                        kB=1.9872041e-3, T=298.15
                    )
    β, n = inv(kB * T), size(X, 1)
    E = E_quasiharmonic(X, C, kB=kB, T=T)
    A = sqrt(det(C)) / (2π)^(n/2)
    return A * exp.(-β*E)
end

function E_quasiharmonic(X::Matrix{Float64}, C::Matrix{Float64}; kB=1.9872041e-3, T=298.15)
    β = inv(kB * T)
    Feff = inv(β) * inv(C)
    return 0.5 * X' * Feff * X
end

function pmf_pca(
                pdbname::String,
                trjname::String;
                selection="protein and name CA",
                mass_weighting=true, quasiharmonic=true, negative=false,
                kB=1.9872041e-3, T=298.15
            )
    
    frames = readcoords(
                    PDBTools.readPDB(pdbname, selection),
                    trjname
                )

    sqrt_M, m = sqrt(massm(pdbname, selection=selection)), mass(pdbname, selection=selection)

    X = if mass_weighting
            dispm(frames, weights=m, method="sqrtmass")
        else
            dispm(frames)
    end
    
    C = if quasiharmonic
            sqrt_M * covarm(X) * sqrt_M
        else
            covarm(X)
    end

    prob = P_quasiharmonic(X, C, kB=kB, T=T)

    pmf = -log.(prob) / β

    for i in eachindex(pmf)
        if pmf[i] == Inf
            pmf[i] = NaN
        end
    end ## dealing with de Inf values from log(0), since it is indefinite

    pmf .-= minimum(pmf[.!isnan.(pmf)]) # normalizing the PMF by the minimum value

    return pmf
end