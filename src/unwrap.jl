function wrap end

# @inline function wrap(x::AbstractVector, xref::AbstractVector, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
#     invu = inv(oneunit(T))
#     unit_cell_matrix = invu * unit_cell_matrix
#     x_f = wrap_cell_fraction(invu * x, unit_cell_matrix)
#     xref_f = wrap_cell_fraction(invu * xref, unit_cell_matrix)
#     xw = wrap(x_f, xref_f, SVector{N,eltype(x_f)}(ntuple(i -> 1, N)))
#     return oneunit(T) * unit_cell_matrix * (xw - xref_f) + xref
# end

@inline function wrap(x::AbstractVector, xref::AbstractVector, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    invu = inv(oneunit(T))
    unit_cell_matrix = invu * unit_cell_matrix
    x_f = invu * x
    xref_f = invu * xref
    xw = wrap(
        wrap_cell_fraction(x_f, unit_cell_matrix),
        wrap_cell_fraction(xref_f, unit_cell_matrix),
        SVector{N,eltype(x_f)}(ntuple(i -> 1, N))
    )
    return oneunit(T) * unit_cell_matrix * (xw - xref_f) + xref + (x_f - xref_f)
end

function wrap_cell_fraction(x::AbstractVector{T}, unit_cell_matrix::AbstractMatrix) where {T}
    # Division by `oneunit` is to support Unitful quantities. 
    # this workaround works here because the units cancel.
    # see: https://github.com/PainterQubits/Unitful.jl/issues/46
    x_stripped = x ./ oneunit(T)
    m_stripped = unit_cell_matrix ./ oneunit(T)
    p = fastmod1.(m_stripped \ x_stripped)
    # Boundary coordinates belong to the lower boundary
    p = ifelse.(p .== one(eltype(x_stripped)), zero(eltype(x_stripped)), p)
    return p
end

@inline fastmod1(x) = x - floor(x)

function unrwap(positions, positions_last, unitcell)
    for i in eachindex(positions, positions_last)
        positions[i] = wrap(positions[i], positions_last[i], unitcell)
    end
    return positions
end

function msd(sim::Simulation, inds)
    d = zeros(length(sim))
    first_frame!(sim)
    p0 = positions(current_frame(sim))[inds]
    plast = copy(p0)
    for iframe in 2:length(sim)
        f = next_frame!(sim)
        p = positions(f)[inds]
        for i in 1:length(inds)
            pat = wrap(p[i], plast[i], unitcell(f))
            d[iframe] += sum(abs2, pat - p0[i])
        end
        d[iframe] /= length(inds)
        plast .= p  
    end
    return d
end