# @inline function wrap(x::AbstractVector, xref::AbstractVector, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
#     invu = inv(oneunit(T))
#     unit_cell_matrix = invu * unit_cell_matrix
#     x_f = wrap_cell_fraction(invu * x, unit_cell_matrix)
#     xref_f = wrap_cell_fraction(invu * xref, unit_cell_matrix)
#     xw = wrap(x_f, xref_f, SVector{N,eltype(x_f)}(ntuple(i -> 1, N)))
#     return oneunit(T) * unit_cell_matrix * (xw - xref_f) + xref
# end

@inline function unwrap(x_w::AbstractVector, xref_u::AbstractVector, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    invu = inv(oneunit(T))
    x_w, xref_u = invu * [ i for i in x_w ], invu * [ i for i in xref_u ]
    unit_cell_matrix = invu * unit_cell_matrix
    xref_w = MolSimToolkit.wrap_to_first(
        xref_u,
        unit_cell_matrix
    )
    ucfloor = floor.( (unit_cell_matrix \ (x_w - xref_w)) .+  0.5 )
    return oneunit(T) * (unit_cell_matrix * ucfloor + xref_u + (x_w - xref_w))

end

# @inline function wrap_to_first(x::AbstractVector{T}, unit_cell_matrix) where {T<:AbstractFloat}
#     p = wrap_cell_fraction(x, unit_cell_matrix)
#     return typeof(x)(unit_cell_matrix * p)
# end

# @inline function wrap_cell_fraction(x::AbstractVector{T}, unit_cell_matrix::AbstractMatrix) where {T}
#     fastmod1(y) = y - floor(y)
#     x_stripped = x ./ oneunit(T)
#     m_stripped = unit_cell_matrix ./ oneunit(T)
#     p = fastmod1.(m_stripped \ x_stripped)
#     # Boundary coordinates belong to the lower boundary
#     p = ifelse.(p .== one(eltype(x_stripped)), zero(eltype(x_stripped)), p)
#     return p
# end

function displace(sim::MolSimToolkit.Simulation, inds; iswrapped=false)
    d = zeros(length(sim))
    MolSimToolkit.first_frame!(sim)
    p0 = MolSimToolkit.positions(MolSimToolkit.current_frame(sim))[inds]
    plast = copy(p0)
    for iframe in 2:length(sim)
        f = MolSimToolkit.next_frame!(sim)
        p = MolSimToolkit.positions(f)[inds]
        for i in 1:length(inds)
            pat = if iswrapped
                    unwrap(p[i], plast[i], MolSimToolkit.unitcell(f))
                else
                    MolSimToolkit.wrap(p[i], plast[i], MolSimToolkit.unitcell(f))
            end
            d[iframe] += sum(abs2, pat - p0[i])
        end
        d[iframe] /= length(inds)
        plast .= p
    end
    return d
end

function originaldisplace(sim::MolSimToolkit.Simulation, inds)
    d = zeros(length(sim))
    MolSimToolkit.first_frame!(sim)
    p0 = MolSimToolkit.positions(MolSimToolkit.current_frame(sim))[inds]
    plast = copy(p0)
    for iframe in 2:length(sim)
        f = MolSimToolkit.next_frame!(sim)
        p = MolSimToolkit.positions(f)[inds]
        for i in 1:length(inds)
            pat = p[i]
            d[iframe] += sum(abs2, pat - p0[i])
        end
        d[iframe] /= length(inds)
        plast .= p
    end
    return d
end


function unrwap(positions, positions_last, unitcell)
    for i in eachindex(positions, positions_last)
        positions[i] = wrap(positions[i], positions_last[i], unitcell)
    end
    return positions
end

# function msd(sim::Simulation, inds)
#     d = zeros(length(sim))
#     first_frame!(sim)
#     p0 = positions(current_frame(sim))[inds]
#     plast = copy(p0)
#     for iframe in 2:length(sim)
#         f = next_frame!(sim)
#         p = positions(f)[inds]
#         for i in 1:length(inds)
#             pat = wrap(p[i], plast[i], unitcell(f))
#             d[iframe] += sum(abs2, pat - p0[i])
#         end
#         d[iframe] /= length(inds)
#         plast .= p  
#     end
#     return d
# end

# @testitem "displacement" begin

#     dtest = zeros(length(sim1))
#     for (iframe, frame) in enumerate(sim1)
#         p = positions(current_frame(sim1))[22211]
#         dtest[iframe] = sum(abs2, p - p0)
#     end
    
# end