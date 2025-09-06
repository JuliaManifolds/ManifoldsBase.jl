function parallel_transport_direction!(M::AbstractManifold, Y, p, X, d; kwargs...)
    return parallel_transport_to!(M, Y, p, X, exp(M, p, d); kwargs...)
end

@doc raw"""
    parallel_transport_direction(M::AbstractManifold, p, X, d)

Compute the parallel transport of ``X`` along the curve ``c(t) = γ_{p,X}(t)`` to ``c(1)=q``,
where ``c(t)=γ_{p,X}(t)`` is the the unique geodesic starting from ``γ_{p,d}(0)=p``
into direction ``̇\dot γ_{p,d}(0)=d``.

By default this function calls [`parallel_transport_to`](@ref)`(M, p, X, q)`, where ``q=\exp_pX``.
"""
function parallel_transport_direction(M::AbstractManifold, p, X, d; kwargs...)
    return parallel_transport_to(M, p, X, exp(M, p, d); kwargs...)
end

function parallel_transport_to! end

@doc raw"""
    parallel_transport_to(M::AbstractManifold, p, X, q)

Compute the parallel transport of ``X`` along the curve ``c(t) = γ_{p,q}(t)``,
i.e. the (assumed to be unique) [`geodesic`](@ref) ``γ_{p,q}`` connecting `p` and `q`.
"""
function parallel_transport_to(M::AbstractManifold, p, X, q; kwargs...)
    Y = allocate_result(M, vector_transport_to, X, p, q)
    return parallel_transport_to!(M, Y, p, X, q; kwargs...)
end
