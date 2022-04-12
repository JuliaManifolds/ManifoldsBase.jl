function parallel_transport_along! end

@doc raw"""
    Y = parallel_transport_along(M::AbstractManifold, p, X, c)

Compute the parallel transport of the vector `X` from the tangent space at `p`
along the curve `c`.

To be precise let ``c(t)`` be a curve ``c(0)=p`` for [`vector_transport_along`](@ref) ``\mathcal P^cY``

THen In the result ``Y\in T_p\mathcal M`` is the vector ``X`` from the tangent space at ``p=c(0)``
to the tangent space at ``c(1)``.

Let ``Z\colon [0,1] \to T\mathcal M``, ``Z(t)\in T_{c(t)}\mathcal M`` be a smooth vector field
along the curve ``c`` with ``Z(0) = Y``, such that ``Z`` is _parallel_, i.e.
its covariant derivative ``\frac{\mathrm{D}}{\mathrm{d}t}Z`` is zero. Note that such a ``Z`` always exists and is unique.

Then the parallel transport is given by ``Z(1)``.
"""
function parallel_transport_along(M::AbstractManifold, p, X, c::AbstractVector)
    Y = allocate_result(M, vector_transport_along, X, p)
    return parallel_transport_along!(M, Y, p, X, c)
end

function parallel_transport_direction!(M::AbstractManifold, Y, p, X, d)
    return parallel_transport_to!(M, Y, p, X, exp(M, p, d))
end

@doc raw"""
    parallel_transport_direction(M::AbstractManifold, p, X, d)

Compute the [`parallel_transport_along`](@ref) the curve ``c(t) = γ_{p,q}(t)``,
i.e. the * the unique geodesic ``c(t)=γ_{p,X}(t)`` from ``γ_{p,d}(0)=p`` into direction ``\dot γ_{p,d}(0)=d``, of the tangent vector `X`.

By default this function calls [`parallel_transport_to`](@ref)`(M, p, X, q)`, where ``q=\exp_pX``.
"""
function parallel_transport_direction(M::AbstractManifold, p, X, d)
    return parallel_transport_to(M, p, X, exp(M, p, d))
end

function parallel_transport_to! end

@doc raw"""
    parallel_transport_to(M::AbstractManifold, p, X, q)

Compute the [`parallel_transport_along`](@ref) the curve ``c(t) = γ_{p,q}(t)``,
i.e. the (assumed to be unique) [`geodesic`](@ref) connecting `p` and `q`, of the tangent vector `X`.
"""
function parallel_transport_to(M::AbstractManifold, p, X, q)
    Y = allocate_result(M, vector_transport_to, X, p, q)
    return parallel_transport_to!(M, Y, p, X, q)
end
