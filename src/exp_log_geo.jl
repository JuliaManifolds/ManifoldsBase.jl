@doc raw"""
    exp(M::Manifold, p, X)
    exp(M::Manifold, p, X, t::Real = 1)

Compute the exponential map of tangent vector `X`, optionally scaled by `t`,  at point `p`
from the manifold [`Manifold`](@ref) `M`, i.e.

```math
\exp_pX = γ_{p,X}(1),
```
where ``γ_{p,X}`` is the unique (shortest) geodesic starting in ``γ(0)=p`` and ``\dot γ(0) = X``.

See also [`shortest_geodesic`](@ref), [`retract`](@ref).
"""
function exp(M::Manifold, p, X)
    q = allocate_result(M, exp, p, X)
    exp!(M, q, p, X)
    return q
end
exp(M::Manifold, p, X, t::Real) = exp(M, p, t * X)

"""
    exp!(M::Manifold, q, p, X)
    exp!(M::Manifold, q, p, X, t::Real = 1)

Compute the exponential map of tangent vector `X`, optionally scaled by `t`,  at point `p`
from the manifold [`Manifold`](@ref) `M`.
The result is saved to `q`.

See also [`exp`](@ref).
"""
function exp!(M::Manifold, q, p, X)
    return error(manifold_function_not_implemented_message(M, exp!, q, p, X))
end
exp!(M::Manifold, q, p, X, t::Real) = exp!(M, q, p, t * X)

@doc raw"""
    geodesic(M::Manifold, p, X) -> Function

Get the geodesic with initial point `p` and velocity `X` on the [`Manifold`](@ref) `M`.
A geodesic is a curve of zero acceleration. That is for the curve ``γ_{p,X}: I → \mathcal M``,
with ``γ_{p,X}(0) = p`` and ``\dot γ_{p,X}(0) = X`` a geodesic further fulfills

```math
∇_{\dot γ_{p,X}(t)} \dot γ_{p,X}(t) = 0,
```

i.e. the curve is acceleration free with respect to the Riemannian metric.
This yields, that the curve has constant velocity that is locally distance-minimizing.

This function returns a function of (time) `t`.

    geodesic(M::Manifold, p, X, t::Real)
    geodesic(M::Manifold, p, X, T::AbstractVector) -> AbstractVector

Return the point at time `t` or points at times `t` in `T` along the geodesic.
"""
geodesic(M::Manifold, p, X) = t -> exp(M, p, X, t)
geodesic(M::Manifold, p, X, t::Real) = exp(M, p, X, t)
geodesic(M::Manifold, p, X, T::AbstractVector) = map(t -> exp(M, p, X, t), T)

"""
    log(M::Manifold, p, q)

Compute the logarithmic map of point `q` at base point `p` on the [`Manifold`](@ref) `M`.
The logarithmic map is the inverse of the [`exp`](@ref)onential map.
Note that the logarithmic map might not be globally defined.

See also [`inverse_retract`](@ref).
"""
function log(M::Manifold, p, q)
    X = allocate_result(M, log, p, q)
    log!(M, X, p, q)
    return X
end

"""
    log!(M::Manifold, X, p, q)

Compute the logarithmic map of point `q` at base point `p` on the [`Manifold`](@ref) `M`.
The result is saved to `X`.
The logarithmic map is the inverse of the [`exp!`](@ref)onential map.
Note that the logarithmic map might not be globally defined.
"""
function log!(M::Manifold, X, p, q)
    return error(manifold_function_not_implemented_message(M, log!, X, p, q))
end

@doc doc"""
    shortest_geodesic(M::Manifold, p, q) -> Function

Get a [`geodesic`](@ref) $γ_{p,q}(t)$ whose length is the shortest path between the
points `p`and `q`, where $γ_{p,q}(0)=p$ and $γ_{p,q}(1)=q$. When there are
multiple shortest geodesics, a deterministic choice will be returned.

This function returns a function of time, which may be a `Real` or an `AbstractVector`.

    shortest_geodesic(M::Manifold, p, q, t::Real)
    shortest_geodesic(M::Manifold, p, q, T::AbstractVector) -> AbstractVector

Return the point at time `t` or points at times `t` in `T` along the shortest [`geodesic`](@ref).
"""
shortest_geodesic(M::Manifold, p, q) = geodesic(M, p, log(M, p, q))
shortest_geodesic(M::Manifold, p, q, t::Real) = geodesic(M, p, log(M, p, q), t)
shortest_geodesic(M::Manifold, p, q, T::AbstractVector) = geodesic(M, p, log(M, p, q), T)
