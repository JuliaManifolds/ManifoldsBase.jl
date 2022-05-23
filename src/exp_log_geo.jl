@doc raw"""
    exp(M::AbstractManifold, p, X)
    exp(M::AbstractManifold, p, X, t::Real = 1)

Compute the exponential map of tangent vector `X`, optionally scaled by `t`,  at point `p`
from the manifold [`AbstractManifold`](@ref) `M`, i.e.

```math
\exp_p X = γ_{p,X}(1),
```
where ``γ_{p,X}`` is the unique geodesic starting in ``γ(0)=p`` such that ``\dot γ(0) = X``.

See also [`shortest_geodesic`](@ref), [`retract`](@ref).
"""
function exp(M::AbstractManifold, p, X)
    q = allocate_result(M, exp, p, X)
    exp!(M, q, p, X)
    return q
end
exp(M::AbstractManifold, p, X, t::Real) = exp(M, p, t * X)

"""
    exp!(M::AbstractManifold, q, p, X)
    exp!(M::AbstractManifold, q, p, X, t::Real = 1)

Compute the exponential map of tangent vector `X`, optionally scaled by `t`,  at point `p`
from the manifold [`AbstractManifold`](@ref) `M`.
The result is saved to `q`.

See also [`exp`](@ref).
"""
exp!(M::AbstractManifold, q, p, X)
exp!(M::AbstractManifold, q, p, X, t::Real) = exp!(M, q, p, t * X)

@doc raw"""
    geodesic(M::AbstractManifold, p, X) -> Function

Get the geodesic with initial point `p` and velocity `X` on the [`AbstractManifold`](@ref) `M`.
A geodesic is a curve of zero acceleration. That is for the curve ``γ_{p,X}: I → \mathcal M``,
with ``γ_{p,X}(0) = p`` and ``\dot γ_{p,X}(0) = X`` a geodesic further fulfills

```math
∇_{\dot γ_{p,X}(t)} \dot γ_{p,X}(t) = 0,
```

i.e. the curve is acceleration free with respect to the Riemannian metric.
This yields, that the curve has constant velocity that is locally distance-minimizing.

This function returns a function of (time) `t`.

    geodesic(M::AbstractManifold, p, X, t::Real)
    geodesic(M::AbstractManifold, p, X, T::AbstractVector) -> AbstractVector

Return the point at time `t` or points at times `t` in `T` along the geodesic.
"""
geodesic(M::AbstractManifold, p, X) = t -> exp(M, p, X, t)
geodesic(M::AbstractManifold, p, X, t::Real) = exp(M, p, X, t)
geodesic(M::AbstractManifold, p, X, T::AbstractVector) = map(t -> exp(M, p, X, t), T)

@doc raw"""
    geodesic!(M::AbstractManifold, p, X) -> Function

Get the geodesic with initial point `p` and velocity `X` on the [`AbstractManifold`](@ref) `M`.
A geodesic is a curve of zero acceleration. That is for the curve ``γ_{p,X}: I → \mathcal M``,
with ``γ_{p,X}(0) = p`` and ``\dot γ_{p,X}(0) = X`` a geodesic further fulfills

```math
∇_{\dot γ_{p,X}(t)} \dot γ_{p,X}(t) = 0,
```

i.e. the curve is acceleration free with respect to the Riemannian metric.
This yields, that the curve has constant velocity that is locally distance-minimizing.

This function returns a function `(q,t)` of (time) `t` that mutates `q``.

    geodesic!(M::AbstractManifold, q, p, X, t::Real)
    geodesic!(M::AbstractManifold, Q, p, X, T::AbstractVector) -> AbstractVector

Return the point at time `t` or points at times `t` in `T` along the geodesic and mutate `q`and `Q`, respectively.
"""
geodesic!(M::AbstractManifold, p, X) = (q, t) -> exp!(M, q, p, X, t)
geodesic!(M::AbstractManifold, q, p, X, t::Real) = exp!(M, q, p, X, t)
function geodesic!(M::AbstractManifold, Q, p, X, T::AbstractVector)
    for (q, t) in zip(Q, T)
        exp!(M, q, p, X, t)
    end
    return Q
end

"""
    log(M::AbstractManifold, p, q)

Compute the logarithmic map of point `q` at base point `p` on the [`AbstractManifold`](@ref) `M`.
The logarithmic map is the inverse of the [`exp`](@ref)onential map.
Note that the logarithmic map might not be globally defined.

See also [`inverse_retract`](@ref).
"""
function log(M::AbstractManifold, p, q)
    X = allocate_result(M, log, p, q)
    log!(M, X, p, q)
    return X
end

"""
    log!(M::AbstractManifold, X, p, q)

Compute the logarithmic map of point `q` at base point `p` on the [`AbstractManifold`](@ref) `M`.
The result is saved to `X`.
The logarithmic map is the inverse of the [`exp!`](@ref)onential map.
Note that the logarithmic map might not be globally defined.

see also [`log`](@ref) and [`inverse_retract!`](@ref),
"""
log!(M::AbstractManifold, X, p, q)

@doc raw"""
    shortest_geodesic(M::AbstractManifold, p, q) -> Function

Get a [`geodesic`](@ref) $γ_{p,q}(t)$ whose length is the shortest path between the
points `p`and `q`, where $γ_{p,q}(0)=p$ and $γ_{p,q}(1)=q$. When there are
multiple shortest geodesics, a deterministic choice will be returned.

This function returns a function of time, which may be a `Real` or an `AbstractVector`.

    shortest_geodesic(M::AabstractManifold, p, q, t::Real)
    shortest_geodesic(M::AbstractManifold, p, q, T::AbstractVector) -> AbstractVector

Return the point at time `t` or points at times `t` in `T` along the shortest [`geodesic`](@ref).
"""
shortest_geodesic(M::AbstractManifold, p, q) = geodesic(M, p, log(M, p, q))
shortest_geodesic(M::AbstractManifold, p, q, t::Real) = geodesic(M, p, log(M, p, q), t)
function shortest_geodesic(M::AbstractManifold, p, q, T::AbstractVector)
    return geodesic(M, p, log(M, p, q), T)
end


@doc raw"""
    shortest_geodesic!(M::AbstractManifold, p, q) -> Function

Get a [`geodesic`](@ref) $γ_{p,q}(t)$ whose length is the shortest path between the
points `p`and `q`, where $γ_{p,q}(0)=p$ and $γ_{p,q}(1)=q$. When there are
multiple shortest geodesics, a deterministic choice will be returned.

This function returns a function `(r,t) -> ... ` of time `t` which mutates `r`.

Further variants

    shortest_geodesic!(M::AabstractManifold, r, p, q, t::Real)
    shortest_geodesic!(M::AbstractManifold, R, p, q, T::AbstractVector) -> AbstractVector

mutate (and return) the point `r` and the vector of points `R`, respectively,
returning the point at time `t` or points at times `t` in `T` along the shortest [`geodesic`](@ref).
"""
shortest_geodesic!(M::AbstractManifold, p, q) = geodesic!(M, p, log(M, p, q))
function shortest_geodesic!(M::AbstractManifold, r, p, q, t::Real)
    return geodesic!(M, r, p, log(M, p, q), t)
end
function shortest_geodesic!(M::AbstractManifold, R, p, q, T::AbstractVector)
    return geodesic!(M, R, p, log(M, p, q), T)
end
