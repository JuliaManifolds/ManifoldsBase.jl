@doc raw"""
    AbstractMetric

Abstract type for the pseudo-Riemannian metric tensor ``g``, a family of smoothly
varying inner products on the tangent space. See [`inner`](@ref).

# Functor

    (metric::Metric)(M::AbstractManifold)
    (metric::Metric)(M::MetricManifold)

Generate the `MetricManifold` that wraps the manifold `M` with given `metric`.
This works for both a variable containing the metric as well as a subtype `T<:AbstractMetric`,
where a zero parameter constructor `T()` is availabe.
If `M` is already a metric manifold, the inner manifold with the new `metric` is returned.
"""
abstract type AbstractMetric end

@doc raw"""
    RiemannianMetric <: AbstractMetric

Abstract type for Riemannian metrics, a family of positive definite inner
products. The positive definite property means that for ``X  ∈ T_p \mathcal M``, the
inner product ``g(X, X) > 0`` whenever ``X`` is not the zero vector.
"""
abstract type RiemannianMetric <: AbstractMetric end

"""
    EuclideanMetric <: RiemannianMetric

A general type for any manifold that employs the Euclidean Metric, for example
the [`Euclidean`](@ref) manifold itself, or the [`Sphere`](@ref), where every
tangent space (as a plane in the embedding) uses this metric (in the embedding).

Since the metric is independent of the field type, this metric is also used for
the Hermitian metrics, i.e. metrics that are analogous to the `EuclideanMetric`
but where the field type of the manifold is `ℂ`.

This metric is the default metric for example for the [`Euclidean`](@ref) manifold.
"""
struct EuclideanMetric <: RiemannianMetric end


@doc raw"""
    change_metric(M::AbstractcManifold, G2::AbstractMetric, p, X)

On the [`AbstractManifold`](@ref) `M` with implicitly given metric ``g_1``
and a second [`AbstractMetric`](@ref)
``g_2`` this function performs a change of metric in the
sense that it returns the tangent vector ``Z=BX`` such that the linear map ``B`` fulfills

```math
g_2(Y_1,Y_2) = g_1(BY_1,BY_2) \quad \text{for all } Y_1, Y_2 ∈ T_p\mathcal M.
```
"""
function change_metric(M::AbstractManifold, G::AbstractMetric, p, X)
    Y = allocate_result(M, change_metric, X, p) # this way we allocate a tangent
    return change_metric!(M, Y, G, p, X)
end

@doc raw"""
    change_metric!(M::AbstractcManifold, Y, G2::AbstractMetric, p, X)

Compute the [`change_metric`](@ref) in place of `Y`.
"""
change_metric!(M::AbstractManifold, Y, G::AbstractMetric, p, X)

@doc raw"""
    change_representer(M::AbstractManifold, G2::AbstractMetric, p, X)

Convert the representer `X` of a linear function (in other words a cotangent vector at `p`)
in the tangent space at `p` on the [`AbstractManifold`](@ref) `M` given with respect to the
[`AbstractMetric`](@ref) `G2` into the representer with respect to the (implicit) metric of `M`.

In order to convert `X` into the representer with respect to the (implicitly given) metric ``g_1`` of `M`,
we have to find the conversion function ``c: T_p\mathcal M \to T_p\mathcal M`` such that

```math
    g_2(X,Y) = g_1(c(X),Y)
```
"""
function change_representer(M::AbstractManifold, G::AbstractMetric, p, X)
    Y = allocate_result(M, change_representer, X, p) # this way we allocate a tangent
    return change_representer!(M, Y, G, p, X)
end

@doc raw"""
    change_representer!(M::AbstractcManifold, Y, G2::AbstractMetric, p, X)

Compute the [`change_metric`](@ref) in place of `Y`.
"""
change_representer!(M::AbstractManifold, Y, G::AbstractMetric, p, X)
