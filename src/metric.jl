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
