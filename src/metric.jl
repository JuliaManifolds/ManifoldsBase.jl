@doc raw"""
    AbstractMetric

Abstract type for the pseudo-Riemannian metric tensor ``g``, a family of smoothly
varying inner products on the tangent space. See [`inner`](@ref).

# Functor

    (metric::Metric)(M::AbstractManifold)
    (metric::Metric)(M::MetricManifold)

Generate the `MetricManifold` that wraps the manifold `M` with given `metric`.
This works for both a variable containing the metric as well as a subtype `T<:AbstractMetric`,
where a zero parameter constructor `T()` is available.
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
the [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) manifold itself, or the [`Sphere`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html), where every
tangent space (as a plane in the embedding) uses this metric (in the embedding).

Since the metric is independent of the field type, this metric is also used for
the Hermitian metrics, i.e. metrics that are analogous to the `EuclideanMetric`
but where the field type of the manifold is `ℂ`.

This metric is the default metric for example for the [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) manifold.
"""
struct EuclideanMetric <: RiemannianMetric end


@doc raw"""
    change_metric(M::AbstractManifold, G2::AbstractMetric, p, X)

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

function change_metric! end

@doc raw"""
    change_metric!(M::AbstractManifold, Y, G2::AbstractMetric, p, X)

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

function change_representer! end

@doc raw"""
    change_representer!(M::AbstractManifold, Y, G2::AbstractMetric, p, X)

Compute the [`change_metric`](@ref) in place of `Y`.
"""
change_representer!(M::AbstractManifold, Y, G::AbstractMetric, p, X)

# piping syntax for decoration
(metric::AbstractMetric)(M::AbstractManifold) = MetricManifold(M, metric)
(::Type{T})(M::AbstractManifold) where {T <: AbstractMetric} = MetricManifold(M, T())

"""
    DefaultMetric <: AbstractMetric

Indicating that a manifold uses the default metric, that one has implicitly assumed
when defining the manifold
"""
struct DefaultMetric <: AbstractMetric end
metric(::AbstractManifold) = DefaultMetric()

"""
    MetricManifold{𝔽,M<:AbstractManifold{𝔽},G<:AbstractMetric} <: AbstractDecoratorManifold{𝔽}

Equip a [`AbstractManifold`](@ref) explicitly with an
[`AbstractMetric`](@ref) `G`.

If the `AbstractMetric` `G` yields closed form formulae for the exponential map or another
function, you can implement it directly. Otherwise, you can use chart-based solvers, see
for example [`solve_chart_exp_ode`](@extref `Manifolds.solve_chart_exp_ode`).

# Constructor

    MetricManifold(M, G)

Generate the [`AbstractManifold`](@ref) `M` as a manifold with the `AbstractMetric` `G`.
"""
struct MetricManifold{𝔽, M <: AbstractManifold{𝔽}, G <: AbstractMetric} <:
    AbstractDecoratorManifold{𝔽}
    manifold::M
    metric::G
end

# remetricise instead of double-decorating
(metric::AbstractMetric)(M::MetricManifold) = MetricManifold(M.manifold, metric)
(::Type{T})(M::MetricManifold) where {T <: AbstractMetric} = MetricManifold(M.manifold, T())

function Base.convert(::Type{MetricManifold{𝔽, MT, GT}}, M::MT) where {𝔽, MT, GT}
    return _convert_with_default(M, GT, Val(is_default_metric(M, GT())))
end

function _convert_with_default(
        M::MT,
        T::Type{<:AbstractMetric},
        ::Val{true},
    ) where {MT <: AbstractManifold}
    return MetricManifold(M, T())
end
function _convert_with_default(
        M::MT,
        T::Type{<:AbstractMetric},
        ::Val{false},
    ) where {MT <: AbstractManifold}
    return error(
        "Can not convert $(M) to a MetricManifold{$(MT),$(T)}, since $(T) is not the default metric.",
    )
end

decorated_manifold(M::MetricManifold) = M.manifold

default_retraction_method(M::MetricManifold) = default_retraction_method(M.manifold)
function default_retraction_method(M::MetricManifold, t::Type)
    return default_retraction_method(M.manifold, t)
end
function default_inverse_retraction_method(M::MetricManifold)
    return default_inverse_retraction_method(M.manifold)
end
function default_inverse_retraction_method(M::MetricManifold, t::Type)
    return default_inverse_retraction_method(M.manifold, t)
end
function default_vector_transport_method(M::MetricManifold)
    return default_vector_transport_method(M.manifold)
end
function default_vector_transport_method(M::MetricManifold, t::Type)
    return default_vector_transport_method(M.manifold, t)
end

get_embedding(M::MetricManifold) = get_embedding(M.manifold)
get_embedding(M::MetricManifold, T::Type) = get_embedding(M.manifold, T)


function get_basis(M::MetricManifold, p, B::AbstractBasis)
    (metric(M.manifold) === M.metric) && (return get_basis(M.manifold, p, B))
    return invoke(get_basis, Tuple{AbstractManifold, Any, AbstractBasis}, M, p, B)
end

function get_coordinates(M::MetricManifold, p, X, B::AbstractBasis)
    (metric(M.manifold) === M.metric) && (return get_coordinates(M.manifold, p, X, B))
    return invoke(
        get_coordinates,
        Tuple{AbstractManifold, Any, Any, AbstractBasis},
        M,
        p,
        X,
        B,
    )
end
function get_coordinates!(M::MetricManifold, Y, p, X, B::AbstractBasis)
    (metric(M.manifold) === M.metric) && (return get_coordinates!(M.manifold, Y, p, X, B))
    return invoke(
        get_coordinates!,
        Tuple{AbstractManifold, Any, Any, Any, AbstractBasis},
        M,
        Y,
        p,
        X,
        B,
    )
end

function get_forwarding_type(::MetricManifold, ::typeof(embed), P::Type)
    return SimpleForwardingType()
end
function get_forwarding_type(::MetricManifold, ::typeof(embed!), P::Type)
    return SimpleForwardingType()
end
function get_forwarding_type(::MetricManifold, ::typeof(rand))
    return SimpleForwardingType()
end
function get_forwarding_type(::MetricManifold, ::typeof(rand), P::Type)
    return SimpleForwardingType()
end
function get_forwarding_type(::MetricManifold, ::typeof(rand!))
    return SimpleForwardingType()
end
function get_forwarding_type(::MetricManifold, ::typeof(rand!), P::Type)
    return SimpleForwardingType()
end
function get_forwarding_type(M::MetricManifold, f::typeof(default_approximation_method))
    is_default_metric(M) && (return SimpleForwardingType())
    return invoke(get_forwarding_type, Tuple{AbstractManifold, typeof(f)}, M, f)
end

function get_vector(M::MetricManifold, p, c, B::AbstractBasis)
    (metric(M.manifold) === M.metric) && (return get_vector(M.manifold, p, c, B))
    return invoke(get_vector, Tuple{AbstractManifold, Any, Any, AbstractBasis}, M, p, c, B)
end
function get_vector!(M::MetricManifold, Y, p, c, B::AbstractBasis)
    (metric(M.manifold) === M.metric) && (return get_vector!(M.manifold, Y, p, c, B))
    return invoke(
        get_vector!,
        Tuple{AbstractManifold, Any, Any, Any, AbstractBasis},
        M,
        Y,
        p,
        c,
        B,
    )
end

function exp(M::MetricManifold, p, X)
    (metric(M.manifold) === M.metric) && (return exp(M.manifold, p, X))
    return invoke(exp, Tuple{AbstractManifold, Any, Any}, M, p, X)
end
function ManifoldsBase.exp_fused(M::MetricManifold, p, X, t::Number)
    (metric(M.manifold) === M.metric) && (return exp_fused(M.manifold, p, X, t))
    return invoke(exp_fused, Tuple{AbstractManifold, Any, Any, Number}, M, p, X, t)
end
function exp!(M::MetricManifold, q, p, X)
    (metric(M.manifold) === M.metric) && (return exp!(M.manifold, q, p, X))
    throw(MethodError(exp!, (M, q, p, X)))
end
function ManifoldsBase.exp_fused!(M::MetricManifold, q, p, X, t::Number)
    (metric(M.manifold) === M.metric) && (return exp_fused!(M.manifold, q, p, X, t))
    return invoke(exp_fused!, Tuple{AbstractManifold, Any, Any, Any, Number}, M, q, p, X, t)
end

injectivity_radius(M::MetricManifold) = injectivity_radius(M.manifold)
function injectivity_radius(M::MetricManifold, m::AbstractRetractionMethod)
    return injectivity_radius(M.manifold, m)
end

@doc raw"""
    inner(N::MetricManifold{M,G}, p, X, Y)

Compute the inner product of `X` and `Y` from the tangent space at `p` on the
[`AbstractManifold`](@ref) `M` using the [`AbstractMetric`](@ref) `G`.

````math
g_p(X, Y) = ⟨X, G_p Y⟩,
````
where ``G_p`` is the local matrix representation of the `AbstractMetric` `G`.
"""
inner(::MetricManifold, ::Any, ::Any, ::Any)

function inner(
        M::MetricManifold{𝔽, TM, G},
        p,
        X,
        Y,
    ) where {𝔽, G <: AbstractMetric, TM <: AbstractManifold}
    (metric(M.manifold) === M.metric) && (return inner(M.manifold, p, X, Y))
    throw(MethodError(inner, (M, p, X, Y)))
end

# For the inverse retraction, we distinguish two cases: if we have the LogarithmicInverseRetraction
# we do not pass to the inner manifold, since this falls to log by default,
# otherwise we do pass to the inner manifold

@doc raw"""
    inverse_retract(M::MetricManifold, p, q)
    inverse_retract!(M::MetricManifold, X, p, q)

Compute the inverse retraction on the [`MetricManifold`](@ref) `M`.
Since every inverse retraction is an inverse retraction with respect to any logarithmic map (induced by the metric),
this method falls back to calling [`inverse_retract`](@ref) on the base manifold.
The two exceptions are the [`LogarithmicInverseRetraction`](@ref) and [`ShootingInverseRetraction`](@ref),
in which case the method falls back to the default, that is to calling, respectively, [`log(::AbstractManifold, ::Any, ::Any)`](@ref) and
`inverse_retract_shooting!`.
"""
inverse_retract(::MetricManifold, ::Any, ::Any)

function inverse_retract(M::MetricManifold, p, q, m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)))
    (metric(M.manifold) === M.metric) && (return inverse_retract(M.manifold, p, q, m))
    return invoke(inverse_retract, Tuple{AbstractManifold, Any, Any, AbstractInverseRetractionMethod}, M, p, q, m)
end

# note that if the default inverse retraction is the Logarithmic or the shooting one, this indeed still dispatches correctly to the next case
function inverse_retract!(M::MetricManifold, X, p, q, m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)))
    return inverse_retract!(M.manifold, X, p, q, m)
end
function inverse_retract!(M::MetricManifold, X, p, q, ::LogarithmicInverseRetraction)
    (metric(M.manifold) === M.metric) && (return log!(M.manifold, X, p, q))
    return log!(M, X, p, q)
end

"""
    is_default_metric(M::AbstractManifold, G::AbstractMetric)

Return whether an [`AbstractMetric`](@ref)
is the default metric on the manifold `M` or not.

If `M` is a |`MetricManifold`](@ref) this indicates whether the metric now used is the same as the
default one on the wrapped manifold.
"""
is_default_metric(M::AbstractManifold, G::AbstractMetric)

is_default_metric(M::MetricManifold) = metric(M.manifold) === M.metric
is_default_metric(M::AbstractManifold, G::AbstractMetric) = metric(M) == G

function is_flat(M::MetricManifold{𝔽, TM, G}) where {𝔽, G <: AbstractMetric, TM <: AbstractManifold}
    is_default_metric(M) && (return is_flat(M.manifold))
    return invoke(is_flat, Tuple{AbstractManifold}, M)
end

is_point(M::MetricManifold, p; kwargs...) = is_point(M.manifold, p; kwargs...)

function is_vector(M::MetricManifold, p, X, cbp::Bool = true; kwargs...)
    return is_vector(M.manifold, p, X, cbp; kwargs...)
end

@doc raw"""
    log(N::MetricManifold{M,G}, p, q)

Compute the logarithmic map on the [`AbstractManifold`](@ref) `M` equipped with the
[`AbstractMetric`](@ref) `G`.

If the metric was declared the default metric, this method falls back to `log(M, p, q)`.
Otherwise, you have to provide an implementation for the non-default `AbstractMetric` `G`
metric within its [`MetricManifold`](@ref)`{M,G}`.
"""
log(::MetricManifold, ::Any...)

function log(M::MetricManifold, p, q)
    (metric(M.manifold) === M.metric) && (return log(M.manifold, p, q))
    return invoke(log, Tuple{AbstractManifold, Any, Any}, M, p, q)
end
function log!(M::MetricManifold, X, p, q)
    (metric(M.manifold) === M.metric) && (return log!(M.manifold, X, p, q))
    throw(MethodError(log!, (M, X, p, q)))
end

manifold_dimension(M::MetricManifold) = manifold_dimension(M.manifold)

@doc raw"""
    metric(M::MetricManifold)

Get the metric ``g`` of the [`AbstractManifold`](@ref)`(M)`.
"""
metric(::AbstractManifold)

function metric(M::MetricManifold)
    return M.metric
end
function parallel_transport_to(M::MetricManifold, p, X, q)
    (metric(M.manifold) === M.metric) && (return parallel_transport_to(M.manifold, p, X, q))
    return invoke(parallel_transport_to, Tuple{AbstractManifold, Any, Any, Any}, M, p, X, q)
end
function parallel_transport_to!(M::MetricManifold, Y, p, X, q)
    (metric(M.manifold) === M.metric) &&
        (return parallel_transport_to!(M.manifold, Y, p, X, q))
    throw(MethodError(parallel_transport_to!, (M, Y, p, X, q)))
end

function project(M::MetricManifold, p)
    (metric(M.manifold) === M.metric) && (return project(M.manifold, p))
    return invoke(project, Tuple{AbstractManifold, Any}, M, p)
end
function project!(M::MetricManifold, q, p)
    (metric(M.manifold) === M.metric) && (return project!(M.manifold, q, p))
    return project!(M.manifold, q, p)
end
function project(M::MetricManifold, p, X)
    (metric(M.manifold) === M.metric) && (return project(M.manifold, p, X))
    return invoke(project, Tuple{AbstractManifold, Any, Any}, M, p, X)
end
function project!(M::MetricManifold, Y, p, X)
    (metric(M.manifold) === M.metric) && (return project!(M.manifold, Y, p, X))
    return project!(M.manifold, Y, p, X)
end

representation_size(M::MetricManifold) = representation_size(M.manifold)

@doc raw"""
    retract(M::MetricManifold, p, X)
    retract!(M::MetricManifold, q, p, X)

Compute the retraction on the [`MetricManifold`](@ref) `M`.
Since every retraction is a retraction with respect to any exponential map (here induced by the metric),
this method falls back to calling [`retract`](@ref) on the inner manifold.
The one exception is the [`ExponentialRetraction`](@ref), in which case the method falls back to
the default, i.e. to calling [`exp(::AbstractManifold, ::Any, ::Any)`](@ref) but still on `M`.
"""
retract(::MetricManifold, ::Any, ::Any)

function retract(M::MetricManifold, p, X, m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)); kwargs...)
    (metric(M.manifold) === M.metric) && (return retract(M.manifold, p, X, m; kwargs...))
    return invoke(retract, Tuple{AbstractManifold, Any, Any, AbstractRetractionMethod}, M, p, X, m; kwargs...)
end
function retract_fused(M::MetricManifold, p, X, t::Number, m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)))
    (metric(M.manifold) === M.metric) && (return retract_fused(M.manifold, p, X, t, m))
    return invoke(retract_fused, Tuple{AbstractManifold, Any, Any, Number, AbstractRetractionMethod}, M, p, X, t, m)
end

# note that if the default retraction is the Exponential, this indeed still dispatches correctly to the next case
function retract!(M::MetricManifold, q, p, X, m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)); kwargs...)
    return retract!(M.manifold, q, p, X, m; kwargs...)
end
function retract!(M::MetricManifold, q, p, X, ::ExponentialRetraction)
    (metric(M.manifold) === M.metric) && (return exp!(M.manifold, q, p, X))
    return exp!(M, q, p, X)
end
# note that if the default retraction is the Exponential, this indeed still dispatches correctly to the next case
function retract_fused!(M::MetricManifold, q, p, X, t::Number, m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)))
    return retract_fused!(M.manifold, q, p, X, t, m)
end
function retract_fused!(M::MetricManifold, q, p, X, t::Number, ::ExponentialRetraction)
    (metric(M.manifold) === M.metric) && (return exp_fused!(M.manifold, q, p, X, t))
    return exp_fused!(M, q, p, X, t)
end

function Base.show(io::IO, M::MetricManifold)
    return print(io, "MetricManifold($(M.manifold), $(M.metric))")
end

@doc raw"""
    vector_transport_direction(M::MetricManifold, p, X, d)
    vector_transport_direction!(M::MetricManifold, Y, p, X, d)

Compute the vector transport of the tangent vector `X` at point `p` in the direction `d`
on the [`MetricManifold`](@ref) `M`.

Since a vector transport is usually defined with respect to a retraction, cf. e.g. [AbsilMahonySepulchre:2008](@cite),
and the vector transport is closely related to an affine connection, it is to some extent metric dependent.
Therefore, this method only falls back to calling its corresponding method on the base manifold, if the metric is the default one.
"""
vector_transport_direction(::MetricManifold, ::Any, ::Any, ::Any)

function vector_transport_direction(
        M::MetricManifold, p, X, d,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    (metric(M.manifold) === M.metric) && (return vector_transport_direction(M.manifold, p, X, d, m))
    return invoke(
        vector_transport_direction,
        Tuple{AbstractManifold, Any, Any, Any, AbstractVectorTransportMethod},
        M, p, X, d, m,
    )
end
function vector_transport_direction!(
        M::MetricManifold, Y, p, X, d,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    (metric(M.manifold) === M.metric) && (return vector_transport_direction!(M.manifold, Y, p, X, d, m))
    return invoke(
        vector_transport_direction!,
        Tuple{AbstractManifold, Any, Any, Any, Any, AbstractVectorTransportMethod},
        M, Y, p, X, d, m,
    )
end

@doc raw"""
    vector_transport_to(M::MetricManifold, p, X, d)
    vector_transport_to!(M::MetricManifold, Y, p, X, d)

Compute the vector transport of the tangent vector `X` at point `p` to a point `q` on the [`MetricManifold`](@ref) `M`.

Since a vector transport is usually defined with respect to a retraction, cf. e.g. [AbsilMahonySepulchre:2008](@cite),
and the vector transport is closely related to an affine connection, it is to some extent metric dependent.
Therefore, this method only falls back to calling its corresponding method on the base manifold, if the metric is the default one.
"""
vector_transport_to(::MetricManifold, ::Any, ::Any, ::Any)

function vector_transport_to(
        M::MetricManifold,
        p,
        X,
        q,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    (metric(M.manifold) === M.metric) && (return vector_transport_to(M.manifold, p, X, q, m))
    return invoke(
        vector_transport_to,
        Tuple{AbstractManifold, Any, Any, Any, AbstractVectorTransportMethod},
        M, p, X, q, m,
    )
end
function vector_transport_to!(
        M::MetricManifold, Y, p, X, q,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    (metric(M.manifold) === M.metric) && (return vector_transport_to!(M.manifold, Y, p, X, q, m))
    return invoke(
        vector_transport_to!,
        Tuple{AbstractManifold, Any, Any, Any, Any, AbstractVectorTransportMethod},
        M, Y, p, X, q, m,
    )
end


function Weingarten(M::MetricManifold, p, X, V)
    (metric(M.manifold) === M.metric) && (return Weingarten(M.manifold, p, X, V))
    return invoke(Weingarten, Tuple{AbstractManifold, Any, Any, Any}, M, p, X, V)
end
function Weingarten!(M::MetricManifold, Y, p, X, V)
    (metric(M.manifold) === M.metric) && (return Weingarten!(M.manifold, Y, p, X, V))
    throw(MethodError(Weingarten!, (M, Y, p, X, V)))
end

zero_vector(M::MetricManifold, p) = zero_vector(M.manifold, p)
zero_vector!(M::MetricManifold, X, p) = zero_vector!(M.manifold, X, p)

is_metric_function(::Any) = false

for mf in [
        change_metric,
        change_metric!,
        change_representer,
        change_representer!,
        exp,
        exp!,
        exp_fused,
        exp_fused!,
        get_basis,
        get_coordinates,
        get_coordinates!,
        get_vector,
        get_vector!,
        get_vectors,
        inner,
        inverse_retract,
        inverse_retract!,
        log,
        log!,
        mid_point,
        norm,
        parallel_transport_direction,
        parallel_transport_direction!,
        parallel_transport_to,
        parallel_transport_to!,
        retract,
        retract!,
        retract_fused,
        retract_fused!,
        vector_transport_direction,
        vector_transport_direction!,
        vector_transport_to,
        vector_transport_to!,
        Weingarten,
        Weingarten!,
    ]
    @eval is_metric_function(::typeof($mf)) = true
end
