"""
    AbstractInverseRetractionMethod

Abstract type for methods for inverting a retraction (see [`inverse_retract`](@ref)).
"""
abstract type AbstractInverseRetractionMethod end

"""
    AbstractRetractionMethod

Abstract type for methods for [`retract`](@ref)ing a tangent vector to a manifold.
"""
abstract type AbstractRetractionMethod end

"""
    ApproximateInverseRetraction <: AbstractInverseRetractionMethod

An abstract type for representing approximate inverse retraction methods.

"""
abstract type ApproximateInverseRetraction <: AbstractInverseRetractionMethod end

"""
    ApproximateRetraction <: AbstractInverseRetractionMethod

An abstract type for representing approximate retraction methods.
"""
abstract type ApproximateRetraction <: AbstractRetractionMethod end

@doc raw"""
    EmbeddedRetraction{T<:AbstractRetractionMethod} <: AbstractRetractionMethod

Compute a retraction by using the retraction of type `T` in the embedding and projecting the result.

# Constructor

    EmbeddedRetraction(r::AbstractRetractionMethod)

Generate the retraction with retraction `r` to use in the embedding.
"""
struct EmbeddedRetraction{T<:AbstractRetractionMethod} <: AbstractRetractionMethod
    retraction::T
end

"""
    ExponentialRetraction <: AbstractRetractionMethod

Retraction using the exponential map.
"""
struct ExponentialRetraction <: AbstractRetractionMethod end

@doc raw"""
    ODEExponentialRetraction{T<:AbstractRetractionMethod, B<:AbstractBasis} <: AbstractRetractionMethod

Approximate the exponential map on the manifold by evaluating the ODE descripting the geodesic at 1,
assuming the default connection of the given manifold by solving the ordinary differential
equation

```math
\frac{d^2}{dt^2} p^k + Γ^k_{ij} \frac{d}{dt} p_i \frac{d}{dt} p_j = 0,
```

where ``Γ^k_{ij}`` are the Christoffel symbols of the second kind, and
the Einstein summation convention is assumed.

# Constructor

    ODEExponentialRetraction(
        r::AbstractRetractionMethod,
        b::AbstractBasis=DefaultOrthogonalBasis(),
    )

Generate the retraction with a retraction to use internally (for some approaches)
and a basis for the tangent space(s).
"""
struct ODEExponentialRetraction{T<:AbstractRetractionMethod,B<:AbstractBasis} <:
       AbstractRetractionMethod
    retraction::T
    basis::B
end
function ODEExponentialRetraction(r::T) where {T<:AbstractRetractionMethod}
    return ODEExponentialRetraction(r, DefaultOrthonormalBasis())
end
function ODEExponentialRetraction(::T, b::CachedBasis) where {T<:AbstractRetractionMethod}
    return throw(
        DomainError(
            b,
            "Cached Bases are currently not supported, since the basis has to be implemented in a surrounding of the start point as well.",
        ),
    )
end
function ODEExponentialRetraction(r::ExponentialRetraction, ::AbstractBasis)
    return throw(
        DomainError(
            r,
            "You can not use the exponential map as an inner method to solve the ode for the exponential map.",
        ),
    )
end
function ODEExponentialRetraction(r::ExponentialRetraction, ::CachedBasis)
    return throw(
        DomainError(
            r,
            "Neither the exponential map nor a Cached Basis can be used with this retraction type.",
        ),
    )
end

"""
    PolarRetraction <: AbstractRetractionMethod

Retractions that are based on singular value decompositions of the matrix / matrices
for point and tangent vector on a [`AbstractManifold`](@ref)
"""
struct PolarRetraction <: AbstractRetractionMethod end

"""
    ProjectionRetraction <: AbstractRetractionMethod

Retractions that are based on projection and usually addition in the embedding.
"""
struct ProjectionRetraction <: AbstractRetractionMethod end

"""
    QRRetraction <: AbstractRetractionMethod

Retractions that are based on a QR decomposition of the
matrix / matrices for point and tangent vector on a [`AbstractManifold`](@ref)
"""
struct QRRetraction <: AbstractRetractionMethod end

"""
    RetractionWithKeywords{R<:AbstractRetractionMethod,K} <: AbstractRetractionMethod

Since retractions might have keywords, this type is a way to set them as an own type to be
used as a specific retraction.
Another reason for this type is that we dispatch on the retraction first and only the
last layer would be implemented with keywords, so this way they can be passed down.

## Fields

* `retraction` the retraction that is decorated with keywords
* `kwargs` the keyword arguments

Note that you can nest this type. Then the most outer specification of a keyword is used.

## Constructor

    RetractionWithKeywords(m::T; kwargs...) where {T <: AbstractRetractionMethod}

Specify the subtype `T <: `[`AbstractRetractionMethod`](@ref) to have keywords `kwargs...`.
"""
struct RetractionWithKeywords{T<:AbstractRetractionMethod,K} <: AbstractRetractionMethod
    retraction::T
    kwargs::K
end
function RetractionWithKeywords(m::T; kwargs...) where {T<:AbstractRetractionMethod}
    return RetractionWithKeywords{T,typeof(kwargs)}(m, kwargs)
end

"""
    SoftmaxRetraction <: AbstractRetractionMethod

Describes a retraction that is based on the softmax function.
"""
struct SoftmaxRetraction <: AbstractRetractionMethod end


@doc raw"""
    PadeRetraction{m} <: AbstractRetractionMethod

A retraction based on the Padé approximation of order $m$
"""
struct PadeRetraction{m} <: AbstractRetractionMethod end

function PadeRetraction(m::Int)
    (m < 1) && error(
        "The Padé based retraction is only available for positive orders, not for order $m.",
    )
    return PadeRetraction{m}()
end
@doc raw"""
    CayleyRetraction <: AbstractRetractionMethod

A retraction based on the Cayley transform, which is realized by using the
[`PadeRetraction`](@ref)`{1}`.
"""
const CayleyRetraction = PadeRetraction{1}

@doc raw"""
   EmbeddedInverseRetraction{T<:AbstractInverseRetractionMethod} <: AbstractInverseRetractionMethod

Compute an inverse retraction by using the inverse retraction of type `T` in the embedding and projecting the result

# Constructor

    EmbeddedInverseRetraction(r::AbstractInverseRetractionMethod)

Generate the inverse retraction with inverse retraction `r` to use in the embedding.
"""
struct EmbeddedInverseRetraction{T<:AbstractInverseRetractionMethod} <:
       AbstractInverseRetractionMethod
    inverse_retraction::T
end


"""
    LogarithmicInverseRetraction <: AbstractInverseRetractionMethod

Inverse retraction using the [`log`](@ref)arithmic map.
"""
struct LogarithmicInverseRetraction <: AbstractInverseRetractionMethod end


@doc raw"""
    PadeInverseRetraction{m} <: AbstractRetractionMethod

An inverse retraction based on the Padé approximation of order $m$ for the retraction.
"""
struct PadeInverseRetraction{m} <: AbstractInverseRetractionMethod end

function PadeInverseRetraction(m::Int)
    (m < 1) && error(
        "The Padé based inverse retraction is only available for positive orders, not for order $m.",
    )
    return PadeInverseRetraction{m}()
end
@doc raw"""
    CayleyInverseRetraction <: AbstractInverseRetractionMethod

A retraction based on the Cayley transform, which is realized by using the
[`PadeRetraction`](@ref)`{1}`.
"""
const CayleyInverseRetraction = PadeInverseRetraction{1}

"""
    PolarInverseRetraction <: AbstractInverseRetractionMethod

Inverse retractions that are based on a singular value decomposition of the
matrix / matrices for point and tangent vector on a [`AbstractManifold`](@ref)
"""
struct PolarInverseRetraction <: AbstractInverseRetractionMethod end

"""
    ProjectionInverseRetraction <: AbstractInverseRetractionMethod

Inverse retractions that are based on a projection (or its inversion).
"""
struct ProjectionInverseRetraction <: AbstractInverseRetractionMethod end

"""
    QRInverseRetraction <: AbstractInverseRetractionMethod

Inverse retractions that are based on a QR decomposition of the
matrix / matrices for point and tangent vector on a [`AbstractManifold`](@ref)
"""
struct QRInverseRetraction <: AbstractInverseRetractionMethod end

"""
    NLSolveInverseRetraction{T<:AbstractRetractionMethod,TV,TK} <:
        ApproximateInverseRetraction

An inverse retraction method for approximating the inverse of a retraction using `NLsolve`.

# Constructor

    NLSolveInverseRetraction(
        method::AbstractRetractionMethod[, X0];
        project_tangent=false,
        project_point=false,
        nlsolve_kwargs...,
    )

Constructs an approximate inverse retraction for the retraction `method` with initial guess
`X0`, defaulting to the zero vector. If `project_tangent` is `true`, then the tangent
vector is projected before the retraction using `project`. If `project_point` is `true`,
then the resulting point is projected after the retraction. `nlsolve_kwargs` are keyword
arguments passed to `NLsolve.nlsolve`.
"""
struct NLSolveInverseRetraction{TR<:AbstractRetractionMethod,TV,TK} <:
       ApproximateInverseRetraction
    retraction::TR
    X0::TV
    project_tangent::Bool
    project_point::Bool
    nlsolve_kwargs::TK
    function NLSolveInverseRetraction(m, X0, project_point, project_tangent, nlsolve_kwargs)
        return new{typeof(m),typeof(X0),typeof(nlsolve_kwargs)}(
            m,
            X0,
            project_point,
            project_tangent,
            nlsolve_kwargs,
        )
    end
end
function NLSolveInverseRetraction(
    m,
    X0 = nothing;
    project_tangent::Bool = false,
    project_point::Bool = false,
    nlsolve_kwargs...,
)
    return NLSolveInverseRetraction(m, X0, project_point, project_tangent, nlsolve_kwargs)
end


"""
    InverseRetractionWithKeywords{R<:AbstractRetractionMethod,K} <: AbstractRetractionMethod

Since inverse retractions might have keywords, this type is a way to set them as an own type to be
used as a specific inverse retraction.
Another reason for this type is that we dispatch on the inverse retraction first and only the
last layer would be implemented with keywords, so this way they can be passed down.

## Fields

* `inverse_retraction` the inverse retraction that is decorated with keywords
* `kwargs` the keyword arguments

Note that you can nest this type. Then the most outer specification of a keyword is used.

## Constructor

    InverseRetractionWithKeywords(m::T; kwargs...) where {T <: AbstractInverseRetractionMethod}

Specify the subtype `T <: `[`AbstractInverseRetractionMethod`](@ref) to have keywords `kwargs...`.
"""
struct InverseRetractionWithKeywords{T<:AbstractInverseRetractionMethod,K} <:
       AbstractInverseRetractionMethod
    inverse_retraction::T
    kwargs::K
end
function InverseRetractionWithKeywords(
    m::T;
    kwargs...,
) where {T<:AbstractInverseRetractionMethod}
    return InverseRetractionWithKeywords{T,typeof(kwargs)}(m, kwargs)
end

"""
    SoftmaxInverseRetraction <: AbstractInverseRetractionMethod

Describes an inverse retraction that is based on the softmax function.
"""
struct SoftmaxInverseRetraction <: AbstractInverseRetractionMethod end

"""
    default_inverse_retraction_method(M::AbstractManifold)
    default_inverse_retraction_method(M::AbstractManifold, ::Type{T}) where {T}

The [`AbstractInverseRetractionMethod`](@ref) that is used when calling
[`inverse_retract`](@ref) without specifying the inverse retraction method.
By default, this is the [`LogarithmicInverseRetraction`](@ref).

This method can also be specified more precisely with a point type `T`, for the case
that on a `M` there are two different representations of points, which provide
different inverse retraction methods.
"""
default_inverse_retraction_method(::AbstractManifold) = LogarithmicInverseRetraction()
function default_inverse_retraction_method(M::AbstractManifold, ::Type{T}) where {T}
    return default_inverse_retraction_method(M)
end

"""
    default_retraction_method(M::AbstractManifold)
    default_retraction_method(M::AbstractManifold, ::Type{T}) where {T}

The [`AbstractRetractionMethod`](@ref) that is used when calling [`retract`](@ref) without
specifying the retraction method. By default, this is the [`ExponentialRetraction`](@ref).

This method can also be specified more precisely with a point type `T`, for the case
that on a `M` there are two different representations of points, which provide
different retraction methods.
"""
default_retraction_method(::AbstractManifold) = ExponentialRetraction()
function default_retraction_method(M::AbstractManifold, ::Type{T}) where {T}
    return default_retraction_method(M)
end

"""
    inverse_retract!(M::AbstractManifold, X, p, q[, method::AbstractInverseRetractionMethod])

Compute the inverse retraction, a cheaper, approximate version of the
[`log`](@ref)arithmic map), of points `p` and `q` on the [`AbstractManifold`](@ref) `M`.
Result is saved to `X`.

Inverse retraction method can be specified by the last argument, defaulting to
[`default_inverse_retraction_method`](@ref)`(M)`. See the documentation of respective manifolds for
available methods.

See also [`retract!`](@ref).
"""
function inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)
    return _inverse_retract!(M, X, p, q, m)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    ::LogarithmicInverseRetraction;
    kwargs...,
)
    return log!(M, X, p, q; kwargs...)
end

#
# dispatch to lower level
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    ::CayleyInverseRetraction;
    kwargs...,
)
    return inverse_retract_cayley!(M, X, p, q; kwargs...)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    m::EmbeddedInverseRetraction;
    kwargs...,
)
    return inverse_retract_embedded!(M, X, p, q, m.inverse_retraction; kwargs...)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    m::NLSolveInverseRetraction;
    kwargs...,
)
    return inverse_retract_nlsolve!(M, X, p, q, m; kwargs...)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    m::PadeInverseRetraction;
    kwargs...,
)
    return inverse_retract_pade!(M, X, p, q, m; kwargs...)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    ::PolarInverseRetraction;
    kwargs...,
)
    return inverse_retract_polar!(M, X, p, q; kwargs...)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    ::ProjectionInverseRetraction;
    kwargs...,
)
    return inverse_retract_project!(M, X, p, q; kwargs...)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, ::QRInverseRetraction; kwargs...)
    return inverse_retract_qr!(M, X, p, q; kwargs...)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    m::InverseRetractionWithKeywords;
    kwargs...,
)
    return _inverse_retract!(M, X, p, q, m.inverse_retraction; kwargs..., m.kwargs...)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    ::SoftmaxInverseRetraction;
    kwargs...,
)
    return inverse_retract_softmax!(M, X, p, q; kwargs...)
end
"""
    inverse_retract_embedded!(M::AbstractManifold, X, p, q, m::AbstractInverseRetractionMethod)

Compute the in-place variant of the [`EmbeddedInverseRetraction`](@ref) using
the [`AbstractInverseRetractionMethod`](@ref) `m` in the embedding (see [`get_embedding`](@ref))
and projecting the result back.
"""
function inverse_retract_embedded!(
    M::AbstractManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)
    return project!(
        M,
        X,
        p,
        inverse_retract(
            get_embedding(M),
            embed(get_embedding(M), p),
            embed(get_embedding(M), q),
            m,
        ),
    )
end

"""
    inverse_retract_cayley!(M::AbstractManifold, X, p, q)

Compute the in-place variant of the [`CayleyInverseRetraction`](@ref),
which by default calls the first order [`PadeInverseRetraction`§(@ref).
"""
function inverse_retract_cayley!(M::AbstractManifold, X, p, q; kwargs...)
    return inverse_retract_pade!(M, X, p, q, 1; kwargs...)
end

"""
    inverse_retract_pade!(M::AbstractManifold, p, q, n)

Compute the in-place variant of the [`PadeInverseRetraction`](@ref)`(n)`,
"""
inverse_retract_pade!(M::AbstractManifold, p, q, n)

function inverse_retract_pade! end

"""
    inverse_retract_qr!(M::AbstractManifold, X, p, q)

Compute the in-place variant of the [`QRInverseRetraction`](@ref).
"""
inverse_retract_qr!(M::AbstractManifold, X, p, q)

function inverse_retract_qr! end

"""
    inverse_retract_project!(M::AbstractManifold, X, p, q)

Compute the in-place variant of the [`ProjectionInverseRetraction`](@ref).
"""
inverse_retract_project!(M::AbstractManifold, X, p, q)

function inverse_retract_project! end

"""
    inverse_retract_polar!(M::AbstractManifold, X, p, q)

Compute the in-place variant of the [`PolarInverseRetraction`](@ref).
"""
inverse_retract_polar!(M::AbstractManifold, X, p, q)

function inverse_retract_polar! end

"""
    inverse_retract_nlsolve!(M::AbstractManifold, X, p, q, m::NLSolveInverseRetraction)

Compute the in-place variant of the [`NLSolveInverseRetraction`](@ref) `m`.
"""
inverse_retract_nlsolve!(M::AbstractManifold, X, p, q, m::NLSolveInverseRetraction)

function inverse_retract_nlsolve! end

"""
    inverse_retract_softmax!(M::AbstractManifold, X, p, q)

Compute the in-place variant of the [`SoftmaxInverseRetraction`](@ref).
"""
inverse_retract_softmax!(M::AbstractManifold, X, p, q)

function inverse_retract_softmax! end

"""
    inverse_retract(M::AbstractManifold, p, q)
    inverse_retract(M::AbstractManifold, p, q, method::AbstractInverseRetractionMethod

Compute the inverse retraction, a cheaper, approximate version of the
[`log`](@ref)arithmic map), of points `p` and `q` on the [`AbstractManifold`](@ref) `M`.

Inverse retraction method can be specified by the last argument, defaulting to
[`default_inverse_retraction_method`](@ref)`(M)`.
For available inverse retractions on certain manifolds see the documentation on the
corresponding manifold.

See also [`retract`](@ref).
"""
function inverse_retract(
    M::AbstractManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)
    return _inverse_retract(M, p, q, m)
end
function _inverse_retract(M::AbstractManifold, p, q, ::CayleyInverseRetraction; kwargs...)
    return inverse_retract_cayley(M, p, q; kwargs...)
end
function _inverse_retract(
    M::AbstractManifold,
    p,
    q,
    m::EmbeddedInverseRetraction;
    kwargs...,
)
    return inverse_retract_embedded(M, p, q, m.inverse_retraction; kwargs...)
end
function _inverse_retract(
    M::AbstractManifold,
    p,
    q,
    ::LogarithmicInverseRetraction;
    kwargs...,
)
    return log(M, p, q; kwargs...)
end
function _inverse_retract(M::AbstractManifold, p, q, m::NLSolveInverseRetraction; kwargs...)
    return inverse_retract_nlsolve(M, p, q, m; kwargs...)
end
function _inverse_retract(
    M::AbstractManifold,
    p,
    q,
    ::PadeInverseRetraction{n};
    kwargs...,
) where {n}
    return inverse_retract_pade(M, p, q, n; kwargs...)
end
function _inverse_retract(M::AbstractManifold, p, q, ::PolarInverseRetraction; kwargs...)
    return inverse_retract_polar(M, p, q; kwargs...)
end
function _inverse_retract(
    M::AbstractManifold,
    p,
    q,
    ::ProjectionInverseRetraction;
    kwargs...,
)
    return inverse_retract_project(M, p, q; kwargs...)
end
function _inverse_retract(M::AbstractManifold, p, q, ::QRInverseRetraction; kwargs...)
    return inverse_retract_qr(M, p, q; kwargs...)
end
function _inverse_retract(
    M::AbstractManifold,
    p,
    q,
    m::InverseRetractionWithKeywords;
    kwargs...,
)
    return _inverse_retract(M, p, q, m.inverse_retraction; kwargs..., m.kwargs...)
end
function _inverse_retract(M::AbstractManifold, p, q, ::SoftmaxInverseRetraction; kwargs...)
    return inverse_retract_softmax(M, p, q; kwargs...)
end
"""
    inverse_retract_embedded(M::AbstractManifold, p, q, m::AbstractInverseRetractionMethod)

computes the allocating variant of the [`EmbeddedInverseRetraction`](@ref) using
the [`AbstractInverseRetractionMethod`](@ref) `m` in the embedding (see [`get_embedding`](@ref))
and projecting the result back.
"""
function inverse_retract_embedded(
    M::AbstractManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)
    return project(
        M,
        p,
        inverse_retract(
            get_embedding(M),
            embed(get_embedding(M), p),
            embed(get_embedding(M), q),
            m,
        ),
    )
end

"""
    inverse_retract_cayley(M::AbstractManifold, p, q)

computes the allocating variant of the [`CayleyInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_cayley!`](@ref).
"""
function inverse_retract_cayley(M::AbstractManifold, p, q; kwargs...)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_cayley!(M, X, p, q; kwargs...)
end
"""
    inverse_retract_pade(M::AbstractManifold, p, q)

computes the allocating variant of the [`PadeInverseRetraction`](@ref)`(n)`,
which by default allocates and calls [`inverse_retract_pade!`](@ref ManifoldsBase.inverse_retract_pade!).
"""
function inverse_retract_pade(M::AbstractManifold, p, q, n; kwargs...)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_pade!(M, X, p, q, n; kwargs...)
end

"""
    inverse_retract_polar(M::AbstractManifold, p, q)

computes the allocating variant of the [`PolarInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_polar!`](@ref ManifoldsBase.inverse_retract_polar!).
"""
function inverse_retract_polar(M::AbstractManifold, p, q; kwargs...)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_polar!(M, X, p, q; kwargs...)
end
"""
    inverse_retract_project(M::AbstractManifold, p, q)

computes the allocating variant of the [`ProjectionInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_project!`](@ref ManifoldsBase.inverse_retract_project!).
"""
function inverse_retract_project(M::AbstractManifold, p, q; kwargs...)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_project!(M, X, p, q; kwargs...)
end
"""
    inverse_retract_qr(M::AbstractManifold, p, q)

computes the allocating variant of the [`QRInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_qr!`](@ref ManifoldsBase.inverse_retract_qr!).
"""
function inverse_retract_qr(M::AbstractManifold, p, q; kwargs...)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_qr!(M, X, p, q; kwargs...)
end
"""
    inverse_retract_nlsolve(M::AbstractManifold, p, q, m::NLSolveInverseRetraction)

computes the allocating variant of the [`NLSolveInverseRetraction`](@ref) `m`,
which by default allocates and calls [`inverse_retract_nlsolve!`](@ref).
"""
function inverse_retract_nlsolve(
    M::AbstractManifold,
    p,
    q,
    m::NLSolveInverseRetraction;
    kwargs...,
)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_nlsolve!(M, X, p, q, m; kwargs...)
end
"""
    inverse_retract_softmax(M::AbstractManifold, p, q)

computes the allocating variant of the [`SoftmaxInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_softmax!`](@ref ManifoldsBase.inverse_retract_softmax!).
"""
function inverse_retract_softmax(M::AbstractManifold, p, q; kwargs...)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_softmax!(M, X, p, q; kwargs...)
end


"""
    retract!(M::AbstractManifold, q, p, X)
    retract!(M::AbstractManifold, q, p, X, t::Real=1)
    retract!(M::AbstractManifold, q, p, X, method::AbstractRetractionMethod)
    retract!(M::AbstractManifold, q, p, X, t::Real=1, method::AbstractRetractionMethod)

Compute a retraction, a cheaper, approximate version of the [`exp`](@ref)onential map,
from `p` into direction `X`, scaled by `t`, on the [`AbstractManifold`](@ref) manifold `M`.
Result is saved to `q`.

Retraction method can be specified by the last argument, defaulting to
[`default_retraction_method`](@ref)`(M)`. See the documentation of respective manifolds for available
methods.

See [`retract`](@ref) for more details.
"""
function retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return _retract!(M, q, p, X, t, m)
end
function retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    method::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return retract!(M, q, p, X, one(number_eltype(X)), method)
end
# dispatch to lower level
function _retract!(M::AbstractManifold, q, p, X, t::Number, ::CayleyRetraction; kwargs...)
    return retract_cayley!(M, q, p, X, t; kwargs...)
end
function _retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::EmbeddedRetraction;
    kwargs...,
)
    return retract_embedded!(M, q, p, X, t, m.retraction; kwargs...)
end
function _retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    ::ExponentialRetraction;
    kwargs...,
)
    return exp!(M, q, p, X, t; kwargs...)
end
function _retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::ODEExponentialRetraction;
    kwargs...,
)
    return retract_exp_ode!(M, q, p, X, t, m.retraction, m.basis; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, t::Number, ::PolarRetraction; kwargs...)
    return retract_polar!(M, q, p, X, t; kwargs...)
end
function _retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    ::ProjectionRetraction;
    kwargs...,
)
    return retract_project!(M, q, p, X, t; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, t::Number, ::QRRetraction; kwargs...)
    return retract_qr!(M, q, p, X, t; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, t::Number, ::SoftmaxRetraction; kwargs...)
    return retract_softmax!(M, q, p, X, t; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, t::Number, m::PadeRetraction; kwargs...)
    return retract_pade!(M, q, p, X, t, m; kwargs...)
end
function _retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::RetractionWithKeywords;
    kwargs...,
)
    return _retract!(M, q, p, X, t, m.retraction; kwargs..., m.kwargs...)
end

"""
    retract_embedded!(M::AbstractManifold, q, p, X, t::Number, m::AbstractRetractionMethod)

Compute the in-place variant of the [`EmbeddedRetraction`](@ref) using
the [`AbstractRetractionMethod`](@ref) `m` in the embedding (see [`get_embedding`](@ref))
and projecting the result back.
"""
function retract_embedded!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod;
    kwargs...,
)
    return project!(
        M,
        q,
        retract(
            get_embedding(M),
            embed(get_embedding(M), p),
            embed(get_embedding(M), p, X),
            t,
            m;
            kwargs...,
        ),
    )
end

"""
    retract_cayley!(M::AbstractManifold, q, p, X, t::Number)

Compute the in-place variant of the [`CayleyRetraction`](@ref),
which by default falls back to calling the first order [`PadeRetraction`](@ref).
"""
function retract_cayley!(M::AbstractManifold, q, p, X, t::Number; kwargs...)
    return retract_pade!(M, q, p, X, t, PadeRetraction(1); kwargs...)
end

"""
    retract_exp_ode!(M::AbstractManifold, q, p, X, t::Number, m::AbstractRetractionMethod, B::AbstractBasis)

Compute the in-place variant of the [`ODEExponentialRetraction`](@ref)`(m, B)`.
"""
retract_exp_ode!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod,
    B::AbstractBasis,
)

function retract_exp_ode! end

"""
    retract_pade!(M::AbstractManifold, q, p, X, t::Number, m::PadeRetraction)

Compute the in-place variant of the [`PadeRetraction`](@ref) `m`.
"""
retract_pade!(M::AbstractManifold, q, p, X, t::Number, m::PadeRetraction)

function retract_pade! end

"""
    retract_project!(M::AbstractManifold, q, p, X, t::Number)

Compute the in-place variant of the [`ProjectionRetraction`](@ref).
"""
retract_project!(M::AbstractManifold, q, p, X, t::Number)

function retract_project! end

"""
    retract_polar!(M::AbstractManifold, q, p, X, t::Number)

Compute the in-place variant of the [`PolarRetraction`](@ref).
"""
retract_polar!(M::AbstractManifold, q, p, X, t::Number)

function retract_polar! end

"""
    retract_qr!(M::AbstractManifold, q, p, X, t::Number)

Compute the in-place variant of the [`QRRetraction`](@ref).
"""
retract_qr!(M::AbstractManifold, q, p, X, t::Number)

function retract_qr! end

"""
    retract_softmax!(M::AbstractManifold, q, p, X, t::Number)

Compute the in-place variant of the [`SoftmaxRetraction`](@ref).
"""
retract_softmax!(M::AbstractManifold, q, p, X, t::Number)

function retract_softmax! end

@doc raw"""
    retract(M::AbstractManifold, p, X, method::AbstractRetractionMethod=default_retraction_method(M, typeof(p)))
    retract(M::AbstractManifold, p, X, t::Number=1, method::AbstractRetractionMethod=default_retraction_method(M, typeof(p)))

Compute a retraction, a cheaper, approximate version of the [`exp`](@ref)onential map,
from `p` into direction `X`, scaled by `t`, on the [`AbstractManifold`](@ref) `M`.

A retraction ``\operatorname{retr}_p: T_p\mathcal M → \mathcal M`` is a smooth map that fulfils

1. ``\operatorname{retr}_p(0) = p``
2. ``D\operatorname{retr}_p(0): T_p\mathcal M \to T_p\mathcal M`` is the identity map,
i.e. ``D\operatorname{retr}_p(0)[X]=X`` holds for all ``X\in T_p\mathcal M``,

where ``D\operatorname{retr}_p`` denotes the differential of the retraction

The retraction is called of second order if for all ``X`` the curves ``c(t) = R_p(tX)``
have a zero acceleration at ``t=0``, i.e. ``c''(0) = 0``.

Retraction method can be specified by the last argument, defaulting to
[`default_retraction_method`](@ref)`(M)`. For further available retractions see the documentation of respective manifolds.

Locally, the retraction is invertible. For the inverse operation, see [`inverse_retract`](@ref).
"""
function retract(
    M::AbstractManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return retract(M, p, X, one(number_eltype(X)), m)
end
function retract(
    M::AbstractManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return _retract(M, p, X, t, m)
end
function _retract(M::AbstractManifold, p, X, t::Number, m::EmbeddedRetraction; kwargs...)
    return retract_embedded(M, p, X, t, m.retraction; kwargs...)
end
function _retract(M::AbstractManifold, p, X, t::Number, ::ExponentialRetraction; kwargs...)
    return exp(M, p, X, t)
end
function _retract(
    M::AbstractManifold,
    p,
    X,
    t::Number,
    m::ODEExponentialRetraction;
    kwargs...,
)
    return retract_exp_ode(M, p, X, t, m.retraction, m.basis)
end
function _retract(M::AbstractManifold, p, X, t::Number, ::PolarRetraction; kwargs...)
    return retract_polar(M, p, X, t; kwargs...)
end
function _retract(M::AbstractManifold, p, X, t::Number, ::ProjectionRetraction; kwargs...)
    return retract_project(M, p, X, t; kwargs...)
end
function _retract(M::AbstractManifold, p, X, t::Number, ::QRRetraction; kwargs...)
    return retract_qr(M, p, X, t; kwargs...)
end
function _retract(M::AbstractManifold, p, X, t::Number, ::SoftmaxRetraction; kwargs...)
    return retract_softmax(M, p, X, t; kwargs...)
end
function _retract(M::AbstractManifold, p, X, t::Number, ::CayleyRetraction; kwargs...)
    return retract_cayley(M, p, X, t; kwargs...)
end
function _retract(M::AbstractManifold, p, X, t::Number, m::PadeRetraction; kwargs...)
    return retract_pade(M, p, X, t, m; kwargs...)
end
function _retract(
    M::AbstractManifold,
    p,
    X,
    t::Number,
    m::RetractionWithKeywords;
    kwargs...,
)
    return _retract(M, p, X, t, m.retraction; kwargs..., m.kwargs...)
end
"""
    retract_embedded(M::AbstractManifold, p, X, t::Number, m::AbstractRetractionMethod)

computes the allocating variant of the [`EmbeddedRetraction`](@ref) using
the [`AbstractRetractionMethod`](@ref) `m` in the embedding (see [`get_embedding`](@ref))
and projecting the result back.
"""
function retract_embedded(
    M::AbstractManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod;
    kwargs...,
)
    return project(
        M,
        retract(
            get_embedding(M),
            embed(get_embedding(M), p),
            embed(get_embedding(M), p, X),
            t,
            m;
            kwargs...,
        ),
    )
end
"""
    retract_polar(M::AbstractManifold, p, X, t::Number)

computes the allocating variant of the [`PolarRetraction`](@ref),
which by default allocates and calls [`retract_polar!`](@ref ManifoldsBase.retract_polar!).
"""
function retract_polar(M::AbstractManifold, p, X, t::Number; kwargs...)
    q = allocate_result(M, retract, p, X)
    return retract_polar!(M, q, p, X, t; kwargs...)
end
"""
    retract_project(M::AbstractManifold, p, X, t::Number)

Compute the allocating variant of the [`ProjectionRetraction`](@ref),
which by default allocates and calls [`retract_project!`](@ref ManifoldsBase.retract_project!).
"""
function retract_project(M::AbstractManifold, p, X, t::Number; kwargs...)
    q = allocate_result(M, retract, p, X)
    return retract_project!(M, q, p, X, t; kwargs...)
end
"""
    retract_qr(M::AbstractManifold, p, X, t::Number)

Compute the allocating variant of the [`QRRetraction`](@ref),
which by default allocates and calls [`retract_qr!`](@ref ManifoldsBase.retract_qr!).
"""
function retract_qr(M::AbstractManifold, p, X, t::Number; kwargs...)
    q = allocate_result(M, retract, p, X, t::Number)
    return retract_qr!(M, q, p, X, t; kwargs...)
end
"""
    retract_exp_ode(M::AbstractManifold, p, q, m::AbstractRetractionMethod, B::AbstractBasis)

Compute the allocating variant of the [`ODEExponentialRetraction`](@ref)`(m,B)`,
which by default allocates and calls [`retract_exp_ode!`](@ref ManifoldsBase.retract_exp_ode!).
"""
function retract_exp_ode(
    M::AbstractManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod,
    B::AbstractBasis;
    kwargs...,
)
    q = allocate_result(M, retract, p, X)
    return retract_exp_ode!(M, q, p, X, t, m, B; kwargs...)
end
"""
    retract_softmax(M::AbstractManifold, p, X, t::Number)

Compute the allocating variant of the [`SoftmaxRetraction`](@ref),
which by default allocates and calls [`retract_softmax!`](@ref ManifoldsBase.retract_softmax!).
"""
function retract_softmax(M::AbstractManifold, p, X, t::Number; kwargs...)
    q = allocate_result(M, retract, p, X)
    return retract_softmax!(M, q, p, X, t; kwargs...)
end

"""
    retract_cayley(M::AbstractManifold, p, X, t::Number)

Compute the allocating variant of the [`CayleyRetraction`](@ref),
which by default allocates and calls [`retract_cayley!`](@ref ManifoldsBase.retract_cayley!).
"""
function retract_cayley(M::AbstractManifold, p, X, t::Number; kwargs...)
    q = allocate_result(M, retract, p, X)
    return retract_cayley!(M, q, p, X, t; kwargs...)
end
"""
    retract_pade(M::AbstractManifold, p, X, t::Number, m::PadeRetraction)

Compute the allocating variant of the [`PadeRetraction`](@ref) `m`,
which by default allocates and calls [`retract_pade!`](@ref ManifoldsBase.retract_pade!).
"""
function retract_pade(M::AbstractManifold, p, X, t::Number, m::PadeRetraction; kwargs...)
    q = allocate_result(M, retract, p, X)
    return retract_pade!(M, q, p, X, t, m; kwargs...)
end

Base.show(io::IO, ::CayleyRetraction) = print(io, "CayleyRetraction()")
Base.show(io::IO, ::PadeRetraction{m}) where {m} = print(io, "PadeRetraction($(m))")
