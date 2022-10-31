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
used as a specific retraction

## Fields

* `retraction` the retraction that is decorated with keywords
* `kwargs` th ekeyword arguments

## Constructor

    RetractionWithKeywords(m::T; kwargs...) where {T <: AbstractRetractionMethod}


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
    SoftmaxInverseRetraction <: AbstractInverseRetractionMethod

Describes an inverse retraction that is based on the softmax function.
"""
struct SoftmaxInverseRetraction <: AbstractInverseRetractionMethod end

"""
    default_inverse_retraction_method(M::AbstractManifold)

The [`AbstractInverseRetractionMethod`](@ref) that is used when calling
[`inverse_retract`](@ref) without specifying the inverse retraction method.
By default, this is the [`LogarithmicInverseRetraction`](@ref).
"""
default_inverse_retraction_method(::AbstractManifold) = LogarithmicInverseRetraction()

"""
    default_retraction_method(M::AbstractManifold)

The [`AbstractRetractionMethod`](@ref) that is used when calling [`retract`](@ref) without
specifying the retraction method. By default, this is the [`ExponentialRetraction`](@ref).
"""
default_retraction_method(::AbstractManifold) = ExponentialRetraction()

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
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return _inverse_retract!(M, X, p, q, m)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, ::LogarithmicInverseRetraction)
    return log!(M, X, p, q)
end

#
# dispatch to lower level
function _inverse_retract!(M::AbstractManifold, X, p, q, ::CayleyInverseRetraction)
    return inverse_retract_caley!(M, X, p, q)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, m::EmbeddedInverseRetraction)
    return inverse_retract_embedded!(M, X, p, q, m.inverse_retraction)
end
function _inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    ::PadeInverseRetraction{n},
) where {n}
    return inverse_retract_pade!(M, X, p, q, n)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, ::PolarInverseRetraction)
    return inverse_retract_polar!(M, X, p, q)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, ::ProjectionInverseRetraction)
    return inverse_retract_project!(M, X, p, q)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, ::QRInverseRetraction)
    return inverse_retract_qr!(M, X, p, q)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, ::SoftmaxInverseRetraction)
    return inverse_retract_softmax!(M, X, p, q)
end
function _inverse_retract!(M::AbstractManifold, X, p, q, m::NLSolveInverseRetraction)
    return inverse_retract_nlsolve!(M, X, p, q, m)
end
"""
    inverse_retract_embedded!(M::AbstractManifold, X, p, q, m::AbstractInverseRetractionMethod)

computes the mutating variant of the [`EmbeddedInverseRetraction`](@ref) using
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
    inverse_retract_caley!(M::AbstractManifold, X, p, q)

computes the mutating variant of the [`CayleyInverseRetraction`](@ref),
which by default calls the first order [`PadeInverseRetraction`§(@ref).
"""
function inverse_retract_caley!(M::AbstractManifold, X, p, q)
    return inverse_retract_pade!(M, X, p, q, 1)
end

"""
    inverse_retract_pade!(M::AbstractManifold, p, q, n)

computes the mutating variant of the [`PadeInverseRetraction`](@ref)`(n)`,
"""
inverse_retract_pade!(M::AbstractManifold, p, q, n)

function inverse_retract_pade! end

"""
    inverse_retract_qr!(M::AbstractManifold, X, p, q)

computes the mutating variant of the [`QRInverseRetraction`](@ref).
"""
inverse_retract_qr!(M::AbstractManifold, X, p, q)

function inverse_retract_qr! end

"""
    inverse_retract_project!(M::AbstractManifold, X, p, q)

computes the mutating variant of the [`ProjectionInverseRetraction`](@ref).
"""
inverse_retract_project!(M::AbstractManifold, X, p, q)

function inverse_retract_project! end

"""
    inverse_retract_polar!(M::AbstractManifold, X, p, q)

computes the mutating variant of the [`PolarInverseRetraction`](@ref).
"""
inverse_retract_polar!(M::AbstractManifold, X, p, q)

function inverse_retract_polar! end

"""
    inverse_retract_nlsolve!(M::AbstractManifold, X, p, q, m::NLSolveInverseRetraction)

computes the mutating variant of the [`NLSolveInverseRetraction`](@ref) `m`.
"""
inverse_retract_nlsolve!(M::AbstractManifold, X, p, q, m::NLSolveInverseRetraction)

function inverse_retract_nlsolve! end

"""
    inverse_retract_softmax!(M::AbstractManifold, X, p, q)

computes the mutating variant of the [`SoftmaxInverseRetraction`](@ref).
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
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return _inverse_retract(M, p, q, m)
end
function _inverse_retract(M::AbstractManifold, p, q, ::CayleyInverseRetraction)
    return inverse_retract_caley(M, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, m::EmbeddedInverseRetraction)
    return inverse_retract_embedded(M, p, q, m.inverse_retraction)
end
_inverse_retract(M::AbstractManifold, p, q, ::LogarithmicInverseRetraction) = log(M, p, q)
function _inverse_retract(M::AbstractManifold, p, q, m::NLSolveInverseRetraction)
    return inverse_retract_nlsolve(M, p, q, m)
end
function _inverse_retract(M::AbstractManifold, p, q, ::PadeInverseRetraction{n}) where {n}
    return inverse_retract_pade(M, p, q, n)
end
function _inverse_retract(M::AbstractManifold, p, q, ::PolarInverseRetraction)
    return inverse_retract_polar(M, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, ::ProjectionInverseRetraction)
    return inverse_retract_project(M, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, ::QRInverseRetraction)
    return inverse_retract_qr(M, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, ::SoftmaxInverseRetraction)
    return inverse_retract_softmax(M, p, q)
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
    inverse_retract_caley(M::AbstractManifold, p, q)

computes the allocating variant of the [`CayleyInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_caley!`](@ref).
"""
function inverse_retract_caley(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_caley!(M, X, p, q)
end
"""
    inverse_retract_pade(M::AbstractManifold, p, q)

computes the allocating variant of the [`PadeInverseRetraction`](@ref)`(n)`,
which by default allocates and calls [`inverse_retract_pade!`](@ref ManifoldsBase.inverse_retract_pade!).
"""
function inverse_retract_pade(M::AbstractManifold, p, q, n)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_pade!(M, X, p, q, n)
end

"""
    inverse_retract_polar(M::AbstractManifold, p, q)

computes the allocating variant of the [`PolarInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_polar!`](@ref ManifoldsBase.inverse_retract_polar!).
"""
function inverse_retract_polar(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_polar!(M, X, p, q)
end
"""
    inverse_retract_project(M::AbstractManifold, p, q)

computes the allocating variant of the [`ProjectionInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_project!`](@ref ManifoldsBase.inverse_retract_project!).
"""
function inverse_retract_project(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_project!(M, X, p, q)
end
"""
    inverse_retract_qr(M::AbstractManifold, p, q)

computes the allocating variant of the [`QRInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_qr!`](@ref ManifoldsBase.inverse_retract_qr!).
"""
function inverse_retract_qr(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_qr!(M, X, p, q)
end
"""
    inverse_retract_nlsolve(M::AbstractManifold, p, q, m::NLSolveInverseRetraction)

computes the allocating variant of the [`NLSolveInverseRetraction`](@ref) `m`,
which by default allocates and calls [`inverse_retract_nlsolve!`](@ref).
"""
function inverse_retract_nlsolve(M::AbstractManifold, p, q, m::NLSolveInverseRetraction)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_nlsolve!(M, X, p, q, m)
end
"""
    inverse_retract_softmax(M::AbstractManifold, p, q)

computes the allocating variant of the [`SoftmaxInverseRetraction`](@ref),
which by default allocates and calls [`inverse_retract_softmax!`](@ref ManifoldsBase.inverse_retract_softmax!).
"""
function inverse_retract_softmax(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_softmax!(M, X, p, q)
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
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return _retract!(M, q, p, X, m)
end
function retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Real,
    method::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract!(M, q, p, t * X, method)
end
# dispatch to lower level
_retract!(M::AbstractManifold, q, p, X, ::CayleyRetraction) = retract_caley!(M, q, p, X)
function _retract!(M::AbstractManifold, q, p, X, m::EmbeddedRetraction)
    return retract_embedded!(M, q, p, X, m.retraction)
end
_retract!(M::AbstractManifold, q, p, X, ::ExponentialRetraction) = exp!(M, q, p, X)
function _retract!(M::AbstractManifold, q, p, X, m::ODEExponentialRetraction)
    return retract_exp_ode!(M, q, p, X, m.retraction, m.basis)
end
_retract!(M::AbstractManifold, q, p, X, ::PolarRetraction) = retract_polar!(M, q, p, X)
function _retract!(M::AbstractManifold, q, p, X, ::ProjectionRetraction)
    return retract_project!(M, q, p, X)
end
_retract!(M::AbstractManifold, q, p, X, ::QRRetraction) = retract_qr!(M, q, p, X)
_retract!(M::AbstractManifold, q, p, X, ::SoftmaxRetraction) = retract_softmax!(M, q, p, X)
function _retract!(M::AbstractManifold, q, p, X, ::PadeRetraction{n}) where {n}
    return retract_pade!(M, q, p, X, n)
end

"""
    retract_embedded!(M::AbstractManifold, X, p, q, m::AbstractRetractionMethod)

computes the mutating variant of the [`EmbeddedRetraction`](@ref) using
the [`AbstractRetractionMethod`](@ref) `m` in the embedding (see [`get_embedding`](@ref))
and projecting the result back.
"""
function retract_embedded!(M::AbstractManifold, q, p, X, m::AbstractRetractionMethod)
    return project!(
        M,
        q,
        retract(
            get_embedding(M),
            embed(get_embedding(M), p),
            embed(get_embedding(M), p, X),
            m,
        ),
    )
end

"""
    retract_caley!(M::AbstractManifold, X, p, q)

computes the mutating variant of the [`CayleyRetraction`](@ref),
which by default falls back to calling the first order [`PadeRetraction`](@ref).
"""
function retract_caley!(M::AbstractManifold, q, p, X)
    return retract_pade!(M, q, p, X, 1)
end

"""
    retract_exp_ode!(M::AbstractManifold, q, p, X, m::AbstractRetractionMethod, B::AbstractBasis)

computes the mutating variant of the [`ODEExponentialRetraction`](@ref)`(m, B)`.
"""
retract_exp_ode!(
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod,
    B::AbstractBasis,
)

function retract_exp_ode! end

"""
    retract_pade!(M::AbstractManifold, q, p, n)

computes the mutating variant of the [`PadeRetraction`](@ref)`(n)`.
"""
retract_pade!(M::AbstractManifold, q, p, n)

function retract_pade! end

"""
    retract_project!(M::AbstractManifold, q, p, X)

computes the mutating variant of the [`ProjectionRetraction`](@ref).
"""
retract_project!(M::AbstractManifold, q, p, X)

function retract_project! end

"""
    retract_polar!(M::AbstractManifold, q, p, X)

computes the mutating variant of the [`PolarRetraction`](@ref).
"""
retract_polar!(M::AbstractManifold, q, p, X)

function retract_polar! end

"""
    retract_qr!(M::AbstractManifold, q, p, X)

computes the mutating variant of the [`QRRetraction`](@ref).
"""
retract_qr!(M::AbstractManifold, q, p, X)

function retract_qr! end

"""
    retract_softmax!(M::AbstractManifold, q, p, X)

computes the mutating variant of the [`SoftmaxRetraction`](@ref).
"""
retract_softmax!(M::AbstractManifold, q, p, X)

function retract_softmax! end

@doc raw"""
    retract(M::AbstractManifold, p, X, method::AbstractRetractionMethod=default_retraction_method(M))
    retract(M::AbstractManifold, p, X, t::Real=1, method::AbstractRetractionMethod=default_retraction_method(M))

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
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return _retract(M, p, X, m)
end
function retract(
    M::AbstractManifold,
    p,
    X,
    t::Real,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract(M, p, t * X, m)
end
function _retract(M::AbstractManifold, p, X, m::EmbeddedRetraction; kwargs...)
    return retract_embedded(M, p, X, m.retraction; kwargs...)
end
_retract(M::AbstractManifold, p, X, ::ExponentialRetraction; kwargs...) = exp(M, p, X)
function _retract(M::AbstractManifold, p, X, m::ODEExponentialRetraction; kwargs...)
    return retract_exp_ode(M, p, X, m.retraction, m.basis)
end
function _retract(M::AbstractManifold, p, X, ::PolarRetraction; kwargs...)
    return retract_polar(M, p, X; kwargs...)
end
function _retract(M::AbstractManifold, p, X, ::ProjectionRetraction; kwargs...)
    return retract_project(M, p, X; kwargs...)
end
function _retract(M::AbstractManifold, p, X, ::QRRetraction; kwargs...)
    return retract_qr(M, p, X; kwargs...)
end
function _retract(M::AbstractManifold, p, X, ::SoftmaxRetraction; kwargs...)
    return retract_softmax(M, p, X; kwargs...)
end
function _retract(M::AbstractManifold, p, X, ::CayleyRetraction; kwargs...)
    return retract_caley(M, p, X; kwargs...)
end
function _retract(M::AbstractManifold, p, X, ::PadeRetraction{n}; kwargs...) where {n}
    return retract_pade(M, p, X, n; kwargs...)
end
function _retract(M::AbstractManifold, p, X, m::RetractionWithKeywords)
    return _retract(M, p, X, m.retraction; m.kwargs...)
end
"""
    retract_embedded(M::AbstractManifold, p, X, m::AbstractRetractionMethod)

computes the allocating variant of the [`EmbeddedRetraction`](@ref) using
the [`AbstractRetractionMethod`](@ref) `m` in the embedding (see [`get_embedding`](@ref))
and projecting the result back.
"""
function retract_embedded(M::AbstractManifold, p, X, m::AbstractRetractionMethod)
    return project(
        M,
        retract(
            get_embedding(M),
            embed(get_embedding(M), p),
            embed(get_embedding(M), p, X),
            m,
        ),
    )
end
"""
    retract_polar(M::AbstractManifold, p, q)

computes the allocating variant of the [`PolarRetraction`](@ref),
which by default allocates and calls [`retract_polar!`](@ref ManifoldsBase.retract_polar!).
"""
function retract_polar(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_polar!(M, q, p, X)
end
"""
    retract_project(M::AbstractManifold, p, q)

computes the allocating variant of the [`ProjectionRetraction`](@ref),
which by default allocates and calls [`retract_project!`](@ref ManifoldsBase.retract_project!).
"""
function retract_project(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_project!(M, q, p, X)
end
"""
    retract_qr(M::AbstractManifold, p, q)

computes the allocating variant of the [`QRRetraction`](@ref),
which by default allocates and calls [`retract_qr!`](@ref ManifoldsBase.retract_qr!).
"""
function retract_qr(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_qr!(M, q, p, X)
end
"""
    retract_exp_ode(M::AbstractManifold, p, q, m::AbstractRetractionMethod, B::AbstractBasis)

computes the allocating variant of the [`ODEExponentialRetraction`](@ref)`(m,B)`,
which by default allocates and calls [`retract_exp_ode!`](@ref ManifoldsBase.retract_exp_ode!).
"""
function retract_exp_ode(
    M::AbstractManifold,
    p,
    X,
    m::AbstractRetractionMethod,
    B::AbstractBasis,
)
    q = allocate_result(M, retract, p, X)
    return retract_exp_ode!(M, q, p, X, m, B)
end
"""
    retract_softmax(M::AbstractManifold, p, q)

computes the allocating variant of the [`SoftmaxRetraction`](@ref),
which by default allocates and calls [`retract_softmax!`](@ref ManifoldsBase.retract_softmax!).
"""
function retract_softmax(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_softmax!(M, q, p, X)
end

"""
    retract_caley(M::AbstractManifold, p, q)

computes the allocating variant of the [`CayleyRetraction`](@ref),
which by default allocates and calls [`retract_caley!`](@ref ManifoldsBase.retract_caley!).
"""
function retract_caley(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_caley!(M, q, p, X)
end
"""
    retract_pade(M::AbstractManifold, p, q)

computes the allocating variant of the [`PadeRetraction`](@ref)`(n)`,
which by default allocates and calls [`retract_pade!`](@ref ManifoldsBase.retract_pade!).
"""
function retract_pade(M::AbstractManifold, p, X, n)
    q = allocate_result(M, retract, p, X)
    return retract_pade!(M, q, p, X, n)
end

Base.show(io::IO, ::CayleyRetraction) = print(io, "CayleyRetraction()")
Base.show(io::IO, ::PadeRetraction{m}) where {m} = print(io, "PadeRetraction($(m))")
