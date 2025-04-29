"""
    AbstractInverseRetractionMethod <: AbstractApproximationMethod

Abstract type for methods for inverting a retraction (see [`inverse_retract`](@ref)).
"""
abstract type AbstractInverseRetractionMethod <: AbstractApproximationMethod end

"""
    AbstractRetractionMethod <: AbstractApproximationMethod

Abstract type for methods for [`retract`](@ref)ing a tangent vector to a manifold.
"""
abstract type AbstractRetractionMethod <: AbstractApproximationMethod end

"""
    ApproximateInverseRetraction <: AbstractInverseRetractionMethod

An abstract type for representing approximate inverse retraction methods.

"""
abstract type ApproximateInverseRetraction <: AbstractInverseRetractionMethod end

"""
    ApproximateRetraction <: AbstractRetractionMethod

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
for point and tangent vectors.

!!! note "Technical Note"
    Though you would call e.g. [`retract`](@ref)`(M, p, X, PolarRetraction())`,
    to implement a polar retraction, define [`retract_polar!`](@ref)`(M, q, p, X, t)`
    for your manifold `M`.
"""
struct PolarRetraction <: AbstractRetractionMethod end

"""
    ProjectionRetraction <: AbstractRetractionMethod

Retractions that are based on projection and usually addition in the embedding.

!!! note "Technical Note"
    Though you would call e.g. [`retract`](@ref)`(M, p, X, ProjectionRetraction())`,
    to implement a projection retraction, define [`retract_project!`](@ref)`(M, q, p, X, t)` for your manifold `M`.
"""
struct ProjectionRetraction <: AbstractRetractionMethod end

"""
    QRRetraction <: AbstractRetractionMethod

Retractions that are based on a QR decomposition of the
matrix / matrices for point and tangent vector on a [`AbstractManifold`](@ref)

!!! note "Technical Note"
    Though you would call e.g. [`retract`](@ref)`(M, p, X, QRRetraction())`,
    to implement a QR retraction, define [`retract_qr!`](@ref)`(M, q, p, X, t)` for your manifold `M`.
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

@doc raw"""
    struct SasakiRetraction <: AbstractRetractionMethod end

Exponential map on [`TangentBundle`](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/vector_bundle.html#Manifolds.TangentBundle) computed via Euler integration as described
in [MuralidharanFletcher:2012](@cite). The system of equations for ``\gamma : ℝ \to T\mathcal M`` such that
``γ(1) = \exp_{p,X}(X_M, X_F)`` and ``γ(0)=(p, X)`` reads

```math
\dot{γ}(t) = (\dot{p}(t), \dot{X}(t)) = (R(X(t), \dot{X}(t))\dot{p}(t), 0)
```

where ``R`` is the Riemann curvature tensor (see [`riemann_tensor`](@ref)).

# Constructor

    SasakiRetraction(L::Int)

In this constructor `L` is the number of integration steps.
"""
struct SasakiRetraction <: AbstractRetractionMethod
    L::Int
end

"""
    SoftmaxRetraction <: AbstractRetractionMethod

Describes a retraction that is based on the softmax function.

!!! note "Technical Note"
    Though you would call e.g. [`retract`](@ref)`(M, p, X, SoftmaxRetraction())`,
    to implement a softmax retraction, define [`retract_softmax!`](@ref)`(M, q, p, X, t)` for your manifold `M`.
"""
struct SoftmaxRetraction <: AbstractRetractionMethod end


"""
    StabilizedRetraction <: AbstractRetractionMethod

A retraction wraps another retraction and projects the resulting point onto the manifold
for numerical stability.

# Constructor

    StabilizedRetraction(::AbstractRetractionMethod=ExponentialRetraction())
"""
struct StabilizedRetraction{TRM<:AbstractRetractionMethod} <: AbstractRetractionMethod
    retraction::TRM
end
StabilizedRetraction() = StabilizedRetraction(ExponentialRetraction())

@doc raw"""
    PadeRetraction{m} <: AbstractRetractionMethod

A retraction based on the Padé approximation of order ``m``

# Constructor
    PadeRetraction(m::Int)

!!! note "Technical Note"
    Though you would call e.g. [`retract`](@ref)`(M, p, X, PadeRetraction(m))`,
    to implement a Padé retraction, define [`retract_pade!`](@ref)`(M, q, p, X, t, m)` for your manifold `M`.
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

!!! note "Technical Note"
    Though you would call e.g. [`retract`](@ref)`(M, p, X, CayleyRetraction())`,
    to implement a caley retraction, define [`retract_cayley!`](@ref)`(M, q, p, X, t)` for your manifold `M`.
    By default both these functions fall back to calling a [`PadeRetraction`](@ref)`(1)`.
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
    PadeInverseRetraction{m} <: AbstractInverseRetractionMethod

An inverse retraction based on the Padé approximation of order ``m`` for the retraction.

!!! note "Technical Note"
    Though you would call e.g. [`inverse_retract`](@ref)`(M, p, q, PadeInverseRetraction(m))`,
    to implement an inverse Padé retraction, define [`inverse_retract_pade!`](@ref)`(M, X, p, q, m)` for your manifold `M`.
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

!!! note "Technical Note"
    Though you would call e.g. [`inverse_retract`](@ref)`(M, p, q, CayleyInverseRetraction())`,
    to implement an inverse caley retraction, define [`inverse_retract_cayley!`](@ref)`(M, X, p, q)` for your manifold `M`.
    By default both these functions fall back to calling a [`PadeInverseRetraction`](@ref)`(1)`.
"""
const CayleyInverseRetraction = PadeInverseRetraction{1}

"""
    PolarInverseRetraction <: AbstractInverseRetractionMethod

Inverse retractions that are based on a singular value decomposition of the
matrix / matrices for point and tangent vector on a [`AbstractManifold`](@ref)

!!! note "Technical Note"
    Though you would call e.g. [`inverse_retract`](@ref)`(M, p, q, PolarInverseRetraction())`,
    to implement an inverse polar retraction, define [`inverse_retract_polar!`](@ref)`(M, X, p, q)` for your manifold `M`.
"""
struct PolarInverseRetraction <: AbstractInverseRetractionMethod end

"""
    ProjectionInverseRetraction <: AbstractInverseRetractionMethod

Inverse retractions that are based on a projection (or its inversion).

!!! note "Technical Note"
    Though you would call e.g. [`inverse_retract`](@ref)`(M, p, q, ProjectionInverseRetraction())`,
    to implement an inverse projection retraction, define [`inverse_retract_project!`](@ref)`(M, X, p, q)` for your manifold `M`.
"""
struct ProjectionInverseRetraction <: AbstractInverseRetractionMethod end

"""
    QRInverseRetraction <: AbstractInverseRetractionMethod

Inverse retractions that are based on a QR decomposition of the
matrix / matrices for point and tangent vector on a [`AbstractManifold`](@ref)

!!! note "Technical Note"
    Though you would call e.g. [`inverse_retract`](@ref)`(M, p, q, QRInverseRetraction())`,
    to implement an inverse QR retraction, define [`inverse_retract_qr!`](@ref)`(M, X, p, q)` for your manifold `M`.
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
    InverseRetractionWithKeywords{R<:AbstractRetractionMethod,K} <: AbstractInverseRetractionMethod

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

!!! note "Technical Note"
    Though you would call e.g. [`inverse_retract`](@ref)`(M, p, q, SoftmaxInverseRetraction())`,
    to implement an inverse softmax retraction, define [`inverse_retract_softmax!`](@ref)`(M, X, p, q)` for your manifold `M`.
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
function _inverse_retract(M::AbstractManifold, p, q, ::LogarithmicInverseRetraction)
    return log(M, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, m::AbstractInverseRetractionMethod)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract!(M, X, p, q, m)
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

_doc_retract = raw"""
    retract(M::AbstractManifold, p, X, method::AbstractRetractionMethod=default_retraction_method(M, typeof(p)))
    retract!(M::AbstractManifold, q, p, X, method::AbstractRetractionMethod=default_retraction_method(M, typeof(p)))

Compute a retraction, an approximate version of the [`exp`](@ref)onential map,
from `p` into direction `X`, scaled by `t`, on the [`AbstractManifold`](@ref) `M`.
This can be computed in-place of `q`.

A retraction ``\operatorname{retr}_p: T_p\mathcal M → \mathcal M`` is a smooth map that fulfils

1. ``\operatorname{retr}_p(0) = p``
2. ``D\operatorname{retr}_p(0): T_p\mathcal M → T_p\mathcal M`` is the identity map,
i.e. ``D\operatorname{retr}_p(0)[X]=X`` holds for all ``X∈ T_p\mathcal M``,

where ``D\operatorname{retr}_p`` denotes the differential of the retraction

The retraction is called of second order if for all ``X`` the curves ``c(t) = R_p(tX)``
have a zero acceleration at ``t=0``, i.e. ``c''(0) = 0``.

Retraction method can be specified by the last argument, defaulting to
[`default_retraction_method`](@ref)`(M)`. For further available retractions see the documentation of respective manifolds.

Locally, the retraction is invertible. For the inverse operation, see [`inverse_retract`](@ref).
"""

@doc "$(_doc_retract)"
function retract(
    M::AbstractManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return _retract(M, p, X, m)
end
# Layer 2 – dispatch on the method
function _retract(M::AbstractManifold, p, X, ::ExponentialRetraction)
    return exp(M, p, X)
end
function _retract(M::AbstractManifold, p, X, m::AbstractRetractionMethod)
    q = allocate_result(M, retract, p, X)
    return retract!(M, q, p, X, m)
end

@doc "$(_doc_retract)"
function retract!(
    M::AbstractManifold,
    q,
    p,
    X,
    method::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return _retract!(M, q, p, X, method)
end

# Retract Layer 2 – dispatch on the method
function _retract!(M::AbstractManifold, q, p, X, ::CayleyRetraction; kwargs...)
    return retract_cayley!(M, q, p, X; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, m::EmbeddedRetraction; kwargs...)
    return retract_embedded!(M, q, p, X, m.retraction; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, ::ExponentialRetraction; kwargs...)
    return exp!(M, q, p, X; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, m::ODEExponentialRetraction; kwargs...)
    return retract_exp_ode!(M, q, p, X, m.retraction, m.basis; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, ::PolarRetraction; kwargs...)
    return retract_polar!(M, q, p, X; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, ::ProjectionRetraction; kwargs...)
    return retract_project!(M, q, p, X; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, ::QRRetraction; kwargs...)
    return retract_qr!(M, q, p, X; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, m::SasakiRetraction)
    return retract_sasaki!(M, q, p, X, m)
end
function _retract!(M::AbstractManifold, q, p, X, ::SoftmaxRetraction; kwargs...)
    return retract_softmax!(M, q, p, X; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, m::StabilizedRetraction; kwargs...)
    return retract_stabilized!(M, q, p, X, m; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, m::PadeRetraction; kwargs...)
    return retract_pade!(M, q, p, X, m; kwargs...)
end
function _retract!(M::AbstractManifold, q, p, X, m::RetractionWithKeywords; kwargs...)
    return _retract!(M, q, p, X, m.retraction; kwargs..., m.kwargs...)
end

_doc_retract_fused = """
    retract_fused(M::AbstractManifold, p, X, t::Number, method::AbstractRetractionMethod=default_retraction_method(M, typeof(p)))
    retract_fused!(M::AbstractManifold, q, p, X, t::Number, method::AbstractRetractionMethod=default_retraction_method(M, typeof(p)))

A variant of [`retract`](@ref) that performs retraction on the vector `X` scaled by `t`.
This can be faster in some cases compared to multiplying `X` by `t`, especially when
performing this for multiple values of `t`.
This can be computed in-place of `q`.

By default, this falls back to calling [`retract`](@ref) with `t*X`.

!!! note "Technical Note"
    This fallback is happening on the in-place variant in [Layer 3](@ref design-layer3).
    Hence implementing this performant variant requires to implement the corresponding
    third layer fused function, like for example `retract_polar_fused!`.
    The “non-fused” variant always also has to be implemented, but can then be just spefied
    to fallback to the fused variant. for example
    ```
    retract_polar!(M, q, p, X) = retract_polar_fused!(M, q, p, X, one(eltype(p)))
    ```
"""

@doc "$(_doc_retract_fused)"
function retract_fused(
    M::AbstractManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return _retract_fused(M, p, X, t, m)
end

function _retract_fused(M::AbstractManifold, p, X, t::Number, ::ExponentialRetraction)
    return exp_fused(M, p, X, t)
end
function _retract_fused(M::AbstractManifold, p, X, t::Number, m::AbstractRetractionMethod)
    q = allocate_result(M, retract, p, X)
    return retract_fused!(M, q, p, X, t, m)
end

@doc "$(_doc_retract_fused)"
function retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return _retract_fused!(M, q, p, X, t, m)
end
# Retract fused Layer 2
# Retract Layer 2 – dispatch on the method
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    ::CayleyRetraction;
    kwargs...,
)
    return retract_cayley_fused!(M, q, p, X, t; kwargs...)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::EmbeddedRetraction;
    kwargs...,
)
    return retract_embedded_fused!(M, q, p, X, t, m.retraction; kwargs...)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    ::ExponentialRetraction;
    kwargs...,
)
    return exp_fused!(M, q, p, X, t; kwargs...)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::ODEExponentialRetraction;
    kwargs...,
)
    return retract_exp_ode_fused!(M, q, p, X, t, m.retraction, m.basis; kwargs...)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    ::PolarRetraction;
    kwargs...,
)
    return retract_polar_fused!(M, q, p, X, t; kwargs...)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    ::ProjectionRetraction;
    kwargs...,
)
    return retract_project_fused!(M, q, p, X, t; kwargs...)
end
function _retract_fused!(M::AbstractManifold, q, p, X, t::Number, ::QRRetraction; kwargs...)
    return retract_qr_fused!(M, q, p, X, t; kwargs...)
end
function _retract_fused!(M::AbstractManifold, q, p, X, t::Number, m::SasakiRetraction)
    return retract_sasaki_fused!(M, q, p, X, t, m)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    ::SoftmaxRetraction;
    kwargs...,
)
    return retract_softmax_fused!(M, q, p, X, t; kwargs...)
end
function _retract_fused!(M::AbstractManifold, q, p, X, t::Number, m::StabilizedRetraction)
    return retract_stabilized_fused!(M, q, p, X, t, m)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::PadeRetraction;
    kwargs...,
)
    return retract_pade_fused!(M, q, p, X, t, m; kwargs...)
end
function _retract_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::RetractionWithKeywords;
    kwargs...,
)
    return _retract_fused!(M, q, p, X, t, m.retraction; kwargs..., m.kwargs...)
end

#
#
# retract and retract_fused layer 3

"""
    retract_embedded!(M::AbstractManifold, q, p, X, m::AbstractRetractionMethod)

Compute the in-place variant of the [`EmbeddedRetraction`](@ref) using
the [`AbstractRetractionMethod`](@ref) `m` in the embedding (see [`get_embedding`](@ref))
and projecting the result back.
"""
function retract_embedded!(
    M::AbstractManifold,
    q,
    p,
    X,
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
            m;
            kwargs...,
        ),
    )
end

"""
    retract_embedded_fused!(M::AbstractManifold, q, p, X, t::Number, m::AbstractRetractionMethod)

Compute the scaled variant of `retract_embedded!`.
"""
function retract_embedded_fused!(
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
        retract_fused(
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
    retract_cayley!(M::AbstractManifold, q, p, X)

Compute the in-place variant of the [`CayleyRetraction`](@ref),
which by default falls back to calling the first order [`PadeRetraction`](@ref).
"""
function retract_cayley!(M::AbstractManifold, q, p, X; kwargs...)
    return retract_pade!(M, q, p, X, PadeRetraction(1); kwargs...)
end

function retract_cayley_fused!(M::AbstractManifold, q, p, X, t::Number; kwargs...)
    return retract_cayley!(M, q, p, t * X; kwargs...)
end

function retract_exp_ode! end
"""
    retract_exp_ode!(M::AbstractManifold, q, p, X, m::AbstractRetractionMethod, B::AbstractBasis)

Compute the in-place variant of the [`ODEExponentialRetraction`](@ref)`(m, B)`.
"""
retract_exp_ode!(
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod,
    B::AbstractBasis,
)

function retract_exp_ode_fused!(M::AbstractManifold, q, p, X, t::Number, m, B)
    return retract_exp_ode!(M, q, p, t * X, m, B)
end

function retract_pade! end
"""
    retract_pade!(M::AbstractManifold, q, p, X, m::PadeRetraction)

Compute the in-place variant of the [`PadeRetraction`](@ref) `m`.
"""
retract_pade!(M::AbstractManifold, q, p, X, m::PadeRetraction)

function retract_pade_fused!(M::AbstractManifold, q, p, X, t::Number, m::PadeRetraction)
    return retract_pade!(M, q, p, t * X, m)
end


function retract_project! end
"""
    retract_project!(M::AbstractManifold, q, p, X)

Compute the in-place variant of the [`ProjectionRetraction`](@ref).
"""
retract_project!(M::AbstractManifold, q, p, X)

"""
    retract_project_fused!(M::AbstractManifold, q, p, X, t::Number)

Compute the in-place variant of the [`ProjectionRetraction`](@ref).
"""
function retract_project_fused!(M::AbstractManifold, q, p, X, t::Number)
    return retract_project!(M, q, p, t * X)
end


function retract_polar! end
"""
    retract_polar!(M::AbstractManifold, q, p, X)

Compute the in-place variant of the [`PolarRetraction`](@ref).
"""
retract_polar!(M::AbstractManifold, q, p, X)

function retract_polar_fused!(M::AbstractManifold, q, p, X, t::Number)
    return retract_polar!(M, q, p, t * X)
end

function retract_qr! end
"""
    retract_qr!(M::AbstractManifold, q, p, X)

Compute the in-place variant of the [`QRRetraction`](@ref).
"""
retract_qr!(M::AbstractManifold, q, p, X)

function retract_qr_fused!(M::AbstractManifold, q, p, X, t::Number)
    return retract_qr!(M, q, p, t * X)
end

function retract_softmax! end
"""
    retract_softmax!(M::AbstractManifold, q, p, X)

Compute the in-place variant of the [`SoftmaxRetraction`](@ref).
"""
retract_softmax!(M::AbstractManifold, q, p, X)

function retract_softmax_fused!(M::AbstractManifold, q, p, X, t::Number)
    return retract_softmax!(M, q, p, t * X)
end

function retract_sasaki! end
"""
    retract_sasaki!(M::AbstractManifold, q, p, X, m::SasakiRetraction)

Compute the in-place variant of the [`SasakiRetraction`](@ref) `m`.
"""
retract_sasaki!(M::AbstractManifold, q, p, X, m::SasakiRetraction)

function retract_sasaki_fused!(M::AbstractManifold, q, p, X, t::Number, m::SasakiRetraction)
    return retract_sasaki!(M, q, p, t * X, m)
end

function retract_stabilized!(M::AbstractManifold, q, p, X, m::StabilizedRetraction)
    retract!(M, q, p, X, m.retraction)
    return embed_project!(M, q, q)
end

function retract_stabilized_fused!(
    M::AbstractManifold,
    q,
    p,
    X,
    t::Number,
    m::StabilizedRetraction,
)
    retract_fused!(M, q, p, X, t, m.retraction)
    return embed_project!(M, q, q)
end

Base.show(io::IO, ::CayleyRetraction) = print(io, "CayleyRetraction()")
Base.show(io::IO, ::PadeRetraction{m}) where {m} = print(io, "PadeRetraction($m)")

#
# default estimation methods pass down with and without the point type
function default_approximation_method(M::AbstractManifold, ::typeof(inverse_retract))
    return default_inverse_retraction_method(M)
end
function default_approximation_method(M::AbstractManifold, ::typeof(inverse_retract), T)
    return default_inverse_retraction_method(M, T)
end
function default_approximation_method(M::AbstractManifold, ::typeof(retract))
    return default_retraction_method(M)
end
function default_approximation_method(M::AbstractManifold, ::typeof(retract), T)
    return default_retraction_method(M, T)
end
