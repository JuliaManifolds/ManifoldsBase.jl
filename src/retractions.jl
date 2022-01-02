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

Compute a retraction by using the retraction of type T imn the embedding and projecting the result

# Constructor

    EmbeddedRetraction(r::AbstractRetractionMethod,)

Generate the retraction with retraction `r` to use in the embedding.
"""
struct EmbeddedRetraction{T<:AbstractRetractionMethod} <: AbstractRetractionMethod
    retraction::T
end

"""
    ExponentialRetraction

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

See [`solve_exp_ode`](@ref) for further details.

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

Compute an inverse retraction by using the inverse retraction of type T imn the embedding and projecting the result

# Constructor

    EmbeddedInverseRetraction(r::AbstractInverseRetractionMethod)

Generate the inverse retraction with inverse retraction `r` to use in the embedding.
"""
struct EmbeddedInverseRetraction{T<:AbstractInverseRetractionMethod} <:
       AbstractInverseRetractionMethod
    inverse_retraction::T
end


"""
    LogarithmicInverseRetraction

Inverse retraction using the [`log`](@ref)arithmic map.
"""
struct LogarithmicInverseRetraction <: AbstractInverseRetractionMethod end

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
function _inverse_retract!(M::AbstractManifold, X, p, q, m::EmbeddedInverseRetraction)
    return inverse_retract_embedded!(M, X, p, q, m.inverse_retraction)
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
# ToDo docu
function inverse_retract_embedded!(M::AbstractManifold, X, p, q, m)
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
function inverse_retract_softmax! end
function inverse_retract_qr! end
function inverse_retract_project! end
function inverse_retract_polar! end
function inverse_retract_nlsolve! end


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
_inverse_retract(M::AbstractManifold, p, q, ::LogarithmicInverseRetraction) = log(M, p, q)
function _inverse_retract(M::AbstractManifold, p, q, m::EmbeddedInverseRetraction)
    return inverse_retract_embedded(M, p, q, m.inverse_retraction)
end
function inverse_retract_embedded(M::AbstractManifold, p, q, m)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_embedded!(M, X, p, q, m)
end
function _inverse_retract(M::AbstractManifold, p, q, ::PolarInverseRetraction)
    return inverse_retract_polar(M, p, q)
end
function inverse_retract_polar(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_polar!(M, X, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, ::ProjectionInverseRetraction)
    return inverse_retract_project(M, p, q)
end
function inverse_retract_project(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_project!(M, X, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, ::QRInverseRetraction)
    return inverse_retract_qr(M, p, q)
end
function inverse_retract_qr(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_qr!(M, X, p, q)
end
function _inverse_retract(M::AbstractManifold, p, q, m::NLSolveInverseRetraction)
    return inverse_retract_nlsolve(M, p, q, m)
end
function inverse_retract_nlsolve(M::AbstractManifold, p, q, m)
    X = allocate_result(M, inverse_retract, p, q)
    return inverse_retract_nlsolve!(M, X, p, q, m)
end

@doc raw"""
    retract(M::AbstractManifold, p, X, method::AbstractRetractionMethod=default_retraction_method(M))
    retract(M::AbstractManifold, p, X, t::Real=1, method::AbstractRetractionMethod=default_retraction_method(M))

Compute a retraction, a cheaper, approximate version of the [`exp`](@ref)onential map,
from `p` into direction `X`, scaled by `t`, on the [`AbstractManifold`](@ref) `M`.

A retraction ``\operatorname{retr}_p: T_p\mathcal M → \mathcal M`` is a smooth map that fulfills

1. ``\operatorname{retr}_p(0) = p``
2. ``D\operatorname{retr}_p(0): T_p\mathcal M \to T_p\mathcal M`` is the identity map, i.e. ``D\operatorname{retr}_p(0)[X]=X``,

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
    return _retract(M, p, t * X, m)
end
function _retract(M::AbstractManifold, p, X, m::EmbeddedRetraction)
    return retract_embedding(M, p, X, m.retraction)
end
_retract(M::AbstractManifold, p, X, ::ExponentialRetraction) = exp(M, p, X)
function _retract(M::AbstractManifold, p, X, m::ODEExponentialRetraction)
    return retract_exp_ode(M, p, X, m.retraction, m.basis)
end
_retract(M::AbstractManifold, p, X, ::PolarRetraction) = retract_polar(M, p, X)
_retract(M::AbstractManifold, p, X, ::ProjectionRetraction) = retract_project(M, p, X)
_retract(M::AbstractManifold, p, X, ::QRRetraction) = retract_qr(M, p, X)
_retract(M::AbstractManifold, p, X, ::SoftmaxRetraction) = retract_softmax(M, p, X)
_retract(M::AbstractManifold, p, X, ::CayleyRetraction) = retract_pade(M, p, X, 1)
function _retract(M::AbstractManifold, p, X, ::PadeRetraction{n}) where {n}
    return retract_pade(M, p, X, n)
end
function retract_embedding(M::AbstractManifold, p, X, m)
    N = get_embedding(M)
    return project(M, retract(get_embedding(M), embed(M, p), embed(N, p, X), m))
end
function retract_polar(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_polar!(M, q, p, X)
end
function retract_project(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_project!(M, q, p, X)
end
function retract_qr(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_qr!(M, q, p, X)
end
function retract_exp_ode(M::AbstractManifold, p, X, m, B)
    q = allocate_result(M, retract, p, X)
    return retract_exp_ode!(M, q, p, X, m, B)
end
function retract_softmax(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    return retract_softmax!(M, q, p, X)
end
function retract_pade(M::AbstractManifold, p, X, n)
    q = allocate_result(M, retract, p, X)
    return retract_pade!(M, q, p, X, n)
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
    return _retract!(M, q, p, t * X, method)
end
# dispatch to lower level
function _retract!(M::AbstractManifold, q, p, X, m::EmbeddedRetraction)
    return retract_embedding!(M, q, p, X, m.retraction)
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
_retract!(M::AbstractManifold, q, p, X, ::CayleyRetraction) = retract_pade!(M, q, p, X, 1)
function _retract!(M::AbstractManifold, q, p, X, ::PadeRetraction{n}) where {n}
    return retract_pade!(M, q, p, X, n)
end

# ToDo - docu
function retract_embedding!(M::AbstractManifold, q, p, X, m)
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
function retract_pade! end
function retract_project! end
function retract_polar! end
function retract_qr! end
function retract_exp_ode! end
function retract_softmax! end

Base.show(io::IO, ::CayleyRetraction) = print(io, "CayleyRetraction()")
Base.show(io::IO, ::PadeRetraction{m}) where {m} = print(io, "PadeRetraction($(m))")
