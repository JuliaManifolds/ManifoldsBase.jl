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

"""
    ExponentialRetraction

Retraction using the exponential map.
"""
struct ExponentialRetraction <: AbstractRetractionMethod end

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
    NLsolveInverseRetraction{T<:AbstractRetractionMethod,TV,TK} <:
        ApproximateInverseRetraction

An inverse retraction method for approximating the inverse of a retraction using `NLsolve`.

# Constructor

    NLsolveInverseRetraction(
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
struct NLsolveInverseRetraction{TR<:AbstractRetractionMethod,TV,TK} <:
       ApproximateInverseRetraction
    retraction::TR
    X0::TV
    project_tangent::Bool
    project_point::Bool
    nlsolve_kwargs::TK
    function NLsolveInverseRetraction(m, X0, project_point, project_tangent, nlsolve_kwargs)
        isdefined(ManifoldsBase, :NLsolve) ||
            @warn "To use NLsolveInverseRetraction, NLsolve must be loaded using `using NLsolve`."
        return new{typeof(m),typeof(X0),typeof(nlsolve_kwargs)}(
            m,
            X0,
            project_point,
            project_tangent,
            nlsolve_kwargs,
        )
    end
end
function NLsolveInverseRetraction(
    m,
    X0 = nothing;
    project_tangent::Bool = false,
    project_point::Bool = false,
    nlsolve_kwargs...,
)
    return NLsolveInverseRetraction(m, X0, project_point, project_tangent, nlsolve_kwargs)
end


"""
    inverse_retract!(M::AbstractManifold, X, p, q[, method::AbstractInverseRetractionMethod])

Compute the inverse retraction, a cheaper, approximate version of the
[`log`](@ref)arithmic map), of points `p` and `q` on the [`AbstractManifold`](@ref) `M`.
Result is saved to `X`.

Inverse retraction method can be specified by the last argument, defaulting to
[`LogarithmicInverseRetraction`](@ref). See the documentation of respective manifolds for
available methods.

See also [`retract!`](@ref).
"""
function inverse_retract!(M::AbstractManifold, X, p, q)
    return inverse_retract!(M, X, p, q, LogarithmicInverseRetraction())
end
function inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    method::LogarithmicInverseRetraction,
)
    return log!(M, X, p, q)
end
function inverse_retract!(
    M::AbstractManifold,
    X,
    p,
    q,
    method::AbstractInverseRetractionMethod,
)
    return error(
        manifold_function_not_implemented_message(M, inverse_retract!, X, p, q, method),
    )
end

"""
    inverse_retract(M::AbstractManifold, p, q)
    inverse_retract(M::AbstractManifold, p, q, method::AbstractInverseRetractionMethod

Compute the inverse retraction, a cheaper, approximate version of the
[`log`](@ref)arithmic map), of points `p` and `q` on the [`AbstractManifold`](@ref) `M`.

Inverse retraction method can be specified by the last argument, defaulting to
[`LogarithmicInverseRetraction`](@ref), since the [`log`](@ref)arithmic map is the inverse of a
retraction, namely the [`exp`](@ref)onential map.
For available inverse retractions on certain manifolds see the documentation on the
correspnonding manifold.

See also [`retract`](@ref).
"""
function inverse_retract(M::AbstractManifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    inverse_retract!(M, X, p, q)
    return X
end
function inverse_retract(M::AbstractManifold, p, q, method::AbstractInverseRetractionMethod)
    X = allocate_result(M, inverse_retract, p, q)
    inverse_retract!(M, X, p, q, method)
    return X
end

@doc raw"""
    retract(M::AbstractManifold, p, X)
    retract(M::AbstractManifold, p, X, t::Real=1)
    retract(M::AbstractManifold, p, X, method::AbstractRetractionMethod)
    retract(M::AbstractManifold, p, X, t::Real=1, method::AbstractRetractionMethod)

Compute a retraction, a cheaper, approximate version of the [`exp`](@ref)onential map,
from `p` into direction `X`, scaled by `t`, on the [`AbstractManifold`](@ref) `M`.

A retraction ``\operatorname{retr}_p: T_p\mathcal M → \mathcal M`` is a smooth map that fulfills

1. ``\operatorname{retr}_p(0) = p``
2. ``D\operatorname{retr}_p(0): T_p\mathcal M \to T_p\mathcal M`` is the identity map, i.e. ``D\operatorname{retr}_p(0)[X]=X``,

where ``D\operatorname{retr}_p`` denotes the differential of the retraction

The retraction is called of second order if for all ``X`` the curves ``c(t) = R_p(tX)``
have a zero acceleration at ``t=0``, i.e. ``c''(0) = 0``.

Retraction method can be specified by the last argument, defaulting to
[`ExponentialRetraction`](@ref), i.e. to use the [`exp`](@ref)onential map, which itself is
a retraction. For further available retractions see the documentation of respective manifolds.

Locally, the retraction is invertible. For the inverse operation, see [`inverse_retract`](@ref).
"""
function retract(M::AbstractManifold, p, X)
    q = allocate_result(M, retract, p, X)
    retract!(M, q, p, X)
    return q
end
retract(M::AbstractManifold, p, X, t::Real) = retract(M, p, t * X)
function retract(M::AbstractManifold, p, X, method::AbstractRetractionMethod)
    q = allocate_result(M, retract, p, X)
    retract!(M, q, p, X, method)
    return q
end
function retract(M::AbstractManifold, p, X, t::Real, method::AbstractRetractionMethod)
    return retract(M, p, t * X, method)
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
[`ExponentialRetraction`](@ref). See the documentation of respective manifolds for available
methods.

See [`retract`](@ref) for more details.
"""
retract!(M::AbstractManifold, q, p, X) = retract!(M, q, p, X, ExponentialRetraction())
retract!(M::AbstractManifold, q, p, X, t::Real) = retract!(M, q, p, t * X)
retract!(M::AbstractManifold, q, p, X, ::ExponentialRetraction) = exp!(M, q, p, X)
function retract!(M::AbstractManifold, q, p, X, t::Real, method::AbstractRetractionMethod)
    return retract!(M, q, p, t * X, method)
end
function retract!(M::AbstractManifold, q, p, X, method::AbstractRetractionMethod)
    return error(manifold_function_not_implemented_message(M, retract!, q, p, method))
end
