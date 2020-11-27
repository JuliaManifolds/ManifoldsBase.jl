module ManifoldsBase

import Base:
    isapprox,
    exp,
    log,
    convert,
    copyto!,
    angle,
    eltype,
    isempty,
    length,
    similar,
    show,
    +,
    -,
    *
import LinearAlgebra: dot, norm, det, cross, I, UniformScaling, Diagonal

import Markdown: @doc_str
using LinearAlgebra

"""
    Manifold{F}

A manifold type. The `Manifold` is used to dispatch to different functions on a manifold,
usually as the first argument of the function. Examples are the [`exp`](@ref)onential and
[`log`](@ref)arithmic maps as well as more general functions that are built on them like the
[`geodesic`](@ref).

The manifold is parametrized by an [`AbstractNumbers`](@ref) to distinguish for example
real (ℝ) and complex (ℂ) manifolds.

For subtypes the preferred order of parameters is: size and simple value parameters,
followed by the [`AbstractNumbers`](@ref) `field`, followed by data type parameters,
which might depend on the abstract number field type.
"""
abstract type Manifold{𝔽} end

"""
    AbstractEstimationMethod

Abstract type for defining statistical estimation methods.
"""
abstract type AbstractEstimationMethod end

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
    ExponentialRetraction

Retraction using the exponential map.
"""
struct ExponentialRetraction <: AbstractRetractionMethod end

"""
    PolarRetraction <: AbstractRetractionMethod

Retractions that are based on singular value decompositions of the matrix / matrices
for point and tangent vector on a [`Manifold`](@ref)
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
matrix / matrices for point and tangent vector on a [`Manifold`](@ref)
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
matrix / matrices for point and tangent vector on a [`Manifold`](@ref)
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
matrix / matrices for point and tangent vector on a [`Manifold`](@ref)
"""
struct QRInverseRetraction <: AbstractInverseRetractionMethod end

"""
    MPoint

Type for a point on a manifold. While a [`Manifold`](@ref) does not necessarily require this
type, for example when it is implemented for `Vector`s or `Matrix` type elements, this type
can be used for more complicated representations, semantic verification, or even dispatch
for different representations of points on a manifold.
"""
abstract type MPoint end

"""
    OutOfInjectivityRadiusError

An error thrown when a function (for example [`log`](@ref)arithmic map or
[`inverse_retract`](@ref)) is given arguments outside of its [`injectivity_radius`](@ref).
"""
struct OutOfInjectivityRadiusError <: Exception end

"""
    TVector

Type for a tangent vector of a manifold. While a [`Manifold`](@ref) does not necessarily
require this type, for example when it is implemented for `Vector`s or `Matrix` type
elements, this type can be used for more complicated representations, semantic verification,
or even dispatch for different representations of tangent vectors and their types on a
manifold.
"""
abstract type TVector end

"""
    CoTVector

Type for a cotangent vector of a manifold. While a [`Manifold`](@ref) does not necessarily
require this type, for example when it is implemented for `Vector`s or `Matrix` type
elements, this type can be used for more complicated representations, semantic verification,
or even dispatch for different representations of cotangent vectors and their types on a
manifold.
"""
abstract type CoTVector end

"""
    allocate(a)
    allocate(a, dims::Integer...)
    allocate(a, dims::Tuple)
    allocate(a, T::Type)
    allocate(a, T::Type, dims::Integer...)
    allocate(a, T::Type, dims::Tuple)

Allocate an object similar to `a`. It is similar to function `similar`, although
instead of working only on the outermost layer of a nested structure, it maps recursively
through outer layers and calls `similar` on the innermost array-like object only.
Type `T` is the new number element type [`number_eltype`](@ref), if it is not given
the element type of `a` is retained. The `dims` argument can be given for non-nested
allocation and is forwarded to the function `similar`.
"""
allocate(a, args...)
allocate(a) = similar(a)
allocate(a, dims::Integer...) = similar(a, dims...)
allocate(a, dims::Tuple) = similar(a, dims)
allocate(a, T::Type) = similar(a, T)
allocate(a, T::Type, dims::Integer...) = similar(a, T, dims...)
allocate(a, T::Type, dims::Tuple) = similar(a, T, dims)
allocate(a::AbstractArray{<:AbstractArray}) = map(allocate, a)
allocate(a::AbstractArray{<:AbstractArray}, T::Type) = map(t -> allocate(t, T), a)
allocate(a::NTuple{N,AbstractArray} where {N}) = map(allocate, a)
allocate(a::NTuple{N,AbstractArray} where {N}, T::Type) = map(t -> allocate(t, T), a)

"""
    allocate_result(M::Manifold, f, x...)

Allocate an array for the result of function `f` on [`Manifold`](@ref) `M` and arguments
`x...` for implementing the non-modifying operation using the modifying operation.

Usefulness of passing a function is demonstrated by methods that allocate results of musical
isomorphisms.
"""
function allocate_result(M::Manifold, f, x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[1], T)
end

"""
    allocate_result_type(M::Manifold, f, args::NTuple{N,Any}) where N

Return type of element of the array that will represent the result of function `f` and the
[`Manifold`](@ref) `M` on given arguments `args` (passed as a tuple).
"""
function allocate_result_type(M::Manifold, f, args::NTuple{N,Any}) where {N}
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, args))
end

"""
    angle(M::Manifold, p, X, Y)

Compute the angle between tangent vectors `X` and `Y` at point `p` from the
[`Manifold`](@ref) `M` with respect to the inner product from [`inner`](@ref).
"""
angle(M::Manifold, p, X, Y) = acos(real(inner(M, p, X, Y)) / norm(M, p, X) / norm(M, p, Y))

"""
    base_manifold(M::Manifold, depth = Val(-1))

Return the internally stored [`Manifold`](@ref) for decorated manifold `M` and the base
manifold for vector bundles or power manifolds. The optional parameter `depth` can be used
to remove only the first `depth` many decorators and return the [`Manifold`](@ref) from that
level, whether its decorated or not. Any negative value deactivates this depth limit.
"""
base_manifold(M::Manifold, depth = Val(-1)) = M

"""
    check_manifold_point(M::Manifold, p; kwargs...) -> Union{Nothing,String}

Return `nothing` when `p` is a point on the [`Manifold`](@ref) `M`. Otherwise, return an
error with description why the point does not belong to manifold `M`.

By default, `check_manifold_point` returns `nothing`, i.e. if no checks are implemented, the
assumption is to be optimistic for a point not deriving from the [`MPoint`](@ref) type.
"""
check_manifold_point(M::Manifold, p; kwargs...) = nothing

"""
    check_tangent_vector(M::Manifold, p, X; kwargs...) -> Union{Nothing,String}

Check whether `X` is a valid tangent vector in the tangent space of `p` on the
[`Manifold`](@ref) `M`. An implementation should first call [`check_manifold_point(M, p;
kwargs...)`](@ref) and then validate `X`. If it is not a tangent vector, an error string
should be returned.

By default, `check_tangent_vector` returns `nothing`, i.e. if no checks are implemented, the
assumption is to be optimistic for tangent vectors not deriving from the [`TVector`](@ref)
type.
"""
check_tangent_vector(M::Manifold, p, X; kwargs...) = nothing

"""
    check_size(M::Manifold, p)
    check_size(M::Manifold, p, X)

Check whether `p` has the right [`representation_size`](@ref) for a [`Manifold`](@ref) `M`.
Additionally if a tangent vector is given, both `p` and `X` are checked to be of
corresponding correct representation sizes for points and tangent vectors on `M`.

By default, `check_size` returns `nothing`, i.e. if no checks are implemented, the
assumption is to be optimistic.
"""
function check_size(M::Manifold, p)
    n = size(p)
    m = representation_size(M)
    if length(n) != length(m)
        return DomainError(
            length(n),
            "The point $(p) can not belong to the manifold $(M), since its size $(n) is not equal to the manifolds representation size ($(m)).",
        )
    end
    if n != m
        return DomainError(
            n,
            "The point $(p) can not belong to the manifold $(M), since its size $(n) is not equal to the manifolds representation size ($(m)).",
        )
    end
end
function check_size(M::Manifold, p, X)
    mse = check_size(M, p)
    mse === nothing || return mse
    n = size(X)
    m = representation_size(M)
    if length(n) != length(m)
        return DomainError(
            length(n),
            "The tangent vector $(X) can not belong to the manifold $(M), since its size $(n) is not equal to the manifolds representation size ($(m)).",
        )
    end
    if n != m
        return DomainError(
            n,
            "The tangent vector $(X) can not belong to the manifold $(M), since its size $(n) is not equal to the manifodls representation size ($(m)).",
        )
    end
end

"""
    distance(M::Manifold, p, q)

Shortest distance between the points `p` and `q` on the [`Manifold`](@ref) `M`.
"""
distance(M::Manifold, p, q) = norm(M, p, log(M, p, q))

"""
    embed(M::Manifold, p)

Embed point `p` from the [`Manifold`](@ref) `M` into the ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `embed` includes changing data representation, if applicable, i.e.
if the points on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly.

See also: [`EmbeddedManifold`](@ref), [`project`](@ref project(M::Manifold,p))
"""
function embed(M::Manifold, p)
    q = allocate_result(M, embed, p)
    embed!(M, q, p)
    return q
end

"""
    embed!(M::Manifold, q, p)

Embed point `p` from the [`Manifold`](@ref) `M` into an ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Not implementing this function means, there is no proper embedding for your manifold.
Additionally, `embed` might include changing data representation, if applicable, i.e.
if points on `M` are not represented in the same way as their counterparts in the embedding,
the representation is changed accordingly.

If you have more than one embedding, see [`EmbeddedManifold`](@ref) for defining a second
embedding. If your point `p` is already represented in some embedding,
see [`AbstractEmbeddedManifold`](@ref) how you can avoid reimplementing code from the embedded manifold

See also: [`EmbeddedManifold`](@ref), [`project!`](@ref project!(M::Manifold, q, p))
"""
function embed!(M::Manifold, q, p)
    return error(manifold_function_not_implemented_message(M, embed!, q, p))
end

"""
    embed(M::Manifold, p, X)

Embed a tangent vector `X` at a point `p` on the [`Manifold`](@ref) `M` into an ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Not implementing this function means, there is no proper embedding for your tangent space(s).

Additionally, `embed` might include changing data representation, if applicable, i.e.
if tangent vectors on `M` are not represented in the same way as their counterparts in the
embedding, the representation is changed accordingly.

If you have more than one embedding, see [`EmbeddedManifold`](@ref) for defining a second
embedding. If your tangent vector `X` is already represented in some embedding,
see [`AbstractEmbeddedManifold`](@ref) how you can avoid reimplementing code from the embedded manifold

See also: [`EmbeddedManifold`](@ref), [`project`](@ref project(M::Manifold, p, X))
"""
function embed(M::Manifold, p, X)
    # Note that the order is switched,
    # since the allocation by default takes the type of the first.
    Y = allocate_result(M, embed, X, p)
    embed!(M, Y, p, X)
    return Y
end

"""
    embed!(M::Manifold, Y, p, X)

Embed a tangent vector `X` at a point `p` on the [`Manifold`](@ref) `M` into the ambient
space and return the result in `Y`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `embed!` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as tangents on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. The embedded tangents are then in
the tangent spaces of the embedded base points.

See also: [`EmbeddedManifold`](@ref), [`project!`](@ref project!(M::Manifold, Y, p, X))
"""
function embed!(M::Manifold, Y, p, X)
    return error(manifold_function_not_implemented_message(M, embed!, Y, p, X))
end

"""
    exp(M::Manifold, p, X)
    exp(M::Manifold, p, X, t::Real = 1)

Compute the exponential map of tangent vector `X`, optionally scaled by `t`,  at point `p`
from manifold the [`Manifold`](@ref) `M`.
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
from manifold the [`Manifold`](@ref) `M`.
The result is saved to `q`.
"""
function exp!(M::Manifold, q, p, X)
    return error(manifold_function_not_implemented_message(M, exp!, q, p, X))
end
exp!(M::Manifold, q, p, X, t::Real) = exp!(M, q, p, t * X)

"""
    geodesic(M::Manifold, p, X) -> Function

Get the geodesic with initial point `p` and velocity `X` on the [`Manifold`](@ref) `M`.
 The geodesic is the curve of constant velocity that is locally distance-minimizing.
 This function returns a function of (time) `t`.

    geodesic(M::Manifold, x, v, t::Real)
    geodesic(M::Manifold, x, v, T::AbstractVector) -> AbstractVector

Return the point at time `t` or points at times `t` in `T` along the geodesic.
"""
geodesic(M::Manifold, p, X) = t -> exp(M, p, X, t)
geodesic(M::Manifold, p, X, t::Real) = exp(M, p, X, t)
geodesic(M::Manifold, p, X, T::AbstractVector) = map(t -> exp(M, p, X, t), T)

@doc doc"""
    injectivity_radius(M::Manifold, p)

Return the distance $d$ such that [`exp(M, p, X)`](@ref exp(::Manifold, ::Any, ::Any)) is
injective for all tangent vectors shorter than $d$ (i.e. has an inverse).

    injectivity_radius(M::Manifold)

Infimum of the injectivity radius of all manifold points.

    injectivity_radius(M::Manifold[, x], method::AbstractRetractionMethod)
    injectivity_radius(M::Manifold, x, method::AbstractRetractionMethod)

Distance $d$ such that
[`retract(M, p, X, method)`](@ref retract(::Manifold, ::Any, ::Any, ::AbstractRetractionMethod))
is injective for all tangent vectors shorter than $d$ (i.e. has an inverse) for point `p`
if provided or all manifold points otherwise.
"""
function injectivity_radius(M::Manifold)
    return error(manifold_function_not_implemented_message(M, injectivity_radius))
end
injectivity_radius(M::Manifold, p) = injectivity_radius(M)
function injectivity_radius(M::Manifold, p, method::AbstractRetractionMethod)
    return injectivity_radius(M, method)
end
function injectivity_radius(M::Manifold, method::AbstractRetractionMethod)
    return error(manifold_function_not_implemented_message(M, injectivity_radius, method))
end
injectivity_radius(M::Manifold, p, ::ExponentialRetraction) = injectivity_radius(M, p)
injectivity_radius(M::Manifold, ::ExponentialRetraction) = injectivity_radius(M)

"""
    inner(M::Manifold, p, X, Y)

Compute the inner product of tangent vectors `X` and `Y` at point `p` from the
[`Manifold`](@ref) `M`.

See also: [`MetricManifold`](@ref Main.Manifolds.MetricManifold)
"""
function inner(M::Manifold, p, X, Y)
    return error(manifold_function_not_implemented_message(M, inner, p, X, Y))
end

"""
    inverse_retract!(M::Manifold, X, p, q[, method::AbstractInverseRetractionMethod])

Compute the inverse retraction, a cheaper, approximate version of the
[`log`](@ref)arithmic map), of points `p` and `q` on the [`Manifold`](@ref) `M`.
Result is saved to `X`.

Inverse retraction method can be specified by the last argument, defaulting to
[`LogarithmicInverseRetraction`](@ref). See the documentation of respective manifolds for
available methods.
"""
function inverse_retract!(M::Manifold, X, p, q)
    return inverse_retract!(M, X, p, q, LogarithmicInverseRetraction())
end
function inverse_retract!(M::Manifold, X, p, q, method::LogarithmicInverseRetraction)
    return log!(M, X, p, q)
end
function inverse_retract!(M::Manifold, X, p, q, method::AbstractInverseRetractionMethod)
    return error(manifold_function_not_implemented_message(
        M,
        inverse_retract!,
        X,
        p,
        q,
        method,
    ))
end

"""
    inverse_retract(M::Manifold, x, y)
    inverse_retract(M::Manifold, x, y, method::AbstractInverseRetractionMethod

Compute the inverse retraction, a cheaper, approximate version of the
[`log`](@ref)arithmic map), of points `p` and `q` on the [`Manifold`](@ref) `M`.

Inverse retraction method can be specified by the last argument, defaulting to
[`LogarithmicInverseRetraction`](@ref). See the documentation of respective manifolds for
available methods.
"""
function inverse_retract(M::Manifold, p, q)
    X = allocate_result(M, inverse_retract, p, q)
    inverse_retract!(M, X, p, q)
    return X
end
function inverse_retract(M::Manifold, p, q, method::AbstractInverseRetractionMethod)
    X = allocate_result(M, inverse_retract, p, q)
    inverse_retract!(M, X, p, q, method)
    return X
end

"""
    isapprox(M::Manifold, p, q; kwargs...)

Check if points `p` and `q` from [`Manifold`](@ref) `M` are approximately equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(M::Manifold, x, y; kwargs...) = isapprox(x, y; kwargs...)

"""
    isapprox(M::Manifold, p, X, Y; kwargs...)

Check if vectors `X` and `Y` tangent at `p` from [`Manifold`](@ref) `M` are approximately
equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(M::Manifold, p, X, Y; kwargs...) = isapprox(X, Y; kwargs...)


"""
    is_manifold_point(M::Manifold, p, throw_error = false; kwargs...)

Return whether `p` is a valid point on the [`Manifold`](@ref) `M`.

If `throw_error` is `false`, the function returns either `true` or `false`. If `throw_error`
is `true`, the function either returns `true` or throws an error. By default the function
calls [`check_manifold_point(M, p; kwargs...)`](@ref) and checks whether the returned value
is `nothing` or an error.
"""
function is_manifold_point(M::Manifold, p, throw_error = false; kwargs...)
    mpe = check_manifold_point(M, p; kwargs...)
    mpe === nothing && return true
    return throw_error ? throw(mpe) : false
end

"""
    is_tangent_vector(M::Manifold, p, X, throw_error = false; kwargs...)

Return whether `X` is a valid tangent vector at point `p` on the [`Manifold`](@ref) `M`.
Returns either `true` or `false`.

If `throw_error` is `false`, the function returns either `true` or `false`. If `throw_error`
is `true`, the function either returns `true` or throws an error. By default the function
calls [`check_tangent_vector(M, p, X; kwargs...)`](@ref) and checks whether the returned
value is `nothing` or an error.
"""
function is_tangent_vector(M::Manifold, p, X, throw_error = false; kwargs...)
    mtve = check_tangent_vector(M, p, X; kwargs...)
    mtve === nothing && return true
    return throw_error ? throw(mtve) : false
end

"""
    log(M::Manifold, p, q)

Compute the logarithmic map of point `q` at base point `p` on the [`Manifold`](@ref) `M`.
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
"""
function log!(M::Manifold, X, p, q)
    return error(manifold_function_not_implemented_message(M, log!, X, p, q))
end

@doc doc"""
    manifold_dimension(M::Manifold)

The dimension $n=\dim_{\mathcal M}$ of real space $\mathbb R^n$ to which the neighborhood of
each point of the [`Manifold`](@ref) `M` is homeomorphic.
"""
function manifold_dimension(M::Manifold)
    return error(manifold_function_not_implemented_message(M, manifold_dimension))
end

function manifold_function_not_implemented_message(M::Manifold, f, x...)
    s = join(map(string, map(typeof, x)), ", ", " and ")
    a = length(x) > 1 ? "arguments" : "argument"
    m = length(x) > 0 ? " for $(a) $(s)." : "."
    return "$(f) not implemented on $(M)$(m)"
end

"""
    mid_point(M::Manifold, p1, p2)

Calculate the middle between the two point `p1` and `p2` from manifold `M`.
By default uses [`log`](@ref), divides the vector by 2 and uses [`exp`](@ref).
"""
function mid_point(M::Manifold, p1, p2)
    q = allocate(p1)
    return mid_point!(M, q, p1, p2)
end

"""
    mid_point!(M::Manifold, q, p1, p2)

Calculate the middle between the two point `p1` and `p2` from manifold `M`.
By default uses [`log`](@ref), divides the vector by 2 and uses [`exp!`](@ref).
Saves the result in `q`.
"""
function mid_point!(M::Manifold, q, p1, p2)
    X = log(M, p1, p2)
    return exp!(M, q, p1, X / 2)
end

@static if VERSION <= v"1.1"
    function mid_point!(
        M::Manifold,
        q::AbstractArray{T1,0},
        p1::AbstractArray{T2,0},
        p2::AbstractArray{T3,0},
    ) where {T1,T2,T3}
        X = log(M, p1, p2)
        return exp!(M, q, p1, fill(X / 2))
    end
end

"""
    norm(M::Manifold, p, X)

Compute the norm of tangent vector `X` at point `p` from a [`Manifold`](@ref) `M`.
By default this is computed using [`inner`](@ref).
"""
norm(M::Manifold, p, X) = sqrt(real(inner(M, p, X, X)))

"""
    number_eltype(x)

Numeric element type of the a nested representation of a point or a vector.
To be used in conjuntion with [`allocate`](@ref) or [`allocate_result`](@ref).
"""
number_eltype(x) = eltype(x)
function number_eltype(x::AbstractArray)
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, x))
end
function number_eltype(x::Tuple)
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, x))
end

"""
    project(M::Manifold, p)

Project point `p` from the ambient space of the [`Manifold`](@ref) `M` to `M`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Additionally, the projection includes changing data representation, if applicable,
i.e. if the points on `M` are not represented in the same array data, the data is changed
accordingly.

See also: [`EmbeddedManifold`](@ref), [`embed`](@ref embed(M::Manifold, p))
"""
function project(M::Manifold, p)
    q = allocate_result(M, project, p)
    project!(M, q, p)
    return q
end

"""
    project!(M::Manifold, q, p)

Project point `p` from the ambient space onto the [`Manifold`](@ref) `M`. The result is
storedin `q`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Additionally, the projection includes changing data representation, if applicable,
i.e. if the points on `M` are not represented in the same array data, the data is changed
accordingly.

See also: [`EmbeddedManifold`](@ref), [`embed!`](@ref embed!(M::Manifold, q, p))
"""
function project!(M::Manifold, q, p)
    return error(manifold_function_not_implemented_message(M, project!, q, p))
end

"""
    project(M::Manifold, p, X)

Project ambient space representation of a vector `X` to a tangent vector at point `p` on
the [`Manifold`](@ref) `M`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `project` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. after projection the change to the
Lie algebra is perfomed, too.

See also: [`EmbeddedManifold`](@ref), [`embed`](@ref embed(M::Manifold, p, X))
"""
function project(M::Manifold, p, X)
    # Note that the order is switched,
    # since the allocation by default takes the type of the first.
    Y = allocate_result(M, project, X, p)
    project!(M, Y, p, X)
    return Y
end

"""
    project!(M::Manifold, Y, p, X)

Project ambient space representation of a vector `X` to a tangent vector at point `p` on
the [`Manifold`](@ref) `M`. The result is saved in vector `Y`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `project!` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. after projection the change to the
Lie algebra is perfomed, too.

See also: [`EmbeddedManifold`](@ref), [`embed!`](@ref embed!(M::Manifold, Y, p, X))
"""
function project!(M::Manifold, Y, p, X)
    return error(manifold_function_not_implemented_message(M, project!, Y, p, X))
end

@doc doc"""
    representation_size(M::Manifold)

The size of an array representing a point on [`Manifold`](@ref) `M`.
Returns `nothing` by default indicating that points are not represented using an
`AbstractArray`.
"""
function representation_size(M::Manifold)
    return nothing
end

"""
    retract(M::Manifold, p, X)
    retract(M::Manifold, p, X, t::Real=1)
    retract(M::Manifold, p, X, method::AbstractRetractionMethod)
    retract(M::Manifold, p, X, t::Real=1, method::AbstractRetractionMethod)

Compute a retraction, a cheaper, approximate version of the [`exp`](@ref)onential map,
from `p` into direction `X`, scaled by `t`, on the [`Manifold`](@ref) `M`.

Retraction method can be specified by the last argument, defaulting to
[`ExponentialRetraction`](@ref). See the documentation of respective manifolds for available
methods.
"""
function retract(M::Manifold, p, X)
    q = allocate_result(M, retract, p, X)
    retract!(M, q, p, X)
    return q
end
retract(M::Manifold, p, X, t::Real) = retract(M, p, t * X)
function retract(M::Manifold, p, X, method::AbstractRetractionMethod)
    q = allocate_result(M, retract, p, X)
    retract!(M, q, p, X, method)
    return q
end
function retract(M::Manifold, p, X, t::Real, method::AbstractRetractionMethod)
    return retract(M, p, t * X, method)
end

"""
    retract!(M::Manifold, q, p, X)
    retract!(M::Manifold, q, p, X, t::Real=1)
    retract!(M::Manifold, q, p, X, method::AbstractRetractionMethod)
    retract!(M::Manifold, q, p, X, t::Real=1, method::AbstractRetractionMethod)

Compute a retraction, a cheaper, approximate version of the [`exp`](@ref)onential map,
from `p` into direction `X`, scaled by `t`, on the [`Manifold`](@ref) manifold `M`.
Result is saved to `q`.

Retraction method can be specified by the last argument, defaulting to
[`ExponentialRetraction`](@ref). See the documentation of respective manifolds for available
methods.
"""
retract!(M::Manifold, q, p, X) = retract!(M, q, p, X, ExponentialRetraction())
retract!(M::Manifold, q, p, X, t::Real) = retract!(M, q, p, t * X)
retract!(M::Manifold, q, p, X, method::ExponentialRetraction) = exp!(M, q, p, X)
function retract!(M::Manifold, q, p, X, t::Real, method::AbstractRetractionMethod)
    return retract!(M, q, p, t * X, method)
end
function retract!(M::Manifold, q, p, X, method::AbstractRetractionMethod)
    return error(manifold_function_not_implemented_message(M, retract!, q, p, method))
end

@doc doc"""
    shortest_geodesic(M::Manifold, p, q) -> Function

Get a [`geodesic`](@ref) $\gamma_{p,q}(t)$ whose length is the shortest path between the
points `p`and `q`, where $\gamma_{p,q}(0)=p$ and $\gamma_{p,q}(1)=q$. When there are
multiple shortest geodesics, there is no guarantee which will be returned.

This function returns a function of time, which may be a `Real` or an `AbstractVector`.

    shortest_geodesic(M::Manifold, p, q, t::Real)
    shortest_geodesic(M::Manifold, p, q, T::AbstractVector) -> AbstractVector

Return the point at time `t` or points at times `t` in `T` along the shortest geodesic.
"""
shortest_geodesic(M::Manifold, p, q) = geodesic(M, p, log(M, p, q))
shortest_geodesic(M::Manifold, p, q, t::Real) = geodesic(M, p, log(M, p, q), t)
shortest_geodesic(M::Manifold, p, q, T::AbstractVector) = geodesic(M, p, log(M, p, q), T)

"""
    size_to_tuple(::Type{S}) where S<:Tuple

Converts a size given by `Tuple{N, M, ...}` into a tuple `(N, M, ...)`.
"""
Base.@pure size_to_tuple(::Type{S}) where {S<:Tuple} = tuple(S.parameters...)

"""
    zero_tangent_vector!(M::Manifold, X, p)

Save to `X` a vector such that retracting `X` to the [`Manifold`](@ref) `M` at `p`
produces `p`.
"""
zero_tangent_vector!(M::Manifold, X, p) = log!(M, X, p, p)

"""
    zero_tangent_vector(M::Manifold, p)

Return the tangent vector from the tangent space at `p` on the [`Manifold`](@ref) `M`, that
represents the zero vector, i.e. such that a retraction at `p` produces `p`.
"""
function zero_tangent_vector(M::Manifold, p)
    X = allocate_result(M, zero_tangent_vector, p)
    zero_tangent_vector!(M, X, p)
    return X
end
include("errors.jl")
include("numbers.jl")
include("vector_transport.jl")
include("DecoratorManifold.jl")
include("bases.jl")
include("ValidationManifold.jl")
include("EmbeddedManifold.jl")
include("DefaultManifold.jl")
include("PowerManifold.jl")

export Manifold, MPoint, TVector, CoTVector
export AbstractDecoratorManifold
export ValidationManifold, ValidationMPoint, ValidationTVector, ValidationCoTVector
export AbstractEmbeddingType,
    TransparentIsometricEmbedding, DefaultIsometricEmbeddingType, DefaultEmbeddingType
export AbstractEmbeddedManifold, EmbeddedManifold, TransparentIsometricEmbedding
export AbstractPowerManifold, PowerManifold
export AbstractPowerRepresentation, NestedPowerRepresentation

export OutOfInjectivityRadiusError

export AbstractRetractionMethod,
    ExponentialRetraction,
    QRRetraction,
    PolarRetraction,
    ProjectionRetraction,
    PowerRetraction,
    InversePowerRetraction,
    SchildsLadderTransport,
    PoleLadderTransport,
    PowerVectorTransport

export AbstractInverseRetractionMethod,
    LogarithmicInverseRetraction,
    QRInverseRetraction,
    PolarInverseRetraction,
    ProjectionInverseRetraction

export ParallelTransport, ProjectionTransport

export CachedBasis,
    DefaultBasis,
    DefaultOrthogonalBasis,
    DefaultOrthonormalBasis,
    DiagonalizingOrthonormalBasis,
    DefaultOrthonormalBasis,
    GramSchmidtOrthonormalBasis,
    ProjectedOrthonormalBasis

export CompositeManifoldError, ComponentManifoldError

export allocate,
    base_manifold,
    check_manifold_point,
    check_tangent_vector,
    check_size,
    distance,
    exp,
    exp!,
    embed,
    embed!,
    geodesic,
    get_basis,
    get_component,
    get_component!,
    get_coordinates,
    get_coordinates!,
    get_embedding,
    get_iterator,
    get_vector,
    get_vector!,
    get_vectors,
    hat,
    hat!,
    shortest_geodesic,
    injectivity_radius,
    inner,
    inverse_retract,
    inverse_retract!,
    isapprox,
    is_manifold_point,
    is_tangent_vector,
    isempty,
    length,
    log,
    log!,
    manifold_dimension,
    mid_point,
    mid_point!,
    norm,
    number_eltype,
    number_of_coordinates,
    number_system,
    power_dimensions,
    project,
    project!,
    real_dimension,
    representation_size,
    set_component!,
    show,
    retract,
    retract!,
    vector_transport_along,
    vector_transport_along!,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
    vector_transport_to!,
    vee,
    vee!,
    zero_tangent_vector,
    zero_tangent_vector!

end # module
