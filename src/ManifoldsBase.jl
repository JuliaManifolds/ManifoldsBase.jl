module ManifoldsBase

import Base: isapprox, exp, log, convert, copyto!, angle, eltype, similar, +, -, *
import LinearAlgebra: dot, norm, det, cross, I, UniformScaling, Diagonal

import Markdown: @doc_str
using LinearAlgebra

"""
    Manifold

A manifold type. The `Manifold` is used to dispatch to different functions on a
manifold, usually as the first argument of the function. Examples are the
[`exp`](@ref)onential and [`log`](@ref)arithmic maps as well as more general
functions that are built on them like the [`geodesic`](@ref).
"""
abstract type Manifold end

"""
    MPoint

Type for a point on a manifold. While a [`Manifold`](@ref) does not necessarily
require this type, for example when it is implemented for `Vector`s or `Matrix`
type elements, this type can be used for more complicated representations,
semantic verification, or even dispatch for different representations of points
on a manifold.
"""
abstract type MPoint end

"""
    TVector

Type for a tangent vector of a manifold. While a [`Manifold`](@ref) does not
necessarily require this type, for example when it is implemented for `Vector`s
or `Matrix` type elements, this type can be used for more complicated
representations, semantic verification, or even dispatch for different
representations of tangent vectors and their types on a manifold.
"""
abstract type TVector end

"""
    CoTVector

Type for a cotangent vector of a manifold. While a [`Manifold`](@ref) does not
necessarily require this type, for example when it is implemented for `Vector`s
or `Matrix` type elements, this type can be used for more complicated
representations, semantic verification, or even dispatch for different
representations of cotangent vectors and their types on a manifold.
"""
abstract type CoTVector end

"""
    is_decorator_manifold(M::Manifold)

Indicate whether a manifold is a decorator manifold, i.e. whether it
encapsulates a [`Manifold`](@ref) with additional features and stores internally
the original manifold instance. An example is the [`ArrayManifold`](@ref).

Certain functions are just calling themselves on the internal manifold and hence
do not need to be reimplemented for decorators again, for example
[`manifold_dimension`](@ref) and especially [`base_manifold`](@ref).

It is assumed that the undecorated (base) manifold is stored in `M.manifold`.
Alternatively, overload [`base_manifold`](@ref).
"""
is_decorator_manifold(::Manifold) = Val(false)

"""
    base_manifold(M::Manifold)

Return the internally stored manifold for decorated manifolds and the base
manifold for vector bundles or power manifolds.
"""
base_manifold(M::Manifold) = base_manifold(M, is_decorator_manifold(M))
base_manifold(M::Manifold, ::Val{true}) = base_manifold(M.manifold)
base_manifold(M::Manifold, ::Val{false}) = M

@doc doc"""
    representation_size(M::Manifold)

The size of an array representing a point on manifold `M`.
"""
representation_size(M::Manifold) = representation_size(M, is_decorator_manifold(M))
representation_size(M::Manifold, ::Val{true}) = representation_size(base_manifold(M))

function representation_size(M::Manifold, ::Val{false})
    error("representation_size not implemented for manifold $(typeof(M)).")
end

@doc doc"""
    manifold_dimension(M::Manifold)

The dimension $n$ of real space $\mathbb R^n$ to which the neighborhood of each
point of the manifold is homeomorphic.
"""
manifold_dimension(M::Manifold) = manifold_dimension(M, is_decorator_manifold(M))
manifold_dimension(M::Manifold, ::Val{true}) = manifold_dimension(base_manifold(M))

function manifold_dimension(M::Manifold, ::Val{false})
    error("manifold_dimension not implemented for manifold $(typeof(M)).")
end

"""
    isapprox(M::Manifold, x, y; kwargs...)

Check if points `x` and `y` from manifold `M` are approximately equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(M::Manifold, x, y; kwargs...) = isapprox(x, y; kwargs...)

"""
    isapprox(M::Manifold, x, v, w; kwargs...)

Check if vectors `v` and `w` tangent at `x` from manifold `M` are approximately
equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(M::Manifold, x, v, w; kwargs...) = isapprox(v, w; kwargs...)

"""
    OutOfInjectivityRadiusError

An error thrown when a function (for example [`log`](@ref)arithmic map or
[`inverse_retract`](@ref)) is given arguments outside of its
[`injectivity_radius`](@ref).
"""
struct OutOfInjectivityRadiusError <: Exception end

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
    retract!(M::Manifold, y, x, v[, t::Real=1], method::AbstractRetractionMethod=ExponentialRetraction())

Retraction (cheaper, approximate version of [`exp`](@ref)onential map) of
tangent vector `t*v` at point `x` from manifold `M`. Result is saved to `y`.

Retraction method can be specified by the last argument. Please look at the
documentation of respective manifolds for available methods.
"""
retract!(M::Manifold, y, x, v, method::ExponentialRetraction) = exp!(M, y, x, v)
retract!(M::Manifold, y, x, v) = retract!(M, y, x, v, ExponentialRetraction())
retract!(M::Manifold, y, x, v, t::Real) = retract!(M, y, x, t * v)

function retract!(M::Manifold, y, x, v, t::Real, method::AbstractRetractionMethod)
    return retract!(M, y, x, t * v, method)
end

"""
    retract(M::Manifold, x, v[, t::Real=1], method::AbstractRetractionMethod=ExponentialRetraction())

Retraction (cheaper, approximate version of [`exp`](@ref)onential map) of tangent
vector `t*v` at point `x` from manifold `M`.
"""
function retract(M::Manifold, x, v, method::AbstractRetractionMethod)
    xr = similar_result(M, retract, x, v)
    retract!(M, xr, x, v, method)
    return xr
end

function retract(M::Manifold, x, v)
    xr = similar_result(M, retract, x, v)
    retract!(M, xr, x, v)
    return xr
end

retract(M::Manifold, x, v, t::Real) = retract(M, x, t * v)

function retract(M::Manifold, x, v, t::Real, method::AbstractRetractionMethod)
    return retract(M, x, t * v, method)
end

"""
    AbstractInverseRetractionMethod

Abstract type for methods for inverting a retraction
(see [`inverse_retract`](@ref)).
"""
abstract type AbstractInverseRetractionMethod end

"""
    LogarithmicInverseRetraction

Inverse retraction using the [`log`](@ref)arithmic map.
"""
struct LogarithmicInverseRetraction <: AbstractInverseRetractionMethod end

"""
    inverse_retract!(M::Manifold, v, x, y)
    inverse_retract!(M::Manifold, v, x, y, method::AbstractInverseRetractionMethod)

Inverse retraction (cheaper, approximate version of [`log`](@ref)arithmic map)
of points `x` and `y`. Result is saved to `v`.

Inverse retraction method can be specified by the last argument, defaulting to
[`LogarithmicInverseRetraction`](@ref). Please look at the documentation of
respective manifolds for available methods.
"""
function inverse_retract!(M::Manifold, v, x, y, method::LogarithmicInverseRetraction)
    return log!(M, v, x, y)
end

function inverse_retract!(M::Manifold, v, x, y)
    return inverse_retract!(M, v, x, y, LogarithmicInverseRetraction())
end

"""
    inverse_retract(M::Manifold, x, y)
    inverse_retract(M::Manifold, x, y, method::AbstractInverseRetractionMethod

Inverse retraction (cheaper, approximate version of [`log`](@ref)arithmic map)
of points `x` and `y` from manifold `M`.

Inverse retraction method can be specified by the last argument, defaulting to
[`LogarithmicInverseRetraction`](@ref). Please look at the documentation of
respective manifolds for available methods.
"""
function inverse_retract(M::Manifold, x, y, method::AbstractInverseRetractionMethod)
    vr = similar_result(M, inverse_retract, x, y)
    inverse_retract!(M, vr, x, y, method)
    return vr
end

function inverse_retract(M::Manifold, x, y)
    vr = similar_result(M, inverse_retract, x, y)
    inverse_retract!(M, vr, x, y)
    return vr
end

"""
    project_point!(M::Manifold, y, x)

Project point `x` from the ambient space onto the manifold `M`. The point `y` is
overwritten by the projection. The function works only for selected embedded
manifolds and is *not* required to return the closest point.
"""
function project_point!(M::Manifold, y, x)
    error("project_point! not implemented for a $(typeof(M)) and points $(typeof(y)) and $(typeof(x)).")
end

"""
    project_point(M::Manifold, x)

Project point from the ambient space onto the manifold `M`. The point `x` is not
modified. The function works only for selected embedded manifolds and is *not*
required to return the closest point.
"""
function project_point(M::Manifold, x)
    y = similar_result(M, project_point, x)
    project_point!(M, y, x)
    return y
end

"""
    project_tangent!(M::Manifold, w, x, v)

Project ambient space representation of a vector `v` to a tangent vector
at point `x` from the manifold `M`. The result is saved in vector `w`.

The function works only for selected embedded manifolds and
is *not* required to return the closest vector.
"""
function project_tangent!(M::Manifold, w, x, v)
    error("project_tangent! not implemented for a $(typeof(M)) and point $(typeof(x)) with input $(typeof(v)).")
end

"""
    project_tangent(M::Manifold, x, v)

Project ambient space representation of a vector `v` to a tangent vector
at point `x` from the manifold `M`.

The function works only for selected embedded manifolds and
is *not* required to return the closest vector.
"""
function project_tangent(M::Manifold, x, v)
    vt = similar_result(M, project_tangent, v, x)
    project_tangent!(M, vt, x, v)
    return vt
end

"""
    inner(M::Manifold, x, v, w)

Inner product of tangent vectors `v` and `w` at point `x` from manifold `M`.
"""
function inner(M::Manifold, x, v, w)
    error("inner not implemented on a $(typeof(M)) for input point $(typeof(x)) and tangent vectors $(typeof(v)) and $(typeof(w)).")
end

"""
    norm(M::Manifold, x, v)

Norm of tangent vector `v` at point `x` from manifold `M`.
"""
norm(M::Manifold, x, v) = sqrt(inner(M, x, v, v))

"""
    distance(M::Manifold, x, y)

Shortest distance between the points `x` and `y` on manifold `M`.
"""
distance(M::Manifold, x, y) = norm(M, x, log(M, x, y))

"""
    angle(M::Manifold, x, v, w)

Angle between tangent vectors `v` and `w` at point `x` from manifold `M`.
"""
angle(M::Manifold, x, v, w) = acos(inner(M, x, v, w) / norm(M, x, v) / norm(M, x, w))

"""
    exp!(M::Manifold, y, x, v, t::Real = 1)

Exponential map of tangent vector `t*v` at point `x` from manifold `M`.
Result is saved to `y`.
"""
exp!(M::Manifold, y, x, v, t::Real) = exp!(M, y, x, t * v)

function exp!(M::Manifold, y, x, v)
    error("exp! not implemented on a $(typeof(M)) for input point $(x) and tangent vector $(v).")
end

"""
    exp(M::Manifold, x, v, t::Real = 1)

Exponential map of tangent vector `t*v` at point `x` from manifold `M`.
"""
function exp(M::Manifold, x, v)
    y = similar_result(M, exp, x, v)
    exp!(M, y, x, v)
    return y
end

exp(M::Manifold, x, v, t::Real) = exp(M, x, t * v)

"""
    exp(M::Manifold, x, v, T::AbstractVector) -> AbstractVector

Exponential map of tangent vector `t*v` at point `x` from manifold `M` for
each `t` in `T`.
"""
exp(M::Manifold, x, v, T::AbstractVector) = map(geodesic(M, x, v), T)

"""
    log!(M::Manifold, v, x, y)

Logarithmic map of point `y` at base point `x` on Manifold `M`. Result is saved
to `v`.
"""
function log!(M::Manifold, v, x, y)
    error("log! not implemented on $(typeof(M)) for points $(typeof(x)) and $(typeof(y))")
end

"""
    log(M::Manifold, x, y)

Logarithmic map of point `y` at base point `x` on Manifold `M`.
"""
function log(M::Manifold, x, y)
    v = similar_result(M, log, x, y)
    log!(M, v, x, y)
    return v
end

"""
    geodesic(M::Manifold, x, v) -> Function

Get the geodesic with initial point `x` and velocity `v`. The geodesic is the
curve of constant velocity that is locally distance-minimizing. This function
returns a function of time, which may be a `Real` or an `AbstractVector`.
"""
geodesic(M::Manifold, x, v) = t -> exp(M, x, v, t)

"""
    geodesic(M::Manifold, x, v, t::Real)

Get the point at time `t` traveling from `x` along the geodesic with initial
point `x` and velocity `v`.
"""
geodesic(M::Manifold, x, v, t::Real) = exp(M, x, v, t)

"""
    geodesic(M::Manifold, x, v, T::AbstractVector) -> AbstractVector

Get the points for each `t` in `T` traveling from `x` along the geodesic with
initial point `x` and velocity `v`.
"""
geodesic(M::Manifold, x, v, T::AbstractVector) = exp(M, x, v, T)

"""
    shortest_geodesic(M::Manifold, x, y) -> Function

Get a [`geodesic`](@ref) with initial point `x` and point `y` at `t=1` whose
length is the shortest path between the two points. When there are multiple
shortest geodesics, there is no guarantee which will be returned. This function
returns a function of time, which may be a `Real` or an `AbstractVector`.
"""
shortest_geodesic(M::Manifold, x, y) = geodesic(M, x, log(M, x, y))

"""
    shortest_geodesic(M::Manifold, x, y, t::Real)

Get the point at time `t` traveling from `x` along a shortest [`geodesic`](@ref)
connecting `x` and `y`, where `y` is reached at `t=1`.
"""
shortest_geodesic(M::Manifold, x, y, t::Real) = geodesic(M, x, log(M, x, y), t)

"""
    shortest_geodesic(M::Manifold, x, y, T::AbstractVector) -> AbstractVector

Get the points for each `t` in `T` traveling from `x` along a shortest
[`geodesic`](@ref) connecting `x` and `y`, where `y` is reached at `t=1`.
"""
shortest_geodesic(M::Manifold, x, y, T::AbstractVector) = geodesic(M, x, log(M, x, y), T)

"""
    AbstractVectorTransportMethod

Abstract type for methods for transporting vectors.
"""
abstract type AbstractVectorTransportMethod end

"""
    ParallelTransport <: AbstractVectorTransportMethod

Specify to use parallel transport as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref).
"""
struct ParallelTransport <: AbstractVectorTransportMethod end

"""
    ProjectionTransport <: AbstractVectorTransportMethod

Specify to use projection onto tangent space as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref). See [`project_tangent`](@ref) for details.
"""
struct ProjectionTransport <: AbstractVectorTransportMethod end

"""
    vector_transport_to!(M::Manifold, vto, x, v, y)
    vector_transport_to!(M::Manifold, vto, x, v, y, method::AbstractVectorTransportMethod)

Vector transport of vector `v` at point `x` to point `y`. The result is saved to
`vto`. By default, the `method` is [`ParallelTransport`](@ref).
"""
function vector_transport_to!(M::Manifold, vto, x, v, y)
    return vector_transport_to!(M, vto, x, v, y, ParallelTransport())
end

"""
    vector_transport_to!(M::Manifold, vto, x, v, y, method::ProjectionTransport)

Transport a vector `v` in the tangent space at `x` on a [`Manifold`](@ref) `M`
by interpreting it as an element of the embedding and then projecting it onto
the tangent space at `y`.
"""
function vector_transport_to!(M::Manifold, vto, x, v, y, ::ProjectionTransport)
    return project_tangent!(M, vto, y, v)
end

function vector_transport_to!(
    M::Manifold,
    vto,
    x,
    v,
    y,
    method::AbstractVectorTransportMethod,
)
    error("vector_transport_to! not implemented from a point of type $(typeof(x)) to a type $(typeof(y)) on a $(typeof(M)) for a vector of type $(v) and the $(typeof(method)).")
end

"""
    vector_transport_to(M::Manifold, x, v, y)
    vector_transport_to(M::Manifold, x, v, y, method::AbstractVectorTransportMethod)

Transport a vector `v` at point `x` to point `y` using the `method`, which
defaults to [`ParallelTransport`](@ref).
"""
function vector_transport_to(M::Manifold, x, v, y)
    return vector_transport_to(M, x, v, y, ParallelTransport())
end

function vector_transport_to(M::Manifold, x, v, y, method::AbstractVectorTransportMethod)
    vto = similar_result(M, vector_transport_to, v, x, y)
    vector_transport_to!(M, vto, x, v, y, method)
    return vto
end

"""
    vector_transport_direction!(M::Manifold, vto, x, v, vdir)
    vector_transport_direction!(M::Manifold, vto, x, v, vdir, method::AbstractVectorTransportMethod)

Transport a vector `v` at point `x` in the direction indicated by the tangent
vector `vdir` at point `x`. The result is saved to `vto`. By default,
[`exp`](@ref) and [`vector_transport_to!`](@ref) are used with the `method`,
which defaults to [`ParallelTransport`](@ref).
"""
function vector_transport_direction!(M::Manifold, vto, x, v, vdir)
    return vector_transport_direction!(M, vto, x, v, vdir, ParallelTransport())
end

function vector_transport_direction!(
    M::Manifold,
    vto,
    x,
    v,
    vdir,
    method::AbstractVectorTransportMethod,
)
    y = exp(M, x, vdir)
    return vector_transport_to!(M, vto, x, v, y, method)
end

"""
    vector_transport_direction(M::Manifold, x, v, vdir)
    vector_transport_direction(M::Manifold, x, v, vdir, method::AbstractVectorTransportMethod)

Transport a vector `v` at point `x` in the direction indicated by the tangent
vector `vdir` at point `x` using the `method`, which defaults to
[`ParallelTransport`](@ref).
"""
function vector_transport_direction(M::Manifold, x, v, vdir)
    return vector_transport_direction(M, x, v, vdir, ParallelTransport())
end

function vector_transport_direction(
    M::Manifold,
    x,
    v,
    vdir,
    method::AbstractVectorTransportMethod,
)
    vto = similar_result(M, vector_transport_direction, v, x, vdir)
    vector_transport_direction!(M, vto, x, v, vdir, method)
    return vto
end

"""
    vector_transport_along!(M::Manifold, vto, x, v, c)
    vector_transport_along!(M::Manifold, vto, x, v, c, method::AbstractVectorTransportMethod)

Transport a vector `v` at point `x` along the curve `c` such that `c(0)` is
equal to `x` to point `c(1)` using the `method`, which defaults to
[`ParallelTransport`](@ref). The result is saved to `vto`.
"""
function vector_transport_along!(M::Manifold, vto, x, v, c)
    return vector_transport_along!(M, vto, x, v, c, ParallelTransport())
end

function vector_transport_along!(
    M::Manifold,
    vto,
    x,
    v,
    c,
    method::AbstractVectorTransportMethod,
)
    error("vector_transport_along! not implemented for manifold $(typeof(M)), vector $(typeof(vto)), point $(typeof(x)), vector $(typeof(v)) along curve $(typeof(c)) with method $(typeof(method)).")
end

"""
    vector_transport_along(M::Manifold, x, v, c)
    vector_transport_along(M::Manifold, x, v, c, method::AbstractVectorTransportMethod)

Transport a vector `v` at point `x` along the curve `c` such that `c(0)` is
equal to `x` to point `c(1)`. The default `method` used is
[`ParallelTransport`](@ref).
"""
function vector_transport_along(M::Manifold, x, v, c)
    return vector_transport_along(M, x, v, c, ParallelTransport())
end

function vector_transport_along(M::Manifold, x, v, c, m::AbstractVectorTransportMethod)
    vto = similar_result(M, vector_transport_along, v, x)
    vector_transport_along!(M, vto, x, v, c, m)
    return vto
end

@doc doc"""
    injectivity_radius(M::Manifold, x)

Distance $d$ such that [`exp(M, x, v)`](@ref exp(::Manifold, ::Any, ::Any)) is
injective for all tangent vectors shorter than $d$ (i.e. has a left inverse).
"""
injectivity_radius(M::Manifold, x) = injectivity_radius(M)

@doc doc"""
    injectivity_radius(M::Manifold, x, method::AbstractRetractionMethod)

Distance $d$ such that
[`retract(M, x, v, method)`](@ref retract(::Manifold, ::Any, ::Any, ::AbstractRetractionMethod))
is injective for all tangent vectors shorter than $d$ (i.e. has a left inverse).
"""
injectivity_radius(M::Manifold, x, ::AbstractRetractionMethod) = injectivity_radius(M, x)

"""
    injectivity_radius(M::Manifold)

Infimum of the [`injectivity_radius`](@ref injectivity_radius(::Manifold, ::Any))
of all manifold points.
"""
injectivity_radius(M::Manifold) = Inf

"""
    zero_tangent_vector(M::Manifold, x)

Vector `v` such that retracting `v` to manifold `M` at `x` produces `x`.
"""
function zero_tangent_vector(M::Manifold, x)
    v = similar_result(M, zero_tangent_vector, x)
    zero_tangent_vector!(M, v, x)
    return v
end

"""
    zero_tangent_vector!(M::Manifold, v, x)

Save to `v` a vector such that retracting `v` to manifold `M` at `x` produces `x`.
"""
zero_tangent_vector!(M::Manifold, v, x) = log!(M, v, x, x)

"""
    similar_result_type(M::Manifold, f, args::NTuple{N,Any}) where N

Return type of element of the array that will represent the result of function
`f` for manifold `M` on given arguments `args` (passed as a tuple).
"""
function similar_result_type(M::Manifold, f, args::NTuple{N,Any}) where {N}
    T = typeof(reduce(+, one(eltype(eti)) for eti ∈ args))
    return T
end

"""
    similar_result(M::Manifold, f, x...)

Allocate an array for the result of function `f` on manifold `M` and arguments
`x...` for implementing the non-modifying operation using the modifying
operation.

Usefulness of passing a function is demonstrated by methods that allocate
results of musical isomorphisms.
"""
function similar_result(M::Manifold, f, x...)
    T = similar_result_type(M, f, x)
    return similar(x[1], T)
end

"""
    check_manifold_point(M::Manifold, x; kwargs...) -> Union{Nothing,String}

Return `nothing` when `x` is a point on manifold `M`. Otherwise, return a string
with a description why the point does not belong to manifold `M`.

By default, `check_manifold_point` returns `nothing`, i.e. if no checks are
implemented, the assumption is to be optimistic for a point not deriving from
the [`MPoint`](@ref) type.
"""
check_manifold_point(M::Manifold, x; kwargs...) = nothing

function check_manifold_point(M::Manifold, x::MPoint; kwargs...)
    error("check_manifold_point not implemented for manifold $(typeof(M)) and point $(typeof(x)).")
end

"""
    is_manifold_point(M::Manifold, x, throw_error = false; kwargs...)

Return whether `x` is a valid point on the [`Manifold`](@ref) `M`.

If `throw_error` is `false`, the function returns either `true` or `false`. If
`throw_error` is `true`, the function either returns `true` or throws an error.
By default the function calls [`check_manifold_point(M, x; kwargs...)`](@ref)
and checks whether the returned value is `nothing` or an error.
"""
function is_manifold_point(M::Manifold, x, throw_error = false; kwargs...)
    mpe = check_manifold_point(M, x; kwargs...)
    mpe === nothing && return true
    return throw_error ? throw(mpe) : false
end

"""
    check_tangent_vector(M::Manifold, x, v; kwargs...) -> Union{Nothing,String}

Check whether `v` is a valid tangent vector in the tangent plane of `x` on the
[`Manifold`](@ref) `M`. An implementation should first call
[`check_manifold_point(M, x; kwargs...)`](@ref) and then validate `v`. If it is
not a tangent vector, an error string should be returned.

By default, `check_tangent_vector` returns `nothing`, i.e. if no checks are
implemented, the assumption is to be optimistic for tangent vectors not deriving
from the [`TVector`](@ref) type.
"""
check_tangent_vector(M::Manifold, x, v; kwargs...) = nothing

function check_tangent_vector(M::Manifold, x::MPoint, v::TVector; kwargs...)
    error("check_tangent_vector not implemented for manifold $(typeof(M)), point $(typeof(x)) and vector $(typeof(v)).")
end

"""
    is_tangent_vector(M::Manifold, x, v, throw_error = false; kwargs...)

Return whether `v` is a valid tangent vector at point `x` on the
[`Manifold`](@ref) `M`. Returns either `true` or `false`.

The default is to return `true`, i.e. if no checks are implemented, the
assumption is to be optimistic.
"""
function is_tangent_vector(M::Manifold, x, v, throw_error = false; kwargs...)
    mtve = check_tangent_vector(M, x, v; kwargs...)
    mtve === nothing && return true
    return throw_error ? throw(mtve) : false
end

"""
    AbstractEstimationMethod

Abstract type for defining statistical estimation methods.
"""
abstract type AbstractEstimationMethod end

include("ArrayManifold.jl")
include("DefaultManifold.jl")

export Manifold,
       MPoint,
       TVector,
       CoTVector,
       ArrayManifold,
       ArrayMPoint,
       ArrayTVector,
       ArrayCoTVector

export ParallelTransport, ProjectionTransport

export base_manifold,
       check_manifold_point,
       check_tangent_vector,
       distance,
       exp,
       exp!,
       geodesic,
       shortest_geodesic,
       injectivity_radius,
       inner,
       inverse_retract,
       inverse_retract!,
       isapprox,
       is_manifold_point,
       is_tangent_vector,
       is_decorator_manifold,
       log,
       log!,
       manifold_dimension,
       norm,
       project_point,
       project_point!,
       project_tangent,
       project_tangent!,
       representation_size,
       retract,
       retract!,
       vector_transport_along,
       vector_transport_along!,
       vector_transport_direction,
       vector_transport_direction!,
       vector_transport_to,
       vector_transport_to!,
       zero_tangent_vector,
       zero_tangent_vector!

end # module
