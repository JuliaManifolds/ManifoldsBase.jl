module ManifoldsBase

import Base: isapprox,
    exp,
    log,
    convert,
    copyto!,
    angle,
    eltype,
    similar,
    +,
    -,
    *
import LinearAlgebra: dot,
    norm,
    det,
    cross,
    I,
    UniformScaling,
    Diagonal

import Markdown: @doc_str
using LinearAlgebra

"""
    Manifold

A manifold type. The `Manifold` is used to dispatch to different exponential
and logarithmic maps as well as other function on manifold.
"""
abstract type Manifold end

"""
    MPoint

Type for a point on a manifold. While a [`Manifold`](@ref) not necessarily
requires this type, for example when it is implemented for `Vector`s or
`Matrix` type elements, this type can be used for more complicated
representations, semantic verification or even dispatch for different
representations of points on a manifold.
"""
abstract type MPoint end

"""
    TVector

Type for a tangent vector of a manifold. While a [`Manifold`](@ref) not
necessarily requires this type, for example when it is implemented for `Vector`s
or `Matrix` type elements, this type can be used for more complicated
representations, semantic verification or even dispatch for different
representations of tangent vectors and their types on a manifold.
"""
abstract type TVector end

"""
    CoTVector

Type for a cotangent vector of a manifold. While a [`Manifold`](@ref) not
necessarily requires this type, for example when it is implemented for `Vector`s
or `Matrix` type elements, this type can be used for more complicated
representations, semantic verification or even dispatch for different
representations of cotangent vectors and their types on a manifold.
"""
abstract type CoTVector end

"""
    is_decorator_manifold(M)

indicate whether a manifold is a decorator manifold, i.e. whether it encapsulates
a [`Manifold`](@ref) with additional features, and stores internally the original
manifold instance. An example is the [`ArrayManifold`](@ref).

Using Tim Holys Traits Trick (THTT), certain functions are just calling themselves
on the internal manifold, and hence have not to be implemented for decorators
again, for example [`manifold_dimension`](@ref) and especially [`base_manifold`](@ref).
This also assumes, that the undecoratord (base) manifold is stored in `M.manifold`.

"""
is_decorator_manifold(::T) where {T <: Manifold} = Val(false)

"""
    base_manifold(M::Manifold)

returns the internally stored manifold for decorated manifolds and the
base manifold for vector bundles or power manifolds.
For decorator manifolds it returns the dimension of the internally stored manifold.
"""
base_manifold(M::Manifold) = base_manifold(M, is_decorator_manifold(M))
base_manifold(M::Manifold, ::Val{true}) = base_manifold(M.manifold)
base_manifold(M::Manifold, ::Val{false}) = M

@doc doc"""
    representation_size(M::Manifold)

The size of array representing a point on manifold `M`.
For decorator manifolds it returns the dimension of the internally stored manifold.
"""
representation_size(M::Manifold) = representation_size(M, is_decorator_manifold(M))
representation_size(M::Manifold, ::Val{true}) = representation_size(base_manifold(M))
representation_size(M::Manifold, ::Val{false}) = error("representation_size not implemented for manifold $(typeof(M)).")


@doc doc"""
    manifold_dimension(M::Manifold)

The dimension $n$ of real space $\mathbb R^n$ to which the neighborhood
of each point of the manifold is homeomorphic.
For decorator manifolds it returns the dimension of the internally stored manifold.
"""
manifold_dimension(M::Manifold) = manifold_dimension(M, is_decorator_manifold(M))
manifold_dimension(M::Manifold, ::Val{true}) = manifold_dimension(base_manifold(M))
manifold_dimension(M::Manifold, ::Val{false}) = error("manifold_dimension not implemented for manifold $(typeof(M)).")

"""
    isapprox(M::Manifold, x, y; kwargs...)

Check if points `x` and `y` from manifold `M` are approximately equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(M::Manifold, x, y; kwargs...) = isapprox(x, y; kwargs...)

"""
    isapprox(M::Manifold, x, v, w; kwargs...)

Check if vectors `v` and `w` tangent at `x` from manifold `M` are
approximately equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(M::Manifold, x, v, w; kwargs...) = isapprox(v, w; kwargs...)

"""
    OutOfInjectivityRadiusError

An error thrown when a function (for example logarithmic map or inverse
retraction) is given arguments outside of its injectivity radius.
"""
struct OutOfInjectivityRadiusError <: Exception end

abstract type AbstractRetractionMethod end

"""
    ExponentialRetraction

Retraction using the exponential map.
"""
struct ExponentialRetraction <: AbstractRetractionMethod end

"""
    retract!(M::Manifold, y, x, v, [t=1], [method::AbstractRetractionMethod=ExponentialRetraction()])

Retraction (cheaper, approximate version of exponential map) of tangent
vector `t*v` at point `x` from manifold `M`.
Result is saved to `y`.

Retraction method can be specified by the last argument. Please look at the
documentation of respective manifolds for available methods.
"""
retract!(M::Manifold, y, x, v, method::ExponentialRetraction) = exp!(M, y, x, v)

retract!(M::Manifold, y, x, v) = retract!(M, y, x, v, ExponentialRetraction())

retract!(M::Manifold, y, x, v, t::Real) = retract!(M, y, x, t*v)

retract!(M::Manifold, y, x, v, t::Real, method::AbstractRetractionMethod) = retract!(M, y, x, t*v, method)

"""
    retract(M::Manifold, x, v, [t=1], [method::AbstractRetractionMethod])

Retraction (cheaper, approximate version of exponential map) of tangent
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

retract(M::Manifold, x, v, t::Real) = retract(M, x, t*v)

retract(M::Manifold, x, v, t::Real, method::AbstractRetractionMethod) = retract(M, x, t*v, method)

abstract type AbstractInverseRetractionMethod end

"""
    LogarithmicInverseRetraction

Inverse retraction using the logarithmic map.
"""
struct LogarithmicInverseRetraction <: AbstractInverseRetractionMethod end

"""
    inverse_retract!(M::Manifold, v, x, y, [method::AbstractInverseRetractionMethod=LogarithmicInverseRetraction()])

Inverse retraction (cheaper, approximate version of logarithmic map) of points
`x` and `y`.
Result is saved to `y`.

Inverse retraction method can be specified by the last argument. Please look
at the documentation of respective manifolds for available methods.
"""
inverse_retract!(M::Manifold, v, x, y, method::LogarithmicInverseRetraction) = log!(M, v, x, y)

inverse_retract!(M::Manifold, v, x, y) = inverse_retract!(M, v, x, y, LogarithmicInverseRetraction())

"""
    inverse_retract(M::Manifold, x, y, [method::AbstractInverseRetractionMethod])

Inverse retraction (cheaper, approximate version of logarithmic map) of points
`x` and `y` from manifold `M`.

Inverse retraction method can be specified by the last argument. Please look
at the documentation of respective manifolds for available methods.
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

Project point `x` from the ambient space onto the manifold `M`.
The point `y` is overwritten by the projection.
The function works only for selected embedded manifolds and
is *not* required to return the closest point.
"""
project_point!(M::Manifold, y, x) = error("project_point! not implemented for a $(typeof(M)) and points $(typeof(y)) and $(typeof(x)).")

"""
    project_point(M::Manifold, x)

Project point from the ambient space onto the manifold `M`. The point `x`
is not modified. The function works only for selected embedded manifolds and
is *not* required to return the closest point.
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
project_tangent!(M::Manifold, w, x, v) = error("project onto tangent space not implemented for a $(typeof(M)) and point $(typeof(x)) with input $(typeof(v)).")

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
inner(M::Manifold, x, v, w) = error("inner: Inner product not implemented on a $(typeof(M)) for input point $(typeof(x)) and tangent vectors $(typeof(v)) and $(typeof(w)).")

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
    exp!(M::Manifold, y, x, v, t=1)

Exponential map of tangent vector `t*v` at point `x` from manifold `M`.
Result is saved to `y`.
"""
exp!(M::Manifold, y, x, v, t::Real) = exp!(M, y, x, t*v)

exp!(M::Manifold, y, x, v) = error("Exponential map not implemented on a $(typeof(M)) for input point $(x) and tangent vector $(v).")

"""
    exp(M::Manifold, x, v, t=1)

Exponential map of tangent vector `t*v` at point `x` from manifold `M`.
"""
function exp(M::Manifold, x, v)
    x2 = similar_result(M, exp, x, v)
    exp!(M, x2, x, v)
    return x2
end

exp(M::Manifold, x, v, t::Real) = exp(M, x, t*v)

"""
    exp(M::Manifold, x, v, T::AbstractVector)

Exponential map of tangent vector `t*v` at point `x` from manifold `M` for
each `t` in `T`.
"""
exp(M::Manifold, x, v, T::AbstractVector) = map(geodesic(M, x, v), T)

log!(M::Manifold, v, x, y) = error("Logarithmic map not implemented on $(typeof(M)) for points $(typeof(x)) and $(typeof(y))")

function log(M::Manifold, x, y)
    v = similar_result(M, log, x, y)
    log!(M, v, x, y)
    return v
end

"""
    geodesic(M::Manifold, x, v)

Get the geodesic with initial point `x` and velocity `v`. The geodesic
is the curve of constant velocity that is locally distance-minimizing. This
function returns a function of time, which may be a `Real` or an
`AbstractVector`.
"""
geodesic(M::Manifold, x, v) = t -> exp(M, x, v, t)

"""
    geodesic(M::Manifold, x, v, t)

Get the point at time `t` traveling from `x` along the geodesic with initial
point `x` and velocity `v`.
"""
geodesic(M::Manifold, x, v, t::Real) = exp(M, x, v, t)

"""
    geodesic(M::Manifold, x, v, T::AbstractVector)

Get the points for each `t` in `T` traveling from `x` along the geodesic with
initial point `x` and velocity `v`.
"""
geodesic(M::Manifold, x, v, T::AbstractVector) = exp(M, x, v, T)

"""
    shortest_geodesic(M::Manifold, x, y)

Get a geodesic with initial point `x` and point `y` at `t=1` whose length is
the shortest path between the two points. When there are multiple shortest
geodesics, there is no guarantee which will be returned. This function returns
a function of time, which may be a `Real` or an `AbstractVector`.
"""
shortest_geodesic(M::Manifold, x, y) = geodesic(M, x, log(M, x, y))

"""
    shortest_geodesic(M::Manifold, x, y, t)

Get the point at time `t` traveling from `x` along a shortest geodesic
connecting `x` and `y`, where `y` is reached at `t=1`.
"""
shortest_geodesic(M::Manifold, x, y, t::Real) = geodesic(M, x, log(M, x, y), t)

"""
    shortest_geodesic(M::Manifold, x, y, T::AbstractVector)

Get the points for each `t` in `T` traveling from `x` along a shortest geodesic
connecting `x` and `y`, where `y` is reached at `t=1`.
"""
function shortest_geodesic(M::Manifold, x, y, T::AbstractVector)
    return geodesic(M, x, log(M, x, y), T)
end

abstract type AbstractVectorTransportMethod end

"""
    ParallelTransport <: AbstractVectorTransportMethod

Specify to use parallel transport as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref) or
[`vector_transport_along`](@ref).
"""
struct ParallelTransport <: AbstractVectorTransportMethod end

"""
    ProjectionTransport <: AbstractVectorTransportMethod

Specify to use projection onto tangent space as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref) or
[`vector_transport_along`](@ref). See [`project_tangent`](@ref) for details.
"""
struct ProjectionTransport <: AbstractVectorTransportMethod end
"""
    vector_transport_to!(M::Manifold, vto, x, v, y, m::AbstractVectorTransportMethod=ParallelTransport())

Vector transport of vector `v` at point `x` to point `y`. The result is saved
to `vto`. By default, the method `m` is [`ParallelTransport`](@ref).
"""
vector_transport_to!(M::Manifold, vto, x, v, y) = vector_transport_to!(M,vto,x,v,y,ParallelTransport())

"""
    vector_transport_to!(M::Manifold, vto, x, v, y, ProjectionTransport())

Implements a default projection based vector transport, that projects a tangent
vector `v` at `x` on a [`Manifold`](@ref) `M` onto the tangent space at `y` by
interperting `v` as an element of the embedding and projecting back.
"""
vector_transport_to!(M::Manifold, vto, x, v, y, m::ProjectionTransport) = project_tangent!(M, vto, y, v)

function vector_transport_to!(M::Manifold, vto, x, v, y, m::AbstractVectorTransportMethod)
    error("vector transport from a point of type $(typeof(x)) to a type $(typeof(y)) on a $(typeof(M)) for a vector of type $(v) and the $(typeof(m)) not yet implemented.")
end


"""
    vector_transport_to(M::Manifold, x, v, y, m::AbstractVectorTransportMethod=ParallelTransport())

Vector transport of vector `v` at point `x` to point `y` using the method `m`,
which defaults to [`ParallelTransport`](@ref).
"""
vector_transport_to(M::Manifold, x, v, y) = vector_transport_to(M,x,v,y,ParallelTransport())
function vector_transport_to(M::Manifold, x, v, y, m::AbstractVectorTransportMethod)
    vto = similar_result(M, vector_transport_to, v, x, y)
    vector_transport_to!(M, vto, x, v, y, m)
    return vto
end

"""
    vector_transport_direction!(M::Manifold, vto, x, v, vdir, m=::ParallelTransport])

Vector transport of vector `v` at point `x` in the direction indicated
by the tangent vector `vdir` at point `x`. The result is saved to `vto`.
By default, `exp` and `vector_transport_to!` are used with the method `m` which
defaults to [`ParallelTransport`](@ref).
"""
vector_transport_direction!(M::Manifold, vto, x, v, vdir) = vector_transport_direction!(M,vto,x,v,vdir,ParallelTransport())
function vector_transport_direction!(M::Manifold, vto, x, v, vdir,m::AbstractVectorTransportMethod)
    y = exp(M, x, vdir)
    return vector_transport_to!(M, vto, x, v, y, m)
end

"""
    vector_transport_direction(M::Manifold, x, v, vdir[, m=::ParallelTransport])

Vector transport of vector `v` at point `x` in the direction indicated
by the tangent vector `vdir` at point `x` using the method `m`, which defaults to [`ParallelTransport`](@ref).
"""
vector_transport_direction(M::Manifold, x, v, vdir) = vector_transport_direction(M,x,v,vdir,ParallelTransport())
function vector_transport_direction(M::Manifold, x, v, vdir, m::AbstractVectorTransportMethod)
    vto = similar_result(M, vector_transport_direction, v, x, vdir)
    vector_transport_direction!(M, vto, x, v, vdir,m)
    return vto
end

"""
    vector_transport_along!(M::Manifold, vto, x, v, c[,m=::ParallelTransport])

Vector transport of vector `v` at point `x` along the curve `c` such that
`c(0)` is equal to `x` to point `c(1)` using the method `m`, which defaults to
[`ParallelTransport`](@ref). The result is saved to `vto`.
"""
vector_transport_along!(M::Manifold, vto, x, v, c) = vector_transport_along!(M, vto, x, v, c, ParallelTransport())
function vector_transport_along!(M::Manifold, vto, x, v, c, m::AbstractVectorTransportMethod)
    error("vector_transport_along! not implemented for manifold $(typeof(M)), vector $(typeof(vto)), point $(typeof(x)), vector $(typeof(v)) along curve $(typeof(c)) with method $(typeof(m)).")
end

"""
    vector_transport_along(M::Manifold, x, v, c[,m])

Vector transport of vector `v` at point `x` along the curve `c` such that
`c(0)` is equal to `x` to point `c(1)`.
The default method `m` used is [`ParallelTransport`](@ref).
"""
vector_transport_along(M::Manifold, x, v, c) = vector_transport_along(M,x,v,c,ParallelTransport())
function vector_transport_along(M::Manifold, x, v, c, m::AbstractVectorTransportMethod)
    vto = similar_result(M, vector_transport_along, v, x)
    vector_transport_along!(M, vto, x, v, c, m)
    return vto
end

@doc doc"""
    injectivity_radius(M::Manifold, x)

Distance $d$ such that `exp(M, x, v)` is injective for all tangent
vectors shorter than $d$ (has a left inverse).
"""
injectivity_radius(M::Manifold, x) = Inf

@doc doc"""
    injectivity_radius(M::Manifold, x, R::AbstractRetractionMethod)

Distance $d$ such that `retract(M, x, v, R)` is injective for all tangent
vectors shorter than $d$ (has a left inverse).
"""
injectivity_radius(M::Manifold, x, ::AbstractRetractionMethod) = injectivity_radius(M, x)

"""
    injectivity_radius(M::Manifold)

Infimum of the injectivity radii of all manifold points.
"""
injectivity_radius(M::Manifold) = Inf

function zero_tangent_vector(M::Manifold, x)
    v = similar_result(M, zero_tangent_vector, x)
    zero_tangent_vector!(M, v, x)
    return v
end

zero_tangent_vector!(M::Manifold, v, x) = log!(M, v, x, x)

"""
    similar_result_type(M::Manifold, f, args::NTuple{N,Any}) where N

Returns type of element of the array that will represent the result of
function `f` for manifold `M` on given arguments (passed at a tuple).
"""
function similar_result_type(M::Manifold, f, args::NTuple{N,Any}) where N
    T = typeof(reduce(+, one(eltype(eti)) for eti ∈ args))
    return T
end

"""
    similar_result(M::Manifold, f, x...)

Allocates an array for the result of function `f` on manifold `M`
and arguments `x...` for implementing the non-modifying operation
using the modifying operation.

Usefulness of passing a function is demonstrated by methods that allocate
results of musical isomorphisms.
"""
function similar_result(M::Manifold, f, x...)
    T = similar_result_type(M, f, x)
    return similar(x[1], T)
end

"""
    check_manifold_point(M::Manifold, x; kwargs...)

Return `nothing` when `x` is a point on manifold `M`.
Otherwise, return a string with description why the point does not belong
to manifold `M`.

By default, `check_manifold_point` returns `nothing`, i.e. if no checks are implmented,
the assumption is to be optimistic for point not deriving from the [`MPoint`](@ref) type.
"""
function check_manifold_point(M::Manifold, x; kwargs...)
    return nothing
end

function check_manifold_point(M::Manifold, x::MPoint; kwargs...)
    error("check_manifold_point not implemented for manifold $(typeof(M)) and point $(typeof(x)).")
end

"""
    is_manifold_point(M, x, throw_error = false; kwargs...)

check, whether `x` is a valid point on the [`Manifold`](@ref) `M`.

If `throw_error` is false, the function returns either `true` or `false`.
If `throw_error` if true, the function either returns `true` or throws an error.
By default the function calls [`check_manifold_point`](@ref)`(M, x; kwargs...)`
and checks whether the returned value is `nothing` or an error.
"""
function is_manifold_point(M::Manifold, x, throw_error = false; kwargs...)
    mpe = check_manifold_point(M, x; kwargs...)
    if throw_error
        if mpe !== nothing
            throw(mpe)
        end
        return true
    else
        return mpe === nothing
    end
end

"""
    check_tangent_vector(M::Manifold, x, v; kwargs...)

check, whether `v` is a valid tangent vector in the tangent plane of `x` on the
[`Manifold`](@ref) `M`. An implementation should first check
[`manifold_point_error`](@ref)`(M, x; kwargs...)` and then validate `v`.
If it is not a tangent vector error string should be returned.

By default, `check_tangent_vector` returns `nothing`, i.e. if no checks are implmented,
the assumption is to be optimistic for tangent vectors not deriving from the [`TVector`](@ref) type.
"""
function check_tangent_vector(M::Manifold, x, v; kwargs...)
    return nothing
end

function check_tangent_vector(M::Manifold, x::MPoint, v::TVector; kwargs...)
    error("check_tangent_vector not implemented for manifold $(typeof(M)), point $(typeof(x)) and vector $(typeof(v)).")
end

"""
    is_tangent_vector(M, x, v, throw_error = false; kwargs...)

check, whether `v` is a valid tangent vector at point `x` on
the [`Manifold`](@ref) `M`.

If `throw_error` is false, the function returns either `true` or `false`.
If `throw_error` if true, the function either returns `true` or throws an error.
By default the function calls [`check_tangent_vector`](@ref)`(M, x, v; kwargs...)`
and checks whether the returned value is `nothing` or an error.
"""
function is_tangent_vector(M::Manifold, x, v, throw_error = false; kwargs...)
    tve = check_tangent_vector(M, x, v; kwargs...)
    if throw_error
        if tve !== nothing
            throw(tve)
        end
        return true
    else
        return tve === nothing
    end
end

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

export ParallelTransport,
    ProjectionTransport

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
    manifold_point_error,
    norm,
    project_point,
    project_point!,
    project_tangent,
    project_tangent!,
    representation_size,
    retract,
    retract!,
    tangent_vector_error,
    vector_transport_along,
    vector_transport_along!,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
    vector_transport_to!,
    zero_tangent_vector,
    zero_tangent_vector!

end # module
