module ManifoldsBase

import Base:
    isapprox,
    exp,
    log,
    convert,
    copy,
    copyto!,
    angle,
    eltype,
    isempty,
    length,
    similar,
    show,
    +,
    -,
    *,
    ==
import LinearAlgebra: dot, norm, det, cross, I, UniformScaling, Diagonal

import Markdown: @doc_str
using LinearAlgebra

include("maintypes.jl")
include("numbers.jl")
include("bases.jl")
include("retractions.jl")
include("exp_log_geo.jl")
include("projections.jl")

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
    allocate_result(M::AbstractManifold, f, x...)

Allocate an array for the result of function `f` on [`AbstractManifold`](@ref) `M` and arguments
`x...` for implementing the non-modifying operation using the modifying operation.

Usefulness of passing a function is demonstrated by methods that allocate results of musical
isomorphisms.
"""
@inline function allocate_result(M::AbstractManifold, f, x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[1], T)
end
@inline function allocate_result(M::AbstractManifold, f)
    T = allocate_result_type(M, f, ())
    return Array{T}(undef, representation_size(M)...)
end

"""
    allocate_result_type(M::AbstractManifold, f, args::NTuple{N,Any}) where N

Return type of element of the array that will represent the result of function `f` and the
[`AbstractManifold`](@ref) `M` on given arguments `args` (passed as a tuple).
"""
@inline function allocate_result_type(
    ::AbstractManifold,
    f::TF,
    args::NTuple{N,Any},
) where {N,TF}
    @inline eti_to_one(eti) = one(number_eltype(eti))
    return typeof(sum(map(eti_to_one, args)))
end
@inline function allocate_result_type(::AbstractManifold, f::TF, args::Tuple{}) where {TF}
    return Float64
end

"""
    angle(M::AbstractManifold, p, X, Y)

Compute the angle between tangent vectors `X` and `Y` at point `p` from the
[`AbstractManifold`](@ref) `M` with respect to the inner product from [`inner`](@ref).
"""
function angle(M::AbstractManifold, p, X, Y)
    return acos(real(inner(M, p, X, Y)) / norm(M, p, X) / norm(M, p, Y))
end
"""
    base_manifold(M::AbstractManifold, depth = Val(-1))

Return the internally stored [`AbstractManifold`](@ref) for decorated manifold `M` and the base
manifold for vector bundles or power manifolds. The optional parameter `depth` can be used
to remove only the first `depth` many decorators and return the [`AbstractManifold`](@ref) from that
level, whether its decorated or not. Any negative value deactivates this depth limit.
"""
base_manifold(M::AbstractManifold, ::Val = Val(-1)) = M

"""
    check_point(M::AbstractManifold, p; kwargs...) -> Union{Nothing,String}

Return `nothing` when `p` is a point on the [`AbstractManifold`](@ref) `M`. Otherwise, return an
error with description why the point does not belong to manifold `M`.

By default, `check_point` returns `nothing`, i.e. if no checks are implemented, the
assumption is to be optimistic for a point not deriving from the [`AbstractManifoldPoint`](@ref) type.
"""
check_point(M::AbstractManifold, p; kwargs...) = nothing

"""
    check_vector(M::AbstractManifold, p, X; kwargs...) -> Union{Nothing,String}

Check whether `X` is a valid tangent vector in the tangent space of `p` on the
[`AbstractManifold`](@ref) `M`. An implementation does not have to validate the point `p`.
If it is not a tangent vector, an error string should be returned.

By default, `check_vector` returns `nothing`, i.e. if no checks are implemented, the
assumption is to be optimistic for tangent vectors not deriving from the [`TVector`](@ref)
type.
"""
check_vector(M::AbstractManifold, p, X; kwargs...) = nothing

"""
    check_size(M::AbstractManifold, p)
    check_size(M::AbstractManifold, p, X)

Check whether `p` has the right [`representation_size`](@ref) for a [`AbstractManifold`](@ref) `M`.
Additionally if a tangent vector is given, both `p` and `X` are checked to be of
corresponding correct representation sizes for points and tangent vectors on `M`.

By default, `check_size` returns `nothing`, i.e. if no checks are implemented, the
assumption is to be optimistic.
"""
function check_size(M::AbstractManifold, p)
    m = representation_size(M)
    m === nothing && return nothing # nothing reasonable in size to check
    n = size(p)
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
function check_size(M::AbstractManifold, p, X)
    mse = check_size(M, p)
    mse === nothing || return mse
    m = representation_size(M)
    m === nothing && return nothing # without a representation size - nothing to check.
    n = size(X)
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

@doc raw"""
    copy(M, p)

Copy the value(s) from the point `p` on the [`AbstractManifold`](@ref) `M` into a new point.
See [`allocate_result`](@ref) for the allocation of new point memory and [`copyto!`](@ref) for the copying.
"""
function copy(M::AbstractManifold, p)
    q = allocate_result(M, copy, p)
    copyto!(M, q, p)
    return q
end

@doc raw"""
    copy(M, p, X)

Copy the value(s) from the tangent vector `X` at a point `p` on the
[`AbstractManifold`](@ref) `M` into a new tangent vector.
See [`allocate_result`](@ref) for the allocation of new point memory and [`copyto!`](@ref) for the copying.
"""
function copy(M::AbstractManifold, p, X)
    # the order of args switched, since the allocation by default takes the type of the first.
    Y = allocate_result(M, copy, X, p)
    copyto!(M, Y, p, X)
    return Y
end


@doc raw"""
    copyto!(M::AbstractManifold, q, p)

Copy the value(s) from `p` to `q`, where both are points on the [`AbstractManifold`](@ref) `M`.
This function defaults to calling `copyto!(q, p)`, but it might be useful to overwrite the
function at the level, where also information from `M` can be accessed.
"""
copyto!(::AbstractManifold, q, p) = copyto!(q, p)

@doc raw"""
    copyto!(M::AbstractManifold, Y, p, X)

Copy the value(s) from `X` to `Y`, where both are tangent vectors from the tangent space at
`p` on the [`AbstractManifold`](@ref) `M`.
This function defaults to calling `copyto!(Y, X)`, but it might be useful to overwrite the
function at the level, where also information from `p` and `M` can be accessed.
"""
copyto!(::AbstractManifold, Y, p, X) = copyto!(Y, X)

@doc raw"""
    distance(M::AbstractManifold, p, q)

Shortest distance between the points `p` and `q` on the [`AbstractManifold`](@ref) `M`,
i.e.

```math
d(p,q) = \inf_{γ} L(γ),
```
where the infimum is over all piecewise smooth curves ``γ: [a,b] \to \mathcal M``
connecting ``γ(a)=p`` and ``γ(b)=q`` and

```math
L(γ) = \displaystyle\int_{a}^{b} \lVert \dotγ(t)\rVert_{γ(t)} \mathrm{d}t
```
is the length of the curve $γ$.

If ``\mathcal M`` is not connected, i.e. consists of several disjoint components,
the distance between two points from different components should be ``∞``.
"""
distance(M::AbstractManifold, p, q) = norm(M, p, log(M, p, q))

"""
    embed(M::AbstractManifold, p)

Embed point `p` from the [`AbstractManifold`](@ref) `M` into the ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `embed` includes changing data representation, if applicable, i.e.
if the points on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly.

See also: [`EmbeddedManifold`](@ref), [`project`](@ref project(M::AbstractManifold,p))
"""
function embed(M::AbstractManifold, p)
    q = allocate_result(M, embed, p)
    embed!(M, q, p)
    return q
end

"""
    embed!(M::AbstractManifold, q, p)

Embed point `p` from the [`AbstractManifold`](@ref) `M` into an ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Not implementing this function means, there is no proper embedding for your manifold.
Additionally, `embed` might include changing data representation, if applicable, i.e.
if points on `M` are not represented in the same way as their counterparts in the embedding,
the representation is changed accordingly.

If you have more than one embedding, see [`EmbeddedManifold`](@ref) for defining a second
embedding. If your point `p` is already represented in some embedding,
see [`AbstractDecoratorManifold`](@ref) how you can avoid reimplementing code from the embedded manifold

See also: [`EmbeddedManifold`](@ref), [`project!`](@ref project!(M::AbstractManifold, q, p))
"""
embed!(M::AbstractManifold, q, p)

"""
    embed(M::AbstractManifold, p, X)

Embed a tangent vector `X` at a point `p` on the [`AbstractManifold`](@ref) `M` into an ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Not implementing this function means, there is no proper embedding for your tangent space(s).

Additionally, `embed` might include changing data representation, if applicable, i.e.
if tangent vectors on `M` are not represented in the same way as their counterparts in the
embedding, the representation is changed accordingly.

If you have more than one embedding, see [`EmbeddedManifold`](@ref) for defining a second
embedding. If your tangent vector `X` is already represented in some embedding,
see [`AbstractDecoratorManifold`](@ref) how you can avoid reimplementing code from the embedded manifold

See also: [`EmbeddedManifold`](@ref), [`project`](@ref project(M::AbstractManifold, p, X))
"""
function embed(M::AbstractManifold, p, X)
    # the order of args switched, since the allocation by default takes the type of the first.
    Y = allocate_result(M, embed, X, p)
    embed!(M, Y, p, X)
    return Y
end

"""
    embed!(M::AbstractManifold, Y, p, X)

Embed a tangent vector `X` at a point `p` on the [`AbstractManifold`](@ref) `M` into the ambient
space and return the result in `Y`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `embed!` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as tangents on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. The embedded tangents are then in
the tangent spaces of the embedded base points.

See also: [`EmbeddedManifold`](@ref), [`project!`](@ref project!(M::AbstractManifold, Y, p, X))
"""
embed!(M::AbstractManifold, Y, p, X)

@doc raw"""
    injectivity_radius(M::AbstractManifold)

Infimum of the injectivity radii `injectivity_radius(M,p)` of all points `p` on the [`AbstractManifold`](@ref).

    injectivity_radius(M::AbstractManifold, p)

Return the distance $d$ such that [`exp(M, p, X)`](@ref exp(::AbstractManifold, ::Any, ::Any)) is
injective for all tangent vectors shorter than $d$ (i.e. has an inverse).

    injectivity_radius(M::AbstractManifold[, x], method::AbstractRetractionMethod)
    injectivity_radius(M::AbstractManifold, x, method::AbstractRetractionMethod)

Distance ``d`` such that
[`retract(M, p, X, method)`](@ref retract(::AbstractManifold, ::Any, ::Any, ::AbstractRetractionMethod))
is injective for all tangent vectors shorter than ``d`` (i.e. has an inverse) for point `p`
if provided or all manifold points otherwise.

In order to dispatch on different retraction methods, please either implement
`_injectivity_radius(M[, p], m::T)` for your retraction `R` or specifically `injectivity_radius_exp(M[, p])` for the exponential map.
By default the variant with a point `p` assumes that the default (without `p`) can ve called as a lower bound.
"""
injectivity_radius(M::AbstractManifold)
injectivity_radius(M::AbstractManifold, p) = injectivity_radius(M)
function injectivity_radius(M::AbstractManifold, p, m::AbstractRetractionMethod)
    return _injectivity_radius(M, p, m)
end
function injectivity_radius(M::AbstractManifold, m::AbstractRetractionMethod)
    return _injectivity_radius(M, m)
end
function _injectivity_radius(M::AbstractManifold, p, m::AbstractRetractionMethod)
    return _injectivity_radius(M, m)
end
function _injectivity_radius(M::AbstractManifold, m::AbstractRetractionMethod)
    return _injectivity_radius(M)
end
function _injectivity_radius(M::AbstractManifold)
    return injectivity_radius(M)
end
function _injectivity_radius(M::AbstractManifold, p, ::ExponentialRetraction)
    return injectivity_radius_exp(M, p)
end
function _injectivity_radius(M::AbstractManifold, ::ExponentialRetraction)
    return injectivity_radius_exp(M)
end
injectivity_radius_exp(M, p) = injectivity_radius_exp(M)
injectivity_radius_exp(M) = injectivity_radius(M)

"""
    inner(M::AbstractManifold, p, X, Y)

Compute the inner product of tangent vectors `X` and `Y` at point `p` from the
[`AbstractManifold`](@ref) `M`.
"""
inner(M::AbstractManifold, p, X, Y)

"""
    isapprox(M::AbstractManifold, p, q; kwargs...)

Check if points `p` and `q` from [`AbstractManifold`](@ref) `M` are approximately equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(::AbstractManifold, x, y; kwargs...) = isapprox(x, y; kwargs...)

"""
    isapprox(M::AbstractManifold, p, X, Y; kwargs...)

Check if vectors `X` and `Y` tangent at `p` from [`AbstractManifold`](@ref) `M` are approximately
equal.

Keyword arguments can be used to specify tolerances.
"""
isapprox(::AbstractManifold, p, X, Y; kwargs...) = isapprox(X, Y; kwargs...)


"""
    is_point(M::AbstractManifold, p, throw_error = false; kwargs...)

Return whether `p` is a valid point on the [`AbstractManifold`](@ref) `M`.

If `throw_error` is `false`, the function returns either `true` or `false`. If `throw_error`
is `true`, the function either returns `true` or throws an error. By default the function
calls [`check_point`](@ref) and checks whether the returned value
is `nothing` or an error.
"""
function is_point(M::AbstractManifold, p, throw_error = false; kwargs...)
    mps = check_size(M, p)
    if mps !== nothing
        throw_error && throw(mps)
        return false
    end
    mpe = check_point(M, p; kwargs...)
    mpe === nothing && return true
    return throw_error ? throw(mpe) : false
end

"""
    is_vector(M::AbstractManifold, p, X, throw_error = false; check_base_point=true, kwargs...)

Return whether `X` is a valid tangent vector at point `p` on the [`AbstractManifold`](@ref) `M`.
Returns either `true` or `false`.

If `throw_error` is `false`, the function returns either `true` or `false`. If `throw_error`
is `true`, the function either returns `true` or throws an error. By default the function
calls [`check_vector`](@ref) and checks whether the returned
value is `nothing` or an error.

If `check_base_point` is true, then the point `p` will be first checked using the
[`check_point`](@ref) function.
"""
function is_vector(
    M::AbstractManifold,
    p,
    X,
    throw_error = false,
    check_base_point = true;
    kwargs...,
)
    if check_base_point
        s = is_point(M, p, throw_error; kwargs...) # if throw_error, is_point throws,
        !s && return false # otherwise if not a point return false
    end
    mXs = check_size(M, p, X)
    if mXs !== nothing
        throw_error && throw(mXs)
        return false
    end
    mXe = check_vector(M, p, X; kwargs...)
    mXe === nothing && return true
    throw_error && throw(mXe)
    return false
end

@doc raw"""
    manifold_dimension(M::AbstractManifold)

The dimension $n=\dim_{\mathcal M}$ of real space $\mathbb R^n$ to which the neighborhood of
each point of the [`AbstractManifold`](@ref) `M` is homeomorphic.
"""
manifold_dimension(M::AbstractManifold)

"""
    mid_point(M::AbstractManifold, p1, p2)

Calculate the middle between the two point `p1` and `p2` from manifold `M`.
By default uses [`log`](@ref), divides the vector by 2 and uses [`exp`](@ref).
"""
function mid_point(M::AbstractManifold, p1, p2)
    q = allocate(p1)
    return mid_point!(M, q, p1, p2)
end

"""
    mid_point!(M::AbstractManifold, q, p1, p2)

Calculate the middle between the two point `p1` and `p2` from manifold `M`.
By default uses [`log`](@ref), divides the vector by 2 and uses [`exp!`](@ref).
Saves the result in `q`.
"""
function mid_point!(M::AbstractManifold, q, p1, p2)
    X = log(M, p1, p2)
    return exp!(M, q, p1, X / 2)
end

@static if VERSION <= v"1.1"
    function mid_point!(
        M::AbstractManifold,
        q::AbstractArray{T1,0},
        p1::AbstractArray{T2,0},
        p2::AbstractArray{T3,0},
    ) where {T1,T2,T3}
        X = log(M, p1, p2)
        return exp!(M, q, p1, fill(X / 2))
    end
end

"""
    norm(M::AbstractManifold, p, X)

Compute the norm of tangent vector `X` at point `p` from a [`AbstractManifold`](@ref) `M`.
By default this is computed using [`inner`](@ref).
"""
norm(M::AbstractManifold, p, X) = sqrt(max(real(inner(M, p, X, X)), 0))

"""
    number_eltype(x)

Numeric element type of the a nested representation of a point or a vector.
To be used in conjuntion with [`allocate`](@ref) or [`allocate_result`](@ref).
"""
number_eltype(x) = eltype(x)
@inline function number_eltype(x::AbstractArray)
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, x))
end
@inline number_eltype(::AbstractArray{T}) where {T<:Number} = T
@inline function number_eltype(x::Tuple)
    @inline eti_to_one(eti) = one(number_eltype(eti))
    return typeof(sum(map(eti_to_one, x)))
end

@doc raw"""
    representation_size(M::AbstractManifold)

The size of an array representing a point on [`AbstractManifold`](@ref) `M`.
Returns `nothing` by default indicating that points are not represented using an
`AbstractArray`.
"""
function representation_size(::AbstractManifold)
    return nothing
end

"""
    size_to_tuple(::Type{S}) where S<:Tuple

Converts a size given by `Tuple{N, M, ...}` into a tuple `(N, M, ...)`.
"""
Base.@pure size_to_tuple(::Type{S}) where {S<:Tuple} = tuple(S.parameters...)

@doc raw"""
    zero_vector!(M::AbstractManifold, X, p)

Save to `X` the tangent vector from the tangent space ``T_p\mathcal M`` at `p` that
represents the zero vector, i.e. such that retracting `X` to the [`AbstractManifold`](@ref) `M` at
`p` produces `p`.
"""
zero_vector!(M::AbstractManifold, X, p) = log!(M, X, p, p)

@doc raw"""
    zero_vector(M::AbstractManifold, p)

Return the tangent vector from the tangent space ``T_p\mathcal M`` at `p` on the
[`AbstractManifold`](@ref) `M`, that represents the zero vector, i.e. such that a retraction at
`p` produces `p`.
"""
function zero_vector(M::AbstractManifold, p)
    X = allocate_result(M, zero_vector, p)
    zero_vector!(M, X, p)
    return X
end
include("errors.jl")
include("parallel_transport.jl")
include("vector_transport.jl")
include("vector_spaces.jl")
include("point_vector_fallbacks.jl")
include("nested_trait.jl")
include("decorator_trait.jl")
include("ValidationManifold.jl")
include("EmbeddedManifold.jl")
include("DefaultManifold.jl")
include("PowerManifold.jl")

export AbstractManifold, AbstractManifoldPoint, TVector, CoTVector, TFVector, CoTFVector
export AbstractDecoratorManifold
export AbstractTrait, IsEmbeddedManifold, IsEmbeddedSubmanifold, IsIsometricEmbeddedManifold
export IsExplicitDecorator
export ValidationManifold, ValidationMPoint, ValidationTVector, ValidationCoTVector
export EmbeddedManifold
export AbstractPowerManifold, PowerManifold
export AbstractPowerRepresentation,
    NestedPowerRepresentation, NestedReplacingPowerRepresentation

export OutOfInjectivityRadiusError

export AbstractRetractionMethod,
    ApproximateInverseRetraction,
    CayleyRetraction,
    EmbeddedRetraction,
    ExponentialRetraction,
    NLSolveInverseRetraction,
    ODEExponentialRetraction,
    QRRetraction,
    PadeRetraction,
    PolarRetraction,
    ProjectionRetraction,
    SoftmaxRetraction

export AbstractInverseRetractionMethod,
    ApproximateInverseRetraction,
    CayleyInverseRetraction,
    EmbeddedInverseRetraction,
    LogarithmicInverseRetraction,
    NLSolveInverseRetraction,
    QRInverseRetraction,
    PadeInverseRetraction,
    PolarInverseRetraction,
    ProjectionInverseRetraction,
    SoftmaxInverseRetraction

export AbstractVectorTransportMethod,
    DifferentiatedRetractionVectorTransport,
    ParallelTransport,
    PoleLadderTransport,
    ProjectionTransport,
    ScaledVectorTransport,
    SchildsLadderTransport,
    VectorTransportDirection,
    VectorTransportTo

export CachedBasis,
    DefaultBasis,
    DefaultOrthogonalBasis,
    DefaultOrthonormalBasis,
    DiagonalizingOrthonormalBasis,
    DefaultOrthonormalBasis,
    GramSchmidtOrthonormalBasis,
    ProjectedOrthonormalBasis,
    VeeOrthogonalBasis

export CompositeManifoldError, ComponentManifoldError

export allocate,
    angle,
    base_manifold,
    copy,
    copyto!,
    default_inverse_retraction_method,
    default_retraction_method,
    default_vector_transport_method,
    distance,
    exp,
    exp!,
    embed,
    embed!,
    geodesic,
    get_basis,
    get_component,
    get_coordinates,
    get_coordinates!,
    get_embedding,
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
    is_point,
    is_vector,
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
    parallel_transport_along,
    parallel_transport_along!,
    parallel_transport_direction,
    parallel_transport_direction!,
    parallel_transport_to,
    parallel_transport_to!,
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
    zero_vector,
    zero_vector!

end # module
