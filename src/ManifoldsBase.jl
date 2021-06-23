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

include("maintypes.jl")
include("retractions.jl")
include("exp_log_geo.jl")
include("projections.jl")


"""
    OutOfInjectivityRadiusError

An error thrown when a function (for example [`log`](@ref)arithmic map or
[`inverse_retract`](@ref)) is given arguments outside of its [`injectivity_radius`](@ref).
"""
struct OutOfInjectivityRadiusError <: Exception end

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
function allocate_result(M::AbstractManifold, f, x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[1], T)
end

"""
    allocate_result_type(M::AbstractManifold, f, args::NTuple{N,Any}) where N

Return type of element of the array that will represent the result of function `f` and the
[`AbstractManifold`](@ref) `M` on given arguments `args` (passed as a tuple).
"""
function allocate_result_type(M::AbstractManifold, f, args::NTuple{N,Any}) where {N}
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, args))
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
base_manifold(M::AbstractManifold, depth = Val(-1)) = M

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
function check_size(M::AbstractManifold, p, X)
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
see [`AbstractEmbeddedManifold`](@ref) how you can avoid reimplementing code from the embedded manifold

See also: [`EmbeddedManifold`](@ref), [`project!`](@ref project!(M::AbstractManifold, q, p))
"""
function embed!(M::AbstractManifold, q, p)
    return error(manifold_function_not_implemented_message(M, embed!, q, p))
end

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
see [`AbstractEmbeddedManifold`](@ref) how you can avoid reimplementing code from the embedded manifold

See also: [`EmbeddedManifold`](@ref), [`project`](@ref project(M::AbstractManifold, p, X))
"""
function embed(M::AbstractManifold, p, X)
    # Note that the order is switched,
    # since the allocation by default takes the type of the first.
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
function embed!(M::AbstractManifold, Y, p, X)
    return error(manifold_function_not_implemented_message(M, embed!, Y, p, X))
end

@doc raw"""
    injectivity_radius(M::AbstractManifold, p)

Return the distance $d$ such that [`exp(M, p, X)`](@ref exp(::AbstractManifold, ::Any, ::Any)) is
injective for all tangent vectors shorter than $d$ (i.e. has an inverse).

    injectivity_radius(M::AbstractManifold)

Infimum of the injectivity radius of all manifold points.

    injectivity_radius(M::AbstractManifold[, x], method::AbstractRetractionMethod)
    injectivity_radius(M::AbstractManifold, x, method::AbstractRetractionMethod)

Distance ``d`` such that
[`retract(M, p, X, method)`](@ref retract(::AbstractManifold, ::Any, ::Any, ::AbstractRetractionMethod))
is injective for all tangent vectors shorter than ``d`` (i.e. has an inverse) for point `p`
if provided or all manifold points otherwise.
"""
function injectivity_radius(M::AbstractManifold)
    return error(manifold_function_not_implemented_message(M, injectivity_radius))
end
injectivity_radius(M::AbstractManifold, p) = injectivity_radius(M)
function injectivity_radius(M::AbstractManifold, p, method::AbstractRetractionMethod)
    return injectivity_radius(M, method)
end
function injectivity_radius(M::AbstractManifold, method::AbstractRetractionMethod)
    return error(manifold_function_not_implemented_message(M, injectivity_radius, method))
end
function injectivity_radius(M::AbstractManifold, p, ::ExponentialRetraction)
    return injectivity_radius(M, p)
end
injectivity_radius(M::AbstractManifold, ::ExponentialRetraction) = injectivity_radius(M)

"""
    inner(M::AbstractManifold, p, X, Y)

Compute the inner product of tangent vectors `X` and `Y` at point `p` from the
[`AbstractManifold`](@ref) `M`.

See also: [`MetricManifold`](@ref Main.Manifolds.MetricManifold)
"""
function inner(M::AbstractManifold, p, X, Y)
    return error(manifold_function_not_implemented_message(M, inner, p, X, Y))
end

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
calls [`check_point(M, p; kwargs...)`](@ref) and checks whether the returned value
is `nothing` or an error.
"""
function is_point(M::AbstractManifold, p, throw_error = false; kwargs...)
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
calls [`check_vector(M, p, X; kwargs...)`](@ref) and checks whether the returned
value is `nothing` or an error.

If `check_base_point` is true, then the point `p` will be first checked using the
[`check_point`](@ref) function.
"""
function is_vector(
    M::AbstractManifold,
    p,
    X,
    throw_error = false;
    check_base_point = true,
    kwargs...,
)
    if check_base_point
        mpe = check_point(M, p; kwargs...)
        if mpe !== nothing
            if throw_error
                throw(mpe)
            else
                return false
            end
        end
    end
    mtve = check_vector(M, p, X; kwargs...)
    mtve === nothing && return true
    return throw_error ? throw(mtve) : false
end

@doc raw"""
    manifold_dimension(M::AbstractManifold)

The dimension $n=\dim_{\mathcal M}$ of real space $\mathbb R^n$ to which the neighborhood of
each point of the [`AbstractManifold`](@ref) `M` is homeomorphic.
"""
function manifold_dimension(M::AbstractManifold)
    return error(manifold_function_not_implemented_message(M, manifold_dimension))
end

function manifold_function_not_implemented_message(M::AbstractManifold, f, x...)
    s = join(map(string, map(typeof, x)), ", ", " and ")
    a = length(x) > 1 ? "arguments" : "argument"
    m = length(x) > 0 ? " for $(a) $(s)." : "."
    return "$(f) not implemented on $(M)$(m)"
end

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
function number_eltype(x::AbstractArray)
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, x))
end
function number_eltype(x::Tuple)
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, x))
end

@doc raw"""
    representation_size(M::AbstractManifold)

The size of an array representing a point on [`AbstractManifold`](@ref) `M`.
Returns `nothing` by default indicating that points are not represented using an
`AbstractArray`.
"""
function representation_size(M::AbstractManifold)
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
include("numbers.jl")
include("vector_transport.jl")
include("DecoratorManifold.jl")
include("bases.jl")
include("vector_spaces.jl")
include("ValidationManifold.jl")
include("EmbeddedManifold.jl")
include("DefaultManifold.jl")
include("PowerManifold.jl")

export AbstractManifold, AbstractManifoldPoint, TVector, CoTVector, TFVector, CoTFVector
export AbstractDecoratorManifold
export ValidationManifold, ValidationMPoint, ValidationTVector, ValidationCoTVector
export AbstractEmbeddingType,
    TransparentIsometricEmbedding, DefaultIsometricEmbeddingType, DefaultEmbeddingType
export AbstractEmbeddedManifold, EmbeddedManifold, TransparentIsometricEmbedding
export AbstractPowerManifold, PowerManifold
export AbstractPowerRepresentation,
    NestedPowerRepresentation, NestedReplacingPowerRepresentation

export AbstractDecoratorType, DefaultDecoratorType

export OutOfInjectivityRadiusError

export AbstractRetractionMethod,
    ApproximateInverseRetraction,
    NLsolveInverseRetraction,
    ExponentialRetraction,
    QRRetraction,
    PolarRetraction,
    ProjectionRetraction,
    PowerRetraction,
    InversePowerRetraction

export AbstractInverseRetractionMethod,
    ApproximateInverseRetraction,
    LogarithmicInverseRetraction,
    QRInverseRetraction,
    PolarInverseRetraction,
    ProjectionInverseRetraction

export AbstractVectorTransportMethod,
    DifferentiatedRetractionVectorTransport,
    ParallelTransport,
    PoleLadderTransport,
    PowerVectorTransport,
    ProjectionTransport,
    ScaledVectorTransport,
    SchildsLadderTransport

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
    check_point,
    check_vector,
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
