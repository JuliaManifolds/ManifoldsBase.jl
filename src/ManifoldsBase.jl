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
    fill,
    fill!,
    isempty,
    length,
    similar,
    show,
    +,
    -,
    *,
    ==
import LinearAlgebra: ×, dot, norm, det, cross, I, UniformScaling, Diagonal
import Random: rand, rand!

using LinearAlgebra
using Markdown: @doc_str
using Printf: @sprintf
using Random

include("maintypes.jl")
include("numbers.jl")
include("Fiber.jl")
include("bases.jl")
include("approximation_methods.jl")
include("retractions.jl")
include("exp_log_geo.jl")
include("projections.jl")
include("metric.jl")

if isdefined(Base, Symbol("@constprop"))
    macro aggressive_constprop(ex)
        return esc(:(Base.@constprop :aggressive $ex))
    end
else
    macro aggressive_constprop(ex)
        return esc(ex)
    end
end

"""
    allocate(a)
    allocate(a, dims::Integer...)
    allocate(a, dims::Tuple)
    allocate(a, T::Type)
    allocate(a, T::Type, dims::Integer...)
    allocate(a, T::Type, dims::Tuple)
    allocate(M::AbstractManifold, a)
    allocate(M::AbstractManifold, a, dims::Integer...)
    allocate(M::AbstractManifold, a, dims::Tuple)
    allocate(M::AbstractManifold, a, T::Type)
    allocate(M::AbstractManifold, a, T::Type, dims::Integer...)
    allocate(M::AbstractManifold, a, T::Type, dims::Tuple)

Allocate an object similar to `a`. It is similar to function `similar`, although
instead of working only on the outermost layer of a nested structure, it maps recursively
through outer layers and calls `similar` on the innermost array-like object only.
Type `T` is the new number element type [`number_eltype`](@ref), if it is not given
the element type of `a` is retained. The `dims` argument can be given for non-nested
allocation and is forwarded to the function `similar`.

It's behavior can be overridden by a specific manifold, for example power manifold with
nested replacing representation can decide that `allocate` for `Array{<:SArray}` returns
another `Array{<:SArray}` instead of `Array{<:MArray}`, as would be done by default.
"""
allocate(a, args...)
allocate(a) = similar(a)
allocate(a, dim1::Integer, dims::Integer...) = similar(a, dim1, dims...)
allocate(a, dims::Tuple) = similar(a, dims)
allocate(a, T::Type) = similar(a, T)
allocate(::Number, T::Type) = fill(zero(T))
allocate(a, T::Type, dim1::Integer, dims::Integer...) = similar(a, T, dim1, dims...)
allocate(a, T::Type, dims::Tuple) = similar(a, T, dims)
allocate(a::AbstractArray{<:AbstractArray}) = map(allocate, a)
allocate(a::AbstractArray{<:AbstractArray}, T::Type) = map(t -> allocate(t, T), a)
allocate(a::NTuple{N,AbstractArray} where {N}) = map(allocate, a)
allocate(a::NTuple{N,AbstractArray} where {N}, T::Type) = map(t -> allocate(t, T), a)

allocate(::AbstractManifold, a) = allocate(a)
function allocate(::AbstractManifold, a, dim1::Integer, dims::Integer...)
    return allocate(a, dim1, dims...)
end
allocate(::AbstractManifold, a, dims::Tuple) = allocate(a, dims)
allocate(::AbstractManifold, a, T::Type) = allocate(a, T)
function allocate(::AbstractManifold, a, T::Type, dim1::Integer, dims::Integer...)
    return allocate(a, T, dim1, dims...)
end
allocate(::AbstractManifold, a, T::Type, dims::Tuple) = allocate(a, T, dims)

"""
    allocate_on(M::AbstractManifold, [T:::Type])
    allocate_on(M::AbstractManifold, F::FiberType, [T:::Type])

Allocate a new point on manifold `M` with optional type given by `T`. Note that `T` is not
number element type as in [`allocate`](@ref) but rather the type of the entire point to be
returned.

If `F` is provided, then an element of the corresponding fiber is allocated, assuming it is
independent of the base point.

To allocate a tangent vector, use ``

# Example

```julia-repl
julia> using ManifoldsBase

julia> M = ManifoldsBase.DefaultManifold(4)
DefaultManifold(4; field = ℝ)

julia> allocate_on(M)
4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0

julia> allocate_on(M, Array{Float64})
4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0

julia> allocate_on(M, TangentSpaceType())
4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0

julia> allocate_on(M, TangentSpaceType(), Array{Float64})
4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0

```
"""
allocate_on(M::AbstractManifold) = similar(Array{Float64}, representation_size(M))
function allocate_on(M::AbstractManifold, T::Type{<:AbstractArray})
    return similar(T, representation_size(M))
end

"""
    _pick_basic_allocation_argument(::AbstractManifold, f, x...)

Pick which one of elements of `x` should be used as a basis for allocation in the
`allocate_result(M::AbstractManifold, f, x...)` method. This can be specialized to, for
example, skip `Identity` arguments in Manifolds.jl group-related functions.
"""
function _pick_basic_allocation_argument(::AbstractManifold, f, x...)
    return x[1]
end

"""
    allocate_result(M::AbstractManifold, f, x...)

Allocate an array for the result of function `f` on [`AbstractManifold`](@ref) `M` and arguments
`x...` for implementing the non-modifying operation using the modifying operation.

Usefulness of passing a function is demonstrated by methods that allocate results of musical
isomorphisms.
"""
@inline function allocate_result(M::AbstractManifold, f, x...)
    T = allocate_result_type(M, f, x)
    picked = _pick_basic_allocation_argument(M, f, x...)
    return allocate(M, picked, T)
end
@inline function allocate_result(M::AbstractManifold, f)
    T = allocate_result_type(M, f, ())
    rs = representation_size(M)
    return allocate_result_array(M, f, T, rs)
end

function allocate_result_array(M::AbstractManifold, f, ::Type, ::Nothing)
    msg = "Could not allocate result of function $f on manifold $M."
    if base_manifold(M) isa ProductManifold
        msg *= " This error could be resolved by importing RecursiveArrayTools.jl. If this is not the case, please open report an issue."
    end
    return error(msg)
end
function allocate_result_array(::AbstractManifold, f, T::Type, rs::Tuple)
    return Array{T}(undef, rs...)
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
@inline function allocate_result_type(M::AbstractManifold, f::TF, args::Tuple{}) where {TF}
    return allocation_promotion_function(M, f, ())(Float64)
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
    are_linearly_independent(M::AbstractManifold, p, X, Y)

Check is vectors `X`, `Y` tangent at `p` to `M` are linearly independent.
"""
function are_linearly_independent(
    M::AbstractManifold,
    p,
    X,
    Y;
    atol::Real = sqrt(eps(number_eltype(X))),
)
    norm_X = norm(M, p, X)
    norm_Y = norm(M, p, Y)
    innerXY = inner(M, p, X, Y)
    return norm_X > atol && norm_Y > atol && !isapprox(abs(innerXY), norm_X * norm_Y)
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
    check_approx(M::AbstractManifold, p, q; kwargs...)
    check_approx(M::AbstractManifold, p, X, Y; kwargs...)

Check whether two elements are approximately equal, either `p`, `q` on the [`AbstractManifold`](@ref)
or the two tangent vectors `X`, `Y` in the tangent space at `p` are approximately the same.
The keyword arguments `kwargs` can be used to set tolerances, similar to Julia's `isapprox`.

This function might use `isapprox` from Julia internally and is similar to [`isapprox`](@ref),
with the difference that is returns an [`ApproximatelyError`](@ref) if the two elements are
not approximately equal, containting a more detailed description/reason.
If the two elements are approximalely equal, this method returns `nothing`.

This method is an internal function and is called by `isapprox` whenever the user specifies
an `error=` keyword therein. [`_isapprox`](@ref) is another related internal function. It is
supposed to provide a fast true/false decision whether points or vectors are equal or not,
while `check_approx` also provides a textual explanation. If no additional explanation
is needed, a manifold may just implement a method of [`_isapprox`](@ref), while it should also
implement `check_approx` if a more detailed explanation could be helpful.
"""
function check_approx(M::AbstractManifold, p, q; kwargs...)
    # fall back to classical approx mode - just that we do not have a reason then
    res = _isapprox(M, p, q; kwargs...)
    # since we can not assume distance to be implemented, we can not provide a default value
    res && return nothing
    s = "The two points $p and $q on $M are not (approximately) equal."
    return ApproximatelyError(s)
end
function check_approx(M::AbstractManifold, p, X, Y; kwargs...)
    # fall back to classical mode - just that we do not have a reason then
    res = _isapprox(M, p, X, Y; kwargs...)
    res && return nothing
    s = "The two tangent vectors $X and $Y in the tangent space at $p on $M are not (approximately) equal."
    v = try
        norm(M, p, X - Y)
    catch e
        NaN
    end
    return ApproximatelyError(v, s)
end

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
assumption is to be optimistic for tangent vectors not deriving from the [`AbstractTangentVector`](@ref)
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

"""
    convert(T::Type, M::AbstractManifold, p)

Convert point `p` from manifold `M` to type `T`.
"""
convert(T::Type, ::AbstractManifold, p) = convert(T, p)
"""
    convert(T::Type, M::AbstractManifold, p, X)

Convert vector `X` tangent at point `p` from manifold `M` to type `T`.
"""
convert(T::Type, ::AbstractManifold, p, X) = convert(T, p, X)

@doc raw"""
    copy(M::AbstractManifold, p)

Copy the value(s) from the point `p` on the [`AbstractManifold`](@ref) `M` into a new point.
See [`allocate_result`](@ref) for the allocation of new point memory and [`copyto!`](@ref) for the copying.
"""
function copy(M::AbstractManifold, p)
    q = allocate_result(M, copy, p)
    copyto!(M, q, p)
    return q
end
function copy(::AbstractManifold, p::Number)
    return copy(p)
end

@doc raw"""
    copy(M::AbstractManifold, p, X)

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
function copy(::AbstractManifold, p, X::Number)
    return copy(X)
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

"""
    default_type(M::AbstractManifold)

Get the default type of points on manifold `M`.
"""
default_type(M::AbstractManifold) = typeof(allocate_on(M))
"""
    default_type(M::AbstractManifold, ft::FiberType)

Get the default type of points from the fiber `ft` of the fiber bundle based on manifold `M`.
For example, call `default_type(MyManifold(), TangentSpaceType())` to get the default type
of a tangent vector.
"""
default_type(M::AbstractManifold, ft::FiberType) = typeof(allocate_on(M, ft))

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
    distance(M::AbstractManifold, p, q, m::AbstractInverseRetractionMethod)

Approximate distance between points `p` and `q` on manifold `M` using
[`AbstractInverseRetractionMethod`](@ref) `m`.
"""
function distance(M::AbstractManifold, p, q, m::AbstractInverseRetractionMethod)
    return norm(M, p, inverse_retract(M, p, q, m))
end
distance(M::AbstractManifold, p, q, ::LogarithmicInverseRetraction) = distance(M, p, q)

"""
    embed(M::AbstractManifold, p)

Embed point `p` from the [`AbstractManifold`](@ref) `M` into the ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `embed` includes changing data representation, if applicable, i.e.
if the points on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly.

The default is set in such a way that memory is allocated and `embed!(M, q, p)` is called.

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

The default is set in such a way that it assumes that the points on `M` are represented in
their embedding (for example like the unit vectors in a space to represent the sphere) and
hence embedding in the identity by default.

If you have more than one embedding, see [`EmbeddedManifold`](@ref) for defining a second
embedding. If your point `p` is already represented in some embedding,
see [`AbstractDecoratorManifold`](@ref) how you can avoid reimplementing code from the embedded manifold

See also: [`EmbeddedManifold`](@ref), [`project!`](@ref project!(M::AbstractManifold, q, p))
"""
embed!(M::AbstractManifold, q, p) = copyto!(M, q, p)

"""
    embed(M::AbstractManifold, p, X)

Embed a tangent vector `X` at a point `p` on the [`AbstractManifold`](@ref) `M` into an ambient space.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Not implementing this function means, there is no proper embedding for your tangent space(s).

Additionally, `embed` might include changing data representation, if applicable, i.e.
if tangent vectors on `M` are not represented in the same way as their counterparts in the
embedding, the representation is changed accordingly.

The default is set in such a way that memory is allocated and `embed!(M, Y, p. X)` is called.

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

The default is set in such a way that it assumes that the points on `M` are represented in
their embedding (for example like the unit vectors in a space to represent the sphere) and
hence embedding also for tangent vectors is the identity by default.

See also: [`EmbeddedManifold`](@ref), [`project!`](@ref project!(M::AbstractManifold, Y, p, X))
"""
embed!(M::AbstractManifold, Y, p, X) = copyto!(M, Y, p, X)

"""
    embed_project(M::AbstractManifold, p)

Embed `p` from manifold `M` an project it back to `M`. For points from `M` this is identity
but in case embedding is defined for points outside of `M`, this can serve as a way
to for example remove numerical inaccuracies caused by some algorithms.
"""
function embed_project(M::AbstractManifold, p)
    return project(M, embed(M, p))
end
"""
    embed_project(M::AbstractManifold, p, X)

Embed vector `X` tangent at `p` from manifold `M` an project it back to tangent space
at `p`. For points from that tangent space this is identity but in case embedding is
defined for tangent vectors from outside of it, this can serve as a way to for example remove
numerical inaccuracies caused by some algorithms.
"""
function embed_project(M::AbstractManifold, p, X)
    return project(M, p, embed(M, p, X))
end

function embed_project!(M::AbstractManifold, q, p)
    return project!(M, q, embed(M, p))
end
function embed_project!(M::AbstractManifold, Y, p, X)
    return project!(M, Y, p, embed(M, p, X))
end

"""
    has_components(M::AbstractManifold)

Return whether the [`AbstractManifold`](@ref)`(M)` consists of components,
like the [`PowerManifold`](@ref) or the [`ProductManifold`](@ref), that one can iterate over.
By default, this function returns `false`.
"""
has_components(M::AbstractManifold) = false

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
function _injectivity_radius(M::AbstractManifold, ::AbstractRetractionMethod)
    return _injectivity_radius(M, ExponentialRetraction())
end
function _injectivity_radius(M::AbstractManifold, p, ::ExponentialRetraction)
    return injectivity_radius_exp(M, p)
end
function _injectivity_radius(M::AbstractManifold, ::ExponentialRetraction)
    return injectivity_radius_exp(M)
end
injectivity_radius_exp(M::AbstractManifold, p) = injectivity_radius_exp(M)
injectivity_radius_exp(M::AbstractManifold) = injectivity_radius(M)

"""
    inner(M::AbstractManifold, p, X, Y)

Compute the inner product of tangent vectors `X` and `Y` at point `p` from the
[`AbstractManifold`](@ref) `M`.
"""
inner(M::AbstractManifold, p, X, Y)

"""
    is_flat(M::AbstractManifold)

Return true if the [`AbstractManifold`](@ref) `M` is flat, i.e. if its Riemann curvature
tensor is everywhere zero.
"""
is_flat(M::AbstractManifold)

"""
    isapprox(M::AbstractManifold, p, q; error::Symbol=:none, kwargs...)

Check if points `p` and `q` from [`AbstractManifold`](@ref) `M` are approximately equal.

The keyword argument can be used to get more information for the case that
the result is false, if the concrete manifold provides such information.
Currently the following are supported

* `:error` - throws an error if `isapprox` evaluates to false, providing possibly a more detailed error.
  Note that this turns `isapprox` basically to an `@assert`.
* `:info` – prints the information in an `@info`
* `:warn` – prints the information in an `@warn`
* `:none` (default) – the function just returns `true`/`false`

Keyword arguments can be used to specify tolerances.
"""
function isapprox(M::AbstractManifold, p, q; error::Symbol = :none, kwargs...)
    if error === :none
        return _isapprox(M, p, q; kwargs...)
    else
        ma = check_approx(M, p, q; kwargs...)
        if ma !== nothing
            (error === :error) && throw(ma)
            # else: collect and info showerror
            io = IOBuffer()
            showerror(io, ma)
            s = String(take!(io))
            (error === :info) && @info s
            (error === :warn) && @warn s
            return false
        end
        return true
    end
end

"""
    isapprox(M::AbstractManifold, p, X, Y; error:Symbol=:none; kwargs...)

Check if vectors `X` and `Y` tangent at `p` from [`AbstractManifold`](@ref) `M` are approximately
equal.

The optional positional argument can be used to get more information for the case that
the result is false, if the concrete manifold provides such information.
Currently the following are supported

* `:error` - throws an error if `isapprox` evaluates to false, providing possibly a more detailed error.
  Note that this turns `isapprox` basically to an `@assert`.
* `:info` – prints the information in an `@info`
* `:warn` – prints the information in an `@warn`
* `:none` (default) – the function just returns `true`/`false`

By default these informations are collected by calling [`check_approx`](@ref).

Keyword arguments can be used to specify tolerances.
"""
function isapprox(M::AbstractManifold, p, X, Y; error::Symbol = :none, kwargs...)
    if error === :none
        return _isapprox(M, p, X, Y; kwargs...)::Bool
    else
        mat = check_approx(M, p, X, Y; kwargs...)
        if mat !== nothing
            (error === :error) && throw(mat)
            # else: collect and info showerror
            io = IOBuffer()
            showerror(io, mat)
            s = String(take!(io))
            (error === :info) && @info s
            (error === :warn) && @warn s
            return false
        end
        return true
    end
end

"""
    _isapprox(M::AbstractManifold, p, q; kwargs...)

An internal function for testing whether points `p` and `q` from manifold `M` are
approximately equal. Returns either `true` or `false` and does not support errors like
[`isapprox`](@ref).

For more details see documentation of [`check_approx`](@ref).
"""
function _isapprox(M::AbstractManifold, p, q; kwargs...)
    return isapprox(p, q; kwargs...)
end
"""
    _isapprox(M::AbstractManifold, p, X, Y; kwargs...)

An internal function for testing whether tangent vectors `X` and `Y` from tangent space
at point `p` from manifold `M` are approximately equal. Returns either `true` or `false`
and does not support errors like [`isapprox`](@ref).

For more details see documentation of [`check_approx`](@ref).
"""
function _isapprox(M::AbstractManifold, p, X, Y; kwargs...)
    return isapprox(X, Y; kwargs...)
end

"""
    is_point(M::AbstractManifold, p; error::Symbol = :none, kwargs...)
    is_point(M::AbstractManifold, p, throw_error::Bool; kwargs...)

Return whether `p` is a valid point on the [`AbstractManifold`](@ref) `M`.
By default the function calls [`check_point`](@ref), which returns an `ErrorException` or `nothing`.

How to report a potential error can be set using the `error=` keyword

* `:error`          - throws an error if `p` is not a point
* `:info`           - displays the error message as an `@info`
* `:warn`           - displays the error message as a `@warning`
* `:none` (default) – the function just returns `true`/`false`

all other symbols are equivalent to `error=:none`.

The second signature is a shorthand, where the boolean is used for `error=:error` (`true`)
and `error=:none` (default, `false`). This case ignores the `error=` keyword
"""
function is_point(
    M::AbstractManifold,
    p,
    throw_error::Bool;
    error::Symbol = :none,
    kwargs...,
)
    return is_point(M, p; error = throw_error ? :error : :none, kwargs...)
end
function is_point(M::AbstractManifold, p; error::Symbol = :none, kwargs...)
    mps = check_size(M, p)
    if mps !== nothing
        (error === :error) && throw(mps)
        if (error === :info) || (error === :warn)
            # else: collect and info showerror
            io = IOBuffer()
            showerror(io, mps)
            s = String(take!(io))
            (error === :info) && @info s
            (error === :warn) && @warn s
        end
        return false
    end
    mpe = check_point(M, p; kwargs...)
    if mpe !== nothing
        (error === :error) && throw(mpe)
        if (error === :info) || (error === :warn)
            # else: collect and info showerror
            io = IOBuffer()
            showerror(io, mpe)
            s = String(take!(io))
            (error === :info) && @info s
            (error === :warn) && @warn s
        end
        return false
    end
    return true
end


"""
    is_vector(M::AbstractManifold, p, X, check_base_point::Bool=true; error::Symbol=:none, kwargs...)
    is_vector(M::AbstractManifold, p, X, check_base_point::Bool=true, throw_error::Boolean; kwargs...)

Return whether `X` is a valid tangent vector at point `p` on the [`AbstractManifold`](@ref) `M`.
Returns either `true` or `false`.

If `check_base_point` is set to true, this function also (first) calls [`is_point`](@ref)
on `p`.
Then, the function calls [`check_vector`](@ref) and checks whether the returned
value is `nothing` or an error.

How to report a potential error can be set using the `error=` keyword

* `:error`          - throws an error if `X` is not a tangent vector and/or `p` is not point
^ `:info`           - displays the error message as an `@info`
* `:warn`           - displays the error message as a `@warn`ing.
* `:none`           - (default) the function just returns `true`/`false`

all other symbols are equivalent to `error=:none`

The second signature is a shorthand, where `throw_error` is used for `error=:error` (`true`)
and `error=:none` (default, `false`). This case ignores the `error=` keyword.
"""
function is_vector(
    M::AbstractManifold,
    p,
    X,
    check_base_point::Bool,
    throw_error::Bool;
    error::Symbol = :none,
    kwargs...,
)
    return is_vector(
        M,
        p,
        X,
        check_base_point;
        error = throw_error ? :error : :none,
        kwargs...,
    )
end
function is_vector(
    M::AbstractManifold,
    p,
    X,
    check_base_point::Bool = true;
    error::Symbol = :none,
    kwargs...,
)
    if check_base_point
        # if error, is_point throws, otherwise if not a point return false
        !is_point(M, p; error = error, kwargs...) && return false
    end
    mXs = check_size(M, p, X)
    if mXs !== nothing
        (error === :error) && throw(mXs)
        if (error === :info) || (error === :warn)
            # else: collect and info showerror
            io = IOBuffer()
            showerror(io, mXs)
            s = String(take!(io))
            (error === :info) && @info s
            (error === :warn) && @warn s
        end
        return false
    end
    mXe = check_vector(M, p, X; kwargs...)
    if mXe !== nothing
        (error === :error) && throw(mXe)
        if (error === :info) || (error === :warn)
            # else: collect and info showerror
            io = IOBuffer()
            showerror(io, mXe)
            s = String(take!(io))
            (error === :info) && @info s
            (error === :warn) && @warn s
        end
        return false
    end
    return true
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

"""
    norm(M::AbstractManifold, p, X)

Compute the norm of tangent vector `X` at point `p` from a [`AbstractManifold`](@ref) `M`.
By default this is computed using [`inner`](@ref).
"""
norm(M::AbstractManifold, p, X) = sqrt(max(real(inner(M, p, X, X)), 0))

"""
    number_eltype(x)

Numeric element type of the a nested representation of a point or a vector.
To be used in conjunction with [`allocate`](@ref) or [`allocate_result`](@ref).
"""
number_eltype(x) = eltype(x)
@inline function number_eltype(x::AbstractArray)
    return typeof(mapreduce(eti -> one(number_eltype(eti)), +, x))
end
@inline number_eltype(::AbstractArray{T}) where {T<:Number} = T
@inline function number_eltype(::Type{<:AbstractArray{T}}) where {T}
    return number_eltype(T)
end
@inline function number_eltype(::Type{<:AbstractArray{T}}) where {T<:Number}
    return T
end
@inline function number_eltype(x::Tuple)
    @inline eti_to_one(eti) = one(number_eltype(eti))
    return typeof(sum(map(eti_to_one, x)))
end

"""
    Random.rand(M::AbstractManifold, [d::Integer]; vector_at=nothing)
    Random.rand(rng::AbstractRNG, M::AbstractManifold, [d::Integer]; vector_at=nothing)


Generate a random point on manifold `M` (when `vector_at` is `nothing`) or a tangent
vector at point `vector_at` (when it is not `nothing`).

Optionally a random number generator `rng` to be used can be specified. An optional integer
`d` indicates that a vector of `d` points or tangent vectors is to be generated.

!!! note

    Usually a uniform distribution should be expected for compact manifolds and a
    Gaussian-like distribution for non-compact manifolds and tangent vectors, although it is
    not guaranteed. The distribution may change between releases.

    `rand` methods for specific manifolds may take additional keyword arguments.

"""
Random.rand(M::AbstractManifold)

function Random.rand(M::AbstractManifold, d::Integer; kwargs...)
    return [rand(M; kwargs...) for _ in 1:d]
end
function Random.rand(rng::AbstractRNG, M::AbstractManifold, d::Integer; kwargs...)
    return [rand(rng, M; kwargs...) for _ in 1:d]
end
function Random.rand(M::AbstractManifold; vector_at = nothing, kwargs...)
    if vector_at === nothing
        pX = allocate_result(M, rand)
    else
        pX = allocate_result(M, rand, vector_at)
    end
    rand!(M, pX; vector_at = vector_at, kwargs...)
    return pX
end
function Random.rand(rng::AbstractRNG, M::AbstractManifold; vector_at = nothing, kwargs...)
    if vector_at === nothing
        pX = allocate_result(M, rand)
    else
        pX = allocate_result(M, rand, vector_at)
    end
    rand!(rng, M, pX; vector_at = vector_at, kwargs...)
    return pX
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

@doc raw"""
    riemann_tensor(M::AbstractManifold, p, X, Y, Z)

Compute the value of the Riemann tensor ``R(X_f,Y_f)Z_f`` at point `p`, where
``X_f``, ``Y_f`` and ``Z_f`` are vector fields defined by parallel transport of,
respectively, `X`, `Y` and `Z` to the desired point. All computations are performed
using the connection associated to manifold `M`.

The formula reads ``R(X_f,Y_f)Z_f = \nabla_X\nabla_Y Z - \nabla_Y\nabla_X Z - \nabla_{[X, Y]}Z``,
where ``[X, Y]`` is the Lie bracket of vector fields.

Note that some authors define this quantity with inverse sign.
"""
function riemann_tensor(M::AbstractManifold, p, X, Y, Z)
    Xresult = allocate_result(M, riemann_tensor, X)
    return riemann_tensor!(M, Xresult, p, X, Y, Z)
end

@doc raw"""
    sectional_curvature(M::AbstractManifold, p, X, Y)

Compute the sectional curvature of a manifold ``\mathcal M`` at a point ``p \in \mathcal M``
on two linearly independent tangent vectors at ``p``.
The formula reads
```math

    \kappa_p(X, Y) = \frac{⟨R(X, Y, Y), X⟩_p}{\lVert X \rVert^2_p \lVert Y \rVert^2_p - ⟨X, Y⟩^2_p}

```
where ``R(X, Y, Y)`` is the [`riemann_tensor`](@ref) on ``\mathcal M``.

# Input

* `M`:   a manifold ``\mathcal M``
* `p`:   a point ``p \in \mathcal M``
* `X`:   a tangent vector ``X \in T_p \mathcal M``
* `Y`:   a tangent vector ``Y \in T_p \mathcal M``

"""
function sectional_curvature(M::AbstractManifold, p, X, Y)
    R = riemann_tensor(M, p, X, Y, Y)
    return inner(M, p, R, X) / ((norm(M, p, X)^2 * norm(M, p, Y)^2) - (inner(M, p, X, Y)^2))
end

@doc raw"""
    sectional_curvature_max(M::AbstractManifold)

Upper bound on sectional curvature of manifold `M`. The formula reads
```math
\omega = \operatorname{sup}_{p\in\mathcal M, X\in T_p\mathcal M, Y\in T_p\mathcal M, ⟨X, Y⟩ ≠ 0} \kappa_p(X, Y)
```
"""
sectional_curvature_max(M::AbstractManifold)

@doc raw"""
    sectional_curvature_min(M::AbstractManifold)

Lower bound on sectional curvature of manifold `M`. The formula reads
```math
\omega = \operatorname{inf}_{p\in\mathcal M, X\in T_p\mathcal M, Y\in T_p\mathcal M, ⟨X, Y⟩ ≠ 0} \kappa_p(X, Y)
```
"""
sectional_curvature_min(M::AbstractManifold)

"""
    size_to_tuple(::Type{S}) where S<:Tuple

Converts a size given by `Tuple{N, M, ...}` into a tuple `(N, M, ...)`.
"""
Base.@pure size_to_tuple(::Type{S}) where {S<:Tuple} = tuple(S.parameters...)

@doc raw"""
    Weingarten!(M, Y, p, X, V)

Compute the Weingarten map ``\mathcal W_p\colon T_p\mathcal M × N_p\mathcal M \to T_p\mathcal M``
in place of `Y`, see [`Weingarten`](@ref).
"""
Weingarten!(M::AbstractManifold, Y, p, X, V)

@doc raw"""
    Weingarten(M, p, X, V)

Compute the Weingarten map ``\mathcal W_p\colon T_p\mathcal M × N_p\mathcal M \to T_p\mathcal M``,
where ``N_p\mathcal M`` is the orthogonal complement of the tangent space ``T_p\mathcal M``
of the embedded submanifold ``\mathcal M``, where we denote the embedding by ``\mathcal E``.

The Weingarten map can be defined by restricting the differential of the orthogonal [`project`](@ref)ion
``\operatorname{proj}_{T_p\mathcal M}\colon T_p \mathcal E \to T_p\mathcal M`` with respect to the base point ``p``,
i.e. defining
```math
\mathcal P_X := D_p\operatorname{proj}_{T_p\mathcal M}(Y)[X],
\qquad Y \in T_p \mathcal E, X \in T_p\mathcal M,
```
the Weingarten map can be written as ``\mathcal W_p(X,V) = \mathcal P_X(V)``.

The Weingarten map is named after [Julius Weingarten](https://en.wikipedia.org/wiki/Julius_Weingarten) (1836–1910).
"""
function Weingarten(M::AbstractManifold, p, X, V)
    Y = copy(M, p, X)
    Weingarten!(M, Y, p, X, V)
    return Y
end


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
include("shooting.jl")
include("vector_spaces.jl")
include("point_vector_fallbacks.jl")
include("nested_trait.jl")
include("decorator_trait.jl")
include("numerical_checks.jl")

include("VectorFiber.jl")
include("TangentSpace.jl")
include("ValidationManifold.jl")
include("EmbeddedManifold.jl")
include("DefaultManifold.jl")
include("ProductManifold.jl")
include("PowerManifold.jl")

include("deprecated.jl")
#
#
# Init
# -----
function __init__()
    #
    # Error Hints
    #
    @static if isdefined(Base.Experimental, :register_error_hint) # COV_EXCL_LINE
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f === plot_slope
                print(
                    io,
                    """

                    `plot_slope` has to be implemented using your favourite plotting package.
                    A default is available when Plots.jl is added to the current environment.
                    To then get the plotting functionality activated, do
                    """,
                )
                printstyled(io, "`using Plots`"; color = :cyan)
            end
            if exc.f === find_best_slope_window
                print(
                    io,
                    """

                    `find_best_slope_window` has to be implemented using some statistics package
                    A default is available when Statistics.jl is added to the current environment.
                    To then get the functionality activated, do
                    """,
                )
                printstyled(io, "`using Statistics`"; color = :cyan)
            end
        end
    end
end
#
#
# Export
# ------
#
# (a) Manifolds and general types
export AbstractManifold,
    AbstractManifoldPoint,
    AbstractTangentVector,
    AbstractCotangentVector,
    TFVector,
    CoTFVector
export VectorSpaceFiber
export TangentSpace, TangentSpaceType
export CotangentSpace, CotangentSpaceType
export AbstractDecoratorManifold
export AbstractTrait, IsEmbeddedManifold, IsEmbeddedSubmanifold, IsIsometricEmbeddedManifold
export IsExplicitDecorator
export ValidationManifold,
    ValidationMPoint, ValidationTangentVector, ValidationCotangentVector
export EmbeddedManifold
export AbstractPowerManifold, PowerManifold
export AbstractPowerRepresentation,
    NestedPowerRepresentation, NestedReplacingPowerRepresentation
export ProductManifold

# (b) Generic Estimation Types

export GeodesicInterpolationWithinRadius,
    CyclicProximalPointEstimation,
    ExtrinsicEstimation,
    GradientDescentEstimation,
    WeiszfeldEstimation,
    AbstractApproximationMethod,
    GeodesicInterpolation


# (b) Retraction Types
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
    ProductRetraction,
    ProjectionRetraction,
    RetractionWithKeywords,
    SasakiRetraction,
    ShootingInverseRetraction,
    SoftmaxRetraction

# (c) Inverse Retraction Types
export AbstractInverseRetractionMethod,
    ApproximateInverseRetraction,
    CayleyInverseRetraction,
    EmbeddedInverseRetraction,
    InverseProductRetraction,
    LogarithmicInverseRetraction,
    NLSolveInverseRetraction,
    QRInverseRetraction,
    PadeInverseRetraction,
    PolarInverseRetraction,
    ProjectionInverseRetraction,
    InverseRetractionWithKeywords,
    SoftmaxInverseRetraction

# (d) Vector Transport Types
export AbstractVectorTransportMethod,
    DifferentiatedRetractionVectorTransport,
    EmbeddedVectorTransport,
    ParallelTransport,
    PoleLadderTransport,
    ProductVectorTransport,
    ProjectionTransport,
    ScaledVectorTransport,
    SchildsLadderTransport,
    VectorTransportDirection,
    VectorTransportTo,
    VectorTransportWithKeywords

# (e) Basis Types
export CachedBasis,
    DefaultBasis,
    DefaultOrthogonalBasis,
    DefaultOrthonormalBasis,
    DiagonalizingOrthonormalBasis,
    DefaultOrthonormalBasis,
    GramSchmidtOrthonormalBasis,
    ProjectedOrthonormalBasis,
    VeeOrthogonalBasis

# (f) Error Messages
export OutOfInjectivityRadiusError, ManifoldDomainError
export ApproximatelyError
export CompositeManifoldError, ComponentManifoldError, ManifoldDomainError

# (g) Functions on Manifolds
export ×,
    ℝ,
    ℂ,
    allocate,
    allocate_on,
    angle,
    base_manifold,
    base_point,
    change_basis,
    change_basis!,
    change_metric,
    change_metric!,
    change_representer,
    change_representer!,
    check_inverse_retraction,
    check_retraction,
    check_vector_transport,
    copy,
    copyto!,
    default_approximation_method,
    default_basis,
    default_inverse_retraction_method,
    default_retraction_method,
    default_type,
    default_vector_transport_method,
    distance,
    exp,
    exp!,
    embed,
    embed!,
    embed_project,
    embed_project!,
    fill,
    fill!,
    geodesic,
    geodesic!,
    get_basis,
    get_component,
    get_coordinates,
    get_coordinates!,
    get_embedding,
    get_vector,
    get_vector!,
    get_vectors,
    has_components,
    hat,
    hat!,
    injectivity_radius,
    inner,
    internal_value,
    inverse_retract,
    inverse_retract!,
    isapprox,
    is_flat,
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
    parallel_transport_direction,
    parallel_transport_direction!,
    parallel_transport_to,
    parallel_transport_to!,
    project,
    project!,
    real_dimension,
    representation_size,
    set_component!,
    shortest_geodesic,
    shortest_geodesic!,
    show,
    rand,
    rand!,
    retract,
    retract!,
    riemann_tensor,
    riemann_tensor!,
    sectional_curvature,
    sectional_curvature_max,
    sectional_curvature_min,
    vector_space_dimension,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
    vector_transport_to!,
    vee,
    vee!,
    Weingarten,
    Weingarten!,
    zero_vector,
    zero_vector!
end # module
