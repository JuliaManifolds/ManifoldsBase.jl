"""
    AbstractManifold{F}

A manifold type. The `AbstractManifold` is used to dispatch to different functions on a manifold,
usually as the first argument of the function. Examples are the [`exp`](@ref)onential and
[`log`](@ref)arithmic maps as well as more general functions that are built on them like the
[`geodesic`](@ref).

The manifold is parametrized by an [`AbstractNumbers`](@ref) to distinguish for example
real (ℝ) and complex (ℂ) manifolds.

For subtypes the preferred order of parameters is: size and simple value parameters,
followed by the [`AbstractNumbers`](@ref) `field`, followed by data type parameters,
which might depend on the abstract number field type.

For more details see [interface-types-and-functions](@ref) in the ManifoldsBase.jl documentation at
[https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html#Types-and-functions](https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html#Types-and-functions).
"""
abstract type AbstractManifold{𝔽} end

"""
    AbstractManifoldPoint

Type for a point on a manifold. While a [`AbstractManifold`](@ref) does not necessarily require this
type, for example when it is implemented for `Vector`s or `Matrix` type elements, this type
can be used either
* for more complicated representations,
* semantic verification, or
* even dispatch for different representations of points on a manifold.

Since semantic verification and different representations usually might still only store a
matrix internally, the default just falls back to using `p.value` of an [`AbstractManifold`](@ref).
This way, introducing a type for a default (`Vector` or `Matrix`-based) implementation does not introduce implementation overhead
"""
abstract type AbstractManifoldPoint end

Base.:(==)(p::T, q::T) where {T<:AbstractManifoldPoint} = (p.value == q.value)

Base.axes(p::T) where {T<:AbstractManifoldPoint} = axes(p.value)

function Broadcast.BroadcastStyle(::Type{T}) where {T<:AbstractManifoldPoint}
    return Broadcast.Style{T}()
end
function Broadcast.BroadcastStyle(
    ::Broadcast.AbstractArrayStyle{0},
    b::Broadcast.Style{T},
) where {T<:AbstractManifoldPoint}
    return b
end

@inline Base.copy(p::T) where {T<:AbstractManifoldPoint} = T(copy(p.value))

function Base.copyto!(q::T, p::T) where {T<:AbstractManifoldPoint}
    copyto!(q.value, p.value)
    return q
end

function Broadcast.instantiate(
    bc::Broadcast.Broadcasted{Broadcast.Style{T},Nothing},
) where {T<:AbstractManifoldPoint}
    return bc
end
function Broadcast.instantiate(
    bc::Broadcast.Broadcasted{Broadcast.Style{T}},
) where {T<:AbstractManifoldPoint}
    Broadcast.check_broadcast_axes(bc.axes, bc.args...)
    return bc
end

Broadcast.broadcastable(v::T) where {T<:AbstractManifoldPoint} = v

@inline function Base.copy(
    bc::Broadcast.Broadcasted{Broadcast.Style{T}},
) where {T<:AbstractManifoldPoint}
    return T(Broadcast._broadcast_getindex(bc, 1))
end

Base.@propagate_inbounds function Broadcast._broadcast_getindex(
    v::T,
    I,
) where {T<:AbstractManifoldPoint}
    return v.value
end

Base.similar(p::T) where {T<:AbstractManifoldPoint} = T(similar(p.value))

@inline function Base.copyto!(
    dest::T,
    bc::Broadcast.Broadcasted{Broadcast.Style{T}},
) where {T<:AbstractManifoldPoint}
    axes(dest) == axes(bc) || Broadcast.throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{T} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Broadcast.preprocess(dest, bc)
    # Performance may vary depending on whether `@inbounds` is placed outside the
    # for loop or not. (cf. https://github.com/JuliaLang/julia/issues/38086)
    copyto!(dest.value, bc′[1])
    return dest
end
