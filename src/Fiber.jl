"""
    abstract type FiberType end

An abstract type for fiber types that can be used within [`Fiber`](@ref).
"""
abstract type FiberType end

@doc raw"""
    Fiber{𝔽,TFiber<:FiberType,TM<:AbstractManifold,TX} <: AbstractManifold{𝔽}

A fiber of a fiber bundle at a point `p` on the manifold.

This fiber itself is also a `manifold`. For vector fibers it's by default flat and hence
isometric to the [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) manifold.

# Fields

* `manifold`    – base space of the fiber bundle
* `point`       – a point ``p`` from the base space; the fiber corresponds to the preimage
                  by bundle projection ``\pi^{-1}(\{p\})``.
* `fiber_type`  – the [`FiberType`](@ref) of the fiber


# Constructor

    Fiber(M::AbstractManifold, p, fiber_type::FiberType; field::AbstractNumbers = ℝ)

A fiber of type `fiber_type` at point `p` from the manifold `M`.
The number system `field` of the fiber defaults to `ℝ`.
"""
struct Fiber{𝔽, TFiber <: FiberType, TM <: AbstractManifold, TX} <: AbstractManifold{𝔽}
    manifold::TM
    point::TX
    fiber_type::TFiber
end

function Fiber(
        manifold::TM,
        point::TX,
        fiber_type::TFiber;
        field::AbstractNumbers = ℝ,
    ) where {TM <: AbstractManifold, TX, TFiber <: FiberType}
    return Fiber{field, TFiber, TM, TX}(manifold, point, fiber_type)
end

base_manifold(B::Fiber) = B.manifold

function Base.show(io::IO, ::MIME"text/plain", vs::Fiber)
    summary(io, vs)
    println(io, "\nFiber:")
    pre = " "
    sf = sprint(show, "text/plain", vs.fiber_type; context = io, sizehint = 0)
    sf = replace(sf, '\n' => "\n$(pre)")
    sm = sprint(show, "text/plain", vs.manifold; context = io, sizehint = 0)
    sm = replace(sm, '\n' => "\n$(pre)")
    println(io, pre, sf)
    println(io, pre, sm)
    println(io, "Base point:")
    sp = sprint(show, "text/plain", vs.point; context = io, sizehint = 0)
    sp = replace(sp, '\n' => "\n$(pre)")
    return print(io, pre, sp)
end
