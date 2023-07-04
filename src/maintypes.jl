"""
    AbstractManifold{𝔽}

A type to represent a (Riemannian) manifold.
The [`AbstractManifold`](@ref) is a central type of this interface.
It allows to distinguish different implementations of functions like the [`exp`](@ref)onential and
[`log`](@ref)arithmic map for different manifolds.
Usually, the manifold is the first parameter in any of these functions within `ManifoldsBase.jl`.
Based on these, say “elementary” functions, as the two mentioned above, more general functions are built,
for example the [`shortest_geodesic`](@ref) and the [`geodesic`](@ref).
These should only be overwritten (reimplemented) if for a certain manifold specific, more efficient implementations are possible, that do not just call the elementary functions.

The [`AbstractManifold`] is parametrized by [`AbstractNumbers`](@ref) to distinguish for example
real (ℝ) and complex (ℂ) manifolds.

For subtypes the preferred order of parameters is: size and simple value parameters,
followed by the [`AbstractNumbers`](@ref) `field`, followed by data type parameters,
which might depend on the abstract number field type.
"""
abstract type AbstractManifold{𝔽} end

"""
    AbstractManifoldPoint

Type for a point on a manifold.
While an [`AbstractManifold`](@ref) does not necessarily require this
type, for example when it is implemented for `Vector`s or `Matrix` type elements, this type
can be used either

* for more complicated representations,
* semantic verification, or
* when dispatching on different representations of points on a manifold.

Since semantic verification and different representations usually might still only store a
matrix internally, it is possible to use [`@manifold_element_forwards`](@ref) and
[`@default_manifold_fallbacks`](@ref) to reduce implementation overhead.
"""
abstract type AbstractManifoldPoint end

"""
    abstract type AbstractManifoldParameter end

Abstract representation of numeric parameters for a manifold type. Can be either
[`TypeParameter`](@ref) or [`FieldParameter`](@ref).
"""
abstract type AbstractManifoldParameter end

"""
    TypeParameter{T}

Represents numeric parameters of a manifold type as type parameters, allowing for static
specialization of methods.
"""
struct TypeParameter{T} <: AbstractManifoldParameter end
TypeParameter(t::NTuple) = TypeParameter{t}()

"""
    FieldParameter{TS<:NTuple{N,Int} where N}

Represents numeric parameters of a manifold type as values in a field, allowing for
less static specialization of methods and faster TTFX.
"""
struct FieldParameter{TS<:NTuple{N,Int} where {N}} <: AbstractManifoldParameter
    parameter::TS
end

get_parameter(::TypeParameter{T}) where {T} = T
get_parameter(P::FieldParameter) = P.parameter
