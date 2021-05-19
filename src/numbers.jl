"""
    AbstractNumbers

An abstract type to represent the number system on which a manifold is built.

This provides concrete number types for dispatch. The two most common number types are
the fields [`RealNumbers`](@ref) (`ℝ` for short) and [`ComplexNumbers`](@ref) (`ℂ`).
"""
abstract type AbstractNumbers end

"""
    RealNumbers <: AbstractNumbers
    ℝ = RealNumbers()

The field of real numbers.
"""
struct RealNumbers <: AbstractNumbers end

"""
    ComplexNumbers <: AbstractNumbers
    ℂ = ComplexNumbers()

The field of complex numbers.
"""
struct ComplexNumbers <: AbstractNumbers end

"""
    QuaternionNumbers <: AbstractNumbers
    ℍ = QuaternionNumbers()

The division algebra of quaternions.
"""
struct QuaternionNumbers <: AbstractNumbers end

const ℝ = RealNumbers()
const ℂ = ComplexNumbers()
const ℍ = QuaternionNumbers()

"""
    _unify_number_systems(𝔽s::AbstractNumbers...)

Compute a number system that includes all given number systems (as sub-systems) and is
closed under addition and multiplication.
"""
function _unify_number_systems(a::AbstractNumbers, rest::AbstractNumbers...)
    return _unify_number_systems(a, _unify_number_systems(rest...))
end
_unify_number_systems(𝔽::AbstractNumbers) = 𝔽
_unify_number_systems(r::RealNumbers, ::RealNumbers) = r
_unify_number_systems(::RealNumbers, c::ComplexNumbers) = c
_unify_number_systems(::RealNumbers, q::QuaternionNumbers) = q
_unify_number_systems(c::ComplexNumbers, ::RealNumbers) = c
_unify_number_systems(c::ComplexNumbers, ::ComplexNumbers) = c
_unify_number_systems(::ComplexNumbers, q::QuaternionNumbers) = q
_unify_number_systems(q::QuaternionNumbers, ::RealNumbers) = q
_unify_number_systems(q::QuaternionNumbers, ::ComplexNumbers) = q
_unify_number_systems(q::QuaternionNumbers, ::QuaternionNumbers) = q

Base.show(io::IO, ::RealNumbers) = print(io, "ℝ")
Base.show(io::IO, ::ComplexNumbers) = print(io, "ℂ")
Base.show(io::IO, ::QuaternionNumbers) = print(io, "ℍ")

@doc raw"""
    real_dimension(𝔽::AbstractNumbers)

Return the real dimension $\dim_ℝ 𝔽$ of the [`AbstractNumbers`](@ref) system `𝔽`.
The real dimension is the dimension of a real vector space with which a number in `𝔽` can be
identified.
For example, [`ComplexNumbers`](@ref) have a real dimension of 2, and
[`QuaternionNumbers`](@ref) have a real dimension of 4.
"""
function real_dimension(𝔽::AbstractNumbers)
    return error("real_dimension not defined for number system $(𝔽)")
end
real_dimension(::RealNumbers) = 1
real_dimension(::ComplexNumbers) = 2
real_dimension(::QuaternionNumbers) = 4

@doc raw"""
    number_system(M::AbstractManifold{𝔽})

Return the number system the manifold `M` is based on, i.e. the parameter `𝔽`.
"""
number_system(M::AbstractManifold{𝔽}) where {𝔽} = 𝔽
