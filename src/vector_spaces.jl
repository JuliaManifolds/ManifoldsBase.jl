
"""
    FVector(type::VectorSpaceType, data, basis::AbstractBasis)

Decorator indicating that the vector `data` contains coordinates of a vector from a fiber
of a vector bundle of type `type`. `basis` is an object describing the basis of that space
in which the coordinates are given.

Conversion between `FVector` representation and the default representation of an object
(for example a tangent vector) for a manifold should be done using [`get_coordinates`](@ref)
and [`get_vector`](@ref).

# Examples

```julia-repl
julia> using Manifolds

julia> M = Sphere(2)
Sphere(2, ℝ)

julia> p = [1.0, 0.0, 0.0]
3-element Vector{Float64}:
 1.0
 0.0
 0.0

julia> X = [0.0, 2.0, -1.0]
3-element Vector{Float64}:
  0.0
  2.0
 -1.0

julia> B = DefaultOrthonormalBasis()
DefaultOrthonormalBasis(ℝ)

julia> fX = TFVector(get_coordinates(M, p, X, B), B)
TFVector([2.0, -1.0], DefaultOrthonormalBasis(ℝ))

julia> X_back = get_vector(M, p, fX.data, fX.basis)
3-element Vector{Float64}:
 -0.0
  2.0
 -1.0
```
"""
struct FVector{TType<:VectorSpaceType,TData,TBasis<:AbstractBasis}
    type::TType
    data::TData
    basis::TBasis
end

const TFVector = FVector{TangentSpaceType}
const CoTFVector = FVector{CotangentSpaceType}

function TFVector(data, basis::AbstractBasis)
    return TFVector{typeof(data),typeof(basis)}(TangentSpace, data, basis)
end
function CoTFVector(data, basis::AbstractBasis)
    return CoTFVector{typeof(data),typeof(basis)}(CotangentSpace, data, basis)
end

function Base.show(io::IO, fX::TFVector)
    return print(io, "TFVector(", fX.data, ", ", fX.basis, ")")
end
function Base.show(io::IO, fX::CoTFVector)
    return print(io, "CoTFVector(", fX.data, ", ", fX.basis, ")")
end


"""
    AbstractFibreVector{TType<:VectorSpaceType}

Type for a vector from a vector space (fibre of a vector bundle) of type `TType` of a manifold.
While a [`AbstractManifold`](@ref) does not necessarily require this type, for example when it is
implemented for `Vector`s or `Matrix` type elements, this type can be used for more
complicated representations, semantic verification, or even dispatch for different
representations of tangent vectors and their types on a manifold.
"""
abstract type AbstractFibreVector{TType<:VectorSpaceType} end

"""
    TVector = AbstractFibreVector{TangentSpaceType}

Type for a tangent vector of a manifold. While a [`AbstractManifold`](@ref) does not necessarily
require this type, for example when it is implemented for `Vector`s or `Matrix` type
elements, this type can be used for more complicated representations, semantic verification,
or even dispatch for different representations of tangent vectors and their types on a
manifold.
"""
const TVector = AbstractFibreVector{TangentSpaceType}

"""
    CoTVector = AbstractFibreVector{CotangentSpaceType}

Type for a cotangent vector of a manifold. While a [`AbstractManifold`](@ref) does not necessarily
require this type, for example when it is implemented for `Vector`s or `Matrix` type
elements, this type can be used for more complicated representations, semantic verification,
or even dispatch for different representations of cotangent vectors and their types on a
manifold.
"""
const CoTVector = AbstractFibreVector{CotangentSpaceType}

Base.:+(X::FVector, Y::FVector) = FVector(X.type, X.data + Y.data, X.basis)

Base.:-(X::FVector, Y::FVector) = FVector(X.type, X.data - Y.data, X.basis)
Base.:-(X::FVector) = FVector(X.type, -X.data, X.basis)

Base.:*(a::Number, X::FVector) = FVector(X.type, a * X.data, X.basis)

allocate(x::FVector) = FVector(x.type, allocate(x.data), x.basis)
allocate(x::FVector, ::Type{T}) where {T} = FVector(x.type, allocate(x.data, T), x.basis)

function Base.copyto!(X::FVector, Y::FVector)
    copyto!(X.data, Y.data)
    return X
end

function number_eltype(
    ::Type{FVector{TType,TData,TBasis}},
) where {TType<:VectorSpaceType,TData,TBasis}
    return number_eltype(TData)
end
number_eltype(v::FVector) = number_eltype(v.data)

Base.show(io::IO, ::TangentSpaceType) = print(io, "TangentSpace")
Base.show(io::IO, ::CotangentSpaceType) = print(io, "CotangentSpace")

"""
    vector_space_dimension(M::AbstractManifold, V::VectorSpaceType)

Dimension of the vector space of type `V` on manifold `M`.
"""
vector_space_dimension(::AbstractManifold, ::VectorSpaceType)

function vector_space_dimension(M::AbstractManifold, ::TCoTSpaceType)
    return manifold_dimension(M)
end
