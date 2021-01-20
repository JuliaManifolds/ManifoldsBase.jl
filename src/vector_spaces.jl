
"""
    FVector(type::VectorSpaceType, data, basis::AbstractBasis)

Decorator indicating that the vector `data` is from a fiber of a vector bundle
of type `type`. `basis` is an object describing the basis of that space.
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

Base.:+(X::FVector, Y::FVector) = FVector(X.type, X.data + Y.data, X.basis)

Base.:-(X::FVector, Y::FVector) = FVector(X.type, X.data - Y.data, X.basis)
Base.:-(X::FVector) = FVector(X.type, -X.data, X.basis)

Base.:*(a::Number, X::FVector) = FVector(X.type, a * X.data, X.basis)

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
