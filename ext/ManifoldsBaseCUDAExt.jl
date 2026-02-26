"""
    ManifoldsBaseCUDAExt

CUDA extension for ManifoldsBase.jl providing GPU-aware allocation via
`allocate_on(M, CuArray{T})` for explicit GPU point allocation.

Base `allocate` methods already preserve CuArray types via `similar()`,
so no `allocate` overrides are needed.
"""
module ManifoldsBaseCUDAExt

using ManifoldsBase
using ManifoldsBase:
    AbstractManifold,
    PowerManifold,
    TangentSpaceType,
    representation_size,
    get_iterator
import ManifoldsBase: allocate_on

using CUDA
using LinearAlgebra

# Type union covering CuArray and common wrapped variants for dispatch.
const AnyCuArray{T,N} = Union{
    CuArray{T,N},
    LinearAlgebra.Transpose{T,<:CuArray},
    LinearAlgebra.Adjoint{T,<:CuArray},
    SubArray{T,N,<:CuArray},
}

# GPU-aware allocate_on: base uses `similar(Array{Float64}, dims)` which
# doesn't handle CuArray{T} type specs. These let users request GPU
# allocation explicitly via `allocate_on(M, CuArray{Float64})`.

function ManifoldsBase.allocate_on(M::AbstractManifold, ::Type{CuArray{T}}) where {T}
    return CUDA.zeros(T, representation_size(M)...)
end
function ManifoldsBase.allocate_on(
    M::AbstractManifold, ::Type{CuArray{T,N}}
) where {T,N}
    return CUDA.zeros(T, representation_size(M)...)
end

# Tangent space variants.
function ManifoldsBase.allocate_on(
    M::AbstractManifold, ::TangentSpaceType, ::Type{CuArray{T}}
) where {T}
    return CUDA.zeros(T, representation_size(M)...)
end
function ManifoldsBase.allocate_on(
    M::AbstractManifold, ::TangentSpaceType, ::Type{CuArray{T,N}}
) where {T,N}
    return CUDA.zeros(T, representation_size(M)...)
end

# PowerManifold variants: return Vector{CuArray} for nested representation.
function ManifoldsBase.allocate_on(M::PowerManifold, ::Type{CuArray{T}}) where {T}
    return [allocate_on(M.manifold, CuArray{T}) for _ in get_iterator(M)]
end
function ManifoldsBase.allocate_on(M::PowerManifold, ::Type{CuArray{T,N}}) where {T,N}
    return [allocate_on(M.manifold, CuArray{T}) for _ in get_iterator(M)]
end

end # module
