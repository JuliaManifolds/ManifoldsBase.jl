"""
    ManifoldsBaseCUDAExt

CUDA extension for ManifoldsBase.jl, enabling GPU-accelerated manifold operations
with `CuArray`-backed points and tangent vectors.

## What this extension does

The core allocation system in ManifoldsBase has two paths:

1. **With point reference**: `allocate_result(M, f, p, X)` → `allocate(M, p, T)` → `similar(p, T)`
   This already preserves `CuArray` types — no fix needed.

2. **Without point reference**: `allocate_result(M, f)` → `allocate_result_array(M, f, T, rs)` → `Array{T}(undef, rs...)`
   This always creates CPU arrays. It is called when constructing solver state
   (e.g. `ArmijoLinesearchStepsize.candidate_point`) before any iterate is available.

Path 2 cannot be fixed purely at this layer because there is no `CuArray` in the
call signature to dispatch on. It must be fixed at the solver level (see ManoptCUDAExt).

This extension instead provides:
- `allocate_on(M, ::Type{<:CuArray})` for explicit GPU point allocation
- `allocate` overrides ensuring correct behavior for `CuArray` edge cases
- `copyto!` methods for mixed CPU/GPU manifold data transfers
"""
module ManifoldsBaseCUDAExt

using ManifoldsBase
using ManifoldsBase:
    AbstractManifold,
    PowerManifold,
    ProductManifold,
    TangentSpaceType,
    number_eltype,
    representation_size,
    get_iterator
import ManifoldsBase: allocate, allocate_on

using CUDA
using LinearAlgebra

# === CuArray type union (following OMEinsum.jl pattern) ===

"""
    AnyCuArray{T,N}

Type union covering `CuArray` and common wrapped variants (`Transpose`, `Adjoint`,
`SubArray`) for robust dispatch in GPU-aware methods.
"""
const AnyCuArray{T,N} = Union{
    CuArray{T,N},
    LinearAlgebra.Transpose{T,<:CuArray},
    LinearAlgebra.Adjoint{T,<:CuArray},
    SubArray{T,N,<:CuArray},
}

# === GPU-aware allocate_on ===
#
# `allocate_on(M)` returns `similar(Array{Float64}, representation_size(M))` by default.
# These overrides let users request GPU allocation explicitly:
#   p = allocate_on(M, CuArray{Float32})

function ManifoldsBase.allocate_on(M::AbstractManifold, ::Type{CuArray{T}}) where {T}
    rs = representation_size(M)
    return CUDA.zeros(T, rs...)
end
function ManifoldsBase.allocate_on(
    M::AbstractManifold, ::Type{CuArray{T,N}}
) where {T,N}
    rs = representation_size(M)
    return CUDA.zeros(T, rs...)
end

# Tangent space variants
function ManifoldsBase.allocate_on(
    M::AbstractManifold, ::TangentSpaceType, ::Type{CuArray{T}}
) where {T}
    rs = representation_size(M)
    return CUDA.zeros(T, rs...)
end
function ManifoldsBase.allocate_on(
    M::AbstractManifold, ::TangentSpaceType, ::Type{CuArray{T,N}}
) where {T,N}
    rs = representation_size(M)
    return CUDA.zeros(T, rs...)
end

# PowerManifold variants — return Vector{CuArray} (nested representation)
function ManifoldsBase.allocate_on(M::PowerManifold, ::Type{CuArray{T}}) where {T}
    return [allocate_on(M.manifold, CuArray{T}) for _ in get_iterator(M)]
end
function ManifoldsBase.allocate_on(M::PowerManifold, ::Type{CuArray{T,N}}) where {T,N}
    return [allocate_on(M.manifold, CuArray{T}) for _ in get_iterator(M)]
end

# === GPU-aware allocate ===
#
# The base `allocate(a) = similar(a)` already preserves CuArray.
# These handle edge cases where the manifold or type arguments matter.

# Allocate on manifold from CuArray prototype — ensure CuArray output
function ManifoldsBase.allocate(::AbstractManifold, p::CuArray{T,N}) where {T,N}
    return CuArray{T,N}(undef, size(p))
end
function ManifoldsBase.allocate(
    ::AbstractManifold, p::CuArray{T,N}, ::Type{S}
) where {T,N,S}
    return CuArray{S,N}(undef, size(p))
end

# Without manifold argument (used by some ManifoldsBase internal paths)
function ManifoldsBase.allocate(p::CuArray{T,N}, ::Type{S}) where {T,N,S}
    return CuArray{S,N}(undef, size(p))
end

# Nested arrays of CuArrays (e.g. PowerManifold with NestedReplacing)
function ManifoldsBase.allocate(a::AbstractArray{<:CuArray})
    return map(allocate, a)
end
function ManifoldsBase.allocate(a::AbstractArray{<:CuArray}, ::Type{S}) where {S}
    return map(t -> allocate(t, S), a)
end

# === copyto! for mixed CPU/GPU transfers ===
#
# When manifold operations produce CPU results that need to be stored in GPU arrays
# (or vice versa), these methods handle the transfer transparently.

function Base.copyto!(
    M::AbstractManifold, dest::CuArray, src::Array
)
    copyto!(dest, src)
    return dest
end
function Base.copyto!(
    M::AbstractManifold, dest::Array, src::CuArray
)
    copyto!(dest, src)
    return dest
end

end # module
