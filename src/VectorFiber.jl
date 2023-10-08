
"""
    VectorSpaceFiberType{TVS<:VectorSpaceType}

`FiberType` of a [`FiberBundle`](@ref) corresponding to [`VectorSpaceType`](@ref)
`TVS`.
"""
struct VectorSpaceFiberType{TVS<:VectorSpaceType} <: FiberType
    fiber::TVS
end

function Base.show(io::IO, vsf::VectorSpaceFiberType)
    return print(io, "VectorSpaceFiberType($(vsf.fiber))")
end

const TangentFiberType = VectorSpaceFiberType{TangentSpaceType}

TangentFiberType() = VectorSpaceFiberType(TangentSpaceType())

"""
    VectorSpaceFiber{𝔽,M,TFiber}

Alias for [`Fiber`](@ref) when the fiber is a vector space.
"""
const VectorSpaceFiber{𝔽,M,TSpaceType} = Fiber{
    𝔽,
    VectorSpaceFiberType{TSpaceType},
    M,
} where {𝔽,M<:AbstractManifold{𝔽},TSpaceType<:VectorSpaceType}

function VectorSpaceFiber(M::AbstractManifold, fiber::VectorSpaceFiberType, p)
    return Fiber(fiber, M, p)
end
function VectorSpaceFiber(M::AbstractManifold, vs::VectorSpaceType, p)
    return VectorSpaceFiber(M, VectorSpaceFiberType(vs), p)
end

"""
    TangentSpace{𝔽,M}

Alias for [`VectorSpaceFiber`](@ref) for the tangent space at a point.
"""
const TangentSpace{𝔽,M} =
    VectorSpaceFiber{𝔽,M,TangentSpaceType} where {𝔽,M<:AbstractManifold{𝔽}}

"""
    TangentSpace(M::AbstractManifold, p)

Return an object of type [`VectorSpaceFiber`](@ref) representing tangent
space at `p` on the [`AbstractManifold`](@ref) `M`.
"""
TangentSpace(M::AbstractManifold, p) = VectorSpaceFiber(M, TangentFiberType(), p)

function allocate_result(M::TangentSpace, ::typeof(rand))
    return zero_vector(M.manifold, M.point)
end

"""
    distance(M::TangentSpace, p, q)

Distance between vectors `p` and `q` from the vector space `M`. It is calculated as the norm
of their difference.
"""
function distance(M::TangentSpace, p, q)
    return norm(M.manifold, M.point, q - p)
end

function embed!(M::TangentSpace, q, p)
    return embed!(M.manifold, q, M.point, p)
end
function embed!(M::TangentSpace, Y, p, X)
    return embed!(M.manifold, Y, M.point, X)
end

@doc raw"""
    exp(M::TangentSpace, p, X)

Exponential map of tangent vectors `X` and `p` from the tangent space `M`. It is
calculated as their sum.
"""
exp(::TangentSpace, ::Any, ::Any)

function exp!(M::TangentSpace, q, p, X)
    copyto!(M.manifold, q, p + X)
    return q
end

fiber_dimension(M::AbstractManifold, V::VectorSpaceFiberType) = fiber_dimension(M, V.fiber)
fiber_dimension(M::AbstractManifold, ::CotangentSpaceType) = manifold_dimension(M)
fiber_dimension(M::AbstractManifold, ::TangentSpaceType) = manifold_dimension(M)

function get_basis(M::TangentSpace, p, B::CachedBasis)
    return invoke(
        get_basis,
        Tuple{AbstractManifold,Any,CachedBasis},
        M.manifold,
        M.point,
        B,
    )
end
function get_basis(M::TangentSpace, p, B::AbstractBasis{<:Any,TangentSpaceType})
    return get_basis(M.manifold, M.point, B)
end

function get_coordinates(M::TangentSpace, p, X, B::AbstractBasis)
    return get_coordinates(M.manifold, M.point, X, B)
end

function get_coordinates!(M::TangentSpace, Y, p, X, B::AbstractBasis)
    return get_coordinates!(M.manifold, Y, M.point, X, B)
end

function get_vector(M::TangentSpace, p, X, B::AbstractBasis)
    return get_vector(M.manifold, M.point, X, B)
end

function get_vector!(M::TangentSpace, Y, p, X, B::AbstractBasis)
    return get_vector!(M.manifold, Y, M.point, X, B)
end

function get_vectors(M::TangentSpace, p, B::CachedBasis)
    return get_vectors(M.manifold, M.point, B)
end

@doc raw"""
    injectivity_radius(M::TangentSpace)

Return the injectivity radius on the [`TangentSpace`](@ref) `M`, which is $∞$.
"""
injectivity_radius(::TangentSpace) = Inf

"""
    inner(M::TangentSpace, p, X, Y)

Inner product of vectors `X` and `Y` from the tangent space at `M`.
"""
function inner(M::TangentSpace, p, X, Y)
    return inner(M.manifold, M.point, X, Y)
end

"""
    is_flat(::TangentSpace)

Return true. [`TangentSpace`](@ref) is a flat manifold.
"""
is_flat(::TangentSpace) = true

function _isapprox(M::TangentSpace, X, Y; kwargs...)
    return isapprox(M.manifold, M.point, X, Y; kwargs...)
end

"""
    log(M::TangentSpace, p, q)

Logarithmic map on the tangent space manifold `M`, calculated as the difference of tangent
vectors `q` and `p` from `M`.
"""
log(::TangentSpace, ::Any...)
function log!(::TangentSpace, X, p, q)
    copyto!(X, q - p)
    return X
end

function manifold_dimension(M::TangentSpace)
    return manifold_dimension(M.manifold)
end

LinearAlgebra.norm(M::VectorSpaceFiber, p, X) = norm(M.manifold, M.point, X)

function parallel_transport_to!(M::TangentSpace, Y, p, X, q)
    return copyto!(M.manifold, Y, p, X)
end

@doc raw"""
    project(M::TangentSpace, p)

Project the point `p` from the tangent space `M`, that is project the vector `p`
tangent at `M.point`.
"""
project(::TangentSpace, ::Any)

function project!(M::TangentSpace, q, p)
    return project!(M.manifold, q, M.point, p)
end

@doc raw"""
    project(M::TangentSpace, p, X)

Project the vector `X` from the tangent space `M`, that is project the vector `X`
tangent at `M.point`.
"""
project(::TangentSpace, ::Any, ::Any)

function project!(M::TangentSpace, Y, p, X)
    return project!(M.manifold, Y, M.point, X)
end

function Random.rand!(M::TangentSpace, X; vector_at = nothing)
    rand!(M.manifold, X; vector_at = M.point)
    return X
end
function Random.rand!(rng::AbstractRNG, M::TangentSpace, X; vector_at = nothing)
    rand!(rng, M.manifold, X; vector_at = M.point)
    return X
end

function representation_size(B::TangentSpace)
    return representation_size(B.manifold)
end

function Base.show(io::IO, ::MIME"text/plain", TpM::TangentSpace)
    println(io, "Tangent space to the manifold $(base_manifold(TpM)) at point:")
    pre = " "
    sp = sprint(show, "text/plain", TpM.point; context = io, sizehint = 0)
    sp = replace(sp, '\n' => "\n$(pre)")
    return print(io, pre, sp)
end

function vector_transport_to!(M::TangentSpace, Y, p, X, q, ::AbstractVectorTransportMethod)
    return copyto!(M.manifold, Y, p, X)
end

@doc raw"""
    Y = Weingarten(M::TangentSpace, p, X, V)
    Weingarten!(M::TangentSpace, Y, p, X, V)

Compute the Weingarten map ``\mathcal W_p`` at `p` on the [`TangentSpace`](@ref) `M` with respect to the
tangent vector ``X \in T_p\mathcal M`` and the normal vector ``V \in N_p\mathcal M``.

Since this a flat space by itself, the result is always the zero tangent vector.
"""
Weingarten(::TangentSpace, p, X, V)

Weingarten!(::TangentSpace, Y, p, X, V) = fill!(Y, 0)

@doc raw"""
    zero_vector(M::TangentSpace, p)

Zero tangent vector at point `p` from the tangent space `M`, that is the zero tangent vector
at point `M.point`.
"""
zero_vector(::TangentSpace, ::Any...)

function zero_vector!(M::TangentSpace, X, p)
    return zero_vector!(M.manifold, X, M.point)
end
