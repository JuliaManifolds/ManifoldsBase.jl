"""
    VectorSpaceFiber{𝔽,M,TSpaceType} = Fiber{𝔽,TSpaceType,M}
        where {𝔽,M<:AbstractManifold,TSpaceType<:VectorSpaceType}

Alias for a [`Fiber`](@ref) when the fiber is a vector space.
"""
const VectorSpaceFiber{𝔽, M, TSpaceType} =
    Fiber{𝔽, TSpaceType, M} where {𝔽, M <: AbstractManifold, TSpaceType <: VectorSpaceType}

is_flat(::VectorSpaceFiber) = true

LinearAlgebra.norm(M::VectorSpaceFiber, X) = norm(M.manifold, M.point, X)
# disambiguation
LinearAlgebra.norm(M::VectorSpaceFiber, X::Real) = norm(M.manifold, M.point, X)
