
"""
    VectorSpaceFiber{𝔽,M,TFiber}

Alias for [`Fiber`](@ref) when the fiber is a vector space.
"""
const VectorSpaceFiber{𝔽,M,TSpaceType} =
    Fiber{𝔽,TSpaceType,M} where {𝔽,M<:AbstractManifold{𝔽},TSpaceType<:VectorSpaceType}

LinearAlgebra.norm(M::VectorSpaceFiber, p, X) = norm(M.manifold, M.point, X)
