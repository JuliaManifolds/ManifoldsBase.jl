module ManifoldsBaseRecursiveArrayToolsExt

if isdefined(Base, :get_extension)
    using ManifoldsBase
    using RecursiveArrayTools
    using ManifoldsBase: AbstractBasis, VectorBundleBasisData

    import ManifoldsBase: get_vectors, _vector_transport_to
    import Base: getindex, setindex!, view
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..ManifoldsBase
    using ..RecursiveArrayTools
    using ..ManifoldsBase: AbstractBasis, VectorBundleBasisData

    import ..ManifoldsBase: get_vectors, _vector_transport_to
    import ..Base: getindex, setindex!, view
end

function get_vectors(
    M::VectorBundle,
    p::ArrayPartition,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:VectorBundleBasisData},
) where {𝔽}
    xp1 = submanifold_component(p, Val(1))
    zero_m = zero_vector(M.manifold, xp1)
    zero_f = zero_vector(M.fiber, xp1)
    vs = typeof(ArrayPartition(zero_m, zero_f))[]
    for bv in get_vectors(M.manifold, xp1, B.data.base_basis)
        push!(vs, ArrayPartition(bv, zero_f))
    end
    for bv in get_vectors(M.fiber, xp1, B.data.vec_basis)
        push!(vs, ArrayPartition(zero_m, bv))
    end
    return vs
end

"""
    getindex(p::ArrayPartition, M::VectorBundle, s::Symbol)
    p[M::VectorBundle, s]

Access the element(s) at index `s` of a point `p` on a [`VectorBundle`](@ref) `M` by
using the symbols `:point` and `:vector` for the base and vector component, respectively.
"""
@inline function Base.getindex(p::ArrayPartition, M::VectorBundle, s::Symbol)
    (s === :point) && return p.x[1]
    (s === :vector) && return p.x[2]
    return throw(DomainError(s, "unknown component $s on $M."))
end

"""
    setindex!(p::ArrayPartition, val, M::VectorBundle, s::Symbol)
    p[M::VectorBundle, s] = val

Set the element(s) at index `s` of a point `p` on a [`VectorBundle`](@ref) `M` to `val` by
using the symbols `:point` and `:vector` for the base and vector component, respectively.

!!! note

    The *content* of element of `p` is replaced, not the element itself.
"""
@inline function Base.setindex!(x::ArrayPartition, val, M::VectorBundle, s::Symbol)
    if s === :point
        return copyto!(x.x[1], val)
    elseif s === :vector
        return copyto!(x.x[2], val)
    else
        throw(DomainError(s, "unknown component $s on $M."))
    end
end

function _vector_transport_to(
    M::VectorBundle,
    p,
    X,
    q,
    m::VectorBundleProductVectorTransport,
)
    px, pVx = submanifold_components(M.manifold, p)
    VXM, VXF = submanifold_components(M.manifold, X)
    qx, qVx = submanifold_components(M.manifold, q)
    return ArrayPartition(
        vector_transport_to(M.manifold, px, VXM, qx, m.method_point),
        vector_transport_to(M.manifold, px, VXF, qx, m.method_vector),
    )
end

@inline function Base.view(x::ArrayPartition, M::VectorBundle, s::Symbol)
    (s === :point) && return x.x[1]
    (s === :vector) && return x.x[2]
    throw(DomainError(s, "unknown component $s on $M."))
end
end
