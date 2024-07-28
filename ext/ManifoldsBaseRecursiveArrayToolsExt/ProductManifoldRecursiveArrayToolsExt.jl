
allocate(a::AbstractArray{<:ArrayPartition}) = map(allocate, a)
allocate(x::ArrayPartition) = ArrayPartition(map(allocate, x.x)...)
function allocate(x::ArrayPartition, T::Type)
    return ArrayPartition(map(t -> allocate(t, T), submanifold_components(x))...)
end

allocate_on(M::ProductManifold) = ArrayPartition(map(N -> allocate_on(N), M.manifolds)...)
function allocate_on(M::ProductManifold, ::Type{ArrayPartition{T,U}}) where {T,U}
    return ArrayPartition(map((N, V) -> allocate_on(N, V), M.manifolds, U.parameters)...)
end
function allocate_on(M::ProductManifold, ft::TangentSpaceType)
    return ArrayPartition(map(N -> allocate_on(N, ft), M.manifolds)...)
end
function allocate_on(
    M::ProductManifold,
    ft::TangentSpaceType,
    ::Type{ArrayPartition{T,U}},
) where {T,U}
    return ArrayPartition(
        map((N, V) -> allocate_on(N, ft, V), M.manifolds, U.parameters)...,
    )
end

@inline function allocate_result(M::ProductManifold, f)
    return ArrayPartition(map(N -> allocate_result(N, f), M.manifolds))
end

function copyto!(M::ProductManifold, q::ArrayPartition, p::ArrayPartition)
    map(copyto!, M.manifolds, submanifold_components(q), submanifold_components(p))
    return q
end
function copyto!(
    M::ProductManifold,
    Y::ArrayPartition,
    p::ArrayPartition,
    X::ArrayPartition,
)
    map(
        copyto!,
        M.manifolds,
        submanifold_components(Y),
        submanifold_components(p),
        submanifold_components(X),
    )
    return Y
end


function default_retraction_method(M::ProductManifold, ::Type{T}) where {T<:ArrayPartition}
    return ProductRetraction(
        map(default_retraction_method, M.manifolds, T.parameters[2].parameters)...,
    )
end

function default_inverse_retraction_method(
    M::ProductManifold,
    ::Type{T},
) where {T<:ArrayPartition}
    return InverseProductRetraction(
        map(default_inverse_retraction_method, M.manifolds, T.parameters[2].parameters)...,
    )
end

function default_vector_transport_method(
    M::ProductManifold,
    ::Type{T},
) where {T<:ArrayPartition}
    return ProductVectorTransport(
        map(default_vector_transport_method, M.manifolds, T.parameters[2].parameters)...,
    )
end

function Base.exp(M::ProductManifold, p::ArrayPartition, X::ArrayPartition)
    return ArrayPartition(
        map(
            exp,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
        )...,
    )
end
function Base.exp(M::ProductManifold, p::ArrayPartition, X::ArrayPartition, t::Number)
    return ArrayPartition(
        map(
            (N, pc, Xc) -> exp(N, pc, Xc, t),
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
        )...,
    )
end

function get_vector(
    M::ProductManifold,
    p::ArrayPartition,
    Xⁱ,
    B::AbstractBasis{𝔽,TangentSpaceType},
) where {𝔽}
    dims = map(manifold_dimension, M.manifolds)
    @assert length(Xⁱ) == sum(dims)
    dim_ranges = _get_dim_ranges(dims)
    tXⁱ = map(dr -> (@inbounds view(Xⁱ, dr)), dim_ranges)
    ts = ziptuples(M.manifolds, submanifold_components(M, p), tXⁱ)
    return ArrayPartition(map((@inline t -> get_vector(t..., B)), ts))
end
function get_vector(
    M::ProductManifold,
    p::ArrayPartition,
    Xⁱ,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:ProductBasisData},
) where {𝔽}
    dims = map(manifold_dimension, M.manifolds)
    @assert length(Xⁱ) == sum(dims)
    dim_ranges = _get_dim_ranges(dims)
    tXⁱ = map(dr -> (@inbounds view(Xⁱ, dr)), dim_ranges)
    ts = ziptuples(M.manifolds, submanifold_components(M, p), tXⁱ, B.data.parts)
    return ArrayPartition(map((@inline t -> get_vector(t...)), ts))
end

function get_vectors(
    M::ProductManifold,
    p::ArrayPartition,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:ProductBasisData},
) where {𝔽}
    N = number_of_components(M)
    xparts = submanifold_components(p)
    BVs = map(t -> get_vectors(t...), ziptuples(M.manifolds, xparts, B.data.parts))
    zero_tvs = map(t -> zero_vector(t...), ziptuples(M.manifolds, xparts))
    vs = typeof(ArrayPartition(zero_tvs...))[]
    for i in 1:N, k in 1:length(BVs[i])
        push!(
            vs,
            ArrayPartition(zero_tvs[1:(i - 1)]..., BVs[i][k], zero_tvs[(i + 1):end]...),
        )
    end
    return vs
end

ManifoldsBase._get_vector_cache_broadcast(::ArrayPartition) = Val(false)

"""
    getindex(p, M::ProductManifold, i::Union{Integer,Colon,AbstractVector})
    p[M::ProductManifold, i]

Access the element(s) at index `i` of a point `p` on a [`ProductManifold`](@ref) `M` by
linear indexing.
See also [Array Indexing](https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing-1) in Julia.
"""
Base.@propagate_inbounds function Base.getindex(
    p::ArrayPartition,
    M::ProductManifold,
    i::Union{Integer,Colon,AbstractVector,Val},
)
    return get_component(M, p, i)
end

function inverse_retract(
    M::ProductManifold,
    p::ArrayPartition,
    q::ArrayPartition,
    method::InverseProductRetraction,
)
    return ArrayPartition(
        map(
            inverse_retract,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, q),
            method.inverse_retractions,
        ),
    )
end

function Base.log(M::ProductManifold, p::ArrayPartition, q::ArrayPartition)
    return ArrayPartition(
        map(
            log,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, q),
        )...,
    )
end

function parallel_transport_direction(
    M::ProductManifold,
    p::ArrayPartition,
    X::ArrayPartition,
    d::ArrayPartition,
)
    return ArrayPartition(
        map(
            parallel_transport_direction,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
            submanifold_components(M, d),
        ),
    )
end

function parallel_transport_to(
    M::ProductManifold,
    p::ArrayPartition,
    X::ArrayPartition,
    q::ArrayPartition,
)
    return ArrayPartition(
        map(
            parallel_transport_to,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
            submanifold_components(M, q),
        ),
    )
end

function project(M::ProductManifold, p::ArrayPartition)
    return ArrayPartition(map(project, M.manifolds, submanifold_components(M, p))...)
end
function project(M::ProductManifold, p::ArrayPartition, X::ArrayPartition)
    return ArrayPartition(
        map(
            project,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
        )...,
    )
end

@doc raw"""
    rand(M::ProductManifold; parts_kwargs = map(_ -> (;), M.manifolds))

Return a random point on [`ProductManifold`](@ref)  `M`. `parts_kwargs` is
a tuple of keyword arguments for `rand` on each manifold in `M.manifolds`.
"""
function Random.rand(
    M::ProductManifold;
    vector_at = nothing,
    parts_kwargs = map(_ -> (;), M.manifolds),
)
    if vector_at === nothing
        return ArrayPartition(
            map((N, kwargs) -> rand(N; kwargs...), M.manifolds, parts_kwargs)...,
        )
    else
        return ArrayPartition(
            map(
                (N, p, kwargs) -> rand(N; vector_at = p, kwargs...),
                M.manifolds,
                submanifold_components(M, vector_at),
                parts_kwargs,
            )...,
        )
    end
end
function Random.rand(
    rng::AbstractRNG,
    M::ProductManifold;
    vector_at = nothing,
    parts_kwargs = map(_ -> (;), M.manifolds),
)
    if vector_at === nothing
        return ArrayPartition(
            map((N, kwargs) -> rand(rng, N; kwargs...), M.manifolds, parts_kwargs)...,
        )
    else
        return ArrayPartition(
            map(
                (N, p, kwargs) -> rand(rng, N; vector_at = p, kwargs...),
                M.manifolds,
                submanifold_components(M, vector_at),
                parts_kwargs,
            )...,
        )
    end
end

function riemann_tensor(
    M::ProductManifold,
    p::ArrayPartition,
    X::ArrayPartition,
    Y::ArrayPartition,
    Z::ArrayPartition,
)
    return ArrayPartition(
        map(
            riemann_tensor,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
            submanifold_components(M, Y),
            submanifold_components(M, Z),
        ),
    )
end

"""
    setindex!(q, p, M::ProductManifold, i::Union{Integer,Colon,AbstractVector})
    q[M::ProductManifold,i...] = p

set the element `[i...]` of a point `q` on a [`ProductManifold`](@ref) by linear indexing to `q`.
See also [Array Indexing](https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing-1) in Julia.
"""
Base.@propagate_inbounds function Base.setindex!(
    q::ArrayPartition,
    p,
    M::ProductManifold,
    i::Union{Integer,Colon,AbstractVector,Val},
)
    return set_component!(M, q, p, i)
end

@inline submanifold_component(p::ArrayPartition, ::Val{I}) where {I} = p.x[I]
@inline submanifold_components(p::ArrayPartition) = p.x

function vector_transport_direction(
    M::ProductManifold,
    p::ArrayPartition,
    X::ArrayPartition,
    d::ArrayPartition,
    m::ProductVectorTransport,
)
    return ArrayPartition(
        map(
            vector_transport_direction,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
            submanifold_components(M, d),
            m.methods,
        ),
    )
end

function vector_transport_to(
    M::ProductManifold,
    p::ArrayPartition,
    X::ArrayPartition,
    q::ArrayPartition,
    m::ProductVectorTransport,
)
    return ArrayPartition(
        map(
            vector_transport_to,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
            submanifold_components(M, q),
            m.methods,
        ),
    )
end
function vector_transport_to(
    M::ProductManifold,
    p::ArrayPartition,
    X::ArrayPartition,
    q::ArrayPartition,
    m::ParallelTransport,
)
    return ArrayPartition(
        map(
            (iM, ip, iX, id) -> vector_transport_to(iM, ip, iX, id, m),
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
            submanifold_components(M, q),
        ),
    )
end
