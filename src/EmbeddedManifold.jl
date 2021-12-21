"""
    EmbeddedManifold{𝔽, MT <: AbstractManifold, NT <: AbstractManifold} <: AbstractDecoratorManifold{𝔽}

A type to represent an explicit embedding of a [`AbstractManifold`](@ref) `M` of type `MT` embedded
into a manifold `N` of type `NT`.

!!! note
    This type is not required if a manifold `M` is to be embedded in one specific manifold `N`. One can then just implement
    [`embed!`](@ref) and [`project!`](@ref). Only for a second –maybe considered non-default–
    embedding, this type should be considered in order to dispatch on different embed
    and project methods for different embeddings `N`.

# Fields

* `manifold` the manifold that is an embedded manifold
* `embedding` a second manifold, the first one is embedded into

# Constructor

    EmbeddedManifold(M, N)

Generate the `EmbeddedManifold` of the [`AbstractManifold`](@ref) `M` into the
[`AbstractManifold`](@ref) `N`.
"""
struct EmbeddedManifold{𝔽,MT<:AbstractManifold{𝔽},NT<:AbstractManifold} <:
       AbstractManifold{𝔽}
    manifold::MT
    embedding::NT
end

function allocate_result(M::EmbeddedManifold, f::typeof(embed), x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[1], T, representation_size(get_embedding(M)))
end

function allocate_result(M::EmbeddedManifold, f::typeof(project), x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[1], T, representation_size(base_manifold(M)))
end

"""
    base_manifold(M::EmbeddedManifold, d::Val{N} = Val(-1))

Return the base manifold of `M` that is enhanced with its embedding. For this specific
type the internally stored enhanced manifold `M.manifold` is returned.
"""
base_manifold(M::EmbeddedManifold, ::Val{N} = Val(-1)) where {N} = M.manifold

function embed(M::EmbeddedManifold, p)
    q = allocate_result(M, embed, p)
    embed!(M, q, p)
    return q
end

"""
    get_embedding(M::AbstractEmbeddedManifold)

Return the embedding [`AbstractManifold`](@ref) `N` of `M`, if it exists.
"""
get_embedding(::AbstractManifold)

"""
    get_embedding(M::EmbeddedManifold)

Return the [`AbstractManifold`](@ref) `N` an [`EmbeddedManifold`](@ref) is embedded into.
"""
function get_embedding(M::EmbeddedManifold)
    return M.embedding
end

function project(M::EmbeddedManifold, p)
    q = allocate_result(M, project, p)
    project!(M, q, p)
    return q
end

function show(
    io::IO,
    M::EmbeddedManifold{𝔽,MT,NT},
) where {𝔽,MT<:AbstractManifold{𝔽},NT<:AbstractManifold}
    return print(io, "EmbeddedManifold($(M.manifold), $(M.embedding))")
end
#=
#
# Abstract parent – i.e. pass to embedding
for f in [
    copy,
    copyto!,
    embed,
    get_basis,
    get_coordinates,
    get_coordinates!,
    get_vector,
    get_vector!,
    inverse_retract!,
    mid_point!,
    project,
    retract!,
    vector_transport_along,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
]
    eval(
        quote
            function $f(M::decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractEmbeddedManifold,
                args...,
            )
                return Val(:parent)
            end
        end,
    )
end
# Abstract generic isometric
for f in [inverse_retract!, retract!]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractEmbeddedManifold{𝔽,<:AbstractIsometricEmbeddingType},
                args...,
            ) where {𝔽}
                return Val(:parent)
            end
        end,
    )
end
for f in [norm, inner]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractEmbeddedManifold{𝔽,<:AbstractIsometricEmbeddingType},
                args...,
            ) where {𝔽}
                return Val(:transparent)
            end
        end,
    )
end
#
# Transparent Isometric Embedding – additionally transparent
for f in [
    distance,
    exp,
    exp!,
    inner,
    inverse_retract,
    inverse_retract!,
    log,
    log!,
    mid_point,
    mid_point!,
    project!,
    project,
    retract,
    retract!,
    vector_transport_along,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
                args...,
            ) where {𝔽}
                return Val(:transparent)
            end
        end,
    )
end
#
# For explicit EmbeddingManifolds the following have to be reimplemented (:intransparent)
for f in [embed, project]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::EmbeddedManifold,
                args...,
            )
                return Val(:intransparent)
            end
        end,
    )
end

# unified vector transports for the three already implemented cases,
# where _direction! still has its nice fallback
for f in [vector_transport_along!, vector_transport_to!]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractEmbeddedManifold{𝔽,<:E},
                Y,
                p,
                X,
                q,
                ::T,
            ) where {𝔽,T,E}
                return Val(:intransparent)
            end
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
                Y,
                p,
                X,
                q,
                ::T,
            ) where {𝔽,T}
                return Val(:transparent)
            end
        end,
    )
    for m in [PoleLadderTransport, SchildsLadderTransport, ScaledVectorTransport]
        eval(
            quote
                function decorator_transparent_dispatch(
                    ::typeof($f),
                    ::AbstractEmbeddedManifold{𝔽,<:E},
                    Y,
                    p,
                    X,
                    q,
                    ::$m,
                ) where {𝔽,E}
                    return Val(:parent)
                end
                function decorator_transparent_dispatch(
                    ::typeof($f),
                    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
                    Y,
                    p,
                    X,
                    q,
                    ::$m,
                ) where {𝔽}
                    return Val(:parent)
                end
            end,
        )
    end
end
=#
