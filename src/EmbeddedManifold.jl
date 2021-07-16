"""
    AbstractEmbeddingType <: AbstractDecoratorType

A type used to specify properties of an [`AbstractEmbeddedManifold`](@ref).
"""
abstract type AbstractEmbeddingType <: AbstractDecoratorType end

"""
    AbstractEmbeddedManifold{𝔽,T<:AbstractEmbeddingType,𝔽} <: AbstractDecoratorManifold{𝔽}

This abstract type indicates that a concrete subtype is an embedded manifold with the
additional property, that its points are given in the embedding. This also means, that
the default implementation of [`embed`](@ref) is just the identity, since the points are
already stored in the form suitable for this embedding specified. This also holds true for
tangent vectors.

Furthermore, depending on the [`AbstractEmbeddingType`](@ref) different methods are
transparently used from the embedding, for example the [`inner`](@ref) product or even the
[`distance`](@ref) function. Specifying such an embedding type transparently passes the
compuation onwards to the embedding (note again, that no [`embed`](@ref) is required)
and hence avoids to reimplement these methods in the manifold that is embedded.

This should be used for example for [`check_point`](@ref) or [`check_vector`](@ref),
which should first invoke the test of the embedding and then test further constraints
the representation in the embedding has for these points to be valid.

Technically this is realised by making the [`AbstractEmbeddedManifold`](@ref) is a decorator
for the [`AbstractManifold`](@ref)s that are subtypes.
"""
abstract type AbstractEmbeddedManifold{𝔽,T<:AbstractEmbeddingType} <:
              AbstractDecoratorManifold{𝔽,T} end

"""
    DefaultEmbeddingType <: AbstractEmbeddingType

A type of default embedding that does not have any special properties.
"""
struct DefaultEmbeddingType <: AbstractEmbeddingType end

"""
    AbstractIsometricEmbeddingType <: AbstractEmbeddingType

Characterizes an embedding as isometric. For this case the [`inner`](@ref) product
is passed from the embedded manifold to the embedding.
"""
abstract type AbstractIsometricEmbeddingType <: AbstractEmbeddingType end


"""
    DefaultIsometricEmbeddingType <: AbstractIsometricEmbeddingType

An isometric embedding type that acts as a default, i.e. it has no specific properties
beyond its isometric property.
"""
struct DefaultIsometricEmbeddingType <: AbstractIsometricEmbeddingType end

"""
    TransparentIsometricEmbedding <: AbstractIsometricEmbeddingType

Specify that an embedding is the default isometric embedding. This even inherits
logarithmic and exponential map as well as retraction and inverse retractions from the
embedding.

For an example, see [`SymmetricMatrices`](@ref Main.Manifolds.SymmetricMatrices) which are
isometrically embedded in the Euclidean space of matrices but also inherit exponential
and logarithmic maps.
"""
struct TransparentIsometricEmbedding <: AbstractIsometricEmbeddingType end

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
       AbstractDecoratorManifold{𝔽,AbstractDecoratorType}
    manifold::MT
    embedding::NT
end

for M in [AbstractEmbeddedManifold, EmbeddedManifold], fname in [:allocate_result_point, :allocate_result_vector, :allocate_result_coords_like]
    eval(quote
        function $fname(M::$M, f::typeof(embed), x...)
            T = allocate_result_type(M, f, x)
            return allocate(x[1], T, representation_size(decorated_manifold(M)))
        end
        function $fname(M::$M, f::typeof(project), x...)
            T = allocate_result_type(M, f, x)
            return allocate(x[1], T, representation_size(base_manifold(M)))
        end
    end)
end


"""
    base_manifold(M::AbstractEmbeddedManifold, d::Val{N} = Val(-1))

Return the base manifold of `M` that is enhanced with its embedding.
While functions like `inner` might be overwritten to use the (decorated) manifold
representing the embedding, the base_manifold is the manifold itself in the sense that
detemining e.g. the [`is_default_metric`](@ref) does not fall back to check with
the embedding but with the manifold itself. For this abstract case, just `M` is returned.
"""
base_manifold(M::AbstractEmbeddedManifold, ::Val{N} = Val(-1)) where {N} = M
"""
    base_manifold(M::EmbeddedManifold, d::Val{N} = Val(-1))

Return the base manifold of `M` that is enhanced with its embedding. For this specific
type the internally stored enhanced manifold `M.manifold` is returned.
"""
base_manifold(M::EmbeddedManifold, ::Val{N} = Val(-1)) where {N} = M.manifold


"""
    check_point(M::AbstractEmbeddedManifold, p; kwargs)

check whether a point `p` is a valid point on the [`AbstractEmbeddedManifold`](@ref),
i.e. that `embed(M, p)` is a valid point on the embedded manifold.
"""
function check_point(M::AbstractEmbeddedManifold, p; kwargs...)
    return invoke(
        check_point,
        Tuple{typeof(get_embedding(M)),typeof(p)},
        get_embedding(M),
        p;
        kwargs...,
    )
end

"""
    check_vector(M::AbstractEmbeddedManifold, p, X; kwargs...)

Check that `embed(M, p, X)` is a valid tangent to `embed(M, p)`.
"""
function check_vector(M::AbstractEmbeddedManifold, p, X; kwargs...)
    return invoke(
        check_vector,
        Tuple{typeof(get_embedding(M)),typeof(p),typeof(X)},
        get_embedding(M),
        p,
        X;
        kwargs...,
    )
end

decorated_manifold(M::EmbeddedManifold) = M.embedding

function embed(M::EmbeddedManifold, p)
    q = allocate_result_point(get_embedding(M), embed, p)
    embed!(M, q, p)
    return q
end

embed!(::AbstractEmbeddedManifold, q, p) = copyto!(q, p)
embed!(::AbstractEmbeddedManifold, Y, p, X) = copyto!(Y, X)

"""
    get_embedding(M::AbstractEmbeddedManifold)

Return the [`AbstractManifold`](@ref) `N` an [`AbstractEmbeddedManifold`](@ref) is embedded into.
"""
get_embedding(::AbstractEmbeddedManifold)

@decorator_transparent_function function get_embedding(M::AbstractEmbeddedManifold)
    return decorated_manifold(M)
end

"""
    get_embedding(M::EmbeddedManifold)

Return the [`AbstractManifold`](@ref) `N` an [`EmbeddedManifold`](@ref) is embedded into.
"""
get_embedding(::EmbeddedManifold)

function get_embedding(M::EmbeddedManifold)
    return M.embedding
end

function project(M::EmbeddedManifold, p)
    q = allocate_result_point(M.manifold, project, p)
    project!(M, q, p)
    return q
end

function show(
    io::IO,
    M::EmbeddedManifold{𝔽,MT,NT},
) where {𝔽,MT<:AbstractManifold{𝔽},NT<:AbstractManifold}
    return print(io, "EmbeddedManifold($(M.manifold), $(M.embedding))")
end

function default_decorator_dispatch(::EmbeddedManifold)
    return Val(true)
end

@doc raw"""
    default_embedding_dispatch(M::AbstractEmbeddedManifold)

This method indicates that an [`AbstractEmbeddedManifold`](@ref) is the default
and hence acts completely transparently and passes all functions transparently onwards.
This is used by the [`AbstractDecoratorManifold`](@ref) within
[`default_decorator_dispatch`](@ref).
By default this is set to `Val(false)`.
"""
default_embedding_dispatch(::AbstractEmbeddedManifold) = Val(false)

#
# Abstract intransparent – i.e. new implementations necessary
for f in
    [check_point, check_vector, embed!, exp!, inner, log!, manifold_dimension, project!]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractEmbeddedManifold,
                args...,
            )
                return Val(:intransparent)
            end
        end,
    )
end
#
# Abstract parent – i.e. pass to embedding
for f in [
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
            function decorator_transparent_dispatch(
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
