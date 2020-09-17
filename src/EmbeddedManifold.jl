"""
    AbstractEmbeddingType

A type used to specify properties of an [`AbstractEmbeddedManifold`](@ref).
"""
abstract type AbstractEmbeddingType end

"""
    AbstractEmbeddedManifold{𝔽,T<:AbstractEmbeddingType,𝔽} <: AbstractDecoratorManifold{𝔽}

An abstract type for embedded manifolds, which acts as an [`AbstractDecoratorManifold`](@ref).
The functions of the manifold that is embedded can hence be just passed on to the embedding.
The embedding is further specified by an [`AbstractEmbeddingType`](@ref).

This means, that technically an embedded manifold is a decorator for the embedding, i.e.
functions of this type get, in the semi-transparent way of the
[`AbstractDecoratorManifold`](@ref), passed on to the embedding.

!!! note

    Points on an `AbstractEmbeddedManifold` are represented using representation
    of the embedded manifold. A [`check_manifold_point`](@ref) should first invoke
    the test of the embedding and then test further constraints for the embedded manifold.
"""
abstract type AbstractEmbeddedManifold{𝔽,T<:AbstractEmbeddingType} <:
              AbstractDecoratorManifold{𝔽} end

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

An isometric embedding type that acts as a default, i.e. it has no specifig properties
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
    EmbeddedManifold{𝔽, MT <: Manifold, NT <: Manifold, ET} <: AbstractEmbeddedManifold{𝔽, ET}

A type to represent that a [`Manifold`](@ref) `M` of type `MT` is indeed an emebedded
manifold and embedded into the manifold `N` of type `NT`.
Based on the [`AbstractEmbeddingType`](@ref) `ET`, this introduces methods for `M` by
passing them through to embedding `N`.

# Fields

* `manifold` the manifold that is an embedded manifold
* `embedding` a second manifold, the first one is embedded into

# Constructor

    EmbeddedManifold(M, N, e=TransparentIsometricEmbedding())

Generate the `EmbeddedManifold` of the [`Manifold`](@ref) `M` into the
[`Manifold`](@ref) `N` with [`AbstractEmbeddingType`](@ref) `e` that by default is the most
transparent [`TransparentIsometricEmbedding`](@ref)
"""
struct EmbeddedManifold{𝔽,MT<:Manifold{𝔽},NT<:Manifold,ET} <: AbstractEmbeddedManifold{𝔽,ET}
    manifold::MT
    embedding::NT
end
function EmbeddedManifold(
    M::MT,
    N::NT,
    e::ET = TransparentIsometricEmbedding(),
) where {𝔽,MT<:Manifold{𝔽},NT<:Manifold,ET<:AbstractEmbeddingType}
    return EmbeddedManifold{𝔽,MT,NT,ET}(M, N)
end

function allocate_result(M::AbstractEmbeddedManifold, f::typeof(embed), x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[end], T, representation_size(get_embedding(M)))
end

function allocate_result(M::AbstractEmbeddedManifold, f::typeof(project), x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[end], T, representation_size(base_manifold(M)))
end

"""
    base_manifold(M::AbstractEmbeddedManifold, d::Val{N} = Val(-1))

Return the base manifold of `M` that is enhanced with its embedding.
While functions like `inner` might be overwritten to use the (decorated) manifold
representing the embedding, the base_manifold is the manifold itself in the sense that
detemining e.g. the [`is_default_metric`](@ref) does not fall back to check with
the embedding but with the manifold itself. For this abstract case, just `M` is returned.
"""
base_manifold(M::AbstractEmbeddedManifold, d::Val{N} = Val(-1)) where {N} = M
"""
    base_manifold(M::EmbeddedManifold, d::Val{N} = Val(-1))

Return the base manifold of `M` that is enhanced with its embedding. For this specific
type the internally stored enhanced manifold `M.manifold` is returned.
"""
base_manifold(M::EmbeddedManifold, d::Val{N} = Val(-1)) where {N} = M.manifold

decorated_manifold(M::AbstractEmbeddedManifold) = base_manifold(M)

"""
    get_embedding(M::AbstractEmbeddedManifold)

Return the [`Manifold`](@ref) `N` an [`AbstractEmbeddedManifold`](@ref) is embedded into.
"""
get_embedding(::AbstractEmbeddedManifold)

@decorator_transparent_function function get_embedding(M::AbstractEmbeddedManifold)
    return M.embedding
end

function show(
    io::IO,
    M::EmbeddedManifold{𝔽,MT,NT,ET},
) where {𝔽,MT<:Manifold{𝔽},NT<:Manifold,ET<:AbstractEmbeddingType}
    return print(io, "EmbeddedManifold($(M.manifold), $(M.embedding), $(ET()))")
end

function default_decorator_dispatch(M::EmbeddedManifold)
    return default_embedding_dispatch(M)
end

@doc doc"""
    default_embedding_dispatch(M::AbstractEmbeddedManifold)

This method indicates that an [`AbstractEmbeddedManifold`](@ref) is the default
and hence acts completely transparently and passes all functions transparently onwards.
This is used by the [`AbstractDecoratorManifold`](@ref) within
[`default_decorator_dispatch`](@ref).
By default this is set to `Val(false)`.
"""
default_embedding_dispatch(M::AbstractEmbeddedManifold) = Val(false)

function decorator_transparent_dispatch(
    ::typeof(check_manifold_point),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:transparent)
end

function decorator_transparent_dispatch(
    ::typeof(check_tangent_vector),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:transparent)
end

function decorator_transparent_dispatch(
    ::typeof(embed),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(embed!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(::typeof(exp), ::AbstractEmbeddedManifold, args...)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(exp),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(::typeof(exp!), ::AbstractEmbeddedManifold, args...)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(exp!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(get_basis),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(get_coordinates),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(get_coordinates!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(get_vector),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(get_vector!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(inner),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(inner),
    ::AbstractEmbeddedManifold{𝔽,<:AbstractIsometricEmbeddingType},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(inverse_retract),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(inverse_retract),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(inverse_retract!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(inverse_retract!),
    ::AbstractEmbeddedManifold{𝔽,<:AbstractIsometricEmbeddingType},
    args...,
) where {𝔽}
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(inverse_retract!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end

function decorator_transparent_dispatch(::typeof(log), ::AbstractEmbeddedManifold, args...)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(log),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(::typeof(log!), ::AbstractEmbeddedManifold, args...)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(log!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(mid_point),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(mid_point),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(mid_point!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(mid_point!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(::typeof(norm), ::AbstractEmbeddedManifold, args...)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(norm),
    ::AbstractEmbeddedManifold{𝔽,<:AbstractIsometricEmbeddingType},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(manifold_dimension),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(project!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(project!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(project),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(project),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(retract),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(retract),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(retract!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(retract!),
    ::AbstractEmbeddedManifold{𝔽,<:AbstractIsometricEmbeddingType},
    args...,
) where {𝔽}
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(retract!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_along),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_along),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_along!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_along!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_direction),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_direction),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_direction!),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_direction!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_to),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_to),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    args...,
) where {𝔽}
    return Val(:transparent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_to!),
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
    ::typeof(vector_transport_to!),
    ::AbstractEmbeddedManifold{𝔽,<:E},
    Y,
    p,
    X,
    q,
    ::Union{PoleLadderTransport,SchildsLadderTransport},
) where {𝔽,E}
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_to!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    Y,
    p,
    X,
    q,
    ::Union{PoleLadderTransport,SchildsLadderTransport},
) where {𝔽}
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(vector_transport_to!),
    ::AbstractEmbeddedManifold{𝔽,<:TransparentIsometricEmbedding},
    Y,
    p,
    X,
    q,
    ::T,
) where {𝔽,T}
    return Val(:transparent)
end
