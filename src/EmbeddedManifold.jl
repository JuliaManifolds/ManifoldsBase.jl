"""
    AbstractEmbeddingType

A type used to specify properties of an [`AbstractEmbeddedManifold`](@ref).
"""
abstract type AbstractEmbeddingType end

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

This should be used for example for [`check_manifold_point`](@ref) or [`check_tangent_vector`](@ref),
which should first invoke the test of the embedding and then test further constraints
the representation in the embedding has for these points to be valid.

Technically this is realised by making the [`AbstractEmbeddedManifold`](@ref) is a decorator
for the [`Manifold`](@ref)s that are subtypes.
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
    EmbeddedManifold{𝔽, MT <: Manifold, NT <: Manifold} <: AbstractDecoratorManifold{𝔽}

A type to represent an explicit embedding of a [`Manifold`](@ref) `M` of type `MT` embedded
into a manifold `N` of type `NT`.

!!!note
    This type is not required if a manifold `M` is to be embedded in one specific manifold `N`. One can then just implement
    [`embed!`](@ref) and [`project!`](@ref). Only for a second –maybe considered non-default–
    embedding, this type should be considered in order to dispatch on different embed
    and project methods for different embeddings `N`.

# Fields

* `manifold` the manifold that is an embedded manifold
* `embedding` a second manifold, the first one is embedded into

# Constructor

    EmbeddedManifold(M, N)

Generate the `EmbeddedManifold` of the [`Manifold`](@ref) `M` into the
[`Manifold`](@ref) `N`.
"""
struct EmbeddedManifold{𝔽,MT<:Manifold{𝔽},NT<:Manifold} <: AbstractDecoratorManifold{𝔽}
    manifold::MT
    embedding::NT
end

function allocate_result(M::AbstractEmbeddedManifold, f::typeof(embed), x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[1], T, representation_size(decorated_manifold(M)))
end

function allocate_result(M::AbstractEmbeddedManifold, f::typeof(project), x...)
    T = allocate_result_type(M, f, x)
    return allocate(x[1], T, representation_size(base_manifold(M)))
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


"""
    check_manifold_point(M::AbstractEmbeddedManifold, p; kwargs)

check whether a point `p` is a valid point on the [`AbstractEmbeddedManifold`](@ref),
i.e. that `embed(M, p)` is a valid point on the embedded manifold.
"""
function check_manifold_point(M::AbstractEmbeddedManifold, p; kwargs...)
    return invoke(
        check_manifold_point,
        Tuple{typeof(get_embedding(M)),typeof(p)},
        get_embedding(M),
        p;
        kwargs...,
    )
end

"""
    check_tangent_vector(M::AbstractEmbeddedManifold, p, X; check_base_point = true, kwargs...)

check that `embed(M, p, X)` is a valid tangent to `embed(M, p)`, where `check_base_point`
determines whether the validity of `p` is checked, too.
"""
function check_tangent_vector(
    M::AbstractEmbeddedManifold,
    p,
    X;
    check_base_point = true,
    kwargs...,
)
    if check_base_point
        mpe = check_manifold_point(M, p; kwargs...)
        mpe === nothing || return mpe
    end
    return invoke(
        check_tangent_vector,
        Tuple{typeof(get_embedding(M)),typeof(p),typeof(X)},
        get_embedding(M),
        p,
        X;
        check_base_point = check_base_point,
        kwargs...,
    )
end

decorated_manifold(M::EmbeddedManifold) = M.embedding

function embed(M::EmbeddedManifold, p)
    q = allocate_result(M, embed, p)
    embed!(M, q, p)
    return q
end

embed!(::AbstractEmbeddedManifold, q, p) = copyto!(q, p)
embed!(::AbstractEmbeddedManifold, Y, p, X) = copyto!(Y, X)

"""
    get_embedding(M::AbstractEmbeddedManifold)

Return the [`Manifold`](@ref) `N` an [`AbstractEmbeddedManifold`](@ref) is embedded into.
"""
get_embedding(::AbstractEmbeddedManifold)

@decorator_transparent_function function get_embedding(M::AbstractEmbeddedManifold)
    return decorated_manifold(M)
end

"""
    get_embedding(M::EmbeddedManifold)

Return the [`Manifold`](@ref) `N` an [`EmbeddedManifold`](@ref) is embedded into.
"""
get_embedding(::EmbeddedManifold)

function get_embedding(M::EmbeddedManifold)
    return M.embedding
end

function project(M::EmbeddedManifold, p)
    q = allocate_result(M, project, p)
    project!(M, q, p)
    return q
end

function show(io::IO, M::EmbeddedManifold{𝔽,MT,NT}) where {𝔽,MT<:Manifold{𝔽},NT<:Manifold}
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
default_embedding_dispatch(M::AbstractEmbeddedManifold) = Val(false)

function decorator_transparent_dispatch(
    ::typeof(check_manifold_point),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(check_tangent_vector),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:intransparent)
end
function decorator_transparent_dispatch(
    ::typeof(embed),
    ::AbstractEmbeddedManifold,
    args...,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(::typeof(embed), ::EmbeddedManifold, args...)
    return Val(:intransparent)
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
function decorator_transparent_dispatch(::typeof(project), ::EmbeddedManifold, args...)
    return Val(:intransparent)
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
