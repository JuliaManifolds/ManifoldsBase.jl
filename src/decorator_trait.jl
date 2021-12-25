@doc raw"""
    AbstractDecoratorManifold{𝔽} <: AbstractManifold{𝔽}

Declare a manifold to be an abstract decorator.
A manifold which is a subtype of is a __decorated manifold__, i.e. has

* certain additional properties or
* delegates certain properties to other manifolds.

Most prominently, a manifold might be an embedded manifold, i.e. points on a manifold ``\mathcal M``
are represented by (some, maybe not all) points on another manifold ``\mathcal N``.
Depending on the type of embedding, several functions are dedicated to the embedding.
For example if the embedding is isometric, then the [`inner`](@ref) does not have to be
implemented for ``\mathcal M`` but can be automatically implemented by deligation to ``\mathcal N``.

This is modelled by the `AbstractDecoratorManifold` and traits. These are mapped to functions,
which determine the types of transparencies.

A dault function to implement determines the generic manifold that is added (decorates the manifold),
see [`decorated_manifold`}(@ref).
"""
abstract type AbstractDecoratorManifold{𝔽} <: AbstractManifold{𝔽} end

@doc raw"""
    decorated_manifold(M::AbstractDecoratorManifold)

Return the manifold `N` that is “attached” to the manifold `M`, i.e. decorates it.
This is for example called by [`get_embedding`](@ref) for the embedding decorator traits.
"""
decorated_manifold(M::AbstractDecoratorManifold)

#
# Base passons
#
representation_size(M::AbstractDecoratorManifold) = representation_size(base_manifold(M))
manifold_dimension(M::AbstractDecoratorManifold) = manifold_dimension(base_manifold(M))

#
# Traits - each passed to a function that is properly documented
#
@traitdef IsEmbeddedManifold{M}
@traitimpl IsEmbeddedManifold{M} < -is_embedded_manifold(M)

"""
    IsEmbeddedManifold{M}
    is_embedded_manifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded manifold.
To activate this for your manifold, set `isembedded_manifold` for your manifold type to true.
Manifolds that are [`is_isometic_embedded_manifold`](@ref)s set this to true as well.
"""
function is_embedded_manifold(M::Type{<:AbstractDecoratorManifold})
    return is_isometric_embedded_manifold(M)
end

@traitdef IsIsometricEmbeddedManifold{M}
@traitimpl IsIsometricEmbeddedManifold{M} < -is_isometic_embedded_manifold(M)

"""
    IsIsometricEmbeddedManifold{M}
    is_isometic_embedded_manifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an isometrically embedded manifold.
To activate this for your manifold, set `is_isometric_embedded_manifold` for your manifold type to true.

Here, for example [`inner`](@ref) and [`norm`](@ref) are passed to the embedding

This is automatically set to true, when we have an [`is_embedded_submanifold`](@ref).
"""
function is_isometric_embedded_manifold(::Type{<:AbstractDecoratorManifold})
    return is_embedded_submanifold(M)
end

@traitdef IsEmbeddedSubmanifold{M}
@traitimpl IsEmbeddedSubmanifold{M} < -is_embedded_submanifold(M)

"""
    IsEmbeddedSubmanifold{M}
    is_embedded_submanifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
To activate this for your manifold, set `is_embedded_submanifold` for your manifold type to true.

Here, all retraction, inverse retractions and vectors transports, especially
[`exp`](@ref), [`log`](@ref), and [`parallel_transport`](@ref) are passed to the embedding.
"""
is_geodesic_embedded_manifold(::Type{<:AbstractDecoratorManifold}) = false

"""
    get_embedding(M::AbstractDecoratorManifold)

Specify the embedding of a manifold that has abstract decorators.
By default this returns the [`get_decorator`](@ref) manifold.
"""
get_embedding(M::AbstractDecoratorManifold) = get_decorator(M)

#
# Implemented Traits
#

@traitfn function distance(
    M::Mfld,
    p,
    q,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return distance(get_embedding(M), p, q)
end

@traitfn function exp(
    M::Mfld,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return exp(get_embedding(M), p, X)
end

@traitfn function exp!(
    M::Mfld,
    q,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return exp!(get_embedding(M), q, p, X)
end

@traitfn function inner(
    M::Mfld,
    p,
    X,
    Y,
) where {Mfld <: AbstractDecoratorManifold; IsIsometricEmbeddedManifold{Mfld}}
    return inner(get_embedding(M), p, X, Y)
end

@traitfn function inverse_retract(
    M::Mfld,
    p,
    q,
    m = default_inverse_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return inverse_retract(get_embedding(M), p, q, m)
end

@traitfn function inverse_retract!(
    M::Mfld,
    X,
    p,
    q,
    m = default_inverse_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return inverse_retract!(get_embedding(M), X, p, q, m)
end


@traitfn function norm(
    M::Mfld,
    p,
    X,
    Y,
) where {Mfld <: AbstractDecoratorManifold; IsIsometricEmbeddedManifold{Mfld}}
    return inner(get_embedding(M), p, X, Y)
end

@traitfn function log(
    M::Mfld,
    p,
    q,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return log(get_embedding(M), p, q)
end

@traitfn function log!(
    M::Mfld,
    X,
    p,
    q,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return log!(get_embedding(M), X, p, q)
end

@traitfn function retract(
    M::Mfld,
    p,
    X,
    m = default_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return retract(get_embedding(M), p, X, m)
end

@traitfn function retract!(
    M::Mfld,
    q,
    p,
    X,
    m = default_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return retract!(get_embedding(M), q, p, X, m)
end

@traitfn function vector_transport_along(
    M::Mfld,
    p,
    X,
    c,
    m = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along(get_embedding(M), p, X, c, m)
end

@traitfn function vector_transport_along!(
    M::Mfld,
    Y,
    p,
    X,
    c,
    m = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return retract!(get_embedding(M), Y, p, X, c, m)
end

@traitfn function vector_transport_direction(
    M::Mfld,
    p,
    X,
    d,
    m = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along(get_embedding(M), p, X, d, m)
end

@traitfn function vector_transport_direction!(
    M::Mfld,
    Y,
    p,
    X,
    d,
    m = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return retract!(get_embedding(M), Y, p, X, d, m)
end

@traitfn function vector_transport_to(
    M::Mfld,
    p,
    X,
    q,
    m = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along(get_embedding(M), p, X, q, m)
end

@traitfn function vector_transport_to!(
    M::Mfld,
    Y,
    p,
    X,
    q,
    m = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return retract!(get_embedding(M), Y, p, X, q, m)
end
