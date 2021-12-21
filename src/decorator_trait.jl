#
# Ansatz: All olt `AbstractXManifolds` get replaced by traits
# The concrete ones will stay

abstract type AbstractDecoratorManifold{𝔽} <: AbstractManifold{𝔽} end

#
# EmbeddedManifold Decorator
#
@traitdef IsEmbeddedManifold{M}
@traitimpl IsEmbeddedManifold{M} < -is_embedded_manifold(M)

function is_embedded_manifold(M::Type{AbstractDecoratorManifold})
    return is_isometric_embedded_manifold(M)
end

@traitdef IsIsometricEmbeddedManifold{M}
@traitimpl IsIsometricEmbeddedManifold{M} < -is_isometic_embedded_manifold(M)

function is_isometric_embedded_manifold(::Type{AbstractDecoratorManifold})
    return is_geodesic_embedded_manifold(M)
end

@traitdef IsGeodesicEmbeddedManifold{M}
@traitimpl IsGeodesicEmbeddedManifold{M} < -is_geodesic_embedded_manifold(M)

is_geodesic_embedded_manifold(::Type{<:AbstractDecoratorManifold}) = false

@traitfn function inner(
    M::Mfld,
    p,
    X,
    Y,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedManifold{Mfld}}
    return inner(get_embedding(M), p, X, Y)
end

@traitfn function exp(
    M::Mfld,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; IsGeodesicEmbeddedManifold{Mfld}}
    return exp(get_embedding(M), p, X)
end
