@doc raw"""
    AbstractAffineConnection

Abstract type for affine connections on a manifold.
"""
abstract type AbstractAffineConnection end

"""
    LeviCivitaConnection

The [Levi-Civita connection](https://en.wikipedia.org/wiki/Levi-Civita_connection) of a Riemannian manifold.
"""
struct LeviCivitaConnection <: AbstractAffineConnection end

"""
    ConnectionManifold{𝔽,,M<:AbstractManifold{𝔽},G<:AbstractAffineConnection} <: AbstractDecoratorManifold{𝔽}

# Constructor

    ConnectionManifold(M, C)

Decorate the [`AbstractManifold`](@ref) `M` with [`AbstractAffineConnection`](@ref) `C`.
"""
struct ConnectionManifold{𝔽, M <: AbstractManifold{𝔽}, C <: AbstractAffineConnection} <:
    AbstractDecoratorManifold{𝔽}
    manifold::M
    connection::C
end

"""
    connection(M::AbstractManifold)

Get the connection (an object of a subtype of [`AbstractAffineConnection`](@ref))
of [`AbstractManifold`](@ref) `M`.

The global default connection is the [`LeviCivitaConnection`](@ref).
"""
connection(::AbstractManifold) = LeviCivitaConnection()

"""
    connection(M::ConnectionManifold)

Return the connection associated with [`ConnectionManifold`](@ref) `M`.
"""
connection(M::ConnectionManifold) = M.connection

decorated_manifold(M::ConnectionManifold) = M.manifold

default_retraction_method(M::ConnectionManifold) = default_retraction_method(M.manifold)
function default_retraction_method(M::ConnectionManifold, t::Type)
    return default_retraction_method(M.manifold, t)
end
function default_inverse_retraction_method(M::ConnectionManifold)
    return default_inverse_retraction_method(M.manifold)
end
function default_inverse_retraction_method(M::ConnectionManifold, t::Type)
    return default_inverse_retraction_method(M.manifold, t)
end
function default_vector_transport_method(M::ConnectionManifold)
    return default_vector_transport_method(M.manifold)
end
function default_vector_transport_method(M::ConnectionManifold, t::Type)
    return default_vector_transport_method(M.manifold, t)
end

"""
    is_default_connection(M::AbstractManifold, c::AbstractAffineConnection)

returns whether an [`AbstractAffineConnection`](@ref) is the default metric on the manifold `M` or not.

This function falls back to check whether [`connection`](@ref)`(M) == c`.
"""
is_default_connection(M::AbstractManifold, c::AbstractAffineConnection)
function is_default_connection(M::AbstractManifold, c::AbstractAffineConnection)
    return connection(M) == c
end
is_default_connection(M::ConnectionManifold) = true

manifold_dimension(M::ConnectionManifold) = manifold_dimension(M.manifold)

representation_size(M::ConnectionManifold) = representation_size(M.manifold)
