@doc raw"""
    decorated_manifold(M::AbstractDecoratorManifold)


For a manifold `M` that is decorated with properties (for example an embedding `N`)
this function returns the manifold that is attached (as a decorator).
Hence for the embedding example this is `N`.
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

"""
    IsEmbeddedManifold <: AbstractTrait

A trait to declare an [`AbstractManifold`](@ref) as an embedded manifold.
"""
struct IsEmbeddedManifold <: AbstractTrait end

"""
    IsIsometricManifoldEmbeddedManifold <: AbstractTrait

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an isometrically embedded manifold.
To activate this for your manifold, set `is_isometric_embedded_manifold` for your manifold type to true.

Here, for example [`inner`](@ref) and [`norm`](@ref) are passed to the embedding
"""
struct IsIsometricEmbeddedManifold <: AbstractTrait end

parent_trait(::IsIsometricEmbeddedManifold) = IsEmbeddedManifold()

"""
    IsEmbeddedSubmanifold{M}
    is_embedded_submanifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
It is a special case of an [`IsIsometricEmbeddedManifold`](@ref).

Here, additionally, all retraction, inverse retractions and vectors transports, especially
[`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref) are passed to the embedding.
"""
struct IsEmbeddedSubmanifoldManifold <: AbstractTrait end

parent_trait(::IsEmbeddedSubmanifoldManifold) = IsIsometricEmbeddedManifold()


#
# Generic Decorator functions
decorated_manifold(M::AbstractDecoratorManifold) = M

#
# Implemented Traits
function base_manifold(M::AbstractDecoratorManifold, depth::Val{N} = Val(-1)) where {N}
    # end recursion I: depth is 0
    N == 0 && return M
    # end recursion II: M is equal to its decorated manifold (avoid stack overflow)
    D = decorated_manifold(M)
    M === D && return M
    # indefinite many steps for negative values of M
    N < 0 && return base_manifold(D, depth)
    # reduce depth otherwise
    return base_manifold(D, Val(N - 1))
end

#
# Embedded specifix functions.
"""
    get_embedding(M::AbstractDecoratorManifold)

Specify the embedding of a manifold that has abstract decorators.
"""
get_embedding(M::AbstractDecoratorManifold)

#
# -----------------------------------------------------------------------------------------
# This is one new function

# INtroduction and default fallbacks could become a macro?
# Introduce trait
function allocate_result(M::AbstractDecoratorManifold, f, x...)
    return allocate_result(trait(M, f, x...), M, f, x...)
end
# disambiguation
function allocate_result(
    M::AbstractDecoratorManifold,
    f::typeof(get_coordinates),
    p,
    X,
    B::ManifoldsBase.AbstractBasis,
)
    return invoke(
        allocate_result,
        Tuple{AbstractManifold,typeof(get_coordinates),Any,Any,AbstractBasis},
        M,
        f,
        p,
        X,
        B,
    )
end
# Introduce fallback
@inline function allocate_result(::EmptyTrait, M::AbstractManifold, f, x...)
    return invoke(
        allocate_result,
        Tuple{AbstractManifold,typeof(f),typeof(x).parameters...},
        M,
        f,
        x...,
    )
end
# Introduce automatic forward
@inline function allocate_result(t::NestedTrait, M::AbstractManifold, f, x...)
    return allocate_result(t.tail, M, f, x...)
end

function allocate_result(
    ::NestedTrait{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    f::typeof(embed),
    x...,
)
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(get_embedding(M)))
end

# -----------------------------------------------------------------------------------------
# old sruff for hisory / functions still need to be adapted
@traitdef TIsEmbeddedManifold{M}
@traitimpl TIsEmbeddedManifold{M} < -is_embedded_manifold(M)

"""
    TIsEmbeddedManifold{M}
    is_embedded_manifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded manifold.
To activate this for your manifold, set `isembedded_manifold` for your manifold type to true.
Manifolds that are [`is_isometric_embedded_manifold`](@ref)s set this to true as well.
"""
function is_embedded_manifold(M::Type{<:AbstractDecoratorManifold})
    return is_isometric_embedded_manifold(M)
end

@traitdef TIsIsometricEmbeddedManifold{M}
@traitimpl TIsIsometricEmbeddedManifold{M} < -is_isometric_embedded_manifold(M)

"""
    IsIsometricEmbeddedManifold{M}
    is_isometric_embedded_manifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an isometrically embedded manifold.
To activate this for your manifold, set `is_isometric_embedded_manifold` for your manifold type to true.

Here, for example [`inner`](@ref) and [`norm`](@ref) are passed to the embedding

This is automatically set to true, when we have an [`is_embedded_submanifold`](@ref).
"""
function is_isometric_embedded_manifold(Mfld::Type{<:AbstractDecoratorManifold})
    return is_embedded_submanifold(Mfld)
end

@traitdef TIsEmbeddedSubmanifold{M}
@traitimpl TIsEmbeddedSubmanifold{M} < -is_embedded_submanifold(M)

"""
    IsEmbeddedSubmanifold{M}
    is_embedded_submanifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
To activate this for your manifold, set `is_embedded_submanifold` for your manifold type to true.

Here, all retraction, inverse retractions and vectors transports, especially
[`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref) are passed to the embedding.
"""
is_embedded_submanifold(::Type{<:AbstractDecoratorManifold}) = false

@traitfn function allocate_result(
    M::Mfld,
    f::typeof(project),
    x...,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedManifold{Mfld}}
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(M))
end

@traitfn function check_size(
    M::Mfld,
    p,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedManifold{Mfld}}
    return check_size(get_embedding(M), p)
end

@traitfn function check_size(
    M::Mfld,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedManifold{Mfld}}
    return check_size(get_embedding(M), p, X)
end

@traitfn function distance(
    M::Mfld,
    p,
    q,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return distance(get_embedding(M), p, q)
end

@traitfn function embed!(
    M::Mfld,
    q,
    p,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return copyto!(M, q, p)
end

@traitfn function embed!(
    M::Mfld,
    Y,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return copyto!(M, Y, p, X)
end

@traitfn function exp(
    M::Mfld,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return exp(get_embedding(M), p, X)
end

@traitfn function exp!(
    M::Mfld,
    q,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return exp!(get_embedding(M), q, p, X)
end

@traitfn function inner(
    M::Mfld,
    p,
    X,
    Y,
) where {Mfld <: AbstractDecoratorManifold; TIsIsometricEmbeddedManifold{Mfld}}
    return inner(get_embedding(M), p, X, Y)
end

@traitfn function inverse_retract(
    M::Mfld,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return inverse_retract(get_embedding(M), p, q, m)
end

@traitfn function inverse_retract!(
    M::Mfld,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return inverse_retract!(get_embedding(M), X, p, q, m)
end

@traitfn function is_point(
    M::Mfld,
    p,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedManifold{Mfld}}
    return is_point(get_embedding(M), p, throw_error; kwargs...)
end

@traitfn function is_point(
    M::Mfld,
    p,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; !TIsEmbeddedManifold{Mfld}}
    return invoke(is_point, Tuple{AbstractManifold,Any,Any}, M, p, throw_error; kwargs...)
end

@traitfn function is_vector(
    M::Mfld,
    p,
    X,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedManifold{Mfld}}
    return is_vector(get_embedding(M), p, X, throw_error; kwargs...)
end

@traitfn function is_vector(
    M::Mfld,
    p,
    X,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; !TIsEmbeddedManifold{Mfld}}
    return invoke(
        is_vector,
        Tuple{AbstractManifold,Any,Any,Any},
        M,
        p,
        X,
        throw_error;
        kwargs...,
    )
end

@traitfn function norm(
    M::Mfld,
    p,
    X,
    Y,
) where {Mfld <: AbstractDecoratorManifold; TIsIsometricEmbeddedManifold{Mfld}}
    return inner(get_embedding(M), p, X, Y)
end

@traitfn function log(
    M::Mfld,
    p,
    q,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return log(get_embedding(M), p, q)
end

@traitfn function log!(
    M::Mfld,
    X,
    p,
    q,
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return log!(get_embedding(M), X, p, q)
end

@traitfn function retract(
    M::Mfld,
    p,
    X,
    m::AbstractVectorTransportMethod = default_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return retract(get_embedding(M), p, X, m)
end

@traitfn function retract!(
    M::Mfld,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return retract!(get_embedding(M), q, p, X, m)
end

@traitfn function vector_transport_along(
    M::Mfld,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along(get_embedding(M), p, X, c, m)
end

@traitfn function vector_transport_along!(
    M::Mfld,
    Y,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along!(get_embedding(M), Y, p, X, c, m)
end

@traitfn function vector_transport_along!(
    M::Mfld,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along!(get_embedding(M), Y, p, X, c, m)
end

@traitfn function vector_transport_direction(
    M::Mfld,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return vector_transport_direction(get_embedding(M), p, X, d, m)
end

@traitfn function vector_transport_direction!(
    M::Mfld,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return vector_transport_direction!(get_embedding(M), Y, p, X, d, m)
end

@traitfn function vector_transport_to(
    M::Mfld,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return vector_transport_to(get_embedding(M), p, X, q, m)
end

@traitfn function vector_transport_to!(
    M::Mfld,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; TIsEmbeddedSubmanifold{Mfld}}
    return vector_transport_to!(get_embedding(M), Y, p, X, q, m)
end
