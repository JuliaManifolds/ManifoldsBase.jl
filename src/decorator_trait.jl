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
see [`decorated_manifold`](@ref).
"""
abstract type AbstractDecoratorManifold{𝔽} <: AbstractManifold{𝔽} end

"""
    AbstractTrait

An abstract trait type to build a sequence of traits
"""
abstract type AbstractTrait end

"""
    EmptyTrait <: AbstractTrait

A Trait indicating that no feature is present.
"""
struct EmptyTrait <: AbstractTrait end

"""
    NestedTrait <; AbstractTrait

Combine two traits into a combined trait.  Note that this introduces a preceedence.
the first of the traits takes preceedence if a trait is implemented for both functions.

# Constructor

    NestedTrait(t1::AbstractTrait, t2::AbstractTrait)
"""
struct NestedTrait{T1<:AbstractTrait,T2<:AbstractTrait} <: AbstractTrait
    t1::T1
    t2::T2
end

@inline base_trait(args...) = EmptyTrait()

"""
    merge_traits(t1,t2,...)

Merge two traits into a nested trait. Note that this takes trait preceedence into account,
i.e. t1 takes preceedence over t2 is any operations.

This means that for
* one argument it just returns the trait itself,
* for two arguments it returns either `NestedTrait(t1,t2)` or if both are traits that are not [`EmptyTrait`](@ref)
* if for two arguments one is nested, the nesting is changed to have the form `NEstedTrait(t1,(NestedTrait(t2,...))`

"""
merge_traits()

@inline merge_traits() = EmptyTrait()
@inline merge_traits(t::EmptyTrait) = t
@inline merge_traits(t::NestedTrait) = t
@inline merge_traits(t::AbstractTrait) = t # Maybe NestedTrait(t, EmptyTrait()) ?
@inline merge_traits(t1::EmptyTrait, ::EmptyTrait) = t1
@inline merge_traits(::EmptyTrait, t2::AbstractTrait) = t2
@inline merge_traits(t1::AbstractTrait, ::EmptyTrait) = t1 #was NestedTrait(t1, t2)
@inline merge_traits(t1::AbstractTrait, t2::AbstractTrait) = NestedTrait(t1, t2) # Maybe NestedTrait(t1, NestedTrait(t2, EmptyTrait()))
@inline merge_traits(t1::NestedTrait, ::EmptyTrait) = t1
@inline function merge_traits(t1::NestedTrait, t2::AbstractTrait)
    return NestedTrait(t1.t1, merge_traits(t1.t2, t2))
end
@inline function merge_traits(
    t1::AbstractTrait,
    t2::AbstractTrait,
    t3::AbstractTrait,
    trest::AbstractTrait...,
)
    return merge_traits(merge_traits(t1, t2), t3, trest...)
end

@inline parent_trait(::AbstractTrait) = EmptyTrait()

@inline function trait(args...)
    bt = base_trait(args...)
    return expand_trait(bt)
end
"""
    expand_trait(::AbstractTrait)

Expand a trait or two traits, such that the order is always
`NestedTrait(t1, NestedTrait(t2,...))`.
"""
expand_trait(::AbstractTrait)

@inline expand_trait(e::EmptyTrait) = e
@inline expand_trait(t::AbstractTrait) = _expand_trait(t, parent_trait(t))
@inline _expand_trait(t1::AbstractTrait, ::EmptyTrait) = t1
@inline _expand_trait(t1::NestedTrait, ::EmptyTrait) = t1
@inline function _expand_trait(t1::NestedTrait, t2::AbstractTrait)
    et1 = expand_trait(t1.t1)
    et2 = expand_trait(t1.t2)
    return merge_traits(et1, et2, t2)
end
@inline function expand_trait(t::NestedTrait)
    et1 = expand_trait(t.t1)
    et2 = expand_trait(t.t2)
    return merge_traits(et1, et2)
end


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
@traitdef IsEmbeddedManifold{M}
@traitimpl IsEmbeddedManifold{M} < -is_embedded_manifold(M)

"""
    IsEmbeddedManifold{M}
    is_embedded_manifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded manifold.
To activate this for your manifold, set `isembedded_manifold` for your manifold type to true.
Manifolds that are [`is_isometric_embedded_manifold`](@ref)s set this to true as well.
"""
function is_embedded_manifold(M::Type{<:AbstractDecoratorManifold})
    return is_isometric_embedded_manifold(M)
end

@traitdef IsIsometricEmbeddedManifold{M}
@traitimpl IsIsometricEmbeddedManifold{M} < -is_isometric_embedded_manifold(M)

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

@traitdef IsEmbeddedSubmanifold{M}
@traitimpl IsEmbeddedSubmanifold{M} < -is_embedded_submanifold(M)

"""
    IsEmbeddedSubmanifold{M}
    is_embedded_submanifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
To activate this for your manifold, set `is_embedded_submanifold` for your manifold type to true.

Here, all retraction, inverse retractions and vectors transports, especially
[`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref) are passed to the embedding.
"""
is_embedded_submanifold(::Type{<:AbstractDecoratorManifold}) = false

"""
    get_embedding(M::AbstractDecoratorManifold)

Specify the embedding of a manifold that has abstract decorators.
"""
get_embedding(M::AbstractDecoratorManifold)

"""
    decorated_maniofold(M::AbstractDecoratorManifold)

return the manifold that is actually decorated. For an abstract case this is
still the menifold itself, but for example for the [`EmbeddedManifold`](@ref)
it returns the base manifold without its embedding.
"""
decorated_manifold(M::AbstractDecoratorManifold) = M

#
# Implemented Traits
#
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

@traitfn function allocate_result(
    M::Mfld,
    f::typeof(embed),
    x...,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedManifold{Mfld}}
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(get_embedding(M)))
end

@traitfn function allocate_result(
    M::Mfld,
    f::typeof(project),
    x...,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedManifold{Mfld}}
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(M))
end

@traitfn function check_size(
    M::Mfld,
    p,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedManifold{Mfld}}
    return check_size(get_embedding(M), p)
end

@traitfn function check_size(
    M::Mfld,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedManifold{Mfld}}
    return check_size(get_embedding(M), p, X)
end

@traitfn function distance(
    M::Mfld,
    p,
    q,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return distance(get_embedding(M), p, q)
end

@traitfn function embed!(
    M::Mfld,
    q,
    p,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return copyto!(M, q, p)
end

@traitfn function embed!(
    M::Mfld,
    Y,
    p,
    X,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return copyto!(M, Y, p, X)
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
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return inverse_retract(get_embedding(M), p, q, m)
end

@traitfn function inverse_retract!(
    M::Mfld,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return inverse_retract!(get_embedding(M), X, p, q, m)
end

@traitfn function is_point(
    M::Mfld,
    p,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedManifold{Mfld}}
    return is_point(get_embedding(M), p, throw_error; kwargs...)
end

@traitfn function is_point(
    M::Mfld,
    p,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; !IsEmbeddedManifold{Mfld}}
    return invoke(is_point, Tuple{AbstractManifold,Any,Any}, M, p, throw_error; kwargs...)
end

@traitfn function is_vector(
    M::Mfld,
    p,
    X,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedManifold{Mfld}}
    return is_vector(get_embedding(M), p, X, throw_error; kwargs...)
end

@traitfn function is_vector(
    M::Mfld,
    p,
    X,
    throw_error = false;
    kwargs...,
) where {Mfld <: AbstractDecoratorManifold; !IsEmbeddedManifold{Mfld}}
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
    m::AbstractVectorTransportMethod = default_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return retract(get_embedding(M), p, X, m)
end

@traitfn function retract!(
    M::Mfld,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_retraction_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return retract!(get_embedding(M), q, p, X, m)
end

@traitfn function vector_transport_along(
    M::Mfld,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along(get_embedding(M), p, X, c, m)
end

@traitfn function vector_transport_along!(
    M::Mfld,
    Y,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along!(get_embedding(M), Y, p, X, c, m)
end

@traitfn function vector_transport_along!(
    M::Mfld,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_along!(get_embedding(M), Y, p, X, c, m)
end

@traitfn function vector_transport_direction(
    M::Mfld,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_direction(get_embedding(M), p, X, d, m)
end

@traitfn function vector_transport_direction!(
    M::Mfld,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_direction!(get_embedding(M), Y, p, X, d, m)
end

@traitfn function vector_transport_to(
    M::Mfld,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_to(get_embedding(M), p, X, q, m)
end

@traitfn function vector_transport_to!(
    M::Mfld,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
) where {Mfld <: AbstractDecoratorManifold; IsEmbeddedSubmanifold{Mfld}}
    return vector_transport_to!(get_embedding(M), Y, p, X, q, m)
end
