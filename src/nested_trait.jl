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
    NestedTrait <: AbstractTrait

Combine two traits into a combined trait.  Note that this introduces a preceedence.
the first of the traits takes preceedence if a trait is implemented for both functions.

# Constructor

    NestedTrait(head::AbstractTrait, tail::AbstractTrait)
"""
struct NestedTrait{T1<:AbstractTrait,T2<:AbstractTrait} <: AbstractTrait
    head::T1
    tail::T2
end

function Base.show(io::IO, t::NestedTrait)
    return print(io, "NestedTrait(", t.head, ", ", t.tail, ")")
end

"""
    active_traits(args...)

Return the list of traits applicable to the given function call. This function should be
overloaded for specific function calls.
"""
@inline active_traits(args...) = EmptyTrait()

"""
    merge_traits(t1, t2, trest...)

Merge two traits into a nested list of traits. Note that this takes trait preceedence into account,
i.e. `t1` takes preceedence over `t2` is any operations.
It always returns either ab [`EmptyTrait`](@ref) or a [`NestedTrait`](@ref).

This means that for
* one argument it just returns the trait itself if it is list-like, or wraps the trait in a
    single-element list otherwise,
* two arguments that are list-like, it merges them,
* two arguments of which only the first one is list-like and the second one is not,
    it appends the second argument to the list,
* two arguments of which only the second one is list-like, it prepends the first one to
    the list,
* two arguments of which none is list-like, it creates a two-element list.
* more than two arguments it recursively performs a left-assiciative recursive reduction
    on arguments, that is for example `merge_traits(t1, t2, t3)` is equivalent to
    `merge_traits(merge_traits(t1, t2), t3)`
"""
merge_traits()

@inline merge_traits() = EmptyTrait()
@inline merge_traits(t::EmptyTrait) = t
@inline merge_traits(t::NestedTrait) = t
@inline merge_traits(t::AbstractTrait) = NestedTrait(t, EmptyTrait())
@inline merge_traits(t1::EmptyTrait, ::EmptyTrait) = t1
@inline merge_traits(::EmptyTrait, t2::AbstractTrait) = merge_traits(t2)
@inline merge_traits(t1::AbstractTrait, t2::EmptyTrait) = NestedTrait(t1, t2)
@inline merge_traits(t1::AbstractTrait, t2::NestedTrait) = NestedTrait(t1, t2)
@inline merge_traits(::EmptyTrait, t2::NestedTrait) = t2
@inline function merge_traits(t1::AbstractTrait, t2::AbstractTrait)
    return NestedTrait(t1, NestedTrait(t2, EmptyTrait()))
end
@inline merge_traits(t1::NestedTrait, ::EmptyTrait) = t1
@inline function merge_traits(t1::NestedTrait, t2::AbstractTrait)
    return NestedTrait(t1.head, merge_traits(t1.tail, t2))
end
@inline function merge_traits(t1::NestedTrait, t2::NestedTrait)
    return NestedTrait(t1.head, merge_traits(t1.tail, t2))
end
@inline function merge_traits(
    t1::AbstractTrait,
    t2::AbstractTrait,
    t3::AbstractTrait,
    trest::AbstractTrait...,
)
    return merge_traits(merge_traits(t1, t2), t3, trest...)
end

"""
    parent_trait(t::AbstractTrait)

Return the parent trait for trait `t`, that is the more general trait whose behaviour it
inherits as a fallback.
"""
@inline parent_trait(::AbstractTrait) = EmptyTrait()

@inline function trait(args...)
    bt = active_traits(args...)
    return expand_trait(bt)
end

"""
    expand_trait(t::AbstractTrait)

Expand given trait into an ordered [`NestedTrait`](@ref) list of traits with their parent
traits obtained using [`parent_trait`](@ref).
"""
expand_trait(::AbstractTrait)

@inline expand_trait(e::EmptyTrait) = e
@inline expand_trait(t::AbstractTrait) = merge_traits(t, expand_trait(parent_trait(t)))
@inline function expand_trait(t::NestedTrait)
    et1 = expand_trait(t.head)
    et2 = expand_trait(t.tail)
    return merge_traits(et1, et2)
end
