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
    IsExplicitDecorator <: AbstractTrait

Specify that a certain type should dispatch per default to its [`decorated_manifold`](@ref).

!!! note
    Any decorator _behind_ this decorator might not have any effect, since the function
    dispatch is moved to its field at this point. Therefore this decorator should always be
    _last_ in the [`TraitList`](@ref).
"""
struct IsExplicitDecorator <: AbstractTrait end

"""
    TraitList <: AbstractTrait

Combine two traits into a combined trait.  Note that this introduces a preceedence.
the first of the traits takes preceedence if a trait is implemented for both functions.

# Constructor

    TraitList(head::AbstractTrait, tail::AbstractTrait)
"""
struct TraitList{T1<:AbstractTrait,T2<:AbstractTrait} <: AbstractTrait
    head::T1
    tail::T2
end

function Base.show(io::IO, t::TraitList)
    return print(io, "TraitList(", t.head, ", ", t.tail, ")")
end

"""
    active_traits(f, args...)

Return the list of traits applicable to the given call of function `f``. This function should be
overloaded for specific function calls.
"""
@inline active_traits(f, args...) = EmptyTrait()

"""
    merge_traits(t1, t2, trest...)

Merge two traits into a nested list of traits. Note that this takes trait preceedence into account,
i.e. `t1` takes preceedence over `t2` is any operations.
It always returns either ab [`EmptyTrait`](@ref) or a [`TraitList`](@ref).

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
@inline merge_traits(t::TraitList) = t
@inline merge_traits(t::AbstractTrait) = TraitList(t, EmptyTrait())
@inline merge_traits(t1::EmptyTrait, ::EmptyTrait) = t1
@inline merge_traits(::EmptyTrait, t2::AbstractTrait) = merge_traits(t2)
@inline merge_traits(t1::AbstractTrait, t2::EmptyTrait) = TraitList(t1, t2)
@inline merge_traits(t1::AbstractTrait, t2::TraitList) = TraitList(t1, t2)
@inline merge_traits(::EmptyTrait, t2::TraitList) = t2
@inline function merge_traits(t1::AbstractTrait, t2::AbstractTrait)
    return TraitList(t1, TraitList(t2, EmptyTrait()))
end
@inline merge_traits(t1::TraitList, ::EmptyTrait) = t1
@inline function merge_traits(t1::TraitList, t2::AbstractTrait)
    return TraitList(t1.head, merge_traits(t1.tail, t2))
end
@inline function merge_traits(t1::TraitList, t2::TraitList)
    return TraitList(t1.head, merge_traits(t1.tail, t2))
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

@inline function trait(f::TF, args...) where {TF}
    bt = active_traits(f, args...)
    return expand_trait(bt)
end

"""
    expand_trait(t::AbstractTrait)

Expand given trait into an ordered [`TraitList`](@ref) list of traits with their parent
traits obtained using [`parent_trait`](@ref).
"""
expand_trait(::AbstractTrait)

@inline expand_trait(e::EmptyTrait) = e
@inline expand_trait(t::AbstractTrait) = merge_traits(t, expand_trait(parent_trait(t)))
@inline function expand_trait(t::TraitList)
    et1 = expand_trait(t.head)
    et2 = expand_trait(t.tail)
    return merge_traits(et1, et2)
end

"""
    next_trait(t::AbstractTrait)

Return the next trait to consider, which by default is no following trait (i.e. [`EmptyTrait`](@ref)).

Expecially for a a [`TraitList`](@ref) this function returns the (remaining) tail of
the remaining traits.
"""
next_trait(::AbstractTrait) = EmptyTrait()

@inline next_trait(t::TraitList) = t.tail

#! format: off
# turn formatting for for the following functions
# due to the if with returns inside (formatter puts a return upfront the if)
function _split_signature(sig::Expr)
    if sig.head == :where
        where_exprs = sig.args[2:end]
        call_expr = sig.args[1]
    elseif sig.head == :call
        where_exprs = []
        call_expr = sig
    else
        error("Incorrect syntax in $sig. Expected a :where or :call expression.")
    end
    fname = call_expr.args[1]
    if isa(call_expr.args[2], Expr) && call_expr.args[2].head == :parameters
        # we have keyword arguments
        callargs = call_expr.args[3:end]
        kwargs_list = call_expr.args[2].args
    else
        callargs = call_expr.args[2:end]
        kwargs_list = []
    end
    argnames = map(callargs) do arg
        if isa(arg, Expr) && arg.head === :kw # default val present
            arg = arg.args[1]
        end
        if isa(arg, Expr) && arg.head === :(::) # typed
            return arg.args[1]
        end
        return arg
    end
    argtypes = map(callargs) do arg
        if isa(arg, Expr) && arg.head === :kw # default val present
            arg = arg.args[1]
        end
        if isa(arg, Expr)
            return arg.args[2]
        else
            return Any
        end
    end

    kwargs_call = map(kwargs_list) do kwarg
        if kwarg.head === :...
            return kwarg
        else
            if isa(kwarg.args[1], Symbol)
                kwargname = kwarg.args[1]
            else
                kwargname = kwarg.args[1].args[1]
            end
            return :($kwargname = $kwargname)
        end
    end

    return (;
        fname = fname,
        where_exprs = where_exprs,
        callargs = callargs,
        kwargs_list = kwargs_list,
        argnames = argnames,
        argtypes = argtypes,
        kwargs_call = kwargs_call,
    )
end
#! format: on

macro invoke_maker(argnum, type, sig)
    parts = ManifoldsBase._split_signature(sig)
    kwargs_list = parts[:kwargs_list]
    callargs = parts[:callargs]
    fname = parts[:fname]
    where_exprs = parts[:where_exprs]
    argnames = parts[:argnames]
    argtypes = parts[:argtypes]
    kwargs_call = parts[:kwargs_call]

    return esc(
        quote
            function ($fname)($(callargs...); $(kwargs_list...)) where {$(where_exprs...)}
                return invoke(
                    $fname,
                    Tuple{
                        $(argtypes[1:(argnum - 1)]...),
                        $type,
                        $(argtypes[(argnum + 1):end]...),
                    },
                    $(argnames...);
                    $(kwargs_call...),
                )
            end
        end,
    )
end


macro trait_function(sig, opts = :())
    parts = ManifoldsBase._split_signature(sig)
    kwargs_list = parts[:kwargs_list]
    callargs = parts[:callargs]
    fname = parts[:fname]
    where_exprs = parts[:where_exprs]
    argnames = parts[:argnames]
    argtypes = parts[:argtypes]
    kwargs_call = parts[:kwargs_call]

    argnametype_exprs = [:(typeof($(argname))) for argname in argnames]

    block = quote
        @inline function ($fname)($(callargs...); $(kwargs_list...)) where {$(where_exprs...)}
            return ($fname)(trait($fname, $(argnames...)), $(argnames...); $(kwargs_call...))
        end
        @inline function ($fname)(
            t::TraitList,
            $(callargs...);
            $(kwargs_list...),
        ) where {$(where_exprs...)}
            return ($fname)(next_trait(t), $(argnames...); $(kwargs_call...))
        end
        @inline function ($fname)(
            t::TraitList{IsExplicitDecorator},
            $(callargs...);
            $(kwargs_list...),
        ) where {$(where_exprs...)}
            arg1 = decorated_manifold($(argnames[1]))
            argt1 = typeof(arg1)
            return invoke(
                $fname,
                Tuple{argt1,$(argnametype_exprs[2:end]...)},
                arg1,
                $(argnames[2:end]...);
                $(kwargs_call...),
            )
        end
    end
    if !(:no_empty in opts.args)
        block = quote
            $block
            # See https://discourse.julialang.org/t/extremely-slow-invoke-when-inlined/90665
            # for the reasoning behind @noinline
            @noinline function ($fname)(
                ::EmptyTrait,
                $(callargs...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return invoke(
                    $fname,
                    Tuple{supertype($(argtypes[1])),$(argnametype_exprs[2:end]...)},
                    $(argnames...);
                    $(kwargs_call...),
                )
            end
        end
    end
    return esc(block)
end
