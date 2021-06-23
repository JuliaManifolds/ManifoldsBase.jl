#
# Helper
#
@inline _extract_val(::Val{T}) where {T} = T

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
        error("Incorrect syntax in $ex. Expected a :where or :call expression.")
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
        if isa(arg, Expr)
            return arg.args[1]
        else
            return arg
        end
    end
    argtypes = map(callargs) do arg
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
        fname__parent = Symbol(string(fname) * "__parent"),
        fname__transparent = Symbol(string(fname) * "__transparent"),
        fname__intransparent = Symbol(string(fname) * "__intransparent"),
    )
end
#! format: on
function _split_function(ex::Expr)
    if ex.head == :function
        sig = ex.args[1]
        body = ex.args[2]
    else
        error("Incorrect syntax in $ex. Expected :function.")
    end

    return (; body = body, _split_signature(sig)...)
end

#
# Transparency types
#
"""
    AbstractDecoratorType

Decorator types can be used to specify a basic transparency for an [`AbstractDecoratorManifold`](@ref).
This can be seen as an initial (rough) transparency pattern to start a type with.

Note that for a function `f` and it's mutating variant `f!`
* The function `f` is set to `:parent` to first invoke allocation and call of `f!`
* The mutating function `f!` is set to `transparent`
"""
abstract type AbstractDecoratorType end

"""
    DefaultDecoratorType <: AbstractDecoratorType

A default decorator type, where all new functions are transparent by default.
"""
struct DefaultDecoratorType <: AbstractDecoratorType end

#
# Type
#
"""
    AbstractDecoratorManifold{𝔽,T<:AbstractDecoratorType} <: AbstractManifold{𝔽}

An `AbstractDecoratorManifold` indicates that to some extent a manifold subtype
decorates another [`AbstractManifold`](@ref) in the sense that it either

* it extends the functionality of a manifold with further features
* it defines a new manifold that internally uses functions from the decorated manifold

with the main intent that several or most functions of [`AbstractManifold`](@ref) are transparently
passed through to the manifold that is decorated. This way a function implemented for a
decorator acts transparent on all other decorators, i.e. they just pass them through. If
the decorator the function is implemented for is not among the decorators, an error is
issued. By default all base manifold functions, for example [`exp`](@ref) and [`log`](@ref)
are transparent for all decorators.

Transparency of functions with respect to decorators can be specified using the macros
[`@decorator_transparent_fallback`](@ref), [`@decorator_transparent_function`](@ref) and
[`@decorator_transparent_signature`](@ref).

There are currently three modes given a new `AbstractDecoratorManifold` `M`
* `:intransparent` – this function has to be implmented for the new manifold `M`
* `:transparent` – this function is transparent, in the sense that the function is invoked
  on the decorated `M.manifold`. This is the default, when introducing a function or signature.
* `:parent` specifies that (unless implemented) for this function, the classical inheritance
  is issued, i.e. the function is invoked on `M`s supertype.
"""
abstract type AbstractDecoratorManifold{𝔽,T<:AbstractDecoratorType} <: AbstractManifold{𝔽} end

#
# Macros
#
"""
    @decorator_transparent_fallback(ex)
    @decorator_transparent_fallback(fallback_case = :intransparent, ex)

This macro introduces an additional implementation for a certain additional case.
This can especially be used if for an already transparent function and an abstract
intermediate type a change in the default is required.
For implementing a concrete type, neither this nor any other trick is necessary. One
just implements the function as before. Note that a decorator that [`is_default_decorator`](@ref)
still dispatches to the transparent case.


* `:transparent` states, that the function is transparently passed on to the manifold that
  is decorated by the [`AbstractDecoratorManifold`](@ref) `M`, which is determined using
  the function [`decorated_manifold`](@ref).
* `:intransparent` states that an implementation for this decorator is required, and if
  none of the types provides one, an error is issued. Since this macro provides such an
  implementation, this is the default.
* `:parent` states, that this function passes on to the supertype instead of to the
  decorated manifold.

Inline definitions are not supported. The function signature however may contain
keyword arguments and a where clause. It does not allow for parameters with default values.

# Examples

```julia
@decorator_transparent_fallback function log!(M::AbstractGroupManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
@decorator_transparent_fallback :transparent function log!(M::AbstractGroupManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
```
"""
macro decorator_transparent_fallback(ex)
    return esc(quote
        @decorator_transparent_fallback :intransparent ($ex)
    end)
end
macro decorator_transparent_fallback(fallback_case, input_ex)
    ex = macroexpand(__module__, input_ex)
    parts = _split_function(ex)
    callargs = parts[:callargs]
    where_exprs = parts[:where_exprs]
    fname_fallback = Symbol(string(parts.fname) * "__" * string(fallback_case)[2:end])
    return esc(
        quote
            function ($(fname_fallback))(
                $(callargs...);
                $(parts[:kwargs_list]...),
            ) where {$(where_exprs...)}
                return ($(parts[:body]))
            end
        end,
    )
end

"""
    @decorator_transparent_function(ex)
    @decorator_transparent_function(fallback_case = :intransparent, ex)

Introduce the function specified by `ex` to act transparently with respect to
[`AbstractDecoratorManifold`](@ref)s. This introduces the possibility to modify the kind of
transparency the implementation is done for. This optional first argument, the `Symbol`
within `fallback_case`. This macro can be used to define a function and introduce it as
transparent to other decorators. Note that a decorator that [`is_default_decorator`](@ref)
still dispatches to the transparent case.

The cases of transparency are

* `:transparent` states, that the function is transparently passed on to the manifold that
  is decorated by the [`AbstractDecoratorManifold`](@ref) `M`, which is determined using
  the function [`decorated_manifold`](@ref).
* `:intransparent` states that an implementation for this decorator is required, and if
  none of the types provides one, an error is issued. Since this macro provides such an
  implementation, this is the default.
* `:parent` states, that this function passes on to the supertype instead of to the
  decorated manifold. Passing is performed using the `invoke` function where the type of
  manifold is replaced by its supertype.

Innkoline-definitions are not yet covered – the function signature however may contain
keyword arguments and a where clause.

# Examples

```julia
@decorator_transparent_function log!(M::AbstractDecoratorManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
@decorator_transparent_function :parent log!(M::AbstractDecoratorManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
```
"""
macro decorator_transparent_function(ex)
    return esc(quote
        @decorator_transparent_function :intransparent ($ex)
    end)
end
macro decorator_transparent_function(fallback_case, input_ex)
    ex = macroexpand(__module__, input_ex)
    parts = _split_function(ex)
    kwargs_list = parts[:kwargs_list]
    callargs = parts[:callargs]
    fname = parts[:fname]
    where_exprs = parts[:where_exprs]
    body = parts[:body]
    argnames = parts[:argnames]
    argtypes = parts[:argtypes]
    kwargs_call = parts[:kwargs_call]
    fname_fallback = Symbol(string(parts.fname) * "__" * string(fallback_case)[2:end])

    return esc(
        quote
            function ($fname)(
                $(argnames[1])::AbstractDecoratorManifold,
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                transparency = ManifoldsBase._acts_transparently($fname, $(argnames...))
                return if transparency === Val(:parent)
                    return ($(parts.fname__parent))($(argnames...); $(kwargs_call...))
                elseif transparency === Val(:transparent)
                    return ($(parts.fname__transparent))($(argnames...); $(kwargs_call...))
                elseif transparency === Val(:intransparent)
                    return ($(parts.fname__intransparent))(
                        $(argnames...);
                        $(kwargs_call...),
                    )
                else
                    error("incorrect transparency: $transparency")
                end
            end
            function ($(parts[:fname__transparent]))(
                $(argnames[1])::AbstractDecoratorManifold,
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return ($fname)(
                    ManifoldsBase.decorated_manifold($(argnames[1])),
                    $(argnames[2:end]...);
                    $(kwargs_call...),
                )
            end
            function ($(parts[:fname__intransparent]))(
                $(argnames[1])::AbstractDecoratorManifold,
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                error_msg = ManifoldsBase.manifold_function_not_implemented_message(
                    $(argnames[1]),
                    $fname,
                    $(argnames[2:end]...),
                )
                return error(error_msg)
            end
            function ($fname)(
                $(argnames[1])::AbstractManifold,
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return error(
                    string(
                        ManifoldsBase.manifold_function_not_implemented_message(
                            $(argnames[1]),
                            $fname,
                            $(argnames[2:end]...),
                        ),
                        " Usually this is implemented for a ",
                        $(argtypes[1]),
                        ". Maybe you missed to implement this function for a default?",
                    ),
                )
            end
            function ($(parts[:fname__parent]))(
                $(argnames[1])::AbstractDecoratorManifold,
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return invoke(
                    $fname,
                    Tuple{supertype($(argtypes[1])),$(argtypes[2:end]...)},
                    $(argnames...);
                    $(kwargs_call...),
                )
            end
            function ($fname_fallback)(
                $(callargs[1]),
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return ($body)
            end
            function decorator_transparent_dispatch(
                ::typeof($fname),
                $(callargs...),
            ) where {$(where_exprs...)}
                return Val($fallback_case)
            end
        end,
    )
end
#! format: off
# due to the if with returns inside (formatter puts a return upfront the if)
"""
    @decorator_transparent_signature(ex)

Introduces a given function to be transparent with respect to all decorators.
The function is adressed by its signature in `ex`.

Supports standard, keyword arguments and `where` clauses. Doesn't support parameters with
default values. It introduces a dispatch on several transparency modes

The cases of transparency are

* `:transparent` states, that the function is transparently passed on to the manifold that
  is decorated by the [`AbstractDecoratorManifold`](@ref) `M`, which is determined using
  the function [`decorated_manifold`](@ref). This is the default.
* `:intransparent` states that an implementation for this decorator is required, and if
  none of the types provides one, an error is issued.
* `:parent` states, that this function passes on to the supertype instead of to the
  decorated manifold.

Inline definitions are not supported. The function signature however may contain
keyword arguments and a where clause.

The dispatch kind can later still be set to something different, see [`decorator_transparent_dispatch`](@ref)

# Examples:

```julia
@decorator_transparent_signature log!(M::AbstractDecoratorManifold, X, p, q)
@decorator_transparent_signature log!(M::TD, X, p, q) where {TD<:AbstractDecoratorManifold}
@decorator_transparent_signature isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
```
"""
macro decorator_transparent_signature(ex)
    parts = _split_signature(ex)
    kwargs_list = parts[:kwargs_list]
    callargs = parts[:callargs]
    fname = parts[:fname]
    where_exprs = parts[:where_exprs]
    argnames = parts[:argnames]
    argtypes = parts[:argtypes]
    kwargs_call = parts[:kwargs_call]
    #! format: off
    return esc(
        quote
            function ($fname)($(callargs...); $(kwargs_list...)) where {$(where_exprs...)}
                transparency = ManifoldsBase._acts_transparently($fname, $(argnames...))
                if transparency === Val(:parent)
                    return ($(parts.fname__parent))($(argnames...); $(kwargs_call...))
                elseif transparency === Val(:transparent)
                    return ($(parts.fname__transparent))($(argnames...); $(kwargs_call...))
                elseif transparency === Val(:intransparent)
                    return ($(parts.fname__intransparent))(
                        $(argnames...);
                        $(kwargs_call...),
                    )
                else
                    error("incorrect transparency: $transparency")
                end
            end
            function ($(parts[:fname__transparent]))(
                $(callargs...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return ($fname)(
                    ManifoldsBase.decorated_manifold($(argnames[1])),
                    $(argnames[2:end]...);
                    $(kwargs_call...),
                )
            end
            function ($(parts[:fname__intransparent]))(
                $(callargs...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                error_msg = ManifoldsBase.manifold_function_not_implemented_message(
                    $(argnames[1]),
                    $fname,
                    $(argnames[2:end]...),
                )
                return error(error_msg)
            end
            function ($(parts[:fname__parent]))(
                $(callargs...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return invoke(
                    $fname,
                    Tuple{supertype($(argtypes[1])),$(argtypes[2:end]...)},
                    $(argnames...);
                    $(kwargs_call...),
                )
            end
        end,
    )
end
#! format: on


#
# Functions
#

"""
    is_default_decorator(M) -> Bool

For any manifold that is a subtype of [`AbstractDecoratorManifold`](@ref), this function
indicates whether a certain manifold `M` acts as a default decorator.

This yields that _all_ functions are passed through to the decorated [`AbstractManifold`](@ref)
if `M` is indicated as default. This overwrites all [`is_decorator_transparent`](@ref)
values.

This yields the following advantange: For a manifold one usually implicitly assumes for
example a metric. To avoid reimplementation of this metric when introducing a second metric,
the first metric can be set to be the default, i.e. its implementaion is already given by
the undecorated case.

Value returned by this function is determined by [`default_decorator_dispatch`](@ref),
which returns a `Val`-wrapped boolean for type stability of certain functions.
"""
is_default_decorator(M::AbstractManifold) = _extract_val(default_decorator_dispatch(M))

"""
    default_decorator_dispatch(M) -> Val

Return whether by default to dispatch the the inner manifold of
a decorator (`Val(true)`) or not (`Val(false`). For more details see
[`is_decorator_transparent`](@ref).
"""
default_decorator_dispatch(::AbstractManifold) = Val(false)

"""
    is_decorator_transparent(f, M::AbstractManifold, args...) -> Bool

Given a [`AbstractManifold`](@ref) `M` and a function `f(M, args...)`, indicate, whether an
[`AbstractDecoratorManifold`](@ref) acts transparently for `f`. This means, it
just passes through down to the internally stored manifold.
Transparency is only defined for decorator manifolds and by default all decorators are transparent.
A function that is affected by the decorator indicates this by returning `false`. To change
this behaviour, see [`decorator_transparent_dispatch`](@ref).

If a decorator manifold is not in general transparent, it might still pass down
for the case that a decorator is the default decorator, see [`is_default_decorator`](@ref).
"""
function is_decorator_transparent(f, M::AbstractManifold, args...)
    return decorator_transparent_dispatch(f, M, args...) === Val(:transparent)
end

"""
    decorator_transparent_dispatch(f, M::AbstractManifold, args...) -> Val

Given a [`AbstractManifold`](@ref) `M` and a function `f(M,args...)`, indicate, whether a
function is `Val(:transparent)` or `Val(:intransparent)` for the (decorated)
[`AbstractManifold`](@ref) `M`. Another possibility is, that for `M` and given `args...`
the function `f` should invoke `M`s `Val(:parent)` implementation, see
[`@decorator_transparent_function`](@ref) for details.
"""
decorator_transparent_dispatch(::Any, ::AbstractManifold, args...) = Val(:transparent)

function _acts_transparently(f, M::AbstractManifold, args...)
    return _val_or(
        default_decorator_dispatch(M),
        decorator_transparent_dispatch(f, M, args...),
    )
end

_val_or(::Val{true}, ::Val{T}) where {T} = Val(:transparent)
_val_or(::Val{false}, val::Val) = val

#
# Functions overwritten with decorators
#

function base_manifold(M::AbstractDecoratorManifold, depth::Val{N} = Val(-1)) where {N}
    N == 0 && return M
    N < 0 && return base_manifold(decorated_manifold(M), depth)
    return base_manifold(decorated_manifold(M), Val(N - 1))
end

@decorator_transparent_signature check_point(M::AbstractDecoratorManifold, p; kwargs...)

@decorator_transparent_signature check_vector(M::AbstractDecoratorManifold, p, X; kwargs...)

"""
    decorated_manifold(M::AbstractDecoratorManifold)

Return the manifold decorated by the decorator `M`. Defaults to `M.manifold`.
"""
decorated_manifold(M::AbstractManifold) = M.manifold

@decorator_transparent_signature distance(M::AbstractDecoratorManifold, p, q)

@decorator_transparent_signature embed(M::AbstractDecoratorManifold, p, X)
@decorator_transparent_signature embed(M::AbstractDecoratorManifold, p)

@decorator_transparent_signature embed!(M::AbstractDecoratorManifold, q, p)
@decorator_transparent_signature embed!(M::AbstractDecoratorManifold, Y, p, X)

@decorator_transparent_signature exp(M::AbstractDecoratorManifold, p, X)

@decorator_transparent_signature exp!(M::AbstractDecoratorManifold, q, p, X)

@decorator_transparent_signature injectivity_radius(M::AbstractDecoratorManifold)
@decorator_transparent_signature injectivity_radius(M::AbstractDecoratorManifold, p)
@decorator_transparent_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    m::AbstractRetractionMethod,
)
@decorator_transparent_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    m::ExponentialRetraction,
)
@decorator_transparent_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    p,
    m::AbstractRetractionMethod,
)
@decorator_transparent_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    p,
    m::ExponentialRetraction,
)

@decorator_transparent_signature inner(M::AbstractDecoratorManifold, p, X, Y)

@decorator_transparent_signature inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)

@decorator_transparent_signature inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::LogarithmicInverseRetraction,
)

@decorator_transparent_signature inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)

@decorator_transparent_signature inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::LogarithmicInverseRetraction,
)

@decorator_transparent_signature inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::NLsolveInverseRetraction,
)

@decorator_transparent_signature inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::NLsolveInverseRetraction,
)

@decorator_transparent_signature isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
@decorator_transparent_signature isapprox(M::AbstractDecoratorManifold, p, X, Y; kwargs...)

@decorator_transparent_signature log(M::AbstractDecoratorManifold, p, q)
@decorator_transparent_signature log!(M::AbstractDecoratorManifold, X, p, q)

@decorator_transparent_signature manifold_dimension(M::AbstractDecoratorManifold)

@decorator_transparent_signature mid_point(M::AbstractDecoratorManifold, p1, p2)
@decorator_transparent_signature mid_point!(M::AbstractDecoratorManifold, q, p1, p2)

@decorator_transparent_signature number_system(M::AbstractDecoratorManifold)

@decorator_transparent_signature project(M::AbstractDecoratorManifold, p)
@decorator_transparent_signature project!(M::AbstractDecoratorManifold, q, p)

@decorator_transparent_signature project(M::AbstractDecoratorManifold, p, X)
@decorator_transparent_signature project!(M::AbstractDecoratorManifold, Y, p, X)

@decorator_transparent_signature representation_size(M::AbstractDecoratorManifold)

@decorator_transparent_signature retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod,
)

@decorator_transparent_signature retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::ExponentialRetraction,
)

@decorator_transparent_signature retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod,
)

@decorator_transparent_signature retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::ExponentialRetraction,
)

@decorator_transparent_signature vector_transport_along(
    M::AbstractDecoratorManifold,
    p,
    X,
    c,
)
@decorator_transparent_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c,
)
@decorator_transparent_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod,
)
@decorator_transparent_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::PoleLadderTransport,
)
@decorator_transparent_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::SchildsLadderTransport,
)

@decorator_transparent_signature vector_transport_direction(
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
)
@decorator_transparent_signature vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
)

@decorator_transparent_signature vector_transport_to(
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)
@decorator_transparent_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)
@decorator_transparent_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::ProjectionTransport,
)
@decorator_transparent_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::PoleLadderTransport,
)
@decorator_transparent_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::ScaledVectorTransport,
)
@decorator_transparent_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::SchildsLadderTransport,
)

@decorator_transparent_signature zero_vector(M::AbstractDecoratorManifold, p)
@decorator_transparent_signature zero_vector!(M::AbstractDecoratorManifold, X, p)

#
# Manually patch getindex not using the whole machinery
#
Base.@propagate_inbounds function Base.getindex(
    p::AbstractArray,
    M::AbstractDecoratorManifold{𝔽,<:AbstractDecoratorType},
    I::Union{Integer,Colon,AbstractVector}...,
) where {𝔽}
    return getindex(p, decorated_manifold(M), I...)
end

DEFAULT_PARENT_FUNCTIONS = [
    distance,
    exp,
    inner,
    inverse_retract,
    log,
    mid_point,
    norm,
    retract,
    vector_transport_along,
    vector_transport_direction,
    vector_transport_to,
]

for f in DEFAULT_PARENT_FUNCTIONS
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractDecoratorManifold{𝔽,<:AbstractDecoratorType},
                args...,
            ) where {𝔽}
                return Val(:parent)
            end
        end,
    )
end
