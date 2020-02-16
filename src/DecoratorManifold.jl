#
# Type
#
"""
    AbstractDecoratorManifold <: Manifold

An `AbstractDecoratorManifold` indicates that to some extend a manifold subtype
decorates another manifold in the sense that

* it extends the functionality of a manifold with further features
* it defines a new manifold that internally uses functions from another manifold

with the main intend that several or most functions of [`Manifold`](@ref) are transparently
passed thrugh to the manifold that is decorated. This way a function implemented for a
decorator acts transparent on all other decorators, i.e. they just pass them through. If
the decorator the function is implemented for is not among the decorators, an error is
issued. By default all base manifold functions, for example [`exp`](@ref) and [`log`](@ref)
are transparent for all decorators.
"""
abstract type AbstractDecoratorManifold <: Manifold end

#
# Macros
#
"""
    @decorator_transparent_fallback(ex)
    @decorator_transparent_fallback(fallback_case = :intransparent, ex)

This macro introdcues an additional implementation for a certain additional case.
This can especially be used if for an already transparent function and an abstract
intermediate type a change in the default is required.
For implementing a concrete type, neither this nor any other trick is necessary. One
just implements the function as before. Not that a decorator that [`is_default_decorator`](@ref)
still dispatches to the transparent case.


* `:transparent` states, that the function is transparently passed on to the manifold that
is decorated by the [`AbstractDecoratorManifold`](@ref) `M`, which by default is assumed to
be stored in `M.manifold`.
* `: intransparent` states that an implementation for this decorator is required, and if
none of the types provides one, an error is issued. Since this macro provides such an
implementation, this is the default.
* `:parent` states, that this function passes on to the supertype instead of to the
decorated manifold.

Currently inline-definitions are not yet covered – the function signature however may contain
keyword arguments and a where clause. It does not allow for parameters with default values.

# Examples

```julia
@decorator_transparent_fallback function log!(M::AbstractGroupManifold, X, p, q)
    log!(M.manifold, X, p, Q)
end
@decorator_transparent_fallback :transparent function log!(M::AbstractGroupManifold, X, p, q)
    log!(M.manifold, X, p, Q)
end
```
"""
macro decorator_transparent_fallback(ex)
    return esc(quote @decorator_transparent_fallback :intransparent ($ex) end)
end
macro decorator_transparent_fallback(fallback_case, ex)
    if ex.head == :function || ex.head == :(=) #complete or inline function
        sig = ex.args[1]
        body = ex.args[2]
    else
        error("Incorrect syntax in $ex. Expected :function of :(=).")
    end
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
    return esc(quote
        function ($fname)($(callargs[1]), ::Val{$fallback_case}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            ($body)
        end
    end)
end
"""
    @decorator_transparent_function(ex, fallback_case = :intransparent)


Introduce the function specified by `to act transparent with respect to
[`AbstractDecoratorManifold`](@ref)s. This intoduces the possibility to modify the kind of
transparency the implementation is done for. This optional first argument, the `Symbol`
within `fallback_case`. This macro can be used to define a function and introduce it as
transparent to other deocorators. Not that a decorator that [`is_default_decorator`](@ref)
still dispatches to the transparent case.

The cases of transparency are

* `:transparent` states, that the function is transparently passed on to the manifold that
is decorated by the [`AbstractDecoratorManifold`](@ref) `M`, which by default is assumed to
be stored in `M.manifold`.
* `: intransparent` states that an implementation for this decorator is required, and if
none of the types provides one, an error is issued. Since this macro provides such an
implementation, this is the default.
* `:parent` states, that this function passes on to the supertype instead of to the
decorated manifold.

currently inline-definitions are not yet covered – the function signature however may contain
keyword arguments and a where clause.

# Examples

```julia
@decorator_transparent_function log!(M::AbstractDecoratorManifold, X, p, q)
    log!(M.manifold, X, p, Q)
end
@decorator_transparent_function :parent log!(M::TD, X, p, q) where {TD<:AbstractDecoratorManifold}
    log!(M.manifold, X, p, Q)
end
```
"""
macro decorator_transparent_function(ex)
    return esc(quote @decorator_transparent_function :intransparent ($ex) end)
end
macro decorator_transparent_function(fallback_case, ex)
    if ex.head == :function
        sig = ex.args[1]
        body = ex.args[2]
    else
        error("Incorrect syntax in $ex. Expected :function. It does not yet work for inline functions.")
    end
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
    return esc(quote
        function ($fname)($(argnames[1])::AbstractDecoratorManifold, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            return ($fname)($(argnames[1]), _acts_transparently($fname, $(argnames...)), $(argnames[2:end]...),; $(kwargs_list...))
        end
        function ($fname)($(argnames[1])::AbstractDecoratorManifold, ::Val{:transparent}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            return ($fname)($(argnames[1]).manifold, $(argnames[2:end]...); $(kwargs_list...))
        end
        function ($fname)($(argnames[1])::AbstractDecoratorManifold, ::Val{:intransparent}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            error(manifold_function_not_implemented_message($(argnames[1]), $fname, $(argnames[2:end]...)))
        end
        function ($fname)($(argnames[1])::Manifold, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            error(string(
                manifold_function_not_implemented_message($(argnames[1]), $fname, $(argnames[2:end]...)),
                "Usually this is implemented for a ",
                $(argtypes[1]),
                ". Maybe you missed to implement this function for a default?"
            ))
        end
        function ($fname)($(argnames[1])::AbstractDecoratorManifold, ::Val{:parent}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            return invoke($fname, Tuple{supertype($(argtypes[1])), $(argtypes[2:end]...)}, $(argnames...); $(kwargs_list...))
        end
        function ($fname)($(callargs[1]), ::Val{$fallback_case}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            ($body)
        end
        decorator_transparent_dispatch(::typeof($fname), $(callargs...)) where {$(where_exprs...)} = Val($fallback_case)
    end)
end
"""
    @decorator_transparent_signature(ex)

Introduces a given function to be transparent with respect to all decorators.
The function is adressed by its signature in `ex`.

Supports standard, keyword arguments and `where` clauses. Doesn't support parameters with
default values. It introduces a dispatch on several transparency modes

The cases of transparency are

* `:transparent` states, that the function is transparently passed on to the manifold that
is decorated by the [`AbstractDecoratorManifold`](@ref) `M`, which by default is assumed to
be stored in `M.manifold`. This is the default.
* `: intransparent` states that an implementation for this decorator is required, and if
none of the types provides one, an error is issued.
* `:parent` states, that this function passes on to the supertype instead of to the
decorated manifold.

currently inline-definitions are not yet covered – the function signature however may contain
keyword arguments and a where clause.

The dispatch kind can later still be set to something diffrent, see [`decorator_transparent_dispatch`](@ref)

# Examples:

```julia
@decorator_transparent_signature log!(M::AbstractDecoratorManifold, X, p, q)
@decorator_transparent_signature log!(M::TD, X, p, q) where {TD<:AbstractDecoratorManifold}
@decorator_transparent_signature isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
```
"""
macro decorator_transparent_signature(ex)
    if ex.head == :where
        where_exprs = ex.args[2:end]
        call_expr = ex.args[1]
    elseif ex.head == :call
        where_exprs = []
        call_expr = ex
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
    return esc(quote
        function ($fname)($(callargs...); $(kwargs_list...)) where {$(where_exprs...)}
            return ($fname)($(argnames[1]), _acts_transparently($fname, $(argnames...)), $(argnames[2:end]...),; $(kwargs_list...))
        end
        function ($fname)($(callargs[1]), ::Val{:transparent}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            return ($fname)($(argnames[1]).manifold, $(argnames[2:end]...); $(kwargs_list...))
        end
        function ($fname)($(callargs[1]), ::Val{:intransparent}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            error(manifold_function_not_implemented_message($(argnames[1]), $fname, $(argnames[2:end]...)))
        end
        function ($fname)($(callargs[1]), ::Val{:parent}, $(callargs[2:end]...); $(kwargs_list...)) where {$(where_exprs...)}
            return invoke($fname, Tuple{supertype($(argtypes[1])), $(argtypes[2:end]...)}, $(argnames...); $(kwargs_list...))
        end
    end)
end

#
# Functions
#

"""
    is_default_decorator(M)

For any manifold that is a subtype of [`AbstractDecoratorManifold`](@ref), this function
indicates whether a certain manifold `M` acts as a default decorator.

This yields that _all_ functions are passed through to the decorated [`Manifold`](@ref)
if `M` is indicated as default. This overwrites all [`is_decorator_transparent`](@ref)
values.

This yields the following advantange: For a manifold one usually implicitly assumes for
example a metric. To avoid reimplementation of this metric when introducing a second metric,
the first metric can be set to be the default, i.e. its implementaion is already given by
the undecorated case.

to change this value, see [`default_decorator_dispatch`](@ref).
"""
is_default_decorator(M::Manifold) = _extract_val(default_decorator_dispatch(M))

"""
    default_decorator_dispatch(M)

A function to decide whether by default to dispatch th the inner manifold of
a decorator (`Val(true)`) or not (`Val(false`). For more details see
[`is_decorator_transparent`](@ref).
"""
default_decorator_dispatch(M::Manifold) = Val(false)

"""
    is_decorator_transparent(f, M, args...)

Given a [`Manifold`](@ref) `M` and a function `f(M,arge...)`, indicate, whether a
[`AbstractDecoratorManifold`](@ref) acts transparent for `f`. This means, it
just passes through down to the internally stored manifold.
Only decorator manifolds can be transparent and their default is, to be transparent.
A function that is affected by the decorator indicates this by returning `false`. To change
this behaviour, see [`decorator_transparent_dispatch`](@ref).

If a decorator manifold is not in general transparent, it might still pass down
for the case that a decorator is the default decorator, see [`is_default_decorator`](@ref).
"""
function is_decorator_transparent(f, M::Manifold, args...)
    return decorator_transparent_dispatch(f, M, args...) == Val(:transparent)
end

"""
    decorator_transparent_dispatch(f, M, arge...)

Given a [`Manifold`](@ref) `M` and a function `f(M,arge...)`, indicate, whether a
function is `Val(:transparent)` or `Val(:intransparent)` for the (decorated)
[`Manifold`](@ref) `M`. Another possibility is, that for `M` and given `args...`
the function `f` should invoke `M`s `Val(:parent)` implementation.
"""
decorator_transparent_dispatch(f, M::Manifold, args...) = Val(:transparent)

function _acts_transparently(f, M::Manifold, args...)
    return _val_or(default_decorator_dispatch(M), decorator_transparent_dispatch(f, M, args...))
end

_val_or(::Val{true}, ::Val{T}) where {T} = Val(:transparent)
_val_or(::Val{false}, ::Val{T}) where {T} = Val(T)

#
# Functions overwritten with decorators
#

function base_manifold(M::AbstractDecoratorManifold, depth::Val{N} = Val(-1)) where {N}
    return (N != 0) ? base_manifold(M.manifold, (N > 0) ? Val(N-1) : depth) : M
end

@decorator_transparent_signature check_manifold_point(
    M::AbstractDecoratorManifold,
    p;
    kwargs...,
)

@decorator_transparent_signature check_tangent_vector(
    M::AbstractDecoratorManifold,
    p,
    X;
    kwargs...,
)

@decorator_transparent_signature distance(
    M::AbstractDecoratorManifold,
    p,
    q,
)

@decorator_transparent_signature exp!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
)

@decorator_transparent_signature exp!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    T
)

@decorator_transparent_signature injectivity_radius(M::AbstractDecoratorManifold)
@decorator_transparent_signature injectivity_radius(M::AbstractDecoratorManifold, p)

@decorator_transparent_signature inner(M::AbstractDecoratorManifold, p, X, Y)

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

@decorator_transparent_signature isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
@decorator_transparent_signature isapprox(M::AbstractDecoratorManifold, p, X, Y; kwargs...)

@decorator_transparent_signature log!(M::AbstractDecoratorManifold, X, p, q)

@decorator_transparent_signature manifold_dimension(M::AbstractDecoratorManifold)

@decorator_transparent_signature project_point!(M::AbstractDecoratorManifold, q, p)

@decorator_transparent_signature project_tangent!(M::AbstractDecoratorManifold, Y, p, X)

@decorator_transparent_signature projected_distribution(M::AbstractDecoratorManifold, d, p)

@decorator_transparent_signature projected_distribution(M::AbstractDecoratorManifold, d)

@decorator_transparent_signature representation_size(M::AbstractDecoratorManifold)

@decorator_transparent_signature retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod
)

@decorator_transparent_signature retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::ExponentialRetraction
)

@decorator_transparent_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c,
)

@decorator_transparent_signature vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
)

@decorator_transparent_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)

@decorator_transparent_signature zero_tangent_vector!(M::AbstractDecoratorManifold, X, p)
