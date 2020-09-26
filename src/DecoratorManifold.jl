#
# Helper
#
@inline _extract_val(::Val{T}) where {T} = T

#! format: off
# turn formatting for for the following functions
# due to the if with returns inside (formatter puts a return upfront the if)
@doc """
    _split_signature(sig::Expr)

this method splits a function signature and returns a named tuple containing
- `fname` the function name
- `where_exprs` the expressions in `where`
- `callargs` the call arguments
- `kwargs_list` the list of keywort arguments
- `argnames` - the names of the arguments
- `argtypes` - the types of the arguments
- `kwargs_call` - ?
"""
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
        if isa(arg, Expr)
            return arg.args[1]
        else
            return arg
        end
    end
    argtypes = map(callargs) do arg
        if isa(arg, Expr)
            return arg.args[end] #if we omit the variable name, its [1] else [2]
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
# Macros
#
"""
    @decorate_case(ex)
    @decorate_case(fallback_case = :intransparent, ex)

Inline definitions are not supported. The function signature however may contain
keyword arguments and a where clause. It does not allow for parameters with default values.

# Examples

```julia
@decorate_case function log!(M::AbstractGroupManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
@decorate_case :undecorate function log!(M::AbstractGroupManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
```
"""
macro decorate_case(ex)
    return esc(quote
        @decorate_case :implement ($ex)
    end)
end
macro decorate_case(case, input_ex)
    ex = macroexpand(__module__, input_ex)
    parts = _split_function(ex)
    callargs = parts[:callargs]
    body = parts[:body]
    argnames = parts[:argnames]
    argtypes = parts[:argtypes]
    fname = parts[:fname]
    where_exprs = parts[:where_exprs]
    kwargs_list = parts[:kwargs_list]
    return esc(
        quote
            function ($(fname))(
                ($(argnames[1]),s)::Tuple{$(argtypes[1]),Val{$(case)}},
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return ($body)
            end
        end,
    )
end

"""
    @decorate_signature(ex)

Introduces a given function to be transparent with respect to all decorators.
The function is adressed by its signature in `ex`.

Supports standard, keyword arguments and `where` clauses. Doesn't support parameters with
default values. Inline definitions are not supported. The function signature however may contain
keyword arguments and a where clause.

It introduces a dispatch on several transparency modes

The cases of transparency are

* `:undecorate` states, that the function is transparently passed on to the manifold that
  is decorated by the [`AbstractDecoratorManifold`](@ref) `M`, which is determined using
  the function [`decorated_manifold`](@ref). This is the default.
* `:implement` states that an implementation for this decorator is required, and if
  none of the types provides one, an error is issued.
* `:inherit` states, that this function passes on to the supertype instead of to the
  decorated manifold.

Any other symbol can be used to access a field from the first (manifold) argument


The dispatch kind can later still be set to something different, see [`decorator_transparent_dispatch`](@ref)

# Examples:

```julia
@decorate_signature log!(M::AbstractDecoratorManifold, X, p, q)
@decorate_signature log!(M::TD, X, p, q) where {TD<:AbstractDecoratorManifold}
@decorate_signature isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
```
"""
macro decorate_signature(input_ex)
    ex = macroexpand(__module__, input_ex)
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
            function ($fname)(
                $(argnames[1])::AbstractDecoratorManifold,
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                transparency = ManifoldsBase._acts_transparently($fname, $(argnames...))
                $(parts.fname)(
                    ($(argnames[1]),Val(transparency)),
                    $(argnames[2:end]...);
                    $(kwargs_call...),
                )
            end
            # (a) :inherit (from parent)
            function ($(fname))(
                ($(argnames[1]),s)::Tuple{$(argtypes[1]),Val{:inherit}},
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
            # (b) :implement (this function is not passed on anywhere and needs to be implemented)
            function ($(fname))(
                ($(argnames[1]),s)::Tuple{$(argtypes[1]),Val{:implement}},
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return error(string(
                    ManifoldsBase.manifold_function_not_implemented_message(
                        $(argnames[1]),
                        $fname,
                        $(argnames[2:end]...),
                    ),
                    " Usually this is implemented for a ",
                    $(argtypes[1]),
                    ". Maybe you missed to implement this function for a default?",
                ))
            end
            # (c) :undecorate act transparently and pass to decorator
            function ($(fname))(
                ($(argnames[1]),s)::Tuple{$(argtypes[1]),Val{:undecorate}},
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...)}
                return ($fname)(
                    ManifoldsBase.decorated_manifold($(argnames[1])),
                    $(argnames[2:end]...);
                    $(kwargs_call...),
                )
            end
            # (e) :s dispatch on the field s of the first arg.
            function ($(fname))(
                ($(argnames[1]),dispatch_symbol)::Tuple{$(argtypes[1]),Val{S}},
                $(callargs[2:end]...);
                $(kwargs_list...),
            ) where {$(where_exprs...), S}
                return ($fname)(
                    getproperty($(argnames[1]),_extract_val(S)),
                    $(argnames[2:end]...);
                    $(kwargs_call...),
                )
            end
        end,
    )
end

"""
    @decorate_function(ex)
    @decorate_function(fallback_case = :implement, ex)

Introduce the function specified by `ex` to be a decorated function.

Inline-definitions are not yet covered – the function signature however may contain
keyword arguments and a where clause.

# Examples

```julia
@decorate_function log!(M::AbstractDecoratorManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
@decorate_function :parent log!(M::AbstractDecoratorManifold, X, p, q)
    log!(decorated_manifold(M), X, p, Q)
end
```
"""
macro decorate_function(ex)
    return esc(quote
        @decorate_function :implement ($ex)
    end)
end
macro decorate_function(case, input_ex)
    ex = macroexpand(__module__, input_ex)
    parts = _split_function(ex)
    callargs = parts[:callargs]
    fname = parts[:fname]
    where_exprs = parts[:where_exprs]
    return esc(
        quote
            @decorate_signature ($(ex.args[1]))
            @decorate_case ($case) ($ex)
            function decorator_dispatch(
                ::typeof($fname),
                $(callargs...),
            ) where {$(where_exprs...)}
                return Val($case)
            end
        end,
    )
end

#! format: on

#
# Type
#
"""
    AbstractDecoratorManifold{𝔽} <: Manifold{𝔽}

An `AbstractDecoratorManifold` indicates that to some extent a manifold subtype
decorates another [`Manifold`](@ref) in the sense that it either

* it extends the functionality of a manifold with further features
* it defines a new manifold that internally uses functions from the decorated manifold

with the main intent that several or most functions of [`Manifold`](@ref) are transparently
passed through to the manifold that is decorated. This way a function implemented for a
decorator acts transparent on all other decorators, i.e. they just pass them through. If
the decorator the function is implemented for is not among the decorators, an error is
issued. By default all base manifold functions, for example [`exp`](@ref) and [`log`](@ref)
are transparent for all decorators.

Transparency of functions with respect to decorators can be specified using the macros
[`@decorator_transparent_fallback`](@ref), [`@decorate_function`](@ref) and
[`@decorate_signature`](@ref).
"""
abstract type AbstractDecoratorManifold{𝔽} <: Manifold{𝔽} end

#
# Functions
#

"""
    is_default_decorator(M) -> Bool

For any manifold that is a subtype of [`AbstractDecoratorManifold`](@ref), this function
indicates whether a certain manifold `M` acts as a default decorator.

This yields that _all_ functions are passed through to the decorated [`Manifold`](@ref)
if `M` is indicated as default. This overwrites all [`is_decorator_transparent`](@ref)
values.

This yields the following advantange: For a manifold one usually implicitly assumes for
example a metric. To avoid reimplementation of this metric when introducing a second metric,
the first metric can be set to be the default, i.e. its implementaion is already given by
the undecorated case.

Value returned by this function is determined by [`default_decorator_dispatch`](@ref),
which returns a `Val`-wrapped boolean for type stability of certain functions.
"""
is_default_decorator(M::Manifold) = _extract_val(default_decorator_dispatch(M))

"""
    default_decorator_dispatch(M) -> Val

Return whether by default to dispatch the the inner manifold of
a decorator (`Val(true)`) or not (`Val(false`). For more details see
[`is_decorator_transparent`](@ref).
"""
default_decorator_dispatch(M::Manifold) = Val(false)

"""
    is_decorator_transparent(f, M::Manifold, args...) -> Bool

Given a [`Manifold`](@ref) `M` and a function `f(M, args...)`, indicate, whether an
[`AbstractDecoratorManifold`](@ref) acts transparently for `f`. This means, it
just passes through down to the internally stored manifold.
Transparency is only defined for decorator manifolds and by default all decorators are transparent.
A function that is affected by the decorator indicates this by returning `false`. To change
this behaviour, see [`decorator_transparent_dispatch`](@ref).

If a decorator manifold is not in general transparent, it might still pass down
for the case that a decorator is the default decorator, see [`is_default_decorator`](@ref).
"""
function is_decorator_transparent(f, M::Manifold, args...)
    return decorator_transparent_dispatch(f, M, args...) === Val(:transparent)
end

"""
    decorator_transparent_dispatch(f, M::Manifold, args...) -> Val

Given a [`Manifold`](@ref) `M` and a function `f(M,args...)`, indicate, whether a
function is `Val(:transparent)` or `Val(:intransparent)` for the (decorated)
[`Manifold`](@ref) `M`. Another possibility is, that for `M` and given `args...`
the function `f` should invoke `M`s `Val(:parent)` implementation, see
[`@decorate_function`](@ref) for details.
"""
decorator_transparent_dispatch(f, M::Manifold, args...) = Val(:manifold)

function _acts_transparently(f, M::Manifold, args...)
    return _val_or(
        default_decorator_dispatch(M),
        decorator_transparent_dispatch(f, M, args...),
    )
end

_val_or(::Val{true}, ::Val{T}) where {T} = Val(:manifold)
_val_or(::Val{false}, val::Val) = val

#
# Functions overwritten with decorators
#

function base_manifold(M::AbstractDecoratorManifold, depth::Val{N} = Val(-1)) where {N}
    N == 0 && return M
    N < 0 && return base_manifold(decorated_manifold(M), depth)
    return base_manifold(decorated_manifold(M), Val(N - 1))
end

@decorate_signature check_manifold_point(
    M::AbstractDecoratorManifold,
    p;
    kwargs...,
)

@decorate_signature check_tangent_vector(
    M::AbstractDecoratorManifold,
    p,
    X;
    kwargs...,
)

"""
    decorated_manifold(M::AbstractDecoratorManifold)

Return the manifold decorated by the decorator `M`. Defaults to `M.manifold`.
"""
decorated_manifold(M::Manifold) = M.manifold

@decorate_signature distance(M::AbstractDecoratorManifold, p, q)

@decorate_signature embed(M::AbstractDecoratorManifold, p, X)
@decorate_signature embed(M::AbstractDecoratorManifold, p)

@decorate_signature embed!(M::AbstractDecoratorManifold, q, p)
@decorate_signature embed!(M::AbstractDecoratorManifold, Y, p, X)

@decorate_signature exp(M::AbstractDecoratorManifold, p, X)

@decorate_signature exp!(M::AbstractDecoratorManifold, q, p, X)

@decorate_signature injectivity_radius(M::AbstractDecoratorManifold)
@decorate_signature injectivity_radius(M::AbstractDecoratorManifold, p)
@decorate_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    m::AbstractRetractionMethod,
)
@decorate_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    m::ExponentialRetraction,
)
@decorate_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    p,
    m::AbstractRetractionMethod,
)
@decorate_signature injectivity_radius(
    M::AbstractDecoratorManifold,
    p,
    m::ExponentialRetraction,
)

@decorate_signature inner(M::AbstractDecoratorManifold, p, X, Y)

@decorate_signature inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)

@decorate_signature inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::LogarithmicInverseRetraction,
)

@decorate_signature inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)

@decorate_signature inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::LogarithmicInverseRetraction,
)

@decorate_signature isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
@decorate_signature isapprox(M::AbstractDecoratorManifold, p, X, Y; kwargs...)

@decorate_signature log(M::AbstractDecoratorManifold, p, q)
@decorate_signature log!(M::AbstractDecoratorManifold, X, p, q)

@decorate_signature manifold_dimension(M::AbstractDecoratorManifold)

@decorate_signature mid_point(M::AbstractDecoratorManifold, p1, p2)
@decorate_signature mid_point!(M::AbstractDecoratorManifold, q, p1, p2)

@decorate_signature number_system(M::AbstractDecoratorManifold)

@decorate_signature project(M::AbstractDecoratorManifold, p)
@decorate_signature project!(M::AbstractDecoratorManifold, q, p)

@decorate_signature project(M::AbstractDecoratorManifold, p, X)
@decorate_signature project!(M::AbstractDecoratorManifold, Y, p, X)

@decorate_signature representation_size(M::AbstractDecoratorManifold)

@decorate_signature retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod,
)

@decorate_signature retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::ExponentialRetraction,
)

@decorate_signature retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod,
)

@decorate_signature retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::ExponentialRetraction,
)

@decorate_signature vector_transport_along(
    M::AbstractDecoratorManifold,
    p,
    X,
    c,
)
@decorate_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c,
)
@decorate_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod,
)
@decorate_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::PoleLadderTransport,
)
@decorate_signature vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::SchildsLadderTransport,
)

@decorate_signature vector_transport_direction(
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
)
@decorate_signature vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
)

@decorate_signature vector_transport_to(
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)
@decorate_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)
@decorate_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::ProjectionTransport,
)
@decorate_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::PoleLadderTransport,
)
@decorate_signature vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::SchildsLadderTransport,
)

@decorate_signature zero_tangent_vector(M::AbstractDecoratorManifold, p)
@decorate_signature zero_tangent_vector!(M::AbstractDecoratorManifold, X, p)
