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
                    Tuple{$(argtypes[1:(argnum - 1)]...), $type, $(argtypes[(argnum + 1):end]...)},
                    $(argnames...);
                    $(kwargs_call...),
                )
            end
        end,
    )
end

#TODO: Document thoroughly
macro trait_function(
        sig,
        include_forwards = :(
            (
                EmbeddedForwardingType{DirectEmbedding},
                SimpleForwardingType,
                StopForwardingType,
            )
        ),
        manifold_arg_no = 1,
    )
    parts = ManifoldsBase._split_signature(sig)
    kwargs_list = parts[:kwargs_list]
    callargs = parts[:callargs]
    fname = parts[:fname]
    fname_fwd = Symbol(:_, fname, :_forwarding)
    where_exprs = parts[:where_exprs]
    argnames = parts[:argnames]
    kwargs_call = parts[:kwargs_call]

    argnametype_exprs = [:(typeof($(argname))) for argname in argnames]
    gft = if :p in callargs
        :(ManifoldsBase.get_forwarding_type(M, $fname, p))
    else
        :(ManifoldsBase.get_forwarding_type(M, $fname))
    end

    ge = if :p in callargs
        :(get_embedding(M, p))
    else
        :(get_embedding(M))
    end


    block = quote
        @inline function ($fname)($(callargs...); $(kwargs_list...)) where {$(where_exprs...)}
            M = $(argnames[manifold_arg_no])
            return ($fname_fwd)(($gft), $(argnames...); $(kwargs_call...))
        end
    end
    if :(EmbeddedForwardingType{DirectEmbedding}) in include_forwards.args
        block = quote
            $block
            @inline function ($fname_fwd)(
                    ::ManifoldsBase.EmbeddedForwardingType{ManifoldsBase.DirectEmbedding},
                    $(callargs...);
                    $(kwargs_list...),
                ) where {$(where_exprs...)}
                M = $(argnames[manifold_arg_no])
                return ($fname)(
                    $(argnames[1:(manifold_arg_no - 1)]...),
                    ($ge),
                    $(argnames[(manifold_arg_no + 1):end]...);
                    $(kwargs_call...),
                )
            end
        end
    end
    if :SimpleForwardingType in include_forwards.args
        block = quote
            $block
            @inline function ($fname_fwd)(
                    ::ManifoldsBase.SimpleForwardingType,
                    $(callargs...);
                    $(kwargs_list...),
                ) where {$(where_exprs...)}
                M = $(argnames[manifold_arg_no])
                return ($fname)(
                    $(argnames[1:(manifold_arg_no - 1)]...),
                    decorated_manifold(M),
                    $(argnames[(manifold_arg_no + 1):end]...);
                    $(kwargs_call...),
                )
            end
        end
    end
    if :StopForwardingType in include_forwards.args
        block = quote
            $block
            @inline function ($fname_fwd)(
                    ::ManifoldsBase.StopForwardingType,
                    $(callargs...);
                    $(kwargs_list...),
                ) where {$(where_exprs...)}
                M = $(argnames[manifold_arg_no])
                return invoke(
                    $fname,
                    Tuple{
                        $(argnametype_exprs[1:(manifold_arg_no - 1)]...),
                        ManifoldsBase.AbstractManifold,
                        $(argnametype_exprs[(manifold_arg_no + 1):end]...),
                    },
                    $(argnames[1:(manifold_arg_no - 1)]...),
                    M,
                    $(argnames[(manifold_arg_no + 1):end]...);
                    $(kwargs_call...),
                )
            end
        end
    end

    return esc(block)
end
