#
# Base pass-ons
#
manifold_dimension(M::AbstractDecoratorManifold) = manifold_dimension(base_manifold(M))

#
# Forwarding types

"""
    AbstractForwardingType

An abstract type to specify the forwarding behaviour of a function.
"""
abstract type AbstractForwardingType end

"""
    AbstractEmbeddedForwardingType

An abstract type to specify the forwarding behaviour of a function when it should forward
to the embedding of a manifold.
"""
abstract type AbstractEmbeddedForwardingType <: AbstractForwardingType end

"""
    StopForwardingType <: AbstractForwardingType

A property of an embedded manifold that indicates that `embed` and `project` are *not*
available.
"""
struct StopForwardingType <: AbstractForwardingType end

"""
    SimpleForwardingType <: AbstractForwardingType

A type that indicates forwarding to the wrapped manifold without any changes.
"""
struct SimpleForwardingType <: AbstractForwardingType end

"""
    EmbeddedForwardingType <: AbstractEmbeddedForwardingType

A property of an embedded manifold that indicates that [`embed`](@ref) and [`project`](@ref) are available
and that a function using this trait type forwards to the embedding using these.
"""
struct EmbeddedForwardingType <: AbstractEmbeddedForwardingType end

"""
    EmbeddedSimpleForwardingType <: AbstractEmbeddedForwardingType

A property of an embedded manifold that indicates that a function should forward to the embedding,
and even if [`embed`](@ref) and [`project`](@ref) are available, to not use these,
since they are the identity, but calling them might allocate unnecessarily.
"""
struct EmbeddedSimpleForwardingType <: AbstractEmbeddedForwardingType end


"""
    get_forwarding_type(M::AbstractManifold, f)
    get_forwarding_type(M::AbstractManifold, f, p)

Get the type of forwarding to manifold wrapped by [`AbstractManifold`](@ref) `M`, for function `f`.
The returned value is an object of a subtype of [`AbstractForwardingType`](@ref).

Point `p` can be optionally specified if different point types correspond
to different representations of the manifold and hence possibly different embeddings.
"""
get_forwarding_type(::AbstractManifold, f) = StopForwardingType()
get_forwarding_type(M::AbstractManifold, f, p) = get_forwarding_type(M, f)


abstract type AbstractEmbeddingNeed end

struct NeedsEmbedding <: AbstractEmbeddingNeed end
struct DoesntNeedEmbedding <: AbstractEmbeddingNeed end

"""
    AbstractEmbeddingType

Within all [`AbstractEmbeddedForwardingType`](@ref)s this type is used to indicate different kinds of embeddings,
for example the default fallback that [`NotEmbeddedManifoldType`](@ref) a manifold is not embedded,
that is is embedded using [`EmbeddedManifoldType`](@ref) or even specifying further that it is
isometrically embedded using [`IsometricallyEmbeddedManifoldType`](@ref) or as furthermore
a submanifold using [`EmbeddedSubmanifoldType`](@ref).
"""
abstract type AbstractEmbeddingType end

"""
    NotEmbeddedManifoldType <: AbstractEmbeddingType

A property of an embedded manifold that indicates that `embed` and `project` are *not*
available.
"""
struct NotEmbeddedManifoldType <: AbstractEmbeddingType end

"""
    EmbeddedManifoldType <: AbstractEmbeddingType

A property of an embedded manifold that indicates that `embed` and `project` are available.
"""
struct EmbeddedManifoldType{EN<:AbstractEmbeddingNeed} <: AbstractEmbeddingType end

function EmbeddedManifoldType(en::AbstractEmbeddingNeed = DoesntNeedEmbedding())
    return EmbeddedManifoldType{typeof(en)}()
end

"""
    IsometricallyEmbeddedManifold <: AbstractEmbeddingType

A property to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is
an isometrically embedded manifold.

Here, additionally, metric related functions like [`inner`](@ref) and [`norm`](@ref) are passed to the embedding
"""
struct IsometricallyEmbeddedManifoldType{EN<:AbstractEmbeddingNeed} <: AbstractEmbeddingType end

function IsometricallyEmbeddedManifoldType(
    en::AbstractEmbeddingNeed = DoesntNeedEmbedding(),
)
    return IsometricallyEmbeddedManifoldType{typeof(en)}()
end

"""
    EmbeddedSubmanifoldType <: AbstractEmbeddingType

A property to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
It is a special case of the [`IsometricallyEmbeddedManifoldType`](@ref) property, i.e. it has all properties of
this property.

In this property, additionally to the isometric embedded manifold, all retractions, inverse retractions,
and vectors transports, especially [`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref)
are passed to the embedding.
"""
struct EmbeddedSubmanifoldType{EN<:AbstractEmbeddingNeed} <: AbstractEmbeddingType end

function EmbeddedSubmanifoldType(en::AbstractEmbeddingNeed = DoesntNeedEmbedding())
    return EmbeddedSubmanifoldType{typeof(en)}()
end

"""
    get_embedding_type(M::AbstractManifold)
    get_embedding_type(M::AbstractManifold, p)

Get embedding type of [`AbstractManifold`](@ref) `M`.
The returned value is an object of a subtype of [`AbstractEmbeddingType`](@ref), either of:
* [`NotEmbeddedManifoldType`](@ref) (default),
* [`EmbeddedManifoldType`](@ref),
* [`IsometricallyEmbeddedManifoldType`](@ref),
* [`EmbeddedSubmanifoldType`](@ref).

Point `p` can be optionally specified if different point types correspond to different
embeddings.
"""
get_embedding_type(::AbstractManifold) = NotEmbeddedManifoldType()
get_embedding_type(M::AbstractManifold, p) = get_embedding_type(M)

function is_embedded_manifold(M::AbstractManifold)
    return get_embedding_type(M) !== NotEmbeddedManifoldType()
end


#
# Generic Decorator functions
@doc raw"""
    decorated_manifold(M::AbstractDecoratorManifold)

For a manifold `M` that is decorated with some properties, this function returns
the manifold without that manifold, i.e. the manifold that _was decorated_.
"""
decorated_manifold(M::AbstractDecoratorManifold)
decorated_manifold(M::AbstractManifold) = M

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
# Embedded specific functions.
"""
    get_embedding(M::AbstractDecoratorManifold)
    get_embedding(M::AbstractDecoratorManifold, p)

Specify the embedding of a manifold that has abstract decorators.
the embedding might depend on a point representation, where different point representations
are distinguished as subtypes of [`AbstractManifoldPoint`](@ref).
A unique or default representation might also just be an `AbstractArray`.
"""
get_embedding(M::AbstractDecoratorManifold, p) = get_embedding(M)

@inline function allocate_result(
    M::AbstractDecoratorManifold,
    f::TF,
    x::Vararg{Any,N},
) where {TF,N}
    return allocate_result(trait(allocate_result, M, f, x...), M, f, x...)
end
# disambiguation
@invoke_maker 1 AbstractManifold allocate_result(
    M::AbstractDecoratorManifold,
    f::typeof(get_coordinates),
    p,
    X,
    B::AbstractBasis,
)
@invoke_maker 1 AbstractManifold allocate_result(
    M::AbstractDecoratorManifold,
    f::typeof(get_vector),
    p,
    c,
)

# Introduce fallback
@inline function allocate_result(
    ::EmptyTrait,
    M::AbstractManifold,
    f::TF,
    x::Vararg{Any,N},
) where {TF,N}
    return invoke(
        allocate_result,
        Tuple{AbstractManifold,typeof(f),typeof(x).parameters...},
        M,
        f,
        x...,
    )
end
# Introduce automatic forward
@inline function allocate_result(
    t::TraitList,
    M::AbstractManifold,
    f::TF,
    x::Vararg{Any,N},
) where {TF,N}
    return allocate_result(next_trait(t), M, f, x...)
end
function allocate_result_embedding(
    M::AbstractManifold,
    f::typeof(embed),
    x::Vararg{Any,N},
) where {N}
    T = allocate_result_type(get_embedding(M, x[1]), f, x)
    return allocate(M, x[1], T, representation_size(get_embedding(M, x[1])))
end
function allocate_result_embedding(
    M::AbstractManifold,
    f::typeof(project),
    x::Vararg{Any,N},
) where {N}
    T = allocate_result_type(M, f, x)
    return allocate(M, x[1], T, representation_size(M))
end
@inline function allocate_result(
    ::TraitList{IsExplicitDecorator},
    M::AbstractDecoratorManifold,
    f::TF,
    x::Vararg{Any,N},
) where {TF,N}
    return allocate_result(decorated_manifold(M), f, x...)
end

@new_trait_function change_metric(M::AbstractDecoratorManifold, G::AbstractMetric, X, p)

@new_trait_function change_metric!(M::AbstractDecoratorManifold, Y, G::AbstractMetric, X, p)

@new_trait_function change_representer(
    M::AbstractDecoratorManifold,
    G::AbstractMetric,
    X,
    p,
)
@new_trait_function change_representer!(
    M::AbstractDecoratorManifold,
    Y,
    G::AbstractMetric,
    X,
    p,
)

@new_trait_function check_size(M::AbstractDecoratorManifold, p) (
    SimpleForwardingType,
    StopForwardingType,
)
function _check_size_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p)
    mpe = check_size(get_embedding(M, p), embed(M, p))
    if mpe !== nothing
        return ManifoldDomainError(
            "$p is not a point on $M because it is not a valid point in its embedding.",
            mpe,
        )
    end
    return nothing
end
function _check_size_forwarding(
    ::EmbeddedSimpleForwardingType,
    M::AbstractDecoratorManifold,
    p,
)
    mpe = check_size(get_embedding(M, p), p)
    if mpe !== nothing
        return ManifoldDomainError(
            "$p is not a point on $M because it is not a valid point in its embedding.",
            mpe,
        )
    end
    return nothing
end
@new_trait_function check_size(M::AbstractDecoratorManifold, p, X) (
    SimpleForwardingType,
    StopForwardingType,
)
function _check_size_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
)
    mpe = check_size(get_embedding(M, p), embed(M, p), embed(M, p, X))
    if mpe !== nothing
        return ManifoldDomainError(
            "$X is not a tangent vector at $p on $M because it is not a valid tangent vector in its embedding.",
            mpe,
        )
    end
    return nothing
end
function _check_size_forwarding(
    ::EmbeddedSimpleForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
)
    mpe = check_size(get_embedding(M, p), p, X)
    if mpe !== nothing
        return ManifoldDomainError(
            "$X is not a tangent vector at $p on $M because it is not a valid tangent vector in its embedding.",
            mpe,
        )
    end
    return nothing
end
# Introduce Deco Trait | automatic forward | fallback
@new_trait_function copyto!(M::AbstractDecoratorManifold, q, p)
@new_trait_function copyto!(M::AbstractDecoratorManifold, Y, p, X)

function _copyto!_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, q, p)
    return copyto!(get_embedding(M, p), q, p)
end
function _copyto!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
)
    return copyto!(get_embedding(M, p), Y, p, X)
end

# Introduce Deco Trait | automatic forward | fallback
@new_trait_function embed(M::AbstractDecoratorManifold, p)

# Introduce Deco Trait | automatic forward | fallback
@new_trait_function embed!(M::AbstractDecoratorManifold, q, p)

# Introduce Deco Trait | automatic forward | fallback
@new_trait_function embed(M::AbstractDecoratorManifold, p, X)

# Introduce Deco Trait | automatic forward | fallback
@new_trait_function embed!(M::AbstractDecoratorManifold, Y, p, X)

# Introduce Deco Trait | automatic forward | fallback

@new_trait_function exp(M::AbstractDecoratorManifold, p, X)

function _exp_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X)
    return exp(get_embedding(M, p), embed(M, p), embed(M, p, X))
end

@new_trait_function exp!(M::AbstractDecoratorManifold, q, p, X)

function _exp!_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, q, p, X)
    return exp!(get_embedding(M, p), q, embed(M, p), embed(M, p, X))
end

@new_trait_function exp_fused(M::AbstractDecoratorManifold, p, X, t::Number)

function _exp_fused_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    t::Number,
)
    return exp_fused(get_embedding(M, p), embed(M, p), embed(M, p, X), t)
end

@new_trait_function exp_fused!(M::AbstractDecoratorManifold, q, p, X, t::Number)

function _exp_fused!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    t::Number,
)
    return exp_fused!(get_embedding(M, p), q, embed(M, p), embed(M, p, X), t)
end

@new_trait_function has_components(M::AbstractDecoratorManifold)

function _has_components_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold)
    return has_components(get_embedding(M))
end

@new_trait_function injectivity_radius(M::AbstractDecoratorManifold)

function _injectivity_radius_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
)
    return injectivity_radius(get_embedding(M))
end

@new_trait_function injectivity_radius(M::AbstractDecoratorManifold, p)

function _injectivity_radius_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
)
    return injectivity_radius(get_embedding(M, p), embed(M, p))
end

@new_trait_function injectivity_radius(
    M::AbstractDecoratorManifold,
    m::AbstractRetractionMethod,
)

function _injectivity_radius_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    m::AbstractRetractionMethod,
)
    return injectivity_radius(get_embedding(M), m)
end

@new_trait_function injectivity_radius(
    M::AbstractDecoratorManifold,
    p,
    m::AbstractRetractionMethod,
)

function _injectivity_radius_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    m::AbstractRetractionMethod,
)
    return injectivity_radius(get_embedding(M, p), embed(M, p), m)
end

@new_trait_function inner(M::AbstractDecoratorManifold, p, X, Y)

function _inner_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X, Y)
    return inner(get_embedding(M, p), embed(M, p), embed(M, p, X), embed(M, p, Y))
end

# Introduce Deco Trait | automatic forward | fallback
@new_trait_function inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)

function _inverse_retract_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)
    return inverse_retract(get_embedding(M, p), embed(M, p), embed(M, q), m)
end


@new_trait_function inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)

function _inverse_retract!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)
    return inverse_retract!(get_embedding(M, p), X, embed(M, p), embed(M, q), m)
end

@new_trait_function is_point(M::AbstractDecoratorManifold, p; kwargs...) (
    StopForwardingType,
    SimpleForwardingType,
)

function _is_point_forwarding(
    T::Union{
        EmbeddedForwardingType,
        EmbeddedSimpleForwardingType,
        IsometricallyEmbeddedManifoldType,
    },
    M::AbstractDecoratorManifold,
    p;
    error::Symbol = :none,
    kwargs...,
)
    # to be safe check_size first
    es = check_size(M, p)
    if es !== nothing
        (error === :error) && throw(es)
        s = "$(typeof(es)) with $(es)"
        (error === :info) && @info s
        (error === :warn) && @warn s
        return false
    end
    try
        if T isa EmbeddedForwardingType
            pt = is_point(get_embedding(M, p), embed(M, p); error = error, kwargs...)
        else
            pt = is_point(get_embedding(M, p), p; error = error, kwargs...)
        end
        !pt && return false # no error thrown (deactivated) but returned false -> return false
    catch e
        if e isa DomainError || e isa AbstractManifoldDomainError
            e = ManifoldDomainError(
                "$p is not a point on $M because it is not a valid point in its embedding.",
                e,
            )
        end
        throw(e) #an error occured that we do not handle ourselves -> rethrow.
    end
    mpe = check_point(M, p; kwargs...)
    if mpe !== nothing
        (error === :error) && throw(mpe)
        # else: collect and info showerror
        io = IOBuffer()
        showerror(io, mpe)
        s = String(take!(io))
        (error === :info) && @info s
        (error === :warn) && @warn s
        return false
    end
    return true
end

@new_trait_function is_vector(
    M::AbstractDecoratorManifold,
    p,
    X,
    check_base_point::Bool = true;
    kwargs...,
) (StopForwardingType, SimpleForwardingType)

function _is_vector_forwarding(
    T::Union{
        EmbeddedForwardingType,
        EmbeddedSimpleForwardingType,
        IsometricallyEmbeddedManifoldType,
    },
    M::AbstractDecoratorManifold,
    p,
    X,
    check_base_point::Bool = true;
    error::Symbol = :none,
    kwargs...,
)
    es = check_size(M, p, X)
    if es !== nothing
        (error === :error) && throw(es)
        # else: collect and info showerror
        io = IOBuffer()
        showerror(io, es)
        s = String(take!(io))
        (error === :info) && @info s
        (error === :warn) && @warn s
        return false
    end
    if check_base_point
        try
            ep = is_point(M, p; error = error, kwargs...)
            !ep && return false
        catch e
            if e isa DomainError || e isa AbstractManifoldDomainError
                ManifoldDomainError(
                    "$X is not a tangent vector to $p on $M because $p is not a valid point on $p",
                    e,
                )
            end
            throw(e)
        end
    end
    try
        if T isa EmbeddedForwardingType
            tv = is_vector(
                get_embedding(M, p),
                embed(M, p),
                embed(M, p, X),
                check_base_point;
                error = error,
                kwargs...,
            )
        else
            tv = is_vector(
                get_embedding(M, p),
                p,
                X,
                check_base_point;
                error = error,
                kwargs...,
            )
        end
        !tv && return false # no error thrown (deactivated) but returned false -> return false
    catch e
        if e isa DomainError || e isa AbstractManifoldDomainError
            e = ManifoldDomainError(
                "$X is not a tangent vector to $p on $M because it is not a valid tangent vector in its embedding.",
                e,
            )
        end
        throw(e)
    end
    # Check (additional) local stuff
    mXe = check_vector(M, p, X; kwargs...)
    mXe === nothing && return true
    (error === :error) && throw(mXe)
    # else: collect and info showerror
    io = IOBuffer()
    showerror(io, mXe)
    s = String(take!(io))
    (error === :info) && @info s
    (error === :warn) && @warn s
    return false
end

@new_trait_function isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
@new_trait_function isapprox(M::AbstractDecoratorManifold, p, X, Y; kwargs...)

function _isapprox_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    q;
    kwargs...,
)
    return isapprox(get_embedding(M, p), embed(M, p), embed(M, q); kwargs...)
end
function _isapprox_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    Y;
    kwargs...,
)
    return isapprox(
        get_embedding(M, p),
        embed(M, p),
        embed(M, p, X),
        embed(M, p, Y);
        kwargs...,
    )
end

@new_trait_function is_flat(M::AbstractDecoratorManifold)

@new_trait_function norm(M::AbstractDecoratorManifold, p, X)

function _norm_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X)
    return norm(get_embedding(M, p), embed(M, p), embed(M, p, X))
end

@new_trait_function log(M::AbstractDecoratorManifold, p, q)

function _log_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, q)
    return log(get_embedding(M, p), embed(M, p), embed(M, q))
end

@new_trait_function log!(M::AbstractDecoratorManifold, X, p, q)

function _log!_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, X, p, q)
    return log!(get_embedding(M, p), X, embed(M, p), embed(M, q))
end

# Introduce Deco Trait | automatic forward | fallback

@new_trait_function parallel_transport_direction(M::AbstractDecoratorManifold, p, X, d)

function _parallel_transport_direction_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
)
    return parallel_transport_direction(
        get_embedding(M, p),
        embed(M, p),
        embed(M, p, X),
        embed(M, p, d),
    )
end

@new_trait_function parallel_transport_direction!(M::AbstractDecoratorManifold, Y, p, X, d)

function _parallel_transport_direction!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
)
    return parallel_transport_direction!(
        get_embedding(M, p),
        Y,
        embed(M, p),
        embed(M, p, X),
        embed(M, p, d),
    )
end

# Introduce Deco Trait | automatic forward | fallback

@new_trait_function parallel_transport_to(M::AbstractDecoratorManifold, p, X, q)

function _parallel_transport_to_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
)
    return parallel_transport_to(
        get_embedding(M, p),
        embed(M, p),
        embed(M, p, X),
        embed(M, q),
    )
end

@new_trait_function parallel_transport_to!(M::AbstractDecoratorManifold, Y, p, X, q)

function _parallel_transport_to!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
)
    return parallel_transport_to!(
        get_embedding(M, p),
        Y,
        embed(M, p),
        embed(M, p, X),
        embed(M, q),
    )
end

@new_trait_function Random.rand(M::AbstractDecoratorManifold; kwargs...)

function _rand_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold; kwargs...)
    return rand(get_embedding(M); kwargs...)
end

@new_trait_function Random.rand!(M::AbstractDecoratorManifold, p; kwargs...)

function _rand!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p;
    kwargs...,
)
    return rand(get_embedding(M, p), p; kwargs...)
end

@new_trait_function Random.rand(rng::AbstractRNG, M::AbstractDecoratorManifold; kwargs...) (
    EmbeddedSimpleForwardingType,
    SimpleForwardingType,
    StopForwardingType,
) 2

function _rand_forwarding(
    ::EmbeddedForwardingType,
    rng::AbstractRNG,
    M::AbstractDecoratorManifold;
    kwargs...,
)
    return rand(rng, get_embedding(M); kwargs...)
end

@new_trait_function Random.rand!(
    rng::AbstractRNG,
    M::AbstractDecoratorManifold,
    p;
    kwargs...,
) (EmbeddedSimpleForwardingType, SimpleForwardingType, StopForwardingType) 2

function _rand!_forwarding(
    ::EmbeddedForwardingType,
    rng::AbstractRNG,
    M::AbstractDecoratorManifold,
    p;
    kwargs...,
)
    return rand!(rng, get_embedding(M, p), p; kwargs...)
end

@new_trait_function representation_size(M::AbstractDecoratorManifold)

function _representation_size_forwarding(
    ::Union{EmbeddedForwardingType,EmbeddedSimpleForwardingType},
    M::AbstractDecoratorManifold,
)
    return representation_size(get_embedding(M))
end


# Introduce Deco Trait | automatic forward | fallback
@new_trait_function retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

function _retract_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod,
)
    return retract(get_embedding(M, p), embed(M, p), embed(M, p, X), m)
end

@new_trait_function retract_fused(
    M::AbstractDecoratorManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

function _retract_fused_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod,
)
    return retract_fused(get_embedding(M, p), embed(M, p), embed(M, p, X), t, m)
end

@new_trait_function retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

function _retract!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod,
)
    return retract!(get_embedding(M, p), q, embed(M, p), embed(M, p, X), m)
end

@new_trait_function retract_fused!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

function _retract_fused!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod,
)
    return retract_fused!(get_embedding(M, p), q, embed(M, p), embed(M, p, X), t, m)
end


@new_trait_function vector_transport_direction(
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

function _vector_transport_direction_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod,
)
    return vector_transport_direction(
        get_embedding(M, p),
        embed(M, p),
        embed(M, p, X),
        embed(M, d),
        m,
    )
end

@new_trait_function vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

function _vector_transport_direction!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod,
)
    return vector_transport_direction!(
        get_embedding(M, p),
        Y,
        embed(M, p),
        embed(M, p, X),
        embed(M, p, d),
        m,
    )
end

@new_trait_function vector_transport_to(
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

function _vector_transport_to_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)
    return vector_transport_to(
        get_embedding(M, p),
        embed(M, p),
        embed(M, p, X),
        embed(M, q),
        m,
    )
end

@new_trait_function vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

function _vector_transport_to!_forwarding(
    ::EmbeddedForwardingType,
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)
    return vector_transport_to!(
        get_embedding(M, p),
        Y,
        embed(M, p),
        embed(M, p, X),
        embed(M, q),
        m,
    )
end

@new_trait_function Weingarten(M::AbstractDecoratorManifold, p, X, V)
@new_trait_function Weingarten!(M::AbstractDecoratorManifold, Y, p, X, V)

@new_trait_function zero_vector(M::AbstractDecoratorManifold, p)

@new_trait_function zero_vector!(M::AbstractDecoratorManifold, X, p)


const forward_functions_embedded = [
    copyto!,
    check_size,
    has_components,
    is_point,
    is_vector,
    isapprox,
    representation_size,
    zero_vector,
    zero_vector!,
]

const forward_functions_isometric = [inner, norm]

const forward_functions_submanifold = [
    change_metric,
    change_metric!,
    change_representer,
    change_representer!,
    exp,
    exp!,
    exp_fused,
    exp_fused!,
    get_basis,
    get_coordinates,
    get_vector,
    get_vectors,
    injectivity_radius,
    inverse_retract,
    inverse_retract!,
    is_flat,
    log,
    log!,
    mid_point,
    parallel_transport_direction,
    parallel_transport_direction!,
    parallel_transport_to,
    parallel_transport_to!,
    retract,
    retract!,
    retract_fused,
    retract_fused!,
    riemann_tensor,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
    vector_transport_to!,
    Weingarten,
    Weingarten!,
]


function get_forwarding_type(M::AbstractDecoratorManifold, f)
    return get_forwarding_type_embedding(get_embedding_type(M), M, f)
end
function get_forwarding_type(M::AbstractDecoratorManifold, f, p)
    return get_forwarding_type_embedding(get_embedding_type(M, p), M, f)
end

function get_forwarding_type_embedding(
    ::Union{
        EmbeddedManifoldType,
        IsometricallyEmbeddedManifoldType,
        NotEmbeddedManifoldType,
        EmbeddedSubmanifoldType,
    },
    M::AbstractDecoratorManifold,
    f,
)
    return StopForwardingType()
end

for mf in vcat(
    forward_functions_submanifold,
    forward_functions_isometric,
    forward_functions_embedded,
)
    @eval begin
        function get_forwarding_type_embedding(
            ::EmbeddedSubmanifoldType{NeedsEmbedding},
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedForwardingType()
        end
        function get_forwarding_type_embedding(
            ::EmbeddedSubmanifoldType{DoesntNeedEmbedding},
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedSimpleForwardingType()
        end
    end
end

for mf in vcat(forward_functions_isometric, forward_functions_embedded)
    @eval begin
        function get_forwarding_type_embedding(
            ::IsometricallyEmbeddedManifoldType{NeedsEmbedding},
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedForwardingType()
        end
        function get_forwarding_type_embedding(
            ::IsometricallyEmbeddedManifoldType{DoesntNeedEmbedding},
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedSimpleForwardingType()
        end
    end
end

for mf in forward_functions_embedded
    @eval begin
        function get_forwarding_type_embedding(
            ::EmbeddedManifoldType{NeedsEmbedding},
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedForwardingType()
        end
        function get_forwarding_type_embedding(
            ::EmbeddedManifoldType{DoesntNeedEmbedding},
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedSimpleForwardingType()
        end
    end
end
