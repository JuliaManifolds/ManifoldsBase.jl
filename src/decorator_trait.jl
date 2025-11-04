#
# Base pass-ons
#
manifold_dimension(M::AbstractDecoratorManifold) = manifold_dimension(base_manifold(M))

#
# Forwarding types


"""
    abstract type AbstractEmbeddingDirectness end

Supertype for [`DirectEmbedding`](@ref) and [`IndirectEmbedding`](@ref) that indicate
whether [`embed`](@ref) on a manifold is an identity or not.
"""
abstract type AbstractEmbeddingDirectness end

"""
    struct DirectEmbedding <: AbstractEmbeddingDirectness end

A struct indicating that `embed` *is* an identity function on a manifold.
"""
struct DirectEmbedding <: AbstractEmbeddingDirectness end

"""
    struct IndirectEmbedding <: AbstractEmbeddingDirectness end

A struct indicating that `embed` *is not* an identity function on a manifold.
"""
struct IndirectEmbedding <: AbstractEmbeddingDirectness end

"""
    AbstractForwardingType

An abstract type to specify the forwarding behaviour of a function for a decorator
manifold or a trait within `ManifoldsBase.jl`.
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


A type that indicates that a function should not forward to a certain other manifold, e.g.
and embedding. This means that the user is asked to implement this function themselfes.
"""
struct StopForwardingType <: AbstractForwardingType end

"""
    SimpleForwardingType <: AbstractForwardingType

A type that indicates forwarding to the wrapped manifold without any changes.
"""
struct SimpleForwardingType <: AbstractForwardingType end

"""
    EmbeddedForwardingType{TED<:AbstractEmbeddingDirectness} <: AbstractEmbeddedForwardingType

A property of an embedded manifold that indicates that [`embed`](@ref) and [`project`](@ref)
are available and that a function using this trait type forwards to the embedding.
The type parameter `TED`, a subtype of [`AbstractEmbeddingDirectness`](@ref), indicates
whether `embed` on points and tangent vectors needs to be called or is an identity and can
be skipped.
"""
struct EmbeddedForwardingType{TED <: AbstractEmbeddingDirectness} <: AbstractEmbeddedForwardingType end

EmbeddedForwardingType(ed::AbstractEmbeddingDirectness = IndirectEmbedding()) = EmbeddedForwardingType{typeof(ed)}()

"""
    get_forwarding_type(M::AbstractManifold, f)
    get_forwarding_type(M::AbstractManifold, f, P::Type)

Get the type of forwarding to manifold wrapped by [`AbstractManifold`](@ref) `M`, for function `f`.
The returned value is an object of a subtype of [`AbstractForwardingType`](@ref).

Point type `P` can be optionally specified if different point types correspond
to different representations of the manifold and hence possibly different embeddings.
"""
get_forwarding_type(::AbstractManifold, f) = StopForwardingType()
get_forwarding_type(M::AbstractManifold, f, P::Type) = get_forwarding_type(M, f)

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
struct EmbeddedManifoldType{EN <: AbstractEmbeddingDirectness} <: AbstractEmbeddingType end

function EmbeddedManifoldType(en::AbstractEmbeddingDirectness = IndirectEmbedding())
    return EmbeddedManifoldType{typeof(en)}()
end

"""
    IsometricallyEmbeddedManifold <: AbstractEmbeddingType

A property to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is
an isometrically embedded manifold.

Here, additionally, metric related functions like [`inner`](@ref) and [`norm`](@ref) are passed to the embedding
"""
struct IsometricallyEmbeddedManifoldType{EN <: AbstractEmbeddingDirectness} <: AbstractEmbeddingType end

function IsometricallyEmbeddedManifoldType(
        en::AbstractEmbeddingDirectness = IndirectEmbedding(),
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
struct EmbeddedSubmanifoldType{EN <: AbstractEmbeddingDirectness} <: AbstractEmbeddingType end

function EmbeddedSubmanifoldType(en::AbstractEmbeddingDirectness = IndirectEmbedding())
    return EmbeddedSubmanifoldType{typeof(en)}()
end

"""
    get_embedding_type(M::AbstractManifold)
    get_embedding_type(M::AbstractManifold, P::Type)

Get embedding type of [`AbstractManifold`](@ref) `M`.
The returned value is an object of a subtype of [`AbstractEmbeddingType`](@ref), either of:

* [`NotEmbeddedManifoldType`](@ref) (default),
* [`EmbeddedManifoldType`](@ref),
* [`IsometricallyEmbeddedManifoldType`](@ref),
* [`EmbeddedSubmanifoldType`](@ref).

Point type `P` can be optionally specified if different point types correspond to different
embeddings.
"""
get_embedding_type(::AbstractManifold) = NotEmbeddedManifoldType()
get_embedding_type(M::AbstractManifold, ::Type) = get_embedding_type(M)


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
    get_embedding(M::AbstractDecoratorManifold, P::Type)

Specify the embedding of a manifold that has abstract decorators.
The embedding might depend on a point representation type `P`,
where different point representations
are distinguished as subtypes of [`AbstractManifoldPoint`](@ref).
A unique or default representation might also just be an `AbstractArray`.
"""
get_embedding(M::AbstractDecoratorManifold, ::Type) = get_embedding(M)

@inline function allocate_result(
        M::AbstractDecoratorManifold,
        f::TF,
        x::Vararg{Any, N},
    ) where {TF, N}
    return _allocate_result_forwarding(
        get_forwarding_type(M, f, typeof(x[1])),
        M,
        f,
        x...,
    )
end
@inline function allocate_result(M::AbstractDecoratorManifold, f::TF) where {TF}
    return _allocate_result_forwarding(get_forwarding_type(M, f), M, f)
end

@inline function _allocate_result_forwarding(
        ::EmbeddedForwardingType{DirectEmbedding},
        M::AbstractDecoratorManifold,
        f::TF,
        x::Vararg{Any, N},
    ) where {TF, N}
    return allocate_result(get_embedding(M, typeof(x[1])), f, x...)
end
@inline function _allocate_result_forwarding(
        ::SimpleForwardingType,
        M::AbstractDecoratorManifold,
        f::TF,
        x::Vararg{Any, N},
    ) where {TF, N}
    return allocate_result(decorated_manifold(M), f, x...)
end
@inline function _allocate_result_forwarding(
        ::Union{StopForwardingType, EmbeddedForwardingType},
        M::AbstractDecoratorManifold,
        f::TF,
        x::Vararg{Any, N},
    ) where {TF, N}
    return invoke(
        allocate_result,
        Tuple{AbstractManifold, typeof(f), typeof(x).parameters...},
        M,
        f,
        x...,
    )
end

# disambiguation
@invoke_maker 1 AbstractManifold allocate_result(
    M::AbstractDecoratorManifold, f::typeof(get_coordinates), p, X, B::AbstractBasis,
)
@invoke_maker 1 AbstractManifold allocate_result(
    M::AbstractDecoratorManifold, f::typeof(get_vector), p, c,
)

function allocate_result_embedding(
        M::AbstractManifold,
        f::typeof(embed),
        x::Vararg{Any, N},
    ) where {N}
    T = allocate_result_type(get_embedding(M, typeof(x[1])), f, x)
    return allocate(M, x[1], T, representation_size(get_embedding(M, typeof(x[1]))))
end
function allocate_result_embedding(
        M::AbstractManifold,
        f::typeof(project),
        x::Vararg{Any, N},
    ) where {N}
    T = allocate_result_type(M, f, x)
    return allocate(M, x[1], T, representation_size(M))
end

@trait_function change_metric(M::AbstractDecoratorManifold, G::AbstractMetric, X, p)

@trait_function change_metric!(M::AbstractDecoratorManifold, Y, G::AbstractMetric, X, p)

@trait_function change_representer(
    M::AbstractDecoratorManifold,
    G::AbstractMetric,
    X,
    p,
)
@trait_function change_representer!(
    M::AbstractDecoratorManifold,
    Y,
    G::AbstractMetric,
    X,
    p,
)

@trait_function check_size(M::AbstractDecoratorManifold, p) (
    SimpleForwardingType,
    StopForwardingType,
)
function _check_size_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p)
    mpe = check_size(get_embedding(M, typeof(p)), embed(M, p))
    if mpe !== nothing
        return ManifoldDomainError(
            "$p is not a point on $M because it is not a valid point in its embedding.", mpe,
        )
    end
    return nothing
end
function _check_size_forwarding(
        ::EmbeddedForwardingType{DirectEmbedding},
        M::AbstractDecoratorManifold,
        p,
    )
    mpe = check_size(get_embedding(M, typeof(p)), p)
    if mpe !== nothing
        return ManifoldDomainError(
            "$p is not a point on $M because it is not a valid point in its embedding.",
            mpe,
        )
    end
    return nothing
end
@trait_function check_size(M::AbstractDecoratorManifold, p, X) (
    SimpleForwardingType,
    StopForwardingType,
)
function _check_size_forwarding(
        ::EmbeddedForwardingType,
        M::AbstractDecoratorManifold,
        p,
        X,
    )
    mpe = check_size(get_embedding(M, typeof(p)), embed(M, p), embed(M, p, X))
    if mpe !== nothing
        return ManifoldDomainError(
            "$X is not a tangent vector at $p on $M because it is not a valid tangent vector in its embedding.",
            mpe,
        )
    end
    return nothing
end
function _check_size_forwarding(
        ::EmbeddedForwardingType{DirectEmbedding},
        M::AbstractDecoratorManifold,
        p,
        X,
    )
    mpe = check_size(get_embedding(M, typeof(p)), p, X)
    if mpe !== nothing
        return ManifoldDomainError(
            "$X is not a tangent vector at $p on $M because it is not a valid tangent vector in its embedding.",
            mpe,
        )
    end
    return nothing
end
# Introduce Deco Trait | automatic forward | fallback
@trait_function copyto!(M::AbstractDecoratorManifold, q, p)
@trait_function copyto!(M::AbstractDecoratorManifold, Y, p, X)

function _copyto!_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, q, p)
    return copyto!(get_embedding(M, typeof(p)), q, p)
end
function _copyto!_forwarding(
        ::EmbeddedForwardingType,
        M::AbstractDecoratorManifold,
        Y,
        p,
        X,
    )
    return copyto!(get_embedding(M, typeof(p)), Y, p, X)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed(M::AbstractDecoratorManifold, p)

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed!(M::AbstractDecoratorManifold, q, p)

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed(M::AbstractDecoratorManifold, p, X)

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed!(M::AbstractDecoratorManifold, Y, p, X)

# Introduce Deco Trait | automatic forward | fallback

@trait_function exp(M::AbstractDecoratorManifold, p, X)

@trait_function exp!(M::AbstractDecoratorManifold, q, p, X)

@trait_function exp_fused(M::AbstractDecoratorManifold, p, X, t::Number)

@trait_function exp_fused!(M::AbstractDecoratorManifold, q, p, X, t::Number)

@trait_function get_coordinates(M::AbstractDecoratorManifold, p, X) (
    SimpleForwardingType,
    StopForwardingType,
)

function _get_coordinates_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X)
    return @invoke get_coordinates(M::AbstractManifold, p, X)
end

@trait_function get_coordinates(M::AbstractDecoratorManifold, p, X, B::AbstractBasis) (
    SimpleForwardingType,
    StopForwardingType,
)

function _get_coordinates_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X, B::AbstractBasis)
    return @invoke get_coordinates(M::AbstractManifold, p, X, B)
end

@trait_function get_coordinates!(M::AbstractDecoratorManifold, c, p, X) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function get_coordinates!(M::AbstractDecoratorManifold, c, p, X, B::AbstractBasis) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function get_vector(M::AbstractDecoratorManifold, p, c) (
    SimpleForwardingType,
    StopForwardingType,
)

function _get_vector_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, c)
    return @invoke get_vector(M::AbstractManifold, p, c)
end

@trait_function get_vector(M::AbstractDecoratorManifold, p, c, B::AbstractBasis) (
    SimpleForwardingType,
    StopForwardingType,
)

function _get_vector_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, c, B::AbstractBasis)
    return @invoke get_vector(M::AbstractManifold, p, c, B)
end

@trait_function get_vector!(M::AbstractDecoratorManifold, X, p, c) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function get_vector!(M::AbstractDecoratorManifold, X, p, c, B::AbstractBasis) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function has_components(M::AbstractDecoratorManifold)

function _has_components_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold)
    return has_components(get_embedding(M))
end

@trait_function injectivity_radius(M::AbstractDecoratorManifold)

@trait_function injectivity_radius(M::AbstractDecoratorManifold, p)

@trait_function injectivity_radius(
    M::AbstractDecoratorManifold,
    m::AbstractRetractionMethod,
)

@trait_function injectivity_radius(
    M::AbstractDecoratorManifold,
    p,
    m::AbstractRetractionMethod,
)

@trait_function inner(M::AbstractDecoratorManifold, p, X, Y)

function _inner_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X, Y)
    return inner(get_embedding(M, typeof(p)), embed(M, p), embed(M, p, X), embed(M, p, Y))
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)

@trait_function inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)

@trait_function is_point(M::AbstractDecoratorManifold, p; kwargs...) (
    StopForwardingType,
    SimpleForwardingType,
)

function _is_point_forwarding(
        T::Union{
            EmbeddedForwardingType{D},
            IsometricallyEmbeddedManifoldType{D},
        },
        M::AbstractDecoratorManifold,
        p;
        error::Symbol = :none,
        kwargs...,
    ) where {D <: AbstractEmbeddingDirectness}
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
        if (D <: IndirectEmbedding)
            pt = is_point(get_embedding(M, typeof(p)), embed(M, p); error = error, kwargs...)
        else
            pt = is_point(get_embedding(M, typeof(p)), p; error = error, kwargs...)
        end
        !pt && return false # no error thrown (deactivated) but returned false -> return false
    catch e
        !(e isa DomainError || e isa AbstractManifoldDomainError) && rethrow(e)
        throw(
            ManifoldDomainError(
                "$p is not a point on $M because it is not a valid point in its embedding.",
                e,
            )
        )
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

@trait_function is_vector(
    M::AbstractDecoratorManifold,
    p,
    X,
    check_base_point::Bool = true;
    kwargs...,
) (StopForwardingType, SimpleForwardingType)

function _is_vector_forwarding(
        T::Union{
            EmbeddedForwardingType{D},
            IsometricallyEmbeddedManifoldType{D},
        },
        M::AbstractDecoratorManifold,
        p,
        X,
        check_base_point::Bool = true;
        error::Symbol = :none,
        kwargs...,
    ) where {D <: AbstractEmbeddingDirectness}
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
            !(e isa DomainError || e isa AbstractManifoldDomainError) && rethrow(e)
            throw(
                ManifoldDomainError(
                    "$X is not a tangent vector to $p on $M because $p is not a valid point on $M",
                    e,
                )
            )
        end
    end
    try
        if (D <: IndirectEmbedding)
            tv = is_vector(
                get_embedding(M, typeof(p)),
                embed(M, p),
                embed(M, p, X),
                check_base_point;
                error = error,
                kwargs...,
            )
        else
            tv = is_vector(
                get_embedding(M, typeof(p)),
                p,
                X,
                check_base_point;
                error = error,
                kwargs...,
            )
        end
        !tv && return false # no error thrown (deactivated) but returned false -> return false
    catch e
        !(e isa DomainError || e isa AbstractManifoldDomainError) && rethrow(e)
        throw(
            ManifoldDomainError(
                "$X is not a tangent vector to $p on $M because it is not a valid tangent vector in its embedding.",
                e,
            )
        )
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

@trait_function _isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
@trait_function _isapprox(M::AbstractDecoratorManifold, p, X, Y; kwargs...)

function __isapprox_forwarding(
        ::EmbeddedForwardingType,
        M::AbstractDecoratorManifold,
        p,
        q;
        kwargs...,
    )
    return _isapprox(get_embedding(M, typeof(p)), embed(M, p), embed(M, q); kwargs...)
end
function __isapprox_forwarding(
        ::EmbeddedForwardingType,
        M::AbstractDecoratorManifold,
        p,
        X,
        Y;
        kwargs...,
    )
    return _isapprox(
        get_embedding(M, typeof(p)),
        embed(M, p),
        embed(M, p, X),
        embed(M, p, Y);
        kwargs...,
    )
end

@trait_function is_flat(M::AbstractDecoratorManifold)

@trait_function norm(M::AbstractDecoratorManifold, p, X)

function _norm_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X)
    return norm(get_embedding(M, typeof(p)), embed(M, p), embed(M, p, X))
end

@trait_function log(M::AbstractDecoratorManifold, p, q)

@trait_function log!(M::AbstractDecoratorManifold, X, p, q)

# Introduce Deco Trait | automatic forward | fallback

@trait_function parallel_transport_direction(M::AbstractDecoratorManifold, p, X, d)

@trait_function parallel_transport_direction!(M::AbstractDecoratorManifold, Y, p, X, d)

# Introduce Deco Trait | automatic forward | fallback

@trait_function parallel_transport_to(M::AbstractDecoratorManifold, p, X, q)

@trait_function parallel_transport_to!(M::AbstractDecoratorManifold, Y, p, X, q)

@trait_function project(M::AbstractDecoratorManifold, p) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function project(M::AbstractDecoratorManifold, p, X) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function project!(M::AbstractDecoratorManifold, q, p) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function project!(M::AbstractDecoratorManifold, Y, p, X) (
    SimpleForwardingType,
    StopForwardingType,
)

@trait_function Random.rand(M::AbstractDecoratorManifold; kwargs...)

@trait_function Random.rand!(M::AbstractDecoratorManifold, p; kwargs...)

@trait_function Random.rand(rng::AbstractRNG, M::AbstractDecoratorManifold; kwargs...) (
    EmbeddedForwardingType{DirectEmbedding},
    SimpleForwardingType,
    StopForwardingType,
) 2

@trait_function Random.rand!(
    rng::AbstractRNG,
    M::AbstractDecoratorManifold,
    p;
    kwargs...,
) (EmbeddedForwardingType{DirectEmbedding}, SimpleForwardingType, StopForwardingType) 2

@trait_function representation_size(M::AbstractDecoratorManifold)
@trait_function representation_size(M::AbstractDecoratorManifold, T::Type)

# ... otherwise we fall back to the decorated manifold
function _representation_size_forwarding(
        ::AbstractForwardingType,
        M::AbstractDecoratorManifold,
        T::Type,
    )
    return representation_size(M)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

@trait_function retract_fused(
    M::AbstractDecoratorManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

@trait_function retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

@trait_function retract_fused!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)

@trait_function vector_transport_direction(
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

@trait_function vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

@trait_function vector_transport_to(
    M::AbstractDecoratorManifold, p, X, q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

@trait_function vector_transport_to!(
    M::AbstractDecoratorManifold, Y, p, X, q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)

@trait_function Weingarten(M::AbstractDecoratorManifold, p, X, V)
@trait_function Weingarten!(M::AbstractDecoratorManifold, Y, p, X, V)

@trait_function zero_vector(M::AbstractDecoratorManifold, p)

@trait_function zero_vector!(M::AbstractDecoratorManifold, X, p)

const forward_functions_embedded = [
    copyto!,
    check_size,
    has_components,
    is_point,
    is_vector,
    _isapprox,
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
function get_forwarding_type(M::AbstractDecoratorManifold, f, P::Type)
    return get_forwarding_type_embedding(get_embedding_type(M, P), M, f)
end

function get_forwarding_type_embedding(
        ::Union{
            EmbeddedManifoldType, IsometricallyEmbeddedManifoldType, EmbeddedSubmanifoldType,
            NotEmbeddedManifoldType,
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
                ::EmbeddedSubmanifoldType{DirectEmbedding},
                M::AbstractDecoratorManifold, ::typeof($mf),
            )
            return EmbeddedForwardingType()
        end
        function get_forwarding_type_embedding(
                ::EmbeddedSubmanifoldType{IndirectEmbedding},
                M::AbstractDecoratorManifold, ::typeof($mf),
            )
            return EmbeddedForwardingType(DirectEmbedding())
        end
    end
end

for mf in vcat(forward_functions_isometric, forward_functions_embedded)
    @eval begin
        function get_forwarding_type_embedding(
                ::IsometricallyEmbeddedManifoldType{DirectEmbedding},
                M::AbstractDecoratorManifold, ::typeof($mf),
            )
            return EmbeddedForwardingType()
        end
        function get_forwarding_type_embedding(
                ::IsometricallyEmbeddedManifoldType{IndirectEmbedding},
                M::AbstractDecoratorManifold, ::typeof($mf),
            )
            return EmbeddedForwardingType(DirectEmbedding())
        end
    end
end

for mf in forward_functions_embedded
    @eval begin
        function get_forwarding_type_embedding(
                ::EmbeddedManifoldType{DirectEmbedding},
                M::AbstractDecoratorManifold, ::typeof($mf),
            )
            return EmbeddedForwardingType()
        end
        function get_forwarding_type_embedding(
                ::EmbeddedManifoldType{IndirectEmbedding},
                M::AbstractDecoratorManifold, ::typeof($mf),
            )
            return EmbeddedForwardingType(DirectEmbedding())
        end
    end
end
