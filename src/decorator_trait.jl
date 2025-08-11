#
# Base pass-ons
#
manifold_dimension(M::AbstractDecoratorManifold) = manifold_dimension(base_manifold(M))

#
# Traits - each passed to a function that is properly documented
#

"""
    IsIsometricManifoldEmbeddedManifold <: AbstractTrait

A Trait to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is
an isometrically embedded manifold.

Here, additionally, metric related functions like [`inner`](@ref) and [`norm`](@ref) are passed to the embedding
"""
struct IsIsometricEmbeddedManifold <: AbstractTrait end

"""
    IsEmbeddedSubmanifold <: AbstractTrait

A trait to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
It is a special case of the [`IsIsometricEmbeddedManifold`](@ref) trait, i.e. it has all properties of
this trait.

In this trait, additionally to the isometric embedded manifold, all retractions, inverse retractions,
and vectors transports, especially [`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref)
are passed to the embedding.
"""
struct IsEmbeddedSubmanifold <: AbstractTrait end

parent_trait(::IsEmbeddedSubmanifold) = IsIsometricEmbeddedManifold()


abstract type AbstractForwardingType end

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
    EmbeddedForwardingType <: AbstractForwardingType

A property of an embedded manifold that indicates that `embed` and `project` are available.
"""
struct EmbeddedForwardingType <: AbstractForwardingType end

"""
    EmbeddedForwardingType <: AbstractForwardingType

A property of an embedded manifold that indicates that `embed` and `project` are available,
although not needed to propagate a function to the embedding.
"""
struct EmbeddedSimpleForwardingType <: AbstractForwardingType end


"""
    get_forwarding_type(M::AbstractManifold, f)
    get_forwarding_type(M::AbstractManifold, f, p)

Get the type of forwarding to manifold wrapped by [`AbstractManifold`](@ref) `M`, for function `f`.
The returned value is an object of a subtype of [`AbstractForwardingType`](@ref), either of:
* [`StopForwardingType`](@ref) (default),
* [`SimpleForwardingType`](@ref),
* [`EmbeddedForwardingType`](@ref).

Point `p` can be optionally specified if different point types correspond to different
embeddings.
"""
get_forwarding_type(::AbstractManifold, f) = StopForwardingType()
get_forwarding_type(M::AbstractManifold, f, p) = get_forwarding_type(M, f)

abstract type AbstractEmbeddingType end

"""
    NotEmbeddedManifoldType <: AbstractEmbeddingType

A property of an embedded manifold that indicates that `embed` and `project` are *not*
available.
"""

struct NotEmbeddedManifoldType <: AbstractEmbeddingType end

"""
    EmbeddedManifold <: AbstractEmbeddingType

A property of an embedded manifold that indicates that `embed` and `project` are available.
"""
struct EmbeddedManifoldType <: AbstractEmbeddingType end

"""
    IsometricallyEmbeddedManifold <: AbstractEmbeddingType

A property to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is
an isometrically embedded manifold.

Here, additionally, metric related functions like [`inner`](@ref) and [`norm`](@ref) are passed to the embedding
"""
struct IsometricallyEmbeddedManifoldType <: AbstractEmbeddingType end

"""
    EmbeddedSubmanifoldType <: AbstractEmbeddingType

A property to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
It is a special case of the [`IsIsometricEmbeddedManifold`](@ref) property, i.e. it has all properties of
this property.

In this property, additionally to the isometric embedded manifold, all retractions, inverse retractions,
and vectors transports, especially [`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref)
are passed to the embedding.
"""
struct EmbeddedSubmanifoldType <: AbstractEmbeddingType end

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
@trait_function decorated_manifold(M::AbstractDecoratorManifold)

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

#
# -----------------------------------------------------------------------------------------
# This is one new function

# Introduction and default fallbacks could become a macro?
# Introduce trait
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


@trait_function change_metric(M::AbstractDecoratorManifold, G::AbstractMetric, X, p)
@trait_function change_metric!(M::AbstractDecoratorManifold, Y, G::AbstractMetric, X, p)

@trait_function change_representer(M::AbstractDecoratorManifold, G::AbstractMetric, X, p)
@trait_function change_representer!(
    M::AbstractDecoratorManifold,
    Y,
    G::AbstractMetric,
    X,
    p,
)

# Introduce Deco Trait | automatic forward | fallback
@trait_function check_size(M::AbstractDecoratorManifold, p)
# Embedded
function check_size_embedding(M::AbstractDecoratorManifold, p)
    mpe = check_size(get_embedding(M, p), embed(M, p))
    if mpe !== nothing
        return ManifoldDomainError(
            "$p is not a point on $M because it is not a valid point in its embedding.",
            mpe,
        )
    end
    return nothing
end
# Introduce Deco Trait | automatic forward | fallback
@trait_function check_size(M::AbstractDecoratorManifold, p, X)
# Embedded
function check_size_embedding(M::AbstractDecoratorManifold, p, X)
    mpe = check_size(get_embedding(M, p), embed(M, p), embed(M, p, X))
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

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed(M::AbstractDecoratorManifold, p)
# EmbeddedManifold

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed!(M::AbstractDecoratorManifold, q, p)
# EmbeddedManifold

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed(M::AbstractDecoratorManifold, p, X)
# EmbeddedManifold

# Introduce Deco Trait | automatic forward | fallback
@trait_function embed!(M::AbstractDecoratorManifold, Y, p, X)
# EmbeddedManifold

# Introduce Deco Trait | automatic forward | fallback
@trait_function exp(M::AbstractDecoratorManifold, p, X)
@trait_function exp_fused(M::AbstractDecoratorManifold, p, X, t::Number)
# EmbeddedSubManifold
function exp(::TraitList{IsEmbeddedSubmanifold}, M::AbstractDecoratorManifold, p, X)
    return exp(get_embedding(M, p), p, X)
end
function exp_fused(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    t::Number,
)
    return exp_fused(get_embedding(M, p), p, X, t)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function exp!(M::AbstractDecoratorManifold, q, p, X)
@trait_function exp_fused!(M::AbstractDecoratorManifold, q, p, X, t::Number)
# EmbeddedSubManifold
function exp!(::TraitList{IsEmbeddedSubmanifold}, M::AbstractDecoratorManifold, q, p, X)
    return exp!(get_embedding(M, p), q, p, X)
end
function exp_fused!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    t::Number,
)
    return exp_fused!(get_embedding(M, p), q, p, X, t)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function get_basis(M::AbstractDecoratorManifold, p, B::AbstractBasis)

# Introduce Deco Trait | automatic forward | fallback
@trait_function get_coordinates(M::AbstractDecoratorManifold, p, X, B::AbstractBasis)

# Introduce Deco Trait | automatic forward | fallback
@trait_function get_coordinates!(M::AbstractDecoratorManifold, Y, p, X, B::AbstractBasis)

# Introduce Deco Trait | automatic forward | fallback
@trait_function get_vector(M::AbstractDecoratorManifold, p, c, B::AbstractBasis)

# Introduce Deco Trait | automatic forward | fallback
@trait_function get_vector!(M::AbstractDecoratorManifold, Y, p, c, B::AbstractBasis)

# Introduce Deco Trait | automatic forward | fallback
@trait_function get_vectors(M::AbstractDecoratorManifold, p, B::AbstractBasis)

@trait_function has_components(M::AbstractDecoratorManifold)

@trait_function injectivity_radius(M::AbstractDecoratorManifold)
function injectivity_radius(
    ::TraitList{IsIsometricEmbeddedManifold},
    M::AbstractDecoratorManifold,
)
    return injectivity_radius(get_embedding(M))
end
@trait_function injectivity_radius(M::AbstractDecoratorManifold, p)
function injectivity_radius(
    ::TraitList{IsIsometricEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
)
    return injectivity_radius(get_embedding(M, p), p)
end
@trait_function injectivity_radius(
    M::AbstractDecoratorManifold,
    m::AbstractRetractionMethod,
)
function injectivity_radius(
    ::TraitList{IsIsometricEmbeddedManifold},
    M::AbstractDecoratorManifold,
    m::AbstractRetractionMethod,
)
    return injectivity_radius(get_embedding(M), m)
end
@trait_function injectivity_radius(
    M::AbstractDecoratorManifold,
    p,
    m::AbstractRetractionMethod,
)
function injectivity_radius(
    ::TraitList{IsIsometricEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    m::AbstractRetractionMethod,
)
    return injectivity_radius(get_embedding(M, p), p, m)
end

function inner(M::AbstractDecoratorManifold, p, X, Y)
    return _inner_forwarding(get_forwarding_type(M, inner, p), M, p, X, Y)
end

function _inner_forwarding(::StopForwardingType, M::AbstractDecoratorManifold, p, X, Y)
    return invoke(inner, Tuple{AbstractManifold,typeof(p),typeof(X),typeof(Y)}, M, p, X, Y)
end
function _inner_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, X, Y)
    return inner(get_embedding(M, p), p, X, Y)
end
function _inner_forwarding(
    ::EmbeddedSimpleForwardingType,
    M::AbstractDecoratorManifold,
    p,
    X,
    Y,
)
    return inner(get_embedding(M, p), p, X, Y)
end
function _inner_forwarding(::SimpleForwardingType, M::AbstractDecoratorManifold, p, X, Y)
    return inner(decorated_manifold(M), p, X, Y)
end


# Introduce Deco Trait | automatic forward | fallback
@trait_function inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)
# Transparent for Submanifolds
function inverse_retract(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)
    return inverse_retract(get_embedding(M, p), p, q, m)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function inverse_retract!(M::AbstractDecoratorManifold, X, p, q)
@trait_function inverse_retract!(
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod,
)
function inverse_retract!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M, typeof(p)),
)
    return inverse_retract!(get_embedding(M, p), X, p, q, m)
end

@trait_function isapprox(M::AbstractDecoratorManifold, p, q; kwargs...)
@trait_function isapprox(M::AbstractDecoratorManifold, p, X, Y; kwargs...)

@trait_function is_flat(M::AbstractDecoratorManifold)

# Introduce Deco Trait | automatic forward | fallback
@trait_function is_point(M::AbstractDecoratorManifold, p; kwargs...)

# Embedded
function is_point_embedding(
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
        pt = is_point(get_embedding(M, p), embed(M, p); error = error, kwargs...)
        !pt && return false # no error thrown (deactivated) but returned false -> return false
    catch e
        if e isa DomainError || e isa AbstractManifoldDomainError
            e = ManifoldDomainError(
                "$p is not a point on $M because it is not a valid point in its embedding.",
                e,
            )
        end
        throw(e) #an error occurred that we do not handle ourselves -> rethrow.
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

# Introduce Deco Trait | automatic forward | fallback
@trait_function is_vector(M::AbstractDecoratorManifold, p, X, cbp::Bool = true; kwargs...)
# EmbeddedManifold
# I am not yet sure how to properly document this embedding behaviour here in a docstring.
function is_vector_embedding(
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
        tv = is_vector(
            get_embedding(M, p),
            embed(M, p),
            embed(M, p, X),
            check_base_point;
            error = error,
            kwargs...,
        )
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

@trait_function norm(M::AbstractDecoratorManifold, p, X)
function norm(::TraitList{IsIsometricEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    return norm(get_embedding(M, p), p, X)
end


function log(M::AbstractDecoratorManifold, p, q)
    return _log_forwarding(get_forwarding_type(M, log, p), M, p, q)
end

function _log_forwarding(::StopForwardingType, M::AbstractDecoratorManifold, p, q)
    return invoke(log, Tuple{AbstractManifold,typeof(p),typeof(q)}, M, p, q)
end
function _log_forwarding(::EmbeddedForwardingType, M::AbstractDecoratorManifold, p, q)
    return log(get_embedding(M, p), p, q)
end
function _log_forwarding(::EmbeddedSimpleForwardingType, M::AbstractDecoratorManifold, p, q)
    return log(get_embedding(M, p), p, q)
end
function _log_forwarding(::SimpleForwardingType, M::AbstractDecoratorManifold, p, q)
    return log(decorated_manifold(M), p, q)
end

function log!(M::AbstractDecoratorManifold, X, p, q)
    return _log_forwarding!(get_forwarding_type(M, log, p), M, X, p, q)
end

function _log_forwarding!(::StopForwardingType, M::AbstractDecoratorManifold, X, p, q)
    return invoke(log!, Tuple{AbstractManifold,typeof(X),typeof(p),typeof(q)}, M, X, p, q)
end
function _log_forwarding!(::EmbeddedForwardingType, M::AbstractDecoratorManifold, X, p, q)
    return log!(get_embedding(M, p), X, p, q)
end
function _log_forwarding!(
    ::EmbeddedSimpleForwardingType,
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
)
    return log!(get_embedding(M, p), X, p, q)
end
function _log_forwarding!(::SimpleForwardingType, M::AbstractDecoratorManifold, X, p, q)
    return log!(decorated_manifold(M), X, p, q)
end



# Introduce Deco Trait | automatic forward | fallback
@trait_function parallel_transport_direction(M::AbstractDecoratorManifold, p, X, q)
# EmbeddedSubManifold
function parallel_transport_direction(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
)
    return parallel_transport_direction(get_embedding(M, p), p, X, q)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function parallel_transport_direction!(M::AbstractDecoratorManifold, Y, p, X, q)
# EmbeddedSubManifold
function parallel_transport_direction!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
)
    return parallel_transport_direction!(get_embedding(M, p), Y, p, X, q)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function parallel_transport_to(M::AbstractDecoratorManifold, p, X, q)
# EmbeddedSubManifold
function parallel_transport_to(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
)
    return parallel_transport_to(get_embedding(M, p), p, X, q)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function parallel_transport_to!(M::AbstractDecoratorManifold, Y, p, X, q)
# EmbeddedSubManifold
function parallel_transport_to!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
)
    return parallel_transport_to!(get_embedding(M, p), Y, p, X, q)
end

# Introduce Deco Trait | automatic forward | fallback
@trait_function project(M::AbstractDecoratorManifold, p)

# Introduce Deco Trait | automatic forward | fallback
@trait_function project!(M::AbstractDecoratorManifold, q, p)

# Introduce Deco Trait | automatic forward | fallback
@trait_function project(M::AbstractDecoratorManifold, p, X)

# Introduce Deco Trait | automatic forward | fallback
@trait_function project!(M::AbstractDecoratorManifold, Y, p, X)

@trait_function Random.rand(M::AbstractDecoratorManifold; kwargs...)

@trait_function Random.rand!(M::AbstractDecoratorManifold, p; kwargs...)

@trait_function Random.rand(rng::AbstractRNG, M::AbstractDecoratorManifold; kwargs...) :() 2

@trait_function Random.rand!(rng::AbstractRNG, M::AbstractDecoratorManifold, p; kwargs...) :() 2

# Introduce Deco Trait | automatic forward | fallback
@trait_function representation_size(M::AbstractDecoratorManifold) (no_empty,)
# Isometric Embedded submanifold
function representation_size_embedding(M::AbstractDecoratorManifold)
    return representation_size(get_embedding(M))
end
function representation_size(::EmptyTrait, M::AbstractDecoratorManifold)
    return representation_size(decorated_manifold(M))
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
function retract(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return retract(get_embedding(M, p), p, X, m)
end
function retract_fused(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return retract_fused(get_embedding(M, p), p, X, t, m)
end

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
function retract!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return retract!(get_embedding(M, p), q, p, X, m)
end
function retract_fused!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
)
    return retract_fused!(get_embedding(M, p), q, p, X, t, m)
end

@trait_function vector_transport_direction(
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
function vector_transport_direction(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
    return vector_transport_direction(get_embedding(M, p), p, X, d, m)
end

@trait_function vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
function vector_transport_direction!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
    return vector_transport_direction!(get_embedding(M, p), Y, p, X, d, m)
end

@trait_function vector_transport_to(
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
function vector_transport_to(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
    return vector_transport_to(get_embedding(M, p), p, X, q, m)
end

@trait_function vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
function vector_transport_to!(
    ::TraitList{IsEmbeddedSubmanifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
)
    return vector_transport_to!(get_embedding(M, p), Y, p, X, q, m)
end

@trait_function Weingarten(M::AbstractDecoratorManifold, p, X, V)
@trait_function Weingarten!(M::AbstractDecoratorManifold, Y, p, X, V)

@trait_function zero_vector(M::AbstractDecoratorManifold, p)
function zero_vector_embedding(M::AbstractDecoratorManifold, p)
    return zero_vector(get_embedding(M, p), p)
end

@trait_function zero_vector!(M::AbstractDecoratorManifold, X, p)
function zero_vector_embedding!(M::AbstractDecoratorManifold, X, p)
    return zero_vector!(get_embedding(M, p), X, p)
end

# Trait recursion breaking
# An unfortunate consequence of Julia's method recursion limitations
# Add more traits and functions as needed

for trait_type in [TraitList{IsEmbeddedSubmanifold}]
    @eval begin
        @next_trait_function $trait_type isapprox(
            M::AbstractDecoratorManifold,
            p,
            q;
            kwargs...,
        )
        @next_trait_function $trait_type isapprox(
            M::AbstractDecoratorManifold,
            p,
            X,
            Y;
            kwargs...,
        )
    end
end



const metric_functions = [
    change_metric,
    change_representer,
    exp,
    exp_fused,
    get_basis,
    get_coordinates,
    get_vector,
    get_vectors,
    inner,
    inverse_retract,
    log,
    mid_point,
    norm,
    parallel_transport_direction,
    parallel_transport_to,
    retract,
    retract_fused,
    riemann_tensor,
    vector_transport_direction,
    vector_transport_to,
    Weingarten,
]


function get_forwarding_type(M::AbstractDecoratorManifold, f)
    return get_forwarding_type_embedding(get_embedding_type(M), M, f)
end
function get_forwarding_type(M::AbstractDecoratorManifold, f, p)
    return get_forwarding_type_embedding(get_embedding_type(M, p), M, f)
end

function get_forwarding_type_embedding(
    ::Union{EmbeddedManifoldType,IsometricallyEmbeddedManifoldType,NotEmbeddedManifoldType},
    M::AbstractDecoratorManifold,
    f,
)
    return StopForwardingType()
end

for mf in metric_functions
    @eval begin
        function get_forwarding_type_embedding(
            ::EmbeddedSubmanifoldType,
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedSimpleForwardingType()
        end
        function get_forwarding_type_embedding(
            ::IsometricallyEmbeddedManifoldType,
            M::AbstractDecoratorManifold,
            ::typeof($mf),
        )
            return EmbeddedSimpleForwardingType()
        end
    end
end
