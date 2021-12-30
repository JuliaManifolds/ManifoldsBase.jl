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

"""
    IsEmbeddedManifold <: AbstractTrait

A trait to declare an [`AbstractManifold`](@ref) as an embedded manifold.
"""
struct IsEmbeddedManifold <: AbstractTrait end

"""
    IsIsometricManifoldEmbeddedManifold <: AbstractTrait

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an isometrically embedded manifold.
To activate this for your manifold, set `is_isometric_embedded_manifold` for your manifold type to true.

Here, for example [`inner`](@ref) and [`norm`](@ref) are passed to the embedding
"""
struct IsIsometricEmbeddedManifold <: AbstractTrait end

parent_trait(::IsIsometricEmbeddedManifold) = IsEmbeddedManifold()

"""
    IsEmbeddedSubmanifold{M}
    is_embedded_submanifold(M::Type{<:AbstractDecoratorManifold})

Determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
It is a special case of an [`IsIsometricEmbeddedManifold`](@ref).

Here, additionally, all retraction, inverse retractions and vectors transports, especially
[`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref) are passed to the embedding.
"""
struct IsEmbeddedSubmanifoldManifold <: AbstractTrait end

parent_trait(::IsEmbeddedSubmanifoldManifold) = IsIsometricEmbeddedManifold()


#
# Generic Decorator functions
decorated_manifold(M::AbstractDecoratorManifold) = M

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
# Embedded specifix functions.
"""
    get_embedding(M::AbstractDecoratorManifold)

Specify the embedding of a manifold that has abstract decorators.
"""
get_embedding(M::AbstractDecoratorManifold)

#
# -----------------------------------------------------------------------------------------
# This is one new function

# INtroduction and default fallbacks could become a macro?
# Introduce trait
function allocate_result(M::AbstractDecoratorManifold, f, x...)
    return allocate_result(trait(M, f, x...), M, f, x...)
end
# disambiguation
@invoke_maker 1 AbstractManifold allocate_result(
    M::AbstractDecoratorManifold,
    f::typeof(get_coordinates),
    p,
    X,
    B::AbstractBasis,
)

# Introduce fallback
@inline function allocate_result(::EmptyTrait, M::AbstractManifold, f, x...)
    return invoke(
        allocate_result,
        Tuple{AbstractManifold,typeof(f),typeof(x).parameters...},
        M,
        f,
        x...,
    )
end
# Introduce automatic forward
@inline function allocate_result(t::NestedTrait, M::AbstractManifold, f, x...)
    return allocate_result(t.tail, M, f, x...)
end
function allocate_result(
    ::NestedTrait{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    f::typeof(embed),
    x...,
)
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(get_embedding(M)))
end
function allocate_result(
    ::NestedTrait{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    f::typeof(project),
    x...,
)
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(M))
end


# Introduce Deco Trait | automatic foward | fallback
function check_size(M::AbstractDecoratorManifold, p)
    return check_size(trait(M, p), M, p)
end
function check_size(t::NestedTrait, M::AbstractManifold, p)
    return check_size(t.tail, M, p)
end
function check_size(::EmptyTrait, M::AbstractManifold, p)
    return invoke(check_size, Tuple{AbstractManifold,typeof(p)}, M, p)
end
# Embedded
function check_size(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p)
    return check_size(get_embedding(M), p)
end

# Introduce Deco Trait | automatic foward | fallback
function check_size(M::AbstractDecoratorManifold, p, X)
    return check_size(trait(M, p), M, p, X)
end
function check_size(t::NestedTrait, M::AbstractManifold, p, X)
    return check_size(t.tail, M, p, X)
end
function check_size(::EmptyTrait, M::AbstractManifold, p, X)
    return invoke(check_size, Tuple{AbstractManifold,typeof(p)}, M, p, X)
end
# Embedded
function check_size(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    return check_size(get_embedding(M), p, X)
end

# Introduce Deco Trait | automatic foward | fallback
function embed(M::AbstractDecoratorManifold, p)
    return embed(trait(M, p), M, p)
end
embed(t::NestedTrait, M::AbstractManifold, p) = embed(t.tail, M, p)
function embed(::EmptyTrait, M::AbstractManifold, p)
    return invoke(embed, Tuple{AbstractManifold,typeof(p)}, M, p)
end
# EmbeddedManifold
function embed(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p)
    q = allocate_result(M, embed, p)
    return embed!(M, q, p)
end

# Introduce Deco Trait | automatic foward | fallback
function embed!(M::AbstractDecoratorManifold, q, p)
    return embed!(trait(M, q, p), M, q, p)
end
embed!(t::NestedTrait, M::AbstractManifold, q, p) = embed!(t.tail, M, q, p)
function embed!(::EmptyTrait, M::AbstractManifold, q, p)
    return invoke(embed!, Tuple{AbstractManifold,typeof(q),typeof(p)}, M, q, p)
end
# EmbeddedManifold
function embed!(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, q, p)
    return copyto!(M, q, p)
end

# Introduce Deco Trait | automatic foward | fallback
function embed(M::AbstractDecoratorManifold, p, X)
    return embed(trait(M, p, X), M, p, X)
end
embed(t::NestedTrait, M::AbstractManifold, p, X) = embed(t.tail, M, p, X)
function embed(::EmptyTrait, M::AbstractManifold, p, X)
    return invoke(embed, Tuple{AbstractManifold,typeof(p),typeof(X)}, M, p, X)
end
# EmbeddedManifold
function embed(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    q = allocate_result(M, embed, p, X)
    return embed!(M, q, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
function embed!(M::AbstractDecoratorManifold, Y, p, X)
    return embed!(trait(M, Y, p, X), M, Y, p, X)
end
embed!(t::NestedTrait, M::AbstractManifold, Y, p, X) = embed!(t.tail, M, Y, p, X)
function embed!(::EmptyTrait, M::AbstractManifold, Y, p, X)
    return invoke(embed!, Tuple{AbstractManifold,typeof(Y),typeof(p),typeof(X)}, M, Y, p, X)
end
# EmbeddedManifold
function embed!(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, Y, p, X)
    return copyto!(M, Y, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
function exp(M::AbstractDecoratorManifold, p, X)
    return exp(trait(M, p, X), M, p, X)
end
exp(t::NestedTrait, M::AbstractManifold, p, X) = exp(t.tail, M, p, X)
function exp(::EmptyTrait, M::AbstractManifold, p, X)
    return invoke(exp, Tuple{AbstractManifold,typeof(p),typeof(X)}, M, p, X)
end
# EmbeddedSubManifold
function exp(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    return exp(get_embedding(M), p, X)
end

# Introduce Deco Trait | automatic foward | fallback
function exp!(M::AbstractDecoratorManifold, q, p, X)
    return exp!(trait(M, q, p, X), M, q, p, X)
end
exp!(t::NestedTrait, M::AbstractManifold, q, p, X) = exp!(t.tail, M, q, p, X)
function exp!(::EmptyTrait, M::AbstractManifold, q, p, X)
    return invoke(exp!, Tuple{AbstractManifold,typeof(q),typeof(p),typeof(X)}, M, q, p, X)
end
# EmbeddedSubManifold
function exp!(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, q, p, X)
    return exp!(get_embedding(M), q, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
function inner(M::AbstractDecoratorManifold, p, X, Y)
    return inner(trait(M, p, X, Y), M, p, X, Y)
end
function inner(t::NestedTrait, M::AbstractManifold, p, X, Y)
    return inner(t.tail, M, p, X, Y)
end
function inner(::EmptyTrait, M::AbstractManifold, p, X, Y)
    return invoke(inner, Tuple{AbstractManifold,typeof(p),typeof(X),typeof(Y)}, M, p, X, Y)
end
# Isometric Embedded submanifold
function inner(
    ::NestedTrait{IsIsometricEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    Y,
)
    return inner(get_embedding(M), p, X, Y)
end

# Introduce Deco Trait | automatic foward | fallback
function inverse_retract(M::AbstractDecoratorManifold, p, q)
    return inverse_retract(trait(M, q, X), M, p, q)
end
function inverse_retract(
    t::NestedTrait,
    M::AbstractManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return inverse_retract(t.tail, M, p, q)
end
function inverse_retract(
    ::EmptyTrait,
    M::AbstractManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return invoke(
        inverse_retract,
        Tuple{AbstractManifold,typeof(p),typeof(q),typeof(m)},
        M,
        p,
        q,
        m,
    )
end
# Transparent for Submanifolds
function inverse_retract(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return inverse_retract(get_embedding(M), p, q, m)
end

# Introduce Deco Trait | automatic foward | fallback
function inverse_retract!(M::AbstractDecoratorManifold, X, p, q)
    return inverse_retract!(trait(M, X, p, q), M, p, p, q)
end
function inverse_retract!(
    t::NestedTrait,
    M::AbstractManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return exp!(t.tail, M, X, p, q, m)
end
function inverse_retract!(
    ::EmptyTrait,
    M::AbstractManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return invoke(
        inverse_retract!,
        Tuple{AbstractManifold,typeof(X),typeof(p),typeof(q),typeof(m)},
        M,
        X,
        p,
        q,
        m,
    )
end
function inverse_retract!(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return inverse_retract!(get_embedding(M), X, p, q, m)
end

function is_point(M::AbstractDecoratorManifold, p, te = false; kwargs...)
    return is_point(trait(M, p), M, p, te; kwargs...)
end
function is_point(t::NestedTrait, M::AbstractManifold, p, te = false; kwargs...)
    return is_point(t.tail, M, p, te; kwargs...)
end
function is_point(::EmptyTrait, M::AbstractManifold, p, te = false; kwargs...)
    return invoke(
        is_point,
        Tuple{AbstractManifold,typeof(p),typeof(te)},
        M,
        p,
        te;
        kwargs...,
    )
end
# EmbeddedManifold
function is_point(
    ::NestedTrait{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    throw_error = false;
    kwargs...,
)
    return is_point(get_embedding(M), p, throw_error; kwargs...)
end

function is_vector(M::AbstractDecoratorManifold, p, X, te = false; kwargs...)
    return is_vector(trait(M, p, X), M, p, X, te; kwargs...)
end
function is_vector(t::NestedTrait, M::AbstractManifold, p, X, te = false; kwargs...)
    return is_vector(t.tail, M, p, X, te; kwargs...)
end
function is_vector(::EmptyTrait, M::AbstractManifold, p, X, te = false; kwargs...)
    return invoke(
        is_vector,
        Tuple{AbstractManifold,typeof(p),typeof(X),typeof(te)},
        M,
        p,
        X,
        te;
        kwargs...,
    )
end
# EmbeddedManifold
function is_vector(
    ::NestedTrait{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    throw_error = false;
    kwargs...,
)
    return is_vector(get_embedding(M), p, X, throw_error; kwargs...)
end

function norm(M::AbstractDecoratorManifold, p, X)
    return norm(trait(M, p, X), M, p, X)
end
norm(t::NestedTrait, M::AbstractManifold, p, X) = norm(t.tail, M, p, X)
function norm(::EmptyTrait, M::AbstractManifold, p, X)
    return invoke(norm, Tuple{AbstractManifold,typeof(p),typeof(X)}, M, p, X)
end
function norm(
    ::NestedTrait{IsIsometricEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
)
    return norm(get_embedding(M), p, X)
end

function log(M::AbstractDecoratorManifold, p, q)
    return log(trait(M, p, q), M, p, q)
end
log(t::NestedTrait, M::AbstractManifold, p, q) = log(t.tail, M, p, q)
function log(::EmptyTrait, M::AbstractManifold, p, q)
    return invoke(log, Tuple{AbstractManifold,typeof(p),typeof(q)}, M, p, q)
end
function log(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    q,
)
    return log(get_embedding(M), p, q)
end

# Introduce Deco Trait | automatic foward | fallback
function log!(M::AbstractDecoratorManifold, X, p, q)
    return log!(trait(M, X, p, q), M, X, p, q)
end
log!(t::NestedTrait, M::AbstractManifold, q, p, X) = log!(t.tail, M, q, p, X)
function log!(::EmptyTrait, M::AbstractManifold, q, p, X)
    return invoke(log!, Tuple{AbstractManifold,typeof(X),typeof(p),typeof(q)}, M, X, p, q)
end
function log!(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
)
    return log!(get_embedding(M), X, p, q)
end

# Introduce Deco Trait | automatic foward | fallback
function project(M::AbstractDecoratorManifold, p)
    return project(trait(M, p), M, p)
end
project(t::NestedTrait, M::AbstractManifold, p) = project(t.tail, M, p)
function project(::EmptyTrait, M::AbstractManifold, p)
    return invoke(project, Tuple{AbstractManifold,typeof(p)}, M, p)
end
# EmbeddedManifold
function project(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p)
    q = allocate_result(M, project, p)
    return project!(M, q, p)
end

# Introduce Deco Trait | automatic foward | fallback
function project!(M::AbstractDecoratorManifold, q, p)
    return project!(trait(M, q, p), M, q, p)
end
project!(t::NestedTrait, M::AbstractManifold, q, p) = project!(t.tail, M, q, p)
function project!(::EmptyTrait, M::AbstractManifold, q, p)
    return invoke(project!, Tuple{AbstractManifold,typeof(q),typeof(p)}, M, q, p)
end
# EmbeddedManifold
function project!(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, q, p)
    return copyto!(M, q, p)
end

# Introduce Deco Trait | automatic foward | fallback
function project(M::AbstractDecoratorManifold, p, X)
    return project(trait(M, p, X), M, p, X)
end
project(t::NestedTrait, M::AbstractManifold, p, X) = project(t.tail, M, p, X)
function project(::EmptyTrait, M::AbstractManifold, p, X)
    return invoke(project, Tuple{AbstractManifold,typeof(p),typeof(X)}, M, p, X)
end
# EmbeddedManifold
function project(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    q = allocate_result(M, project, p, X)
    return project!(M, q, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
function project!(M::AbstractDecoratorManifold, Y, p, X)
    return project!(trait(M, Y, p, X), M, Y, p, X)
end
project!(t::NestedTrait, M::AbstractManifold, Y, p, X) = project!(t.tail, M, Y, p, X)
function project!(::EmptyTrait, M::AbstractManifold, Y, p, X)
    return invoke(
        project!,
        Tuple{AbstractManifold,typeof(Y),typeof(p),typeof(X)},
        M,
        Y,
        p,
        X,
    )
end
# EmbeddedManifold
function project!(::NestedTrait{IsEmbeddedManifold}, M::AbstractDecoratorManifold, Y, p, X)
    return copyto!(M, Y, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
function retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractVectorTransportMethod = default_retraction_method(M),
)
    return retract(trait(M, p, X, m), M, p, X, m)
end
function retract(
    t::NestedTrait,
    M::AbstractManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract(t.tail, M, p, X, m)
end
function retract(
    ::EmptyTrait,
    M::AbstractManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return invoke(
        retract,
        Tuple{AbstractManifold,typeof(p),typeof(X),typeof(m)},
        M,
        p,
        X,
        m,
    )
end
function retract(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract(get_embedding(M), p, X, m)
end

function retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract!(trait(M, q, p, X), M, q, p, X, m)
end
function retract!(
    t::NestedTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract!(t.tail, M, q, p, X, m)
end
function retract!(
    ::EmptyTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return invoke(
        retract!,
        Tuple{AbstractManifold,typeof(q),typeof(p),typeof(X),typeof(m)},
        M,
        q,
        p,
        X,
        m,
    )
end
function retract!(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract!(get_embedding(M), q, p, X, m)
end

function vector_transport_along(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along(trait(M, q, p, X), M, q, p, X, m)
end
function vector_transport_along(
    t::NestedTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along(t.tail, M, q, p, X, m)
end
function vector_transport_along(
    ::EmptyTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return invoke(
        vector_transport_along,
        Tuple{AbstractManifold,typeof(q),typeof(p),typeof(X),typeof(m)},
        M,
        q,
        p,
        X,
        m,
    )
end
function vector_transport_along(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along(get_embedding(M), p, X, c, m)
end

function vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along!(trait(M, Y, q, p, X), M, Y, q, p, X, m)
end
function vector_transport_along!(
    t::NestedTrait,
    M::AbstractManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along!(t.tail, M, Y, q, p, X, m)
end
function vector_transport_along!(
    ::EmptyTrait,
    M::AbstractManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return invoke(
        vector_transport_along!,
        Tuple{AbstractManifold,typeof(Y),typeof(q),typeof(p),typeof(X),typeof(m)},
        M,
        Y,
        q,
        p,
        X,
        m,
    )
end
function vector_transport_along!(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along!(get_embedding(M), Y, p, X, c, m)
end

@invoke_maker 1 AbstractManifold vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    method::AbstractVectorTransportMethod,
)

function vector_transport_direction(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction(trait(M, q, p, X), M, q, p, X, m)
end
function vector_transport_direction(
    t::NestedTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction(t.tail, M, q, p, X, m)
end
function vector_transport_direction(
    ::EmptyTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return invoke(
        vector_transport_direction,
        Tuple{AbstractManifold,typeof(q),typeof(p),typeof(X),typeof(m)},
        M,
        q,
        p,
        X,
        m,
    )
end
function vector_transport_direction(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction(get_embedding(M), p, X, d, m)
end

function vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction!(trait(M, Y, q, p, X), M, Y, q, p, X, m)
end
function vector_transport_direction!(
    t::NestedTrait,
    M::AbstractManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction!(t.tail, M, Y, q, p, X, m)
end
function vector_transport_direction!(
    ::EmptyTrait,
    M::AbstractManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return invoke(
        vector_transport_direction!,
        Tuple{AbstractManifold,typeof(Y),typeof(q),typeof(p),typeof(X),typeof(m)},
        M,
        Y,
        q,
        p,
        X,
        m,
    )
end
function vector_transport_direction!(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction!(get_embedding(M), Y, p, X, d, m)
end

function vector_transport_to(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to(trait(M, q, p, X), M, q, p, X, m)
end
function vector_transport_to(
    t::NestedTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to(t.tail, M, q, p, X, m)
end
function vector_transport_to(
    ::EmptyTrait,
    M::AbstractManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return invoke(
        vector_transport_to,
        Tuple{AbstractManifold,typeof(q),typeof(p),typeof(X),typeof(m)},
        M,
        q,
        p,
        X,
        m,
    )
end
function vector_transport_to(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to(get_embedding(M), p, X, q, m)
end

function vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to!(trait(M, Y, q, p, X), M, Y, q, p, X, m)
end
function vector_transport_to!(
    t::NestedTrait,
    M::AbstractManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to!(t.tail, M, Y, q, p, X, m)
end
function vector_transport_to!(
    ::EmptyTrait,
    M::AbstractManifold,
    Y,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return invoke(
        vector_transport_to!,
        Tuple{AbstractManifold,typeof(Y),typeof(q),typeof(p),typeof(X),typeof(m)},
        M,
        Y,
        q,
        p,
        X,
        m,
    )
end
function vector_transport_to!(
    ::NestedTrait{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to!(get_embedding(M), Y, p, X, q, m)
end
