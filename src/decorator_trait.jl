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

A Trait to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is
an isometrically embedded manifold.
It is a special case of the [`IsEmbeddedManifold`](@ref) trait, i.e. it has all properties of this trait.

Here, additionally, netric related functions like [`inner`](@ref) and [`norm`](@ref) are passed to the embedding
"""
struct IsIsometricEmbeddedManifold <: AbstractTrait end

parent_trait(::IsIsometricEmbeddedManifold) = IsEmbeddedManifold()

"""
    IsEmbeddedSubmanifold <: AbstractTrait

A trait to determine whether an [`AbstractDecoratorManifold`](@ref) `M` is an embedded submanifold.
It is a special case of the [`IsIsometricEmbeddedManifold`](@ref) trait, i.e. it has all properties of
this trait.

In this trait, additionally to the isometric embedded manifold, all retractions, inverse retractions,
and vectors transports, especially [`exp`](@ref), [`log`](@ref), and [`parallel_transport_to`](@ref)
are passed to the embedding.
"""
struct IsEmbeddedSubmanifoldManifold <: AbstractTrait end

parent_trait(::IsEmbeddedSubmanifoldManifold) = IsIsometricEmbeddedManifold()


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
@inline function allocate_result(t::TraitList, M::AbstractManifold, f, x...)
    return allocate_result(next_trait(t), M, f, x...)
end
function allocate_result(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    f::typeof(embed),
    x...,
)
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(get_embedding(M)))
end
function allocate_result(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    f::typeof(project),
    x...,
)
    T = allocate_result_type(get_embedding(M), f, x)
    return allocate(x[1], T, representation_size(M))
end


# Introduce Deco Trait | automatic foward | fallback
@trait_function check_size(M::AbstractDecoratorManifold, p)
# Embedded
function check_size(::TraitList{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p)
    return check_size(get_embedding(M), p)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function check_size(M::AbstractDecoratorManifold, p, X)
# Embedded
function check_size(::TraitList{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    return check_size(get_embedding(M), p, X)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function embed(M::AbstractDecoratorManifold, p)
# EmbeddedManifold
function embed(::TraitList{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p)
    q = allocate_result(M, embed, p)
    return embed!(M, q, p)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function embed!(M::AbstractDecoratorManifold, q, p)
# EmbeddedManifold
function embed!(::TraitList{IsEmbeddedManifold}, M::AbstractDecoratorManifold, q, p)
    return copyto!(M, q, p)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function embed(M::AbstractDecoratorManifold, p, X)
# EmbeddedManifold
function embed(::TraitList{IsEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    q = allocate_result(M, embed, p, X)
    return embed!(M, q, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function embed!(M::AbstractDecoratorManifold, Y, p, X)
# EmbeddedManifold
function embed!(::TraitList{IsEmbeddedManifold}, M::AbstractDecoratorManifold, Y, p, X)
    return copyto!(M, Y, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function exp(M::AbstractDecoratorManifold, p, X)
# EmbeddedSubManifold
function exp(::TraitList{IsEmbeddedSubmanifoldManifold}, M::AbstractDecoratorManifold, p, X)
    return exp(get_embedding(M), p, X)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function exp!(M::AbstractDecoratorManifold, q, p, X)
# EmbeddedSubManifold
function exp!(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
)
    return exp!(get_embedding(M), q, p, X)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function inner(M::AbstractDecoratorManifold, p, X, Y)
# Isometric Embedded submanifold
function inner(
    ::TraitList{IsIsometricEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    Y,
)
    return inner(get_embedding(M), p, X, Y)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function inverse_retract(
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
# Transparent for Submanifolds
function inverse_retract(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return inverse_retract(get_embedding(M), p, q, m)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function inverse_retract!(M::AbstractDecoratorManifold, X, p, q)
function inverse_retract!(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
)
    return inverse_retract!(get_embedding(M), X, p, q, m)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function is_point(M::AbstractDecoratorManifold, p, te = false; kwargs...)
# Embedded
function is_point(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    te = false;
    kwargs...,
)
    ep = is_point(get_embedding(M), embed(M, p), te; kwargs...)
    (!ep && !te) && return false # no point in E? above throws error - or we false here
    mpe = check_point(M, p; kwargs...)
    mpe === nothing && return true
    te && throw(mpe)
    return false
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function is_vector(
    M::AbstractDecoratorManifold,
    p,
    X,
    te = false,
    cbp = true;
    kwargs...,
)
# EmbeddedManifold
# I am not yet sure how to properly document this embedding behaviour here in a docstring.
function is_vector(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    te = false,
    cbp = true;
    kwargs...,
)
    ev = is_vector(get_embedding(M), embed(M, p), embed(M, p, X), te, cbp; kwargs...)
    (!ev && !te) && return false # if te, the line before throws an error, otherwuse we end with false early here
    if cbp
        # if we are here p is a valid point in embedding from the first code line
        mpe = check_point(M, p; kwargs...)
        if mpe !== nothing
            te && throw(mpe)
            return false
        end
    end
    mtve = check_vector(M, p, X; kwargs...)
    mtve === nothing && return true
    te && throw(mtve)
    return false
end

@trait_function norm(M::AbstractDecoratorManifold, p, X)
function norm(::TraitList{IsIsometricEmbeddedManifold}, M::AbstractDecoratorManifold, p, X)
    return norm(get_embedding(M), p, X)
end

@trait_function log(M::AbstractDecoratorManifold, p, q)
function log(::TraitList{IsEmbeddedSubmanifoldManifold}, M::AbstractDecoratorManifold, p, q)
    return log(get_embedding(M), p, q)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function log!(M::AbstractDecoratorManifold, X, p, q)
function log!(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
)
    return log!(get_embedding(M), X, p, q)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function parallel_transport_along(M::AbstractDecoratorManifold, p, X, c)
# EmbeddedSubManifold
function parallel_transport_along(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    c,
)
    return parallel_transport_along(get_embedding(M), p, X, c)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function parallel_transport_along!(M::AbstractDecoratorManifold, Y, p, X, c)
# EmbeddedSubManifold
function parallel_transport_along!(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c,
)
    return parallel_transport_along!(get_embedding(M), Y, p, X, c)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function parallel_transport_direction(M::AbstractDecoratorManifold, p, X, q)
# EmbeddedSubManifold
function parallel_transport_direction(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
)
    return parallel_transport_direction(get_embedding(M), p, X, q)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function parallel_transport_direction!(M::AbstractDecoratorManifold, Y, p, X, q)
# EmbeddedSubManifold
function parallel_transport_direction!(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
)
    return parallel_transport_direction!(get_embedding(M), Y, p, X, q)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function parallel_transport_to(M::AbstractDecoratorManifold, p, X, q)
# EmbeddedSubManifold
function parallel_transport_to(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
)
    return parallel_transport_to(get_embedding(M), p, X, q)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function parallel_transport_to!(M::AbstractDecoratorManifold, Y, p, X, q)
# EmbeddedSubManifold
function parallel_transport_to!(
    ::TraitList{IsEmbeddedManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
)
    return parallel_transport_to!(get_embedding(M), Y, p, X, q)
end

# Introduce Deco Trait | automatic foward | fallback
@trait_function project(M::AbstractDecoratorManifold, p)

# Introduce Deco Trait | automatic foward | fallback
@trait_function project!(M::AbstractDecoratorManifold, q, p)

# Introduce Deco Trait | automatic foward | fallback
@trait_function project(M::AbstractDecoratorManifold, p, X)

# Introduce Deco Trait | automatic foward | fallback
@trait_function project!(M::AbstractDecoratorManifold, Y, p, X)

# Introduce Deco Trait | automatic foward | fallback
@trait_function retract(
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
function retract(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract(get_embedding(M), p, X, m)
end

@trait_function retract!(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
function retract!(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = default_retraction_method(M),
)
    return retract!(get_embedding(M), q, p, X, m)
end

@trait_function vector_transport_along(
    M::AbstractDecoratorManifold,
    q,
    p,
    X,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
function vector_transport_along(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along(get_embedding(M), p, X, c, m)
end

@trait_function vector_transport_along!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
function vector_transport_along!(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_along!(get_embedding(M), Y, p, X, c, m)
end

@trait_function vector_transport_direction(
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
function vector_transport_direction(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction(get_embedding(M), p, X, d, m)
end

@trait_function vector_transport_direction!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
function vector_transport_direction!(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_direction!(get_embedding(M), Y, p, X, d, m)
end

@trait_function vector_transport_to(
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
function vector_transport_to(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to(get_embedding(M), p, X, q, m)
end

@trait_function vector_transport_to!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
function vector_transport_to!(
    ::TraitList{IsEmbeddedSubmanifoldManifold},
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return vector_transport_to!(get_embedding(M), Y, p, X, q, m)
end
