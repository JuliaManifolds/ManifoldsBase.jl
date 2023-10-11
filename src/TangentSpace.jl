
@doc raw"""
    TangentSpace{𝔽,M} = Fiber{𝔽,TangentSpaceType,M} where {𝔽,M<:AbstractManifold{𝔽}}

A manifold for the tangent space ``T_p\mathcal M`` at a point ``p\in\mathcal M``.
This is modelled as an alias for [`VectorSpaceFiber`](@ref) corresponding to
[`TangentSpaceType`](@ref).

# Constructor

    TangentSpace(M::AbstractManifold, p)

Return the manifold (vector space) representing the tangent space ``T_p\mathcal M``
at point `p`, ``p\in\mathcal M``.
"""
const TangentSpace{𝔽,M} = Fiber{𝔽,TangentSpaceType,M} where {𝔽,M<:AbstractManifold{𝔽}}

TangentSpace(M::AbstractManifold, p) = Fiber(M, p, TangentSpaceType())

@doc raw"""
    CotangentSpace{𝔽,M} = Fiber{𝔽,CotangentSpaceType,M} where {𝔽,M<:AbstractManifold{𝔽}}

A manifold for the Cotangent space ``T^*_p\mathcal M`` at a point ``p\in\mathcal M``.
This is modelled as an alias for [`VectorSpaceFiber`](@ref) corresponding to
[`CotangentSpaceType`](@ref).

# Constructor

    CotangentSpace(M::AbstractManifold, p)

Return the manifold (vector space) representing the cotangent space ``T^*_p\mathcal M``
at point `p`, ``p\in\mathcal M``.
"""
const CotangentSpace{𝔽,M} = Fiber{𝔽,CotangentSpaceType,M} where {𝔽,M<:AbstractManifold{𝔽}}

CotangentSpace(M::AbstractManifold, p) = Fiber(M, p, CotangentSpaceType())


function allocate_result(M::TangentSpace, ::typeof(rand))
    return zero_vector(M.manifold, M.point)
end

"""
    distance(M::TangentSpace, p, q)

Distance between vectors `p` and `q` from the vector space `M`. It is calculated as the norm
of their difference.
"""
function distance(M::TangentSpace, p, q)
    return norm(M.manifold, M.point, q - p)
end

function embed!(M::TangentSpace, q, p)
    return embed!(M.manifold, q, M.point, p)
end
function embed!(M::TangentSpace, Y, p, X)
    return embed!(M.manifold, Y, M.point, X)
end

@doc raw"""
    exp(M::TangentSpace, p, X)

Exponential map of tangent vectors `X` and `p` from the tangent space `M`. It is
calculated as their sum.
"""
exp(::TangentSpace, ::Any, ::Any)

function exp!(M::TangentSpace, q, p, X)
    copyto!(M.manifold, q, p + X)
    return q
end

fiber_dimension(M::AbstractManifold, ::CotangentSpaceType) = manifold_dimension(M)
fiber_dimension(M::AbstractManifold, ::TangentSpaceType) = manifold_dimension(M)

function get_basis(M::TangentSpace, p, B::CachedBasis)
    return invoke(
        get_basis,
        Tuple{AbstractManifold,Any,CachedBasis},
        M.manifold,
        M.point,
        B,
    )
end
function get_basis(M::TangentSpace, p, B::AbstractBasis{<:Any,TangentSpaceType})
    return get_basis(M.manifold, M.point, B)
end

function get_coordinates(M::TangentSpace, p, X, B::AbstractBasis)
    return get_coordinates(M.manifold, M.point, X, B)
end

function get_coordinates!(M::TangentSpace, Y, p, X, B::AbstractBasis)
    return get_coordinates!(M.manifold, Y, M.point, X, B)
end

function get_vector(M::TangentSpace, p, X, B::AbstractBasis)
    return get_vector(M.manifold, M.point, X, B)
end

function get_vector!(M::TangentSpace, Y, p, X, B::AbstractBasis)
    return get_vector!(M.manifold, Y, M.point, X, B)
end

function get_vectors(M::TangentSpace, p, B::CachedBasis)
    return get_vectors(M.manifold, M.point, B)
end

@doc raw"""
    injectivity_radius(M::TangentSpace)

Return the injectivity radius on the [`TangentSpace`](@ref) `M`, which is $∞$.
"""
injectivity_radius(::TangentSpace) = Inf

@doc raw"""
    inner(M::TangentSpace, X, Y, Z)

For any ``X∈T_p\mathcal M`` we identify the tangent space ``T_X(T_p\mathcal M)``
with ``T_p\mathcal M`` again. Hence an inner product of ``Y,Z`` is just the inner product of
the tangent space itself. ``⟨Y,Z⟩_X = ⟨Y,Z⟩_p``.
"""
function inner(M::TangentSpace, X, Y, Z)
    return inner(M.manifold, M.point, Y, Z)
end

"""
    is_flat(::TangentSpace)

Return true. [`TangentSpace`](@ref) is a flat manifold.
"""
is_flat(::TangentSpace) = true

function _isapprox(M::TangentSpace, X, Y; kwargs...)
    return isapprox(M.manifold, M.point, X, Y; kwargs...)
end

"""
    log(TpM::TangentSpace, p, q)

Logarithmic map on the tangent space manifold `TpM`, calculated as the difference of tangent
vectors `q` and `p` from `TpM`.
"""
log(::TangentSpace, ::Any...)
function log!(::TangentSpace, X, p, q)
    copyto!(X, q - p)
    return X
end

@doc raw"""
    manifold_dimension(TpM::TangentSpace)

Return the dimension of the tangent space ``T_p\mathcal M`` at ``p∈\mathcal M``,
which is the same as the dimension of the manifold ``\mathcal M``.
"""
function manifold_dimension(TpM::TangentSpace)
    return manifold_dimension(TpM.manifold)
end

@doc raw"""
    parallel_transport_to(::TangentSpace, X, Z, Y)

Transport the tangent vector ``Z ∈ T_X(T_p\mathcal M)`` from `X` to `Y`.
Since we identify ``T_X\mathcal M = T_p\mathcal M`` and the tangent space is a vector space,
parallel transport simplifies to the identity, so this function yield ``Z`` as a result.
"""
parallel_transport_to(TpM::TangentSpace, X, Z, Y)

function parallel_transport_to!(TpM::TangentSpace, Y, p, X, q)
    return copyto!(TpM.manifold, Y, p, X)
end

@doc raw"""
    project(M::TangentSpace, p)

Project the point `p` from the tangent space `M`, that is project the vector `p`
tangent at `M.point`.
"""
project(::TangentSpace, ::Any)

function project!(M::TangentSpace, q, p)
    return project!(M.manifold, q, M.point, p)
end

@doc raw"""
    project(M::TangentSpace, p, X)

Project the vector `X` from the tangent space `M`, that is project the vector `X`
tangent at `M.point`.
"""
project(::TangentSpace, ::Any, ::Any)

function project!(M::TangentSpace, Y, p, X)
    return project!(M.manifold, Y, M.point, X)
end

function Random.rand!(M::TangentSpace, X; vector_at = nothing)
    rand!(M.manifold, X; vector_at = M.point)
    return X
end
function Random.rand!(rng::AbstractRNG, M::TangentSpace, X; vector_at = nothing)
    rand!(rng, M.manifold, X; vector_at = M.point)
    return X
end

function representation_size(B::TangentSpace)
    return representation_size(B.manifold)
end

function Base.show(io::IO, ::MIME"text/plain", TpM::TangentSpace)
    println(io, "Tangent space to the manifold $(base_manifold(TpM)) at point:")
    pre = " "
    sp = sprint(show, "text/plain", TpM.point; context = io, sizehint = 0)
    sp = replace(sp, '\n' => "\n$(pre)")
    return print(io, pre, sp)
end


function Base.show(io::IO, ::MIME"text/plain", cTpM::CotangentSpace)
    println(io, "Cotangent space to the manifold $(base_manifold(cTpM)) at point:")
    pre = " "
    sp = sprint(show, "text/plain", cTpM.point; context = io, sizehint = 0)
    sp = replace(sp, '\n' => "\n$(pre)")
    return print(io, pre, sp)
end

@doc raw"""
    Y = Weingarten(M::TangentSpace, p, X, V)
    Weingarten!(M::TangentSpace, Y, p, X, V)

Compute the Weingarten map ``\mathcal W_p`` at `p` on the [`TangentSpace`](@ref) `M` with respect to the
tangent vector ``X \in T_p\mathcal M`` and the normal vector ``V \in N_p\mathcal M``.

Since this a flat space by itself, the result is always the zero tangent vector.
"""
Weingarten(::TangentSpace, p, X, V)

Weingarten!(::TangentSpace, Y, p, X, V) = fill!(Y, 0)

@doc raw"""
    zero_vector(M::TangentSpace, p)

Zero tangent vector at point `p` from the tangent space `M`, that is the zero tangent vector
at point `M.point`.
"""
zero_vector(::TangentSpace, ::Any...)

function zero_vector!(M::TangentSpace, X, p)
    return zero_vector!(M.manifold, X, M.point)
end
