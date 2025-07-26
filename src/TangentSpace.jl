@doc raw"""
    TangentSpace{𝔽,M} = Fiber{𝔽,TangentSpaceType,M} where {𝔽,M<:AbstractManifold}

A manifold for the tangent space ``T_p\mathcal M`` at a point ``p\in\mathcal M``.
This is modelled as an alias for [`VectorSpaceFiber`](@ref) corresponding to
[`TangentSpaceType`](@ref).

# Constructor

    TangentSpace(M::AbstractManifold, p)

Return the manifold (vector space) representing the tangent space ``T_p\mathcal M``
at point `p`, ``p\in\mathcal M``.
"""
const TangentSpace{𝔽, M} = Fiber{𝔽, TangentSpaceType, M} where {𝔽, M <: AbstractManifold}

TangentSpace(M::AbstractManifold, p) = Fiber(M, p, TangentSpaceType())

@doc raw"""
    CotangentSpace{𝔽,M} = Fiber{𝔽,CotangentSpaceType,M} where {𝔽,M<:AbstractManifold}

A manifold for the Cotangent space ``T^*_p\mathcal M`` at a point ``p\in\mathcal M``.
This is modelled as an alias for [`VectorSpaceFiber`](@ref) corresponding to
[`CotangentSpaceType`](@ref).

# Constructor

    CotangentSpace(M::AbstractManifold, p)

Return the manifold (vector space) representing the cotangent space ``T^*_p\mathcal M``
at point `p`, ``p\in\mathcal M``.
"""
const CotangentSpace{𝔽, M} = Fiber{𝔽, CotangentSpaceType, M} where {𝔽, M <: AbstractManifold}

CotangentSpace(M::AbstractManifold, p) = Fiber(M, p, CotangentSpaceType())


function allocate_result(M::TangentSpace, ::typeof(rand))
    return zero_vector(M.manifold, M.point)
end

@doc raw"""
    base_point(TpM::TangentSpace)

Return the base point of the [`TangentSpace`](@ref).
"""
base_point(TpM::TangentSpace) = TpM.point

# forward both point checks to tangent vector checks
function check_point(TpM::TangentSpace, p; kwargs...)
    return check_vector(TpM.manifold, TpM.point, p; kwargs...)
end
function check_size(TpM::TangentSpace, p; kwargs...)
    return check_size(TpM.manifold, TpM.point, p; kwargs...)
end
# fix tangent vector checks to use the right base point
function check_vector(TpM::TangentSpace, p, X; kwargs...)
    return check_vector(TpM.manifold, TpM.point, X; kwargs...)
end
function check_size(TpM::TangentSpace, p, X; kwargs...)
    return check_size(TpM.manifold, TpM.point, X; kwargs...)
end


"""
    distance(M::TangentSpace, X, Y)

Distance between vectors `X` and `Y` from the [`TangentSpace`](@ref) `TpM`.
It is calculated as the [`norm`](@ref) (induced by the metric on `TpM`) of their difference.
"""
function distance(TpM::TangentSpace, X, Y)
    return norm(TpM.manifold, TpM.point, Y - X)
end

function embed!(TpM::TangentSpace, Y, X)
    return embed!(TpM.manifold, Y, TpM.point, X)
end
function embed!(TpM::TangentSpace, W, X, V)
    return embed!(TpM.manifold, W, TpM.point, V)
end

@doc raw"""
    exp(TpM::TangentSpace, X, V)

Exponential map of tangent vectors `X` from `TpM` and a direction `V`,
which is also from the [`TangentSpace`](@ref) `TpM` since we identify the tangent space of `TpM` with `TpM`.
The exponential map then simplifies to the sum `X+V`.
"""
exp(::TangentSpace, ::Any, ::Any)

function exp!(TpM::TangentSpace, Y, X, V)
    copyto!(TpM.manifold, Y, TpM.point, X + V)
    return Y
end

fiber_dimension(M::AbstractManifold, ::CotangentSpaceType) = manifold_dimension(M)
fiber_dimension(M::AbstractManifold, ::TangentSpaceType) = manifold_dimension(M)

function get_basis(TpM::TangentSpace, X, B::CachedBasis)
    return invoke(
        get_basis,
        Tuple{AbstractManifold, Any, CachedBasis},
        TpM.manifold,
        TpM.point,
        B,
    )
end
function get_basis(TpM::TangentSpace, X, B::AbstractBasis{<:Any, TangentSpaceType})
    return get_basis(TpM.manifold, TpM.point, B)
end

function get_coordinates(TpM::TangentSpace, X, V, B::AbstractBasis)
    return get_coordinates(TpM.manifold, TpM.point, V, B)
end

function get_coordinates!(TpM::TangentSpace, c, X, V, B::AbstractBasis)
    return get_coordinates!(TpM.manifold, c, TpM.point, V, B)
end

function get_vector(TpM::TangentSpace, X, c, B::AbstractBasis)
    return get_vector(TpM.manifold, TpM.point, c, B)
end

function get_vector!(TpM::TangentSpace, V, X, c, B::AbstractBasis)
    return get_vector!(TpM.manifold, V, TpM.point, c, B)
end

function get_vectors(TpM::TangentSpace, X, B::CachedBasis)
    return get_vectors(TpM.manifold, TpM.point, B)
end

@doc raw"""
    injectivity_radius(TpM::TangentSpace)

Return the injectivity radius on the [`TangentSpace`](@ref) `TpM`, which is $∞$.
"""
injectivity_radius(::TangentSpace) = Inf

@doc raw"""
    inner(M::TangentSpace, X, V, W)

For any ``X ∈ T_p\mathcal M`` we identify the tangent space ``T_X(T_p\mathcal M)``
with ``T_p\mathcal M`` again. Hence an inner product of ``V,W`` is just the inner product of
the tangent space itself. ``⟨V,W⟩_X = ⟨V,W⟩_p``.
"""
function inner(TpM::TangentSpace, X, V, W)
    return inner(TpM.manifold, TpM.point, V, W)
end

"""
    is_flat(::TangentSpace)

The [`TangentSpace`](@ref) is a flat manifold, so this returns `true`.
"""
is_flat(::TangentSpace) = true

function _isapprox(TpM::TangentSpace, X, Y; kwargs...)
    return isapprox(TpM.manifold, TpM.point, X, Y; kwargs...)
end

function _isapprox(TpM::TangentSpace, X, V, W; kwargs...)
    return isapprox(TpM.manifold, TpM.point, V, W; kwargs...)
end

"""
    log(TpM::TangentSpace, X, Y)

Logarithmic map on the [`TangentSpace`](@ref) `TpM`, calculated as the difference of tangent
vectors `q` and `p` from `TpM`.
"""
log(::TangentSpace, ::Any...)
function log!(TpM::TangentSpace, V, X, Y)
    copyto!(TpM, V, TpM.point, Y - X)
    return V
end

@doc raw"""
    manifold_dimension(TpM::TangentSpace)

Return the dimension of the [`TangentSpace`](@ref) ``T_p\mathcal M`` at ``p∈\mathcal M``,
which is the same as the dimension of the manifold ``\mathcal M``.
"""
function manifold_dimension(TpM::TangentSpace)
    return manifold_dimension(TpM.manifold)
end

@doc raw"""
    parallel_transport_to(::TangentSpace, X, V, Y)

Transport the tangent vector ``Z ∈ T_X(T_p\mathcal M)`` from `X` to `Y`.
Since we identify ``T_X(T_p\mathcal M) = T_p\mathcal M`` and the tangent space is a vector space,
parallel transport simplifies to the identity, so this function yields ``V`` as a result.
"""
parallel_transport_to(TpM::TangentSpace, X, V, Y)

function parallel_transport_to!(TpM::TangentSpace, W, X, V, Y)
    return copyto!(TpM.manifold, W, TpM.point, V)
end

@doc raw"""
    project(TpM::TangentSpace, X)

Project the point `X` from embedding of the [`TangentSpace`](@ref) `TpM` onto `TpM`.
"""
project(::TangentSpace, ::Any)

function project!(TpM::TangentSpace, Y, X)
    return project!(TpM.manifold, Y, TpM.point, X)
end

@doc raw"""
    project(TpM::TangentSpace, X, V)

Project the vector `V` from the embedding of the tangent space `TpM` (identified with ``T_X(T_p\mathcal M)``),
that is project the vector `V` onto the tangent space at `TpM.point`.
"""
project(::TangentSpace, ::Any, ::Any)

function project!(TpM::TangentSpace, W, X, V)
    return project!(TpM.manifold, W, TpM.point, V)
end

function Random.rand!(TpM::TangentSpace, X; vector_at = nothing)
    rand!(TpM.manifold, X; vector_at = TpM.point)
    return X
end
function Random.rand!(rng::AbstractRNG, TpM::TangentSpace, X; vector_at = nothing)
    rand!(rng, TpM.manifold, X; vector_at = TpM.point)
    return X
end

function representation_size(TpM::TangentSpace)
    return representation_size(TpM.manifold)
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
    Y = Weingarten(TpM::TangentSpace, X, V, A)
    Weingarten!(TpM::TangentSpace, Y, p, X, V)

Compute the Weingarten map ``\mathcal W_X`` at `X` on the [`TangentSpace`](@ref) `TpM` with respect to the
tangent vector ``V \in T_p\mathcal M`` and the normal vector ``A \in N_p\mathcal M``.

Since this a flat space by itself, the result is always the zero tangent vector.
"""
Weingarten(::TangentSpace, ::Any, ::Any, ::Any)

Weingarten!(::TangentSpace, W, X, V, A) = fill!(W, 0)

@doc raw"""
    zero_vector(TpM::TangentSpace)

Zero tangent vector in the [`TangentSpace`](@ref) `TpM`,
that is the zero tangent vector at point `TpM.point`.
"""
zero_vector(TpM::TangentSpace) = zero_vector(TpM.manifold, TpM.point)

@doc raw"""
    zero_vector(TpM::TangentSpace, X)

Zero tangent vector at point `X` from the [`TangentSpace`](@ref) `TpM`,
that is the zero tangent vector at point `TpM.point`,
since we identify the tangent space ``T_X(T_p\mathcal M)`` with ``T_p\mathcal M``.
"""
zero_vector(::TangentSpace, ::Any...)

function zero_vector!(M::TangentSpace, V, X)
    return zero_vector!(M.manifold, V, M.point)
end
