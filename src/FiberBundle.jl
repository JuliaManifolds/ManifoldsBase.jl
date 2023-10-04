
@doc raw"""
    FiberBundleProductVectorTransport{
        TMP<:AbstractVectorTransportMethod,
        TMV<:AbstractVectorTransportMethod,
    } <: AbstractVectorTransportMethod

Vector transport type on [`FiberBundle`](@ref). `method_point` is used for vector transport
of the point part and `method_fiber` is used for transport of the fiber part.

The vector transport is derived as a product manifold-style vector transport. The considered
product manifold is the product between the manifold $\mathcal M$ and the topological vector
space isometric to the fiber.

# Constructor

    FiberBundleProductVectorTransport(
        method_point::AbstractVectorTransportMethod,
        method_fiber::AbstractVectorTransportMethod,
    )
    FiberBundleProductVectorTransport()

By default both methods are set to `ParallelTransport`.
"""
struct FiberBundleProductVectorTransport{
    TMP<:AbstractVectorTransportMethod,
    TMV<:AbstractVectorTransportMethod,
} <: AbstractVectorTransportMethod
    method_point::TMP
    method_fiber::TMV
end

function FiberBundleProductVectorTransport()
    return FiberBundleProductVectorTransport(ParallelTransport(), ParallelTransport())
end
function FiberBundleProductVectorTransport(M::AbstractManifold)
    m = default_vector_transport_method(M)
    return FiberBundleProductVectorTransport(m, m)
end

"""
    FiberBundle{𝔽,TVS<:FiberType,TM<:AbstractManifold{𝔽},TVT<:FiberBundleProductVectorTransport} <: AbstractManifold{𝔽}

Fiber bundle on a [`AbstractManifold`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/types.html#ManifoldsBase.AbstractManifold) `M`
of type [`FiberType`](@ref). Examples include vector bundles, principal bundles or unit tangent bundles.

# Constructor

    FiberBundle(M::AbstractManifold, type::FiberType)
"""
struct FiberBundle{
    𝔽,
    TF<:FiberType,
    TM<:AbstractManifold{𝔽},
    TVT<:FiberBundleProductVectorTransport,
} <: AbstractManifold{𝔽}
    type::TF
    manifold::TM
    fiber::BundleFibers{TF,TM}
    vector_transport::TVT
end

vector_bundle_transport(::FiberType, M::AbstractManifold) = ParallelTransport()

function FiberBundle(
    fiber::TVS,
    M::TM,
    vtm::FiberBundleProductVectorTransport,
) where {TVS<:FiberType,TM<:AbstractManifold{𝔽}} where {𝔽}
    return FiberBundle{𝔽,TVS,TM,typeof(vtm)}(fiber, M, BundleFibers(fiber, M), vtm)
end
function FiberBundle(fiber::FiberType, M::AbstractManifold)
    vtmm = vector_bundle_transport(fiber, M)
    vtbm = FiberBundleProductVectorTransport(vtmm, vtmm)
    return FiberBundle(fiber, M, vtbm)
end

struct FiberBundleBasisData{BBasis<:CachedBasis,TBasis<:CachedBasis}
    base_basis::BBasis
    fiber_basis::TBasis
end


base_manifold(B::FiberBundle) = base_manifold(B.manifold)

"""
    bundle_projection(B::FiberBundle, p)

Projection of point `p` from the bundle `M` to the base manifold.
Returns the point on the base manifold `B.manifold` at which the vector part
of `p` is attached.
"""
bundle_projection(B::FiberBundle, p) = submanifold_component(B.manifold, p, Val(1))

function get_basis(M::FiberBundle, p, B::AbstractBasis)
    xp1 = submanifold_component(p, Val(1))
    base_basis = get_basis(M.manifold, xp1, B)
    fiber_basis = get_basis(M.fiber, xp1, B)
    return CachedBasis(B, FiberBundleBasisData(base_basis, fiber_basis))
end
function get_basis(M::FiberBundle, p, B::CachedBasis)
    return invoke(get_basis, Tuple{AbstractManifold,Any,CachedBasis}, M, p, B)
end

function get_basis(M::FiberBundle, p, B::DiagonalizingOrthonormalBasis)
    xp1 = submanifold_component(p, Val(1))
    bv1 = DiagonalizingOrthonormalBasis(submanifold_component(B.frame_direction, Val(1)))
    b1 = get_basis(M.manifold, xp1, bv1)
    bv2 = DiagonalizingOrthonormalBasis(submanifold_component(B.frame_direction, Val(2)))
    b2 = get_basis(M.fiber, xp1, bv2)
    return CachedBasis(B, FiberBundleBasisData(b1, b2))
end

function get_coordinates(M::FiberBundle, p, X, B::AbstractBasis)
    px, Vx = submanifold_components(M.manifold, p)
    VXM, VXF = submanifold_components(M.manifold, X)
    n = manifold_dimension(M.manifold)
    return vcat(
        get_coordinates(M.manifold, px, VXM, B),
        get_coordinates(M.fiber, px, VXF, B),
    )
end

function get_coordinates!(M::FiberBundle, Y, p, X, B::AbstractBasis)
    px, Vx = submanifold_components(M.manifold, p)
    VXM, VXF = submanifold_components(M.manifold, X)
    n = manifold_dimension(M.manifold)
    get_coordinates!(M.manifold, view(Y, 1:n), px, VXM, B)
    get_coordinates!(M.fiber, view(Y, (n + 1):length(Y)), px, VXF, B)
    return Y
end

function get_coordinates(
    M::FiberBundle,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:FiberBundleBasisData},
) where {𝔽}
    px, Vx = submanifold_components(M.manifold, p)
    VXM, VXF = submanifold_components(M.manifold, X)
    return vcat(
        get_coordinates(M.manifold, px, VXM, B.data.base_basis),
        get_coordinates(M.fiber, px, VXF, B.data.fiber_basis),
    )
end

function get_coordinates!(
    M::FiberBundle,
    Y,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:FiberBundleBasisData},
) where {𝔽}
    px, Vx = submanifold_components(M.manifold, p)
    VXM, VXF = submanifold_components(M.manifold, X)
    n = manifold_dimension(M.manifold)
    get_coordinates!(M.manifold, view(Y, 1:n), px, VXM, B.data.base_basis)
    get_coordinates!(M.fiber, view(Y, (n + 1):length(Y)), px, VXF, B.data.fiber_basis)
    return Y
end

function get_vector!(M::FiberBundle, Y, p, X, B::AbstractBasis)
    n = manifold_dimension(M.manifold)
    xp1 = submanifold_component(p, Val(1))
    get_vector!(M.manifold, submanifold_component(Y, Val(1)), xp1, X[1:n], B)
    get_vector!(M.fiber, submanifold_component(Y, Val(2)), xp1, X[(n + 1):end], B)
    return Y
end

function get_vector!(
    M::FiberBundle,
    Y,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:FiberBundleBasisData},
) where {𝔽}
    n = manifold_dimension(M.manifold)
    xp1 = submanifold_component(p, Val(1))
    get_vector!(
        M.manifold,
        submanifold_component(Y, Val(1)),
        xp1,
        X[1:n],
        B.data.base_basis,
    )
    get_vector!(
        M.fiber,
        submanifold_component(Y, Val(2)),
        xp1,
        X[(n + 1):end],
        B.data.fiber_basis,
    )
    return Y
end

function get_vectors(M::BundleFibers, p, B::CachedBasis)
    return get_vectors(M.manifold, p, B)
end

function _isapprox(B::FiberBundle, p, q; kwargs...)
    xp, Vp = submanifold_components(B.manifold, p)
    xq, Vq = submanifold_components(B.manifold, q)
    return isapprox(B.manifold, xp, xq; kwargs...) &&
           isapprox(FiberAtPoint(B.fiber, xp), Vp, Vq; kwargs...)
end
function _isapprox(B::FiberBundle, p, X, Y; kwargs...)
    px, Vx = submanifold_components(B.manifold, p)
    VXM, VXF = submanifold_components(B.manifold, X)
    VYM, VYF = submanifold_components(B.manifold, Y)
    return isapprox(B.manifold, VXM, VYM; kwargs...) &&
           isapprox(FiberAtPoint(B.fiber, px), VXF, VYF; kwargs...)
end

function manifold_dimension(B::FiberBundle)
    return manifold_dimension(B.manifold) + fiber_dimension(B.fiber)
end

function Random.rand!(M::FiberBundle, pX; vector_at = nothing)
    return rand!(Random.default_rng(), M, pX; vector_at = vector_at)
end
function Random.rand!(rng::AbstractRNG, M::FiberBundle, pX; vector_at = nothing)
    pXM, pXF = submanifold_components(M.manifold, pX)
    if vector_at === nothing
        rand!(rng, M.manifold, pXM)
        rand!(rng, FiberAtPoint(M.fiber, pXM), pXF)
    else
        vector_atM, vector_atF = submanifold_components(M.manifold, vector_at)
        rand!(rng, M.manifold, pXM; vector_at = vector_atM)
        rand!(rng, FiberAtPoint(M.fiber, pXM), pXF; vector_at = vector_atF)
    end
    return pX
end

@doc raw"""
    zero_vector(B::FiberBundle, p)

Zero tangent vector at point `p` from the fiber bundle `B`
over manifold `B.fiber` (denoted $\mathcal M$). The zero vector belongs to the space $T_{p}B$

Notation:
  * The point $p = (x_p, V_p)$ where $x_p ∈ \mathcal M$ and $V_p$ belongs to the
    fiber $F=π^{-1}(\{x_p\})$ of the vector bundle $B$ where $π$ is the
    canonical projection of that vector bundle $B$.

The zero vector is calculated as

$\mathbf{0}_{p} = (\mathbf{0}_{x_p}, \mathbf{0}_F)$

where $\mathbf{0}_{x_p}$ is the zero tangent vector from $T_{x_p}\mathcal M$ and
$\mathbf{0}_F$ is the zero element of the vector space $F$.
"""
zero_vector(::FiberBundle, ::Any...)

function zero_vector!(B::FiberBundle, X, p)
    xp, Vp = submanifold_components(B.manifold, p)
    VXM, VXF = submanifold_components(B.manifold, X)
    zero_vector!(B.manifold, VXM, xp)
    zero_vector!(B.fiber, VXF, Vp)
    return X
end
