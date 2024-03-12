using ManifoldsBase: ℝ, ℂ, DefaultManifold, RealNumbers, EuclideanMetric
using LinearAlgebra, Random

# minimal sphere implementation for testing more complicated manifolds
struct TestSphere{N,𝔽} <: AbstractManifold{𝔽} end
TestSphere(N::Int, 𝔽 = ℝ) = TestSphere{N,𝔽}()

ManifoldsBase.representation_size(::TestSphere{N}) where {N} = (N + 1,)

function ManifoldsBase.change_representer!(
    M::TestSphere,
    Y,
    ::ManifoldsBase.EuclideanMetric,
    p,
    X,
)
    return copyto!(M, Y, p, X)
end

function ManifoldsBase.change_metric!(
    M::TestSphere,
    Y,
    ::ManifoldsBase.EuclideanMetric,
    p,
    X,
)
    return copyto!(M, Y, p, X)
end

function ManifoldsBase.check_point(M::TestSphere, p; kwargs...)
    if !isapprox(norm(p), 1.0; kwargs...)
        return DomainError(
            norm(p),
            "The point $(p) does not lie on the $(M) since its norm is not 1.",
        )
    end
    return nothing
end

function ManifoldsBase.check_vector(M::TestSphere, p, X; kwargs...)
    if !isapprox(abs(real(dot(p, X))), 0.0; kwargs...)
        return DomainError(
            abs(dot(p, X)),
            "The vector $(X) is not a tangent vector to $(p) on $(M), since it is not orthogonal in the embedding.",
        )
    end
    return nothing
end

function ManifoldsBase.exp!(M::TestSphere, q, p, X)
    return exp!(M, q, p, X, one(number_eltype(X)))
end
function ManifoldsBase.exp!(::TestSphere, q, p, X, t::Number)
    θ = abs(t) * norm(X)
    if θ == 0
        copyto!(q, p)
    else
        X_scale = t * sin(θ) / θ
        q .= p .* cos(θ) .+ X .* X_scale
    end
    return q
end

function ManifoldsBase.get_basis_diagonalizing(
    M::TestSphere{n},
    p,
    B::DiagonalizingOrthonormalBasis{ℝ},
) where {n}
    A = zeros(n + 1, n + 1)
    A[1, :] = transpose(p)
    A[2, :] = transpose(B.frame_direction)
    V = nullspace(A)
    κ = ones(n)
    if !iszero(B.frame_direction)
        # if we have a nonzero direction for the geodesic, add it and it gets curvature zero from the tensor
        V = hcat(B.frame_direction / norm(M, p, B.frame_direction), V)
        κ[1] = 0 # no curvature along the geodesic direction, if x!=y
    end
    T = typeof(similar(B.frame_direction))
    Ξ = [convert(T, V[:, i]) for i in 1:n]
    return CachedBasis(B, κ, Ξ)
end

@inline function nzsign(z, absz = abs(z))
    psignz = z / absz
    return ifelse(iszero(absz), one(psignz), psignz)
end


function ManifoldsBase.get_coordinates_orthonormal!(M::TestSphere, Y, p, X, ::RealNumbers)
    n = manifold_dimension(M)
    p1 = p[1]
    cosθ = abs(p1)
    λ = nzsign(p1, cosθ)
    pend, Xend = view(p, 2:(n + 1)), view(X, 2:(n + 1))
    factor = λ * X[1] / (1 + cosθ)
    Y .= Xend .- pend .* factor
    return Y
end

function ManifoldsBase.get_vector_orthonormal!(M::TestSphere, Y, p, X, ::RealNumbers)
    n = manifold_dimension(M)
    p1 = p[1]
    cosθ = abs(p1)
    λ = nzsign(p1, cosθ)
    pend = view(p, 2:(n + 1))
    pX = dot(pend, X)
    factor = pX / (1 + cosθ)
    Y[1] = -λ * pX
    Y[2:(n + 1)] .= X .- pend .* factor
    return Y
end

ManifoldsBase.inner(::TestSphere, p, X, Y) = dot(X, Y)

ManifoldsBase.injectivity_radius(::TestSphere) = π

function ManifoldsBase.inverse_retract_project!(M::TestSphere, X, p, q)
    X .= q .- p
    project!(M, X, p, X)
    return X
end

ManifoldsBase.is_flat(M::TestSphere) = manifold_dimension(M) == 1

function ManifoldsBase.log!(::TestSphere, X, p, q)
    cosθ = clamp(real(dot(p, q)), -1, 1)
    X .= q .- p .* cosθ
    nrmX = norm(X)
    if nrmX > 0
        X .*= acos(cosθ) / nrmX
    end
    return X
end

ManifoldsBase.manifold_dimension(::TestSphere{N}) where {N} = N

function ManifoldsBase.parallel_transport_to!(::TestSphere, Y, p, X, q)
    m = p .+ q
    mnorm2 = real(dot(m, m))
    factor = 2 * real(dot(X, q)) / mnorm2
    Y .= X .- m .* factor
    return Y
end

function ManifoldsBase.project!(::TestSphere, q, p)
    q .= p ./ norm(p)
    return q
end

function ManifoldsBase.project!(::TestSphere, Y, p, X)
    Y .= X .- p .* real(dot(p, X))
    return Y
end

function Random.rand!(M::TestSphere, pX; vector_at = nothing, σ = one(eltype(pX)))
    return rand!(Random.default_rng(), M, pX; vector_at = vector_at, σ = σ)
end
function Random.rand!(
    rng::AbstractRNG,
    M::TestSphere,
    pX;
    vector_at = nothing,
    σ = one(eltype(pX)),
)
    if vector_at === nothing
        project!(M, pX, randn(rng, eltype(pX), representation_size(M)))
    else
        n = σ * randn(rng, eltype(pX), size(pX)) # Gaussian in embedding
        project!(M, pX, vector_at, n) #project to TpM (keeps Gaussianness)
    end
    return pX
end

function ManifoldsBase.riemann_tensor!(M::TestSphere, Xresult, p, X, Y, Z)
    innerZX = inner(M, p, Z, X)
    innerZY = inner(M, p, Z, Y)
    Xresult .= innerZY .* X .- innerZX .* Y
    return Xresult
end

function ManifoldsBase.sectional_curvature_max(M::TestSphere)
    return ifelse(manifold_dimension(M) == 1, 0.0, 1.0)
end
function ManifoldsBase.sectional_curvature_min(M::TestSphere)
    return ifelse(manifold_dimension(M) == 1, 0.0, 1.0)
end

# from Absil, Mahony, Trumpf, 2013 https://sites.uclouvain.be/absil/2013-01/Weingarten_07PA_techrep.pdf
function ManifoldsBase.Weingarten!(::TestSphere, Y, p, X, V)
    return Y .= -X * p' * V
end

# minimal implementation of SPD (for its nontrivial metric conversion)

struct TestSPD <: AbstractManifold{ℝ}
    n::Int
end

function ManifoldsBase.change_representer!(::TestSPD, Y, ::EuclideanMetric, p, X)
    Y .= p * X * p
    return Y
end

function ManifoldsBase.change_metric!(::TestSPD, Y, ::EuclideanMetric, p, X)
    Y .= p * X
    return Y
end

function ManifoldsBase.inner(::TestSPD, p, X, Y)
    F = cholesky(Symmetric(convert(AbstractMatrix, p)))
    return dot((F \ Symmetric(X)), (Symmetric(Y) / F))
end

function spd_sqrt_and_sqrt_inv(p::AbstractMatrix)
    e = eigen(Symmetric(p))
    U = e.vectors
    S = max.(e.values, floatmin(eltype(e.values)))
    Ssqrt = Diagonal(sqrt.(S))
    SsqrtInv = Diagonal(1 ./ sqrt.(S))
    return (Symmetric(U * Ssqrt * transpose(U)), Symmetric(U * SsqrtInv * transpose(U)))
end

function ManifoldsBase.log!(::TestSPD, X, p, q)
    (p_sqrt, p_sqrt_inv) = spd_sqrt_and_sqrt_inv(p)
    T = Symmetric(p_sqrt_inv * convert(AbstractMatrix, q) * p_sqrt_inv)
    e2 = eigen(T)
    Se = Diagonal(log.(max.(e2.values, eps())))
    pUe = p_sqrt * e2.vectors
    return mul!(X, pUe, Se * transpose(pUe))
end

ManifoldsBase.representation_size(M::TestSPD) = (M.n, M.n)
