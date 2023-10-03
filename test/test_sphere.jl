using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: ℝ, ℂ, DefaultManifold, RealNumbers
using Test

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

# from Absil, Mahony, Trumpf, 2013 https://sites.uclouvain.be/absil/2013-01/Weingarten_07PA_techrep.pdf
function ManifoldsBase.Weingarten!(::TestSphere, Y, p, X, V)
    return Y .= -X * p' * V
end

@testset "TestSphere" begin
    @testset "ShootingInverseRetraction" begin
        vector_transports =
            [ScaledVectorTransport(ProjectionTransport()), PoleLadderTransport()]
        num_transport_points = [2, 4]
        @testset for M in [TestSphere(2), TestSphere(3, ℂ)]
            T = number_system(M) === ℝ ? Float64 : ComplexF64
            p = project(M, randn(T, representation_size(M)))
            X = project(M, p, randn(T, representation_size(M)))
            X ./= norm(M, p, X)
            q = exp(M, p, X)
            @testset for vector_transport in vector_transports, ntp in num_transport_points
                inverse_retraction = ShootingInverseRetraction(
                    ExponentialRetraction(),
                    ProjectionInverseRetraction(),
                    vector_transport,
                    ntp,
                    1e-9,
                    10_000,
                )
                Y = inverse_retract(M, p, q, inverse_retraction)
                @test isapprox(M, p, Y, X; atol = 1e-3)
                Y2 = similar(Y)
                inverse_retract!(M, Y2, p, q, inverse_retraction)
                @test isapprox(M, p, Y2, Y)
            end

            @testset "isapprox" begin
                @test_logs (:info,) !isapprox(M, p, q; error = :info)
                @test_logs (:info,) !isapprox(M, p, X, zero_vector(M, p); error = :info)
            end
        end
    end
    @testset "Weingarten" begin
        M = TestSphere(2)
        p = [1.0, 0.0, 0.0]
        X = [0.0, 0.2, 0.0]
        V = [0.1, 0.0, 0.0] #orthogonal to TpM -> parallel to p
        @test isapprox(M, p, Weingarten(M, p, X, V), -0.1 * X)
    end
end
