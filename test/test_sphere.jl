using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: ℝ, ℂ, DefaultManifold
using Test

# minimal sphere implementation for testing more complicated manifolds
struct TestSphere{N,𝔽} <: AbstractManifold{𝔽} end
TestSphere(N::Int, 𝔽 = ℝ) = TestSphere{N,𝔽}()

ManifoldsBase.representation_size(::TestSphere{N}) where {N} = (N + 1,)

function ManifoldsBase.exp!(M::TestSphere, q, p, X)
    return exp!(M, q, p, X, one(number_eltype(X)))
end
function ManifoldsBase.exp!(::TestSphere, q, p, X, t::Number)
    θ = norm(X)
    if θ == 0
        copyto!(q, p)
    else
        X_scale = t * sin(θ) / θ
        q .= p .* cos(θ) .+ X .* X_scale
    end
    return q
end

ManifoldsBase.inner(::TestSphere, p, X, Y) = dot(X, Y)

ManifoldsBase.injectivity_radius(::TestSphere) = π

function ManifoldsBase.inverse_retract_project!(M::TestSphere, X, p, q)
    X .= q .- p
    project!(M, X, p, X)
    return X
end

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

function parallel_transport_to!(::TestSphere, Y, p, X, q)
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
