using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: ℝ, ℂ, DefaultManifold
using Test

# minimal sphere implementation for testing more complicated manifolds
struct TestSphere{N,𝔽} <: AbstractManifold{𝔽} end
TestSphere(N::Int, 𝔽 = ℝ) = TestSphere{N,𝔽}()

ManifoldsBase.representation_size(::TestSphere{N}) where {N} = (N + 1,)

function ManifoldsBase.exp!(::TestSphere, q, p, X)
    θ = norm(X)
    if θ == 0
        copyto!(q, p)
    else
        q .= p .* cos(θ) .+ X .* sin(θ) ./ θ
    end
    return q
end

ManifoldsBase.inner(::TestSphere, p, X, Y) = dot(X, Y)

function ManifoldsBase.inverse_retract_project!(::TestSphere, X, p, q)
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

function ManifoldsBase.project!(::TestSphere, q, p)
    q .= p ./ norm(p)
    return q
end

function ManifoldsBase.project!(::TestSphere, Y, p, X)
    Y .= X .- p .* real(dot(p, X))
    return Y
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
                @test isapprox(M, p, Y, X; atol = 1e-6)
            end
        end
    end
end
