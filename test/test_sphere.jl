using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: ℝ, ℂ, DefaultManifold, RealNumbers
using Test

s = @__DIR__
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))

using ManifoldsBaseTestUtils

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
                    1.0e-9,
                    10_000,
                )
                Y = inverse_retract(M, p, q, inverse_retraction)
                @test isapprox(M, p, Y, X; atol = 1.0e-3)
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
    @testset "Tangent Space" begin
        M = TestSphere(2)
        p = [1.0, 0.0, 0.0]
        X = [0.0, 0.2, 0.0]
        TpM = TangentSpace(M, p)
        @test is_point(TpM, X)
        @test !is_point(TpM, p)
        @test is_vector(TpM, X, X)
        @test !is_vector(TpM, X, p)
        @test zero_vector(TpM) == zero_vector(M, p)
    end
end
