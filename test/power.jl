using Test
using ManifoldsBase
using ManifoldsBase: AbstractNumbers, ℝ, ℂ

struct DummyPowerRepresentation <: AbstractPowerRepresentation end

@testset "Power Manifold" begin
    M = ManifoldsBase.DefaultManifold(3)
    N = PowerManifold(M, NestedPowerRepresentation(), 2)
    p = [zeros(3), ones(3)]
    q = [ones(3), zeros(3)]
    @testset "Constructors" begin
        @test repr(N) ==
              "PowerManifold(DefaultManifold(3; field = ℝ), NestedPowerRepresentation(), 2)"
        # add to type
        @test repr(N^3) ==
              "PowerManifold(DefaultManifold(3; field = ℝ), NestedPowerRepresentation(), 2, 3)"
        # add to type
        @test repr(PowerManifold(N, 3)) ==
              "PowerManifold(DefaultManifold(3; field = ℝ), NestedPowerRepresentation(), 2, 3)"
        # switch type
        @test repr(PowerManifold(N, DummyPowerRepresentation(), 3)) ==
              "PowerManifold(DefaultManifold(3; field = ℝ), DummyPowerRepresentation(), 2, 3)"
        # nest
        @test repr(PowerManifold(N, NestedPowerRepresentation(), 3)) ==
              "PowerManifold(PowerManifold(DefaultManifold(3; field = ℝ), NestedPowerRepresentation(), 2), NestedPowerRepresentation(), 3)"
    end
    @testset "checks" begin
        pE1 = [zeros(3), ones(2)] # one component wrong
        pE2 = [zeros(2), ones(2)] # both wrong
        @test is_manifold_point(N, p)
        @test is_manifold_point(N, p, true)
        @test !is_manifold_point(N, pE1)
        @test !is_manifold_point(N, pE2)
        @test_throws ComponentManifoldError is_manifold_point(N, pE1, true)
        @test_throws CompositeManifoldError is_manifold_point(N, pE2, true)
        # tangent - test base
        @test is_tangent_vector(N, p, p)
        @test is_tangent_vector(N, p, p, true)
        @test !is_tangent_vector(N, pE1, p)
        @test !is_tangent_vector(N, pE2, p)
        # tangents - with proper base
        @test is_tangent_vector(N, p, p, true)
        @test !is_tangent_vector(N, p, pE1)
        @test !is_tangent_vector(N, p, pE2)
        @test_throws ComponentManifoldError is_tangent_vector(N, p, pE1, true)
        @test_throws CompositeManifoldError is_tangent_vector(N, p, pE2, true)
    end
    @testset "specific functions" begin
        @test distance(N, p, q) == sqrt(sum(distance.(Ref(M), p, q) .^ 2))
        @test exp(N, p, q) == p .+ q
        @test log(N, p, q) == q .- p
        @test p[N, 1] == p[1]
        p[N, 1] = 2 .* ones(3)
        @test p[N, 1] == 2 .* ones(3)
    end

end
