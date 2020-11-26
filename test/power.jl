using Test
using ManifoldsBase
using ManifoldsBase: AbstractNumbers, ℝ, ℂ
@testset "Power Manifold" begin
    M = ManifoldsBase.DefaultManifold(3)
    N = PowerManifold(M, NestedPowerRepresentation(), 2)
    @test repr(O) ==
          "PowerManifold(DefaultManifold(3; field = ℝ), NestedPowerRepresentation(), 2)"
    O = N^3
    @test repr(O) ==
          "PowerManifold(DefaultManifold(3; field = ℝ), NestedPowerRepresentation(), 2, 3)"
end
