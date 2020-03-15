using Test
using ManifoldsBase
using ManifoldsBase: AbstractNumbers, ℝ, ℂ, ℍ

struct NotImplementedNumbers <: ManifoldsBase.AbstractNumbers end

@testset "Number systems" begin
    @test_throws ErrorException real_dimension(NotImplementedNumbers())

    @test ℝ isa ManifoldsBase.RealNumbers
    @test ManifoldsBase.RealNumbers() === ℝ
    @test real_dimension(ℝ) == 1
    @test repr(ℝ) == "ℝ"

    @test ℂ isa ManifoldsBase.ComplexNumbers
    @test ManifoldsBase.ComplexNumbers() === ℂ
    @test real_dimension(ℂ) == 2
    @test repr(ℂ) == "ℂ"

    @test ℍ isa ManifoldsBase.QuaternionNumbers
    @test ManifoldsBase.QuaternionNumbers() === ℍ
    @test real_dimension(ℍ) == 4
    @test repr(ℍ) == "ℍ"
end
