using Test
using ManifoldsBase

@testset "Metrics" begin
    @test ManifoldsBase.EuclideanMetric() isa ManifoldsBase.AbstractMetric
end
