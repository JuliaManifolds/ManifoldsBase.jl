using Test
using ManifoldsBase
import ManifoldsBase: inner, change_representer, change_metric

@testset "Metrics" begin
    @test ManifoldsBase.EuclideanMetric() isa ManifoldsBase.AbstractMetric
end
