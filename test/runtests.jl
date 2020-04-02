using Test
using ManifoldsBase

@testset "ManifoldsBase" begin
    # This should remain at 0
    @test length(Test.detect_ambiguities(ManifoldsBase)) == 4
    include("allocation.jl")
    include("numbers.jl")
    include("bases.jl")
    include("decorator_manifold.jl")
    include("empty_manifold.jl")
    include("default_manifold.jl")
    include("complex_manifold.jl")
    include("validation_manifold.jl")
    include("embedded_manifold.jl")
    include("domain_errors.jl")
end
