using Test
using ManifoldsBase
@testset "ManifoldsBase" begin

    num_ambiguities = length(Test.detect_ambiguities(ManifoldsBase))
    #num_ambiguities > 0 && @warn "The number of ambiguities in ManifoldsBase is $(num_ambiguities)."
    @test num_ambiguities == 0
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
