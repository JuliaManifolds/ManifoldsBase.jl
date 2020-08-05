using Test
using ManifoldsBase
@testset "ManifoldsBase" begin

    num_ambiguities = length(Test.detect_ambiguities(ManifoldsBase))
    #num_ambiguities > 0 && @warn "The number of ambiguities in ManifoldsBase is $(num_ambiguities)."
    if VERSION >= v"1.1"
        @test num_ambiguities == 0
    end
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
    include("vector_transport_along.jl")
end
