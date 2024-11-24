using Test
using ManifoldsBase
@testset "ManifoldsBase" begin
    bound = 6
    # six ambiguities come from possible incorrectly formed calls to `allocate`
    num_ambiguities = length(Test.detect_ambiguities(ManifoldsBase))
    #num_ambiguities > 0 && @warn "The number of ambiguities in ManifoldsBase is $(num_ambiguities)."
    if VERSION >= v"1.9-DEV"
        @test num_ambiguities <= bound + 6
    elseif VERSION >= v"1.8-DEV"
        @test num_ambiguities <= bound + 5
    elseif VERSION >= v"1.6-DEV"
        # At the time of writing there seem to be two ambiguities regarding `getindex`,
        # one with a method from SparseArrays and one from VSCode's JSON processing
        # that's automatically loaded when running code in VSCode.
        @test num_ambiguities <= bound + 3
    else
        @test num_ambiguities == bound + 1
    end

    include("decorator_traits.jl")
    include("passthrough_decorator.jl")
    include("allocation.jl")
    include("numbers.jl")
    include("bases.jl")
    include("manifold_fallbacks.jl")
    include("empty_manifold.jl")
    include("errors.jl")
    include("default_manifold.jl")
    include("complex_manifold.jl")
    include("validation_manifold.jl")
    include("embedded_manifold.jl")
    include("test_sphere.jl")
    include("product_manifold.jl")
    include("power.jl")
    include("domain_errors.jl")
    include("vector_transport.jl")
    include("metric.jl")
    include("fibers.jl")
    include("numerical_checks.jl")
    include("test_aqua.jl")
end
