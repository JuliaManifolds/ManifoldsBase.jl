using Test

@testset "ManifoldsBase" begin
    include("allocation.jl")
    include("decorator_manifold.jl")
    include("empty_manifold.jl")
    include("default_manifold.jl")
    include("complex_manifold.jl")
    include("array_manifold.jl")
    include("domain_errors.jl")
end
