using ManifoldsBase
using Test

struct AllocManifold <: Manifold end

function ManifoldsBase.exp!(::AllocManifold, v, x, y)
    v[1] .= x[1] .+ y[1]
    v[2] .= x[2] .+ y[2]
    return v
end

@testset "Allocation" begin
    a = [[1.0], [2.0], [3.0]]
    b = [[2.0], [3.0], [-3.0]]
    M = AllocManifold()
    v = exp(M, a, b)
    @test v ≈ [[3.0], [5.0], [0.0]]
    @test allocate(([1.0], [2.0])) isa Tuple{Vector{Float64}, Vector{Float64}}
    @test allocate(([1.0], [2.0]), Int) isa Tuple{Vector{Int}, Vector{Int}}
end
