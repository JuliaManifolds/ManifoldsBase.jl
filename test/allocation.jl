using ManifoldsBase
using Test
using ManifoldsBase: combine_allocation_promotion_functions

struct AllocManifold <: Manifold end

function ManifoldsBase.exp!(::AllocManifold, v, x, y)
    v[1] .= x[1] .+ y[1]
    v[2] .= x[2] .+ y[2]
    v[3] .= x[3] .+ y[3]
    return v
end

@testset "Allocation" begin
    a = [[1.0], [2.0], [3.0]]
    b = [[2.0], [3.0], [-3.0]]
    M = AllocManifold()
    v = exp(M, a, b)
    @test v ≈ [[3.0], [5.0], [0.0]]
    @test allocate(([1.0], [2.0])) isa Tuple{Vector{Float64},Vector{Float64}}
    @test allocate(([1.0], [2.0]), Int) isa Tuple{Vector{Int},Vector{Int}}
    @test allocate([[1.0], [2.0]]) isa Vector{Vector{Float64}}
    @test allocate([[1.0], [2.0]], Int) isa Vector{Vector{Int}}

    a1 = allocate([1], 2, 3)
    @test a1 isa Matrix{Int}
    @test size(a1) == (2, 3)
    a2 = allocate([1], (2, 3))
    @test a2 isa Matrix{Int}
    @test size(a2) == (2, 3)
    a3 = allocate([1], Float64, 2, 3)
    @test a3 isa Matrix{Float64}
    @test size(a3) == (2, 3)
    a4 = allocate([1], Float64, (2, 3))
    @test a4 isa Matrix{Float64}
    @test size(a4) == (2, 3)

    @test combine_allocation_promotion_functions(identity, identity) === identity
    @test combine_allocation_promotion_functions(identity, complex) === complex
    @test combine_allocation_promotion_functions(complex, identity) === complex
    @test combine_allocation_promotion_functions(complex, complex) === complex

    @test number_eltype([2.0]) == Float64
    @test number_eltype([[2.0], [3]]) == Float64
    @test number_eltype([[2], [3.0]]) == Float64
    @test number_eltype([[2], [3]]) == Int
    @test number_eltype(([2.0], [3])) == Float64
    @test number_eltype(([2], [3.0])) == Float64
    @test number_eltype(([2], [3])) == Int
end
