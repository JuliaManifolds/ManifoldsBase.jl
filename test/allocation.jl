using ManifoldsBase
using Test
using ManifoldsBase:
    combine_allocation_promotion_functions, allocation_promotion_function, ℝ, ℂ

struct AllocManifold{𝔽} <: AbstractManifold{𝔽} end
AllocManifold() = AllocManifold{ℝ}()

function ManifoldsBase.exp!(::AllocManifold, v, x, y)
    v[1] .= x[1] .+ y[1]
    v[2] .= x[2] .+ y[2]
    v[3] .= x[3] .+ y[3]
    return v
end

struct AllocManifold2 <: AbstractManifold{ℝ} end
ManifoldsBase.representation_size(::AllocManifold2) = (2, 3)

struct AllocManifold3 <: AbstractManifold{ℂ} end
ManifoldsBase.representation_size(::AllocManifold3) = (2, 3)

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
    @test allocate(M, a, 5) isa Vector{Vector{Float64}}
    @test length(allocate(M, a, 5)) == 5
    @test allocate(M, a, (5,)) isa Vector{Vector{Float64}}
    @test length(allocate(M, a, (5,))) == 5
    @test allocate(M, a, Int, (5,)) isa Vector{Int}
    @test length(allocate(M, a, Int, (5,))) == 5
    @test allocate(M, a, Int, 5) isa Vector{Int}
    @test length(allocate(M, a, Int, 5)) == 5

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

    @test allocation_promotion_function(M, exp, (a, b)) === identity
    @test combine_allocation_promotion_functions(identity, identity) === identity
    @test combine_allocation_promotion_functions(identity, complex) === complex
    @test combine_allocation_promotion_functions(complex, identity) === complex
    @test combine_allocation_promotion_functions(complex, complex) === complex

    @test number_eltype([2.0]) === Float64
    @test number_eltype([[2.0], [3]]) === Float64
    @test number_eltype([[2], [3.0]]) === Float64
    @test number_eltype([[2], [3]]) === Int
    @test number_eltype(([2.0], [3])) === Float64
    @test number_eltype(([2], [3.0])) === Float64
    @test number_eltype(([2], [3])) === Int
    @test number_eltype(Any[[2.0], [3.0]]) === Float64
    @test number_eltype(typeof([[1.0, 2.0]])) === Float64

    M2 = AllocManifold2()
    alloc2 = ManifoldsBase.allocate_result(M2, rand)
    @test alloc2 isa Matrix{Float64}
    @test size(alloc2) == representation_size(M2)

    @test ManifoldsBase.allocate_result(AllocManifold3(), rand) isa Matrix{ComplexF64}

    an = allocate_as(M2)
    @test an isa Matrix{Float64}
    @test size(an) == representation_size(M2)
    an = allocate_as(M2, Array{Float32})
    @test an isa Matrix{Float32}
    @test size(an) == representation_size(M2)

    an = allocate_as(M2, TangentSpaceType())
    @test an isa Matrix{Float64}
    @test size(an) == representation_size(M2)
    an = allocate_as(M2, TangentSpaceType(), Array{Float32})
    @test an isa Matrix{Float32}
    @test size(an) == representation_size(M2)

    @test default_type(M2) === Matrix{Float64}
    @test default_type(M2, TangentSpaceType()) === Matrix{Float64}
end
