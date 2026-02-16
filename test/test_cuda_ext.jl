using ManifoldsBase
using ManifoldsBase: ℝ, ℂ
using Test
using LinearAlgebra

@testset "ManifoldsBaseCUDAExt" begin
    cuda_loaded = false
    try
        using CUDA
        cuda_loaded = CUDA.functional()
    catch
        cuda_loaded = false
    end

    if cuda_loaded
        @eval using CUDA

        @testset "allocate preserves CuArray" begin
            # allocate(p) → similar(p) — should preserve CuArray
            p = CUDA.zeros(Float64, 4)
            q = ManifoldsBase.allocate(p)
            @test q isa CuArray{Float64,1}
            @test size(q) == (4,)

            # allocate(p, T) — change element type, keep on GPU
            q32 = ManifoldsBase.allocate(p, Float32)
            @test q32 isa CuArray{Float32,1}
            @test size(q32) == (4,)

            # allocate(M, p) — with manifold argument
            M = ManifoldsBase.DefaultManifold(4)
            q_m = ManifoldsBase.allocate(M, p)
            @test q_m isa CuArray{Float64,1}
            @test size(q_m) == (4,)

            # allocate(M, p, T) — with manifold and type
            q_m32 = ManifoldsBase.allocate(M, p, Float32)
            @test q_m32 isa CuArray{Float32,1}
            @test size(q_m32) == (4,)
        end

        @testset "allocate preserves CuMatrix" begin
            P = CUDA.zeros(Float64, 3, 3)
            Q = ManifoldsBase.allocate(P)
            @test Q isa CuArray{Float64,2}
            @test size(Q) == (3, 3)

            Q32 = ManifoldsBase.allocate(P, ComplexF64)
            @test Q32 isa CuArray{ComplexF64,2}
            @test size(Q32) == (3, 3)
        end

        @testset "allocate_result with point preserves CuArray" begin
            M = ManifoldsBase.DefaultManifold(4)
            p = CUDA.zeros(Float64, 4)
            # allocate_result(M, f, p) uses allocate(M, p, T) — should stay on GPU
            q = ManifoldsBase.allocate_result(M, exp, p, p)
            @test q isa CuArray{Float64,1}
            @test size(q) == (4,)
        end

        @testset "allocate_result without point returns CPU Array" begin
            # This is the ROOT CAUSE of GPU issues — allocate_result(M, f) has no
            # point reference to infer CuArray from, so it returns Array.
            # This is expected behavior — the fix is at the solver level.
            M = ManifoldsBase.DefaultManifold(4)
            q = ManifoldsBase.allocate_result(M, rand)
            @test q isa Array{Float64}
            @test size(q) == (4,)
        end

        @testset "allocate_on with CuArray type" begin
            M = ManifoldsBase.DefaultManifold(4)
            p = ManifoldsBase.allocate_on(M, CuArray{Float64})
            @test p isa CuArray{Float64,1}
            @test size(p) == (4,)

            p2 = ManifoldsBase.allocate_on(M, CuArray{Float32})
            @test p2 isa CuArray{Float32,1}
            @test size(p2) == (4,)
        end

        @testset "allocate nested CuArrays" begin
            # Vector{CuArray} — used by PowerManifold with NestedReplacing
            a = [CUDA.zeros(Float64, 2, 2) for _ in 1:3]
            b = ManifoldsBase.allocate(a)
            @test b isa Vector{<:CuArray{Float64,2}}
            @test length(b) == 3
            @test size(b[1]) == (2, 2)

            # With type change
            c = ManifoldsBase.allocate(a, Float32)
            @test c isa Vector{<:CuArray{Float32,2}}
            @test length(c) == 3
        end

        @testset "zero_vector preserves CuArray" begin
            M = ManifoldsBase.DefaultManifold(4)
            p = CuArray(ones(Float64, 4))
            X = ManifoldsBase.zero_vector(M, p)
            @test X isa CuArray{Float64,1}
            @test all(Array(X) .== 0.0)
        end

        @testset "copyto! mixed CPU/GPU" begin
            M = ManifoldsBase.DefaultManifold(4)
            p_cpu = ones(Float64, 4)
            p_gpu = CUDA.zeros(Float64, 4)

            # CPU → GPU
            copyto!(p_gpu, p_cpu)
            @test all(Array(p_gpu) .== 1.0)

            # GPU → CPU
            p_cpu2 = zeros(Float64, 4)
            copyto!(p_cpu2, p_gpu)
            @test all(p_cpu2 .== 1.0)
        end

        @testset "basic manifold operations with CuArray" begin
            M = ManifoldsBase.DefaultManifold(4)
            p = CuArray(randn(Float64, 4))
            X = CuArray(randn(Float64, 4))

            # exp
            q = exp(M, p, X)
            @test q isa CuArray{Float64,1}
            @test isapprox(Array(q), Array(p) + Array(X))

            # retract
            q2 = ManifoldsBase.retract(M, p, X)
            @test q2 isa CuArray{Float64,1}

            # inner, norm
            Y = CuArray(randn(Float64, 4))
            ip = ManifoldsBase.inner(M, p, X, Y)
            @test ip isa Real
            n = ManifoldsBase.norm(M, p, X)
            @test n isa Real
            @test n >= 0
        end
    else
        @info "CUDA not functional, skipping ManifoldsBaseCUDAExt tests"
    end
end
