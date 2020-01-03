using ManifoldsBase

using Test
using LinearAlgebra

struct ProjManifold <: Manifold end

ManifoldsBase.inner(::ProjManifold, x, w, v) = dot(w, v)
ManifoldsBase.project_tangent!(S::ProjManifold, w, x, v) = (w .= v .- dot(x, v) .* x)
ManifoldsBase.representation_size(::ProjManifold) = (2,3)
ManifoldsBase.manifold_dimension(::ProjManifold) = 5

@testset "Projected OrthonormalBasis" begin
    M = ProjManifold()
    x = [sqrt(2)/2 0.0 0.0;
         0.0 sqrt(2)/2 0.0]

    pb = basis(M, x, ProjectedOrthonormalBasis())
    N = manifold_dimension(M)
    @test length(pb) == N
    # test orthonormality
    for i in 1:N
        @test norm(M, x, pb[i]) ≈ 1
        for j in i+1:N
            @test inner(M, x, pb[i], pb[j]) ≈ 0 atol = 1e-15
        end
    end
    # check projection idempotency
    for i in 1:N
        @test project_tangent(M, x, pb[i]) ≈ pb[i]
    end
end
