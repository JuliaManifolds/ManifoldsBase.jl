using ManifoldsBase
using Test

@testset "Testing decorator manifold functions" begin
    M = ManifoldsBase.DefaultManifold(3)
    A = ArrayManifold(M)

    @test (@inferred base_manifold(M)) == M
    @test (@inferred base_manifold(M, Val(false))) == M
    @test (@inferred base_manifold(A)) == M
    @test (@inferred base_manifold(A, Val(true))) == M

    @test (@inferred base_manifold(M, Val(1))) == M
    @test (@inferred base_manifold(M, Val(0))) == M
    @test (@inferred base_manifold(A, Val(1))) == M
    @test (@inferred base_manifold(A, Val(0))) == A

    @test representation_size(M) == (3,)
    @test_throws ErrorException representation_size(M, Val(false))
    @test representation_size(A) == (3,)
    @test representation_size(A, Val(true)) == (3,)

    @test manifold_dimension(M) == 3
    @test_throws ErrorException manifold_dimension(M, Val(false))
    @test manifold_dimension(A) == 3
    @test manifold_dimension(A, Val(true)) == 3
end
