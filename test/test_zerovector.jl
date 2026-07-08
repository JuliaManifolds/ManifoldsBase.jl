using ManifoldsBase, Test

@testset "ZeroVector" begin
    X = ZeroVector()
    Y = [1.0, 2.0, 3.0]
    @test X + Y == Y
    @test Y + X == Y
    @test X + X == X
    @test X - Y == -Y
    @test Y - X == Y
    @test X - X == X
    @test 1.0 * X == X
    @test 0.1 * X == X

    M = ManifoldsBase.DefaultManifold(3)
    p = [0.0, 1.0, 2.0]
    @test copy(M, p, X) == X
    zero_vector(M, p, false) == X
end
