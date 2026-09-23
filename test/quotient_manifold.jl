using ManifoldsBase, Test
using ManifoldsBase.Test: DummyQuotientManifold, DummyTotalSpace

@testset "Allocations on a dummy quotient manifold" begin
    M = DummyQuotientManifold()
    p = [1.0, 2.0]
    X = [3.0, 4.0]
    q = canonical_project(M, p)
    @test q == p
    Y = diff_canonical_project(M, p, X)
    @test Y == X
    @test get_total_space(M) == DummyTotalSpace()
    Yh = horizontal_component(M, p, X)
    @test Yh == X
    Yl = horizontal_lift(M, p, X)
    @test Yl == X
    # Since all of the above are the identity, in this dummy v_space is 0
    Yv = vertical_component(M, p, X)
    @test Yv == zeros(2)
    Yv2 = copy(X)
    vertical_component!(M, Yv2, p, X)
    @test Yv2 == zeros(2)
end
