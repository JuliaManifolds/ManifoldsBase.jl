using ManifoldsBase, Test

struct DummyQuotientManifold <: AbstractManifold{ℝ} end
struct DummyTotalSpace <: AbstractManifold{ℝ} end

ManifoldsBase.canonical_project!(M::DummyQuotientManifold, q, p) = copyto!(q, p)
ManifoldsBase.diff_canonical_project!(M::DummyQuotientManifold, Y, p, X) = copyto!(Y, X)
ManifoldsBase.get_total_space(::DummyQuotientManifold) = DummyTotalSpace()
ManifoldsBase.horizontal_component!(N::DummyQuotientManifold, Y, p, X) = copyto!(Y, X)
ManifoldsBase.horizontal_lift!(N::DummyQuotientManifold, Y, q, X) = copyto!(Y, X)
ManifoldsBase.zero_vector(::DummyQuotientManifold, p) = zeros(2)
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
end
