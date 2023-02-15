#
# Test the specific vector trnasport along implementations that do iterative transport,
# also the Schild and pole special cases
#
using ManifoldsBase, Test
import ManifoldsBase: parallel_transport_to!, parallel_transport_along!

struct NonDefaultEuclidean <: AbstractManifold{ManifoldsBase.ℝ} end
ManifoldsBase.log!(::NonDefaultEuclidean, v, x, y) = (v .= y .- x)
ManifoldsBase.exp!(::NonDefaultEuclidean, y, x, v) = (y .= x .+ v)
function ManifoldsBase.parallel_transport_to!(::NonDefaultEuclidean, Y, p, X, q; a = 0)
    return copyto!(Y, X .+ a)
end
function ManifoldsBase.parallel_transport_along!(
    ::NonDefaultEuclidean,
    Y,
    p,
    X,
    c::AbstractVector,
)
    return copyto!(Y, X)
end

@testset "vector_transport_along" begin
    M = NonDefaultEuclidean()
    types = [Vector{Float64}, Vector{Float32}]
    for T in types
        @testset "Type $T" begin
            pts = convert.(Ref(T), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            X2 = log(M, pts[1], pts[3])
            c = [0.5 * (pts[1] + pts[2]), pts[2], 0.5 * (pts[2] + pts[3]), pts[3]]
            @test vector_transport_along(M, pts[1], X2, c, SchildsLadderTransport()) == X2
            @test vector_transport_along(M, pts[1], X2, c, PoleLadderTransport()) == X2
            @test vector_transport_along(M, pts[1], X2, c, ParallelTransport()) == X2
            @test vector_transport_along(M, pts[1], X2, [], SchildsLadderTransport()) == X2
            @test vector_transport_along(M, pts[1], X2, [], PoleLadderTransport()) == X2
            @test vector_transport_along(M, pts[1], X2, [], ParallelTransport()) == X2
            kwP = VectorTransportWithKeywords(ParallelTransport())
            @test vector_transport_along(M, pts[1], X2, [], kwP) == X2
            @test vector_transport_along!(M, X2, pts[1], X2, [], kwP) == X2
            # check in-place ones with defaults
            p = allocate(pts[1])
            ManifoldsBase.pole_ladder!(M, p, pts[1], pts[2], pts[3])
            # -log_p3 p == log_p1 p2
            @test isapprox(M, -log(M, pts[3], p), log(M, pts[1], pts[2]))
            ManifoldsBase.schilds_ladder!(M, p, pts[1], pts[2], pts[3])
            @test isapprox(M, log(M, pts[3], p), log(M, pts[1], pts[2]))
        end
    end
end

@testset "vector-transport fallback types" begin
    VT = VectorTransportDirection()
    M = NonDefaultEuclidean()
    p = [1.0, 0.0, 0.0]
    q = [0.0, 1.0, 0.0]
    X = [0.0, 0.0, 1.0]
    Y = similar(X)
    @test vector_transport_direction(M, p, X, p - q, VT) == X
    @test vector_transport_direction!(M, Y, p, X, p - q, VT) == X
    VT2 = VectorTransportTo()
    @test vector_transport_to(M, p, X, q, VT2) == X
    @test vector_transport_to!(M, Y, p, X, q, VT2) == X
    VT3 = VectorTransportWithKeywords(VT; a = 1)
    @test vector_transport_direction(M, p, X, q, VT3) == (X .+ 1)
    @test vector_transport_direction!(M, Y, p, X, q, VT3) == (X .+ 1)
    VT4 = VectorTransportWithKeywords(VT2; a = 2)
    @test vector_transport_to(M, p, X, q, VT4) == (X .+ 2)
    @test vector_transport_to!(M, Y, p, X, q, VT4) == (X .+ 2)
end
