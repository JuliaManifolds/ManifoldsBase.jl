#
# Test the specific vector trnasport along implementations that do iterative transport,
# also the Schild and pole special cases
#
using ManifoldsBase, Test


struct NonDefaultEuclidean <: Manifold{ManifoldsBase.ℝ} end
ManifoldsBase.log!(::NonDefaultEuclidean, v, x, y) = (v .= y .- x)
ManifoldsBase.exp!(::NonDefaultEuclidean, y, x, v) = (y .= x .+ v)
function ManifoldsBase.vector_transport_to!(
    ::NonDefaultEuclidean,
    vto,
    x,
    v,
    y,
    ::ParallelTransport,
)
    return copyto!(vto, v)
end

@testset "vector_transport_along" begin
    M = NonDefaultEuclidean();
    types = [
        Vector{Float64},
        Vector{Float32},
    ]
    for T in types
        @testset "Type $T" begin
            pts = convert.(Ref(T), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            v2 = log(M, pts[1], pts[3])
            c = [0.5 * (pts[1] + pts[2]), pts[2], 0.5 * (pts[2] + pts[3]), pts[3]]
            @test vector_transport_along(M, pts[1], v2, c, SchildsLadderTransport()) == v2
            @test vector_transport_along(M, pts[1], v2, c, PoleLadderTransport()) == v2
            @test vector_transport_along(M, pts[1], v2, c, ParallelTransport()) == v2
            # check mutating ones with defaults
            p = allocate(pts[1])
            ManifoldsBase.pole_ladder!(M, p, pts[1], pts[2], pts[3])
            # -log_p3 p == log_p1 p2
            @test isapprox(M, -log(M, pts[3], p), log(M, pts[1], pts[2]))
            ManifoldsBase.schilds_ladder!(M, p, pts[1], pts[2], pts[3])
            @test isapprox(M, log(M, pts[3], p), log(M, pts[1], pts[2]))
        end
    end
end