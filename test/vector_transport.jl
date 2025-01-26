#
# Test the specific vector trnasport along implementations that do iterative transport,
# also the Schild and pole special cases
#
using ManifoldsBase, Test
import ManifoldsBase: parallel_transport_to!

struct NonDefaultEuclidean <: AbstractManifold{ManifoldsBase.ℝ} end
ManifoldsBase.log!(::NonDefaultEuclidean, v, x, y) = (v .= y .- x)
ManifoldsBase.exp!(::NonDefaultEuclidean, y, x, v) = (y .= x .+ v)
function ManifoldsBase.parallel_transport_to!(::NonDefaultEuclidean, Y, p, X, q; a = 0)
    return copyto!(Y, X .+ a)
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
