using ManifoldsBase
import ManifoldsBase: representation_size, manifold_dimension, inner
using LinearAlgebra
using Test

struct ComplexEuclidean{N} <: AbstractManifold{ℂ} where {N} end
ComplexEuclidean(n::Int) = ComplexEuclidean{n}()
representation_size(::ComplexEuclidean{N}) where {N} = (N,)
manifold_dimension(::ComplexEuclidean{N}) where {N} = 2N
@inline inner(::ComplexEuclidean, x, v, w) = dot(v, w)

@testset "Complex Euclidean" begin
    M = ComplexEuclidean(3)

    x = complex.([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    v1 = complex.([0.1, 0.2, 0.3], [0.4, 0.5, 0.6])
    v2 = complex.([0.3, 0.2, 0.1], [0.6, 0.5, 0.4])

    @test norm(M, x, v1) isa Real
    @test norm(M, x, v1) ≈ norm(v1)
    @test angle(M, x, v1, v2) isa Real
    @test angle(M, x, v1, v2) ≈ acos(real(dot(v1, v2)) / norm(v1) / norm(v2))

    vv1 = vcat(reim(v1)...)
    vv2 = vcat(reim(v2)...)
    @test angle(M, x, v1, v2) ≈ acos(dot(normalize(vv1), normalize(vv2)))
end
