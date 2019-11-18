using ManifoldsBase
using LinearAlgebra

@testset "Array manifold" begin
    M = ManifoldsBase.DefaultManifold(3)
    A = ArrayManifold(M)
    x = [1., 0., 0.]
    y = 1/sqrt(2)*[1., 1., 0.]
    z = [0., 1., 0.]
    v = log(M,x,y)
    x2 = ArrayMPoint(x)
    y2 = ArrayMPoint(y)
    v2 = log(A,x,y) # auto convert
    y2 = exp(A,x,v2)
    w = log(M,x,z)
    w2 = log(A,x,z; atol=10^(-15))
    @test isapprox(y2.value,y)
    @test distance(A,x,y) == distance(M,x,y)
    @test norm(A,x,v) == norm(M,x,v)
    @test inner(A,x,v2,w2; atol=10^(-15)) == inner(M,x,v,w)
    @test isapprox(A, x2, y2 ) == isapprox(M,x,y)
    @test isapprox(A,x,y) == isapprox(A,x2,y2)
    @test isapprox(A,x, v2,v2 ) == isapprox(M,x,v,v)
    @test isapprox(A, exp(A,x,v),y2)
    @test isapprox(A, zero_tangent_vector(A,x), zero_tangent_vector(M,x))
end
