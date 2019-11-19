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
    @testset "Types and Conversion" begin
        @test convert(typeof(M), A) == M
        @test convert(typeof(A),M) == A
        @test is_decorator_manifold(A) == Val(true)
        @test base_manifold(A) == M
        @test ManifoldsBase.representation_size(A) == ManifoldsBase.representation_size(M)
        @test manifold_dimension(A) == manifold_dimension(M)
        for T in [ArrayMPoint, ArrayTVector, ArrayCoTVector]
            p = T(x)
            @test convert(typeof(x),p) == x
            @test convert(typeof(p),y) == T(y)
            @test eltype(typeof(p)) == eltype(x)
            @test eltype(p) == eltype(x)
            @test typeof(similar(p)) == typeof(p)
            @test typeof(similar(p,eltype(x))) == typeof(p)
            q = similar(p)
            copyto!(q,p)
            @test isapprox(A,q,p)
            @test ManifoldsBase.array_value(p) == x
            @test ManifoldsBase.array_value(x) == x
        end
    end
    @testset "Vector functions" begin
        for T in [ArrayTVector, ArrayCoTVector]
            a = T(v)
            b = T(w)
            @test isapprox(A, a+b, T(v+w))
            @test isapprox(A, (a-b), T(v-w) )
            @test isapprox(A, -b, T(-w) )
            @test isapprox(A, 2*a, T(2 .* v) )
        end
    end
    @testset "Manifold functions" begin
        @test manifold_dimension(A) == manifold_dimension(M)
        @test isapprox(y2.value,y)
        @test distance(A,x,y) == distance(M,x,y)
        @test norm(A,x,v) == norm(M,x,v)
        @test inner(A,x,v2,w2; atol=10^(-15)) == inner(M,x,v,w)
        @test isapprox(A, x2, y2 ) == isapprox(M,x,y)
        @test isapprox(A,x,y) == isapprox(A,x2,y2)
        @test isapprox(A,x, v2,v2 ) == isapprox(M,x,v,v)
        v2s = similar(v2)
        project_tangent!(A,v2s,x2,v2)
        @test isapprox(A, v2, v2s)
        y2s = similar(y2)
        exp!(A,y2s,x2,v2)
        @test isapprox(A,y2s,y2)
        log!(A, v2s, x, y)
        @test isapprox(A, x, v2s, v2)
        @test isapprox(A, exp(A,x,v),y2)
        @test isapprox(A, zero_tangent_vector(A,x), zero_tangent_vector(M,x))
        vector_transport_to!(A, v2s, x2, v2, y2)
        @test isapprox(A, x2, v2, v2s)
        zero_tangent_vector!(A, v2s, x)
        @test isapprox(A, v2s, zero_tangent_vector(M,x))
        @test_throws ErrorException vector_transport_along!(A,v2s,x2,v2,ParallelTransport())
    end
end
