using ManifoldsBase
using Test

import ManifoldsBase.decorator_transparent_dispatch

struct TestDecorator{M<:Manifold} <: AbstractDecoratorManifold
    manifold::M
end

test1(M::Manifold, p; a = 0) = 101 + a
test2(M::Manifold, p; a = 0) = 102 + a
test3(M::Manifold, p; a = 0) = 103 + a
function test4(M::Manifold, p; a = 0)
    error(ManifoldsBase.manifold_function_not_implemented_message(M, test4, p))
end

function test1(M::TestDecorator, p; a = 0)
    return 1 + a
end

decorator_transparent_dispatch(::typeof(test1), M::TestDecorator, args...) = Val(:intransparent)
decorator_transparent_dispatch(::typeof(test2), M::TestDecorator, args...) = Val(:transparent)
decorator_transparent_dispatch(::typeof(test3), M::TestDecorator, args...) = Val(:parent)
decorator_transparent_dispatch(::typeof(test4), M::TestDecorator, args...) = Val(:intransparent)

@testset "Testing decorator manifold functions" begin
    M = ManifoldsBase.DefaultManifold(3)
    A = ArrayManifold(M)

    @test (@inferred base_manifold(M)) == M
    @test (@inferred base_manifold(A)) == M
    @test ManifoldsBase._extract_val(Val(:transparent)) === :transparent

    @test (@inferred base_manifold(M, Val(1))) == M
    @test (@inferred base_manifold(M, Val(0))) == M
    @test (@inferred base_manifold(A, Val(1))) == M
    @test (@inferred base_manifold(A, Val(0))) == A

    @test representation_size(M) == (3,)
    @test representation_size(A) == (3,)

    @test manifold_dimension(M) == 3
    @test manifold_dimension(A) == 3

    p = [1.0, 0.0, 0.0]
    X = [2.0, 1.0, 3.0]
    @test inner(A, p, X, X) ≈ inner(A, Val(:transparent), p, X, X)
    @test_throws ErrorException inner(A, Val(:intransparent), p, X, X)

    TD = TestDecorator(M)

    @test test1(TD, p) == 1
    @test test1(TD, p; a = 1000) == 1001
    @test test2(TD, p) == 102
    @test test2(TD, p; a = 1000) == 1102
    @test test3(TD, p) == 103
    @test test3(TD, p; a = 1000) == 1103
    @test_throws ErrorException test4(TD, p)
    @test_throws ErrorException test4(TD, p; a = 1000)
end
