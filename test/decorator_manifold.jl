using ManifoldsBase
using Test

import ManifoldsBase.decorator_transparent_dispatch
using ManifoldsBase: @decorator_transparent_function,
    @decorator_transparent_fallback,
    @decorator_transparent_signature,
    is_decorator_transparent

struct TestDecorator{M<:Manifold} <: AbstractDecoratorManifold
    manifold::M
end

abstract type AbstractTestDecorator <: AbstractDecoratorManifold end

struct TestDecorator2{M<:Manifold} <: AbstractTestDecorator
    manifold::M
end

struct TestDecorator3{M<:Manifold} <: AbstractTestDecorator
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

@decorator_transparent_function :transparent function test5(M::AbstractDecoratorManifold, p)
    return 5
end

@decorator_transparent_function @inline function test6(M::TestDecorator, p)
    return 6
end

@decorator_transparent_function :parent function test7(M::TestDecorator, p)
    return 7
end

@decorator_transparent_fallback :parent @inline function test7(M::TestDecorator, p)
    return 17
end

test8(M::Manifold, p; a = 0) = 8 + a

@decorator_transparent_function :parent function test9(M::AbstractDecoratorManifold, p; a = 0, kwargs...)
    return 9 + a + (haskey(kwargs, :b) ? kwargs[:b] : 0)
end

@decorator_transparent_fallback :parent @inline function test9(M::AbstractTestDecorator, p::TP; a = 0, kwargs...) where {TP}
    return 19 + a + (haskey(kwargs, :b) ? kwargs[:b] : 0)
end

function test9(M::TestDecorator3, p::TP; a = 0, kwargs...) where {TP}
    return 109 + a + (haskey(kwargs, :b) ? kwargs[:b] : 0)
end

test10(M::AbstractTestDecorator, p::TP; a=0) where {TP} = 10*a
@decorator_transparent_function function test10(M::TestDecorator3, p::TP; a=0) where {TP}
    return 5*a
end
# the following then ignores the previous definition and passes again to the parent above
decorator_transparent_dispatch(::typeof(test10), M::TestDecorator3, args...) = Val(:parent)

@decorator_transparent_function function test11(M::TestDecorator3, p::TP; a::Int=0) where {TP}
    return 15*a
end

@testset "Testing decorator manifold functions" begin
    M = ManifoldsBase.DefaultManifold(3)
    A = ArrayManifold(M)

    @test (@inferred base_manifold(M)) == M
    @test (@inferred base_manifold(A)) == M
    @test (@inferred ManifoldsBase.decorated_manifold(A)) == M
    @test ManifoldsBase._extract_val(Val(:transparent)) === :transparent

    @test (@inferred base_manifold(M, Val(1))) == M
    @test (@inferred base_manifold(M, Val(0))) == M
    @test (@inferred base_manifold(A, Val(1))) == M
    @test (@inferred base_manifold(A, Val(0))) == A

    x = 0
    @test_throws LoadError eval(:(@decorator_transparent_fallback x = x+1))
    @test_throws LoadError eval(:(@decorator_transparent_function x = x+1))
    @test_throws LoadError eval(:(@decorator_transparent_signature x = x+1))

    @test representation_size(M) == (3,)
    @test representation_size(A) == (3,)

    @test manifold_dimension(M) == 3
    @test manifold_dimension(A) == 3

    p = [1.0, 0.0, 0.0]
    X = [2.0, 1.0, 3.0]
    @test inner(A, p, X, X) ≈ ManifoldsBase.inner__transparent(A, p, X, X)
    @test_throws ErrorException ManifoldsBase.inner__intransparent(A, p, X, X)

    TD = TestDecorator(M)

    @test (@inferred ManifoldsBase.default_decorator_dispatch(M)) === Val(false)
    @test ManifoldsBase.is_default_decorator(M) === false

    @test test1(TD, p) == 1
    @test test1(TD, p; a = 1000) == 1001
    @test test2(TD, p) == 102
    @test test2(TD, p; a = 1000) == 1102
    @test test3(TD, p) == 103
    @test test3(TD, p; a = 1000) == 1103
    @test_throws ErrorException test4(TD, p)
    @test_throws ErrorException test4(TD, p; a = 1000)
    @test (@inferred decorator_transparent_dispatch(test5, TD, p)) === Val(:transparent)
    @test is_decorator_transparent(test5, TD, p)
    @test test5(TD, p) == 5
    @test (@inferred decorator_transparent_dispatch(test6, TD, p)) === Val(:intransparent)
    @test_throws ErrorException test7(M, p)
    @test test7(TD, p) == 17
    @test (@inferred decorator_transparent_dispatch(test8, M, p)) === Val(:transparent)
    @test is_decorator_transparent(test8, M, p)
    @test_throws ErrorException test9(M, p; a = 1000)
    @test test9(TD, p; a = 1000) == 1009
    @test test9(TD, p; a = 1000, b = 10000) == 11009
    @test test9(TestDecorator2(TD), p; a = 1000) == 1019
    @test test9(TestDecorator2(TD), p; a = 1000, b = 10000) == 11019
    @test test9(TestDecorator3(TestDecorator2(TD)), p; a = 1000) == 1109
    @test test9(TestDecorator3(TestDecorator2(TD)), p; a = 1000, b = 10000) == 11109
    @test test9(TestDecorator3(TD), p; a = 1000) == 1109
    @test test9(TestDecorator3(TD), p; a = 1000, b = 10000) == 11109
    @test test10(TestDecorator3(TD), p; a = 11) == 110
    @test test11(TestDecorator3(TD), p; a = 12) == 180
end
