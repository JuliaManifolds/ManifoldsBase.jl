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

abstract type AbstractParentDecorator <: AbstractDecoratorManifold end

struct ChildDecorator{M<:Manifold} <: AbstractParentDecorator
    manifold::M
end

struct DefaultDecorator{M<:Manifold} <: AbstractDecoratorManifold
    manifold::M
end
ManifoldsBase.default_decorator_dispatch(::DefaultDecorator) = Val(true)

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

@decorator_transparent_function function test12(M::ManifoldsBase.DefaultManifold, p)
    return 12*p
end
ManifoldsBase._acts_transparently(test12, TestDecorator3, p) = Val(:foo)

@decorator_transparent_function :none function test13(M::TestDecorator3, p)
    return 13.5*p
end
decorator_transparent_dispatch(::typeof(test13), M::TestDecorator, args...) = Val(:intransparent)
decorator_transparent_dispatch(::typeof(test13), M::TestDecorator2, args...) = Val(:transparent)
test13(::ManifoldsBase.DefaultManifold,a) = 13*a

function test14(M::AbstractDecoratorManifold, p)
    return 14.5*p
end
@decorator_transparent_signature test14(M::AbstractDecoratorManifold,p)
decorator_transparent_dispatch(::typeof(test14), M::TestDecorator3, args...) = Val(:none)
decorator_transparent_dispatch(::typeof(test14), M::TestDecorator, args...) = Val(:intransparent)
decorator_transparent_dispatch(::typeof(test14), M::TestDecorator2, args...) = Val(:transparent)
test14(::ManifoldsBase.DefaultManifold,a) = 14*a

test15(::ManifoldsBase.DefaultManifold,a) = 15.5*a
@decorator_transparent_function function test15(M::AbstractDecoratorManifold, p)
    error("Not yet implemented")
end
test15(::AbstractParentDecorator,p) = 15*p
decorator_transparent_dispatch(::typeof(test15), M::ChildDecorator, args...) = Val(:parent)

function test16(::AbstractParentDecorator, p)
    return 16*p
end
test16(::ManifoldsBase.DefaultManifold, a) = 16.5*a
@decorator_transparent_signature test16(M::AbstractDecoratorManifold, p)
decorator_transparent_dispatch(::typeof(test16), M::ChildDecorator, args...) = Val(:parent)

function test17(M::ManifoldsBase.DefaultManifold, p)
    return 17*p
end
@decorator_transparent_signature test17(M::AbstractDecoratorManifold, p)
decorator_transparent_dispatch(::typeof(test17), M::AbstractDecoratorManifold, args...) = Val(:intransparent)
default_decorator_dispatch(::DefaultDecorator) = Val(true)

@decorator_transparent_function function test18(M::AbstractDecoratorManifold, p)
    return 18.25*p
end
decorator_transparent_dispatch(::typeof(test18), M::ChildDecorator, args...) = Val(:parent)

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

    @test injectivity_radius(TD, ManifoldsBase.ExponentialRetraction()) == Inf

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
    @test_throws ErrorException test12(TestDecorator3(TD), p)

    @test_throws ErrorException test13(TestDecorator3(M),1) # :none nonexistent
    @test_throws ErrorException test13(TestDecorator(M),1) # not implemented
    @test test13(TestDecorator2(M),2) == 26 # from parent

    @test_throws ErrorException test14(TestDecorator3(M),1) # :none nonexistent
    @test_throws ErrorException test14(TestDecorator(M),1) # not implemented
    @test test14(TestDecorator2(M),2) == 28 # from parent

    @test test15(ChildDecorator(M),1) == 15
    @test test16(ChildDecorator(M),1) == 16

    @test_throws ErrorException test17(TestDecorator(ManifoldsBase.DefaultManifold(3)),1)
    @test test17(DefaultDecorator(ManifoldsBase.DefaultManifold(3)),1) == 17

    # states that child has to implement at least a parent case
    @test_throws ErrorException test18(ChildDecorator(ManifoldsBase.DefaultManifold(3)), 1)
end
