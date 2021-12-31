using Test
using ManifoldsBase

using ManifoldsBase: AbstractTrait, NestedTrait, EmptyTrait, trait, merge_traits
using ManifoldsBase: expand_trait, next_trait
import ManifoldsBase: active_traits, parent_trait

struct IsCool <: AbstractTrait end
struct IsNice <: AbstractTrait end

abstract type AbstractA end
abstract type DecoA <: AbstractA end
# A few concrete types
struct A <: AbstractA end
struct A1 <: DecoA end
struct A2 <: DecoA end
struct A3 <: DecoA end
struct A4 <: DecoA end
struct A5 <: DecoA end

# just some function
f(::AbstractA, x) = x + 2
g(::DecoA, x, y) = x + y

struct IsGreat <: AbstractTrait end # a special case of IsNice
parent_trait(::IsGreat) = IsNice()

active_traits(::A1, ::Any) = merge_traits(IsNice())
active_traits(::A2, ::Any) = merge_traits(IsCool())
active_traits(::A3, ::Any) = merge_traits(IsCool(), IsNice())
active_traits(::A5, ::Any) = merge_traits(IsGreat())

f(a::DecoA, b) = f(trait(a, b), a, b)

f(::NestedTrait{IsNice}, a, b) = g(a, b, 3)
f(::NestedTrait{IsCool}, a, b) = g(a, b, 5)

# generic forward to the next trait to be looked at
f(t::NestedTrait, a, b) = f(next_trait(t), a, b)
# generic fallback when no traits are defined
f(::EmptyTrait, a, b) = invoke(f, Tuple{AbstractA,typeof(b)}, a, b)

@testset "Decorator trait tests" begin
    t = ManifoldsBase.EmptyTrait()
    t2 = ManifoldsBase.NestedTrait(t, t)
    @test merge_traits() == t
    @test merge_traits(t) == t
    @test merge_traits(t, t) == t
    @test merge_traits(t2) == t2
    @test merge_traits(
        merge_traits(IsGreat(), IsNice()),
        merge_traits(IsGreat(), IsNice()),
    ) === merge_traits(IsGreat(), IsNice(), IsGreat(), IsNice())
    @test expand_trait(merge_traits(IsGreat(), IsCool())) ===
          merge_traits(IsGreat(), IsNice(), IsCool())
    @test expand_trait(merge_traits(IsCool(), IsGreat())) ===
          merge_traits(IsCool(), IsGreat(), IsNice())

    @test string(merge_traits(IsGreat(), IsNice())) ==
          "NestedTrait(IsGreat(), NestedTrait(IsNice(), EmptyTrait()))"

    global f
    @test f(A(), 0) == 2
    @test f(A2(), 0) == 5
    @test f(A3(), 0) == 5
    @test f(A4(), 0) == 2
    @test f(A5(), 890) == 893
    f(::NestedTrait{IsGreat}, a, b) = g(a, b, 54)
    @test f(A5(), 890) == 944
end

#
# A Manifold decorator test - check that EmptyTrait cases call Abstract and those fail with
# MethodError due to ambiguities (between Abstract and Decorator)
struct NonDecoratorManifold <: AbstractDecoratorManifold{ManifoldsBase.ℝ} end
ManifoldsBase.representation_size(::NonDecoratorManifold) = (2,)

@testset "Testing a NonDecoratorManifold - emptytrait fallbacks" begin
    M = NonDecoratorManifold()
    p = [2.0, 1.0]
    q = similar(p)
    X = [1.0, 0.0]
    Y = similar(X)
    @test check_size(M, p) === nothing
    # Ambiguous since not implemented
    @test_throws MethodError embed(M, p)
    @test_throws MethodError embed!(M, q, p)
    @test_throws MethodError embed(M, p, X)
    @test_throws MethodError embed!(M, Y, p, X)
    # the following is implemented but passes to the second and hence fails
    @test_throws MethodError exp(M, p, X)
    @test_throws MethodError exp!(M, q, p, X)
    @test_throws MethodError retract(M, p, X)
    @test_throws MethodError retract!(M, q, p, X)
    @test_throws MethodError log(M, p, q)
    @test_throws MethodError log!(M, Y, p, q)
    @test_throws MethodError inverse_retract(M, p, q)
    @test_throws MethodError inverse_retract!(M, Y, p, q)
    @test_throws MethodError vector_transport_along(M, p, X, :curve)
    @test_throws MethodError vector_transport_along!(M, Y, p, X, :curve)
end

# With even less, check that representation size stacjk overflows
struct NonDecoratorNonManifold <: AbstractDecoratorManifold{ManifoldsBase.ℝ} end
@testset "Testing a NonDecoratorNonManifold - emptytrait fallback Errors" begin
    N = NonDecoratorNonManifold()
    @test_throws StackOverflowError representation_size(N)
end
