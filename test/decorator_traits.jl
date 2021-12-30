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
