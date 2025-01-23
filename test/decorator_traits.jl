using Test
using ManifoldsBase

using ManifoldsBase: AbstractTrait, TraitList, EmptyTrait, trait, merge_traits
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

active_traits(f, ::A1, ::Any) = merge_traits(IsNice())
active_traits(f, ::A2, ::Any) = merge_traits(IsCool())
active_traits(f, ::A3, ::Any) = merge_traits(IsCool(), IsNice())
active_traits(f, ::A5, ::Any) = merge_traits(IsGreat())

f(a::DecoA, b) = f(trait(f, a, b), a, b)

f(::TraitList{IsNice}, a, b) = g(a, b, 3)
f(::TraitList{IsCool}, a, b) = g(a, b, 5)

# generic forward to the next trait to be looked at
f(t::TraitList, a, b) = f(next_trait(t), a, b)
# generic fallback when no traits are defined
f(::EmptyTrait, a, b) = invoke(f, Tuple{AbstractA,typeof(b)}, a, b)

@testset "Decorator trait tests" begin
    t = ManifoldsBase.EmptyTrait()
    t2 = ManifoldsBase.TraitList(t, t)
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
          "TraitList(IsGreat(), TraitList(IsNice(), EmptyTrait()))"

    global f
    @test f(A(), 0) == 2
    @test f(A2(), 0) == 5
    @test f(A3(), 0) == 5
    @test f(A4(), 0) == 2
    @test f(A5(), 890) == 893
    f(::TraitList{IsGreat}, a, b) = g(a, b, 54)
    @test f(A5(), 890) == 944

    @test next_trait(EmptyTrait()) === EmptyTrait()
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
    @test ManifoldsBase.check_size(M, p) === nothing
    @test ManifoldsBase.check_size(M, p, X) === nothing
    # default to identity
    @test embed(M, p) == p
    @test embed!(M, q, p) == p
    @test embed(M, p, X) == X
    @test embed!(M, Y, p, X) == X
    # the following is implemented but passes to the second and hence fails
    @test_throws MethodError exp(M, p, X)
    @test_throws MethodError ManifoldsBase.ManifoldsBase.ManifoldsBase.expt(M, p, X, 2.0)
    @test_throws MethodError exp!(M, q, p, X)
    @test_throws MethodError ManifoldsBase.expt!(M, q, p, X, 2.0)
    @test_throws MethodError retract(M, p, X)
    @test_throws MethodError retract(M, p, X, 2.0)
    @test_throws MethodError retract!(M, q, p, X)
    @test_throws MethodError retract!(M, q, p, X, 2.0)
    @test_throws MethodError log(M, p, q)
    @test_throws MethodError log!(M, Y, p, q)
    @test_throws MethodError inverse_retract(M, p, q)
    @test_throws MethodError inverse_retract!(M, Y, p, q)
    @test_throws MethodError parallel_transport_along(M, p, X, :curve)
    @test_throws MethodError parallel_transport_along!(M, Y, p, X, :curve)
    @test_throws MethodError parallel_transport_direction(M, p, X, X)
    @test_throws MethodError parallel_transport_direction!(M, Y, p, X, X)
    @test_throws MethodError parallel_transport_to(M, p, X, q)
    @test_throws MethodError parallel_transport_to!(M, Y, p, X, q)
    @test_throws MethodError vector_transport_along(M, p, X, :curve)
    @test_throws MethodError vector_transport_along!(M, Y, p, X, :curve)
end

# With even less, check that representation size stack overflows
struct NonDecoratorNonManifold <: AbstractDecoratorManifold{ManifoldsBase.ℝ} end
@testset "Testing a NonDecoratorNonManifold - emptytrait fallback Errors" begin
    N = NonDecoratorNonManifold()
    @test_throws StackOverflowError representation_size(N)
end


h(::AbstractA, x::Float64, y; a = 1) = x + y - a
h(::DecoA, x, y) = x + y
ManifoldsBase.@invoke_maker 1 AbstractA h(A::DecoA, x::Float64, y)

@testset "@invoke_maker" begin
    @test h(A1(), 1.0, 2) == 2
    sig = ManifoldsBase._split_signature(:(fname(x::T; k::Float64 = 1) where {T}))
    @test sig.fname == :fname
    @test sig.where_exprs == Any[:T]
    @test sig.callargs == Any[:(x::T)]
    @test sig.kwargs_list == Any[:($(Expr(:kw, :(k::Float64), 1)))]
    @test sig.argnames == [:x]
    @test sig.argtypes == [:T]
    @test sig.kwargs_call == Expr[:(k = k)]
    @test_throws ErrorException ManifoldsBase._split_signature(:(a = b))

    sig2 = ManifoldsBase._split_signature(:(fname(x; kwargs...)))
    @test sig2.kwargs_call == Expr[:(kwargs...)]

    sig3 = ManifoldsBase._split_signature(:(fname(x::T, y::Int = 10; k1 = 1) where {T}))
    @test sig3.kwargs_call == Expr[:(k1 = k1)]
    @test sig3.callargs[2] == :($(Expr(:kw, :(y::Int), 10)))
    @test sig3.argnames == [:x, :y]
end
