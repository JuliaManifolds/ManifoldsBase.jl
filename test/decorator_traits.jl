using Test
using ManifoldsBase

s = @__DIR__
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using ManifoldsBaseTestUtils

#
# A Manifold decorator test - check that StopForwarding cases call Abstract and those fail with
# MethodError due to ambiguities (between Abstract and Decorator)
struct NonDecoratorManifold <: AbstractDecoratorManifold{ManifoldsBase.ℝ} end
ManifoldsBase.representation_size(::NonDecoratorManifold) = (2,)

@testset "Testing a NonDecoratorManifold - StopForwarding fallbacks" begin
    M = NonDecoratorManifold()
    p = [2.0, 1.0]
    q = similar(p)
    X = [1.0, 0.0]
    Xc = [1.0, 0.0]
    Y = similar(X)
    Yc = similar(Xc)
    @test ManifoldsBase.check_size(M, p) === nothing
    @test ManifoldsBase.check_size(M, p, X) === nothing
    # default to identity
    @test embed(M, p) == p
    @test embed!(M, q, p) == p
    @test embed(M, p, X) == X
    @test embed!(M, Y, p, X) == X
    # the following is implemented but passes to the second and hence fails
    @test_throws MethodError exp(M, p, X)
    @test_throws MethodError ManifoldsBase.exp_fused(M, p, X, 2.0)
    @test_throws MethodError exp!(M, q, p, X)
    @test_throws MethodError ManifoldsBase.exp_fused!(M, q, p, X, 2.0)
    @test_throws MethodError retract(M, p, X)
    @test_throws MethodError ManifoldsBase.retract_fused(M, p, X, 2.0)
    @test_throws MethodError ManifoldsBase.retract_fused(
        M,
        p,
        X,
        2.0,
        ExponentialRetraction(),
    )
    @test_throws MethodError retract!(M, q, p, X)
    @test_throws MethodError ManifoldsBase.retract_fused!(M, q, p, X, 2.0)
    @test_throws MethodError ManifoldsBase.retract_fused!(
        M,
        q,
        p,
        X,
        2.0,
        ExponentialRetraction(),
    )
    @test_throws MethodError log(M, p, q)
    @test_throws MethodError log!(M, Y, p, q)
    @test_throws MethodError inverse_retract(M, p, q)
    @test_throws MethodError inverse_retract!(M, Y, p, q)
    @test_throws MethodError parallel_transport_direction(M, p, X, X)
    @test_throws MethodError parallel_transport_direction!(M, Y, p, X, X)
    @test_throws MethodError parallel_transport_to(M, p, X, q)
    @test_throws MethodError parallel_transport_to!(M, Y, p, X, q)
    # @test_throws StackOverflowError get_coordinates(M, p, X)
    @test_throws MethodError get_vector(M, p, Xc)
    # @test_throws StackOverflowError get_coordinates(M, p, X, DefaultOrthogonalBasis())
    @test_throws MethodError get_vector(M, p, Xc, DefaultOrthogonalBasis())
    @test_throws MethodError get_coordinates!(M, Yc, p, X)
    @test_throws MethodError get_vector!(M, Y, p, Xc)
    @test_throws MethodError get_coordinates!(M, Yc, p, X, DefaultOrthogonalBasis())
    @test_throws MethodError get_vector!(M, Y, p, Xc, DefaultOrthogonalBasis())

    @test_throws MethodError project(M, p)
    @test_throws MethodError project(M, p, X)
    @test_throws MethodError project!(M, q, p)
    @test_throws MethodError project!(M, Y, p, X)
end

# With even less, check that representation size stack overflows
struct NonDecoratorNonManifold <: AbstractDecoratorManifold{ManifoldsBase.ℝ} end
@testset "Testing a NonDecoratorNonManifold - StopForwarding fallback returns nothing" begin
    N = NonDecoratorNonManifold()
    @test representation_size(N) === nothing
end

@testset "Non-decorator manifold trait defaults" begin
    M = ProjManifold()
    p = [
        sqrt(2) / 2 0.0 0.0
        0.0 sqrt(2) / 2 0.0
    ]

    @test ManifoldsBase.get_forwarding_type(M, exp) === StopForwardingType()
    @test ManifoldsBase.get_forwarding_type(M, exp, p) === StopForwardingType()
end

abstract type AbstractA end
abstract type DecoA <: AbstractA end
# A few concrete types
struct A <: AbstractA end
struct A1 <: DecoA end

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
