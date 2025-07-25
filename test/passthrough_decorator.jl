using ManifoldsBase
using Test
using Random

using ManifoldsBase: TraitList, merge_traits

struct PassthoughTrait <: AbstractTrait end

struct PassthroughDecorator{𝔽, MT <: AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}
    manifold::MT
end

function ManifoldsBase.active_traits(f, ::PassthroughDecorator, ::Any...)
    return merge_traits(PassthoughTrait())
end
function ManifoldsBase.active_traits(f, ::AbstractRNG, ::PassthroughDecorator, ::Any...)
    return merge_traits(PassthoughTrait())
end

function ManifoldsBase.log!(
        ::TraitList{PassthoughTrait},
        M::AbstractDecoratorManifold,
        X,
        p,
        q,
    )
    return log!(M.manifold, X, p, q)
end
function ManifoldsBase.exp!(
        ::TraitList{PassthoughTrait},
        M::AbstractDecoratorManifold,
        q,
        p,
        X,
    )
    return exp!(M.manifold, q, p, X)
end

function ManifoldsBase.rand(::TraitList{PassthoughTrait}, M::AbstractDecoratorManifold)
    return rand(M.manifold)
end
function ManifoldsBase.rand!(::TraitList{PassthoughTrait}, M::AbstractDecoratorManifold, p)
    return rand!(M.manifold, p)
end
function ManifoldsBase.rand(
        ::TraitList{PassthoughTrait},
        rng::AbstractRNG,
        M::AbstractDecoratorManifold,
    )
    return rand(rng, M.manifold)
end
function ManifoldsBase.rand!(
        ::TraitList{PassthoughTrait},
        rng::AbstractRNG,
        M::AbstractDecoratorManifold,
        p,
    )
    return rand!(rng, M.manifold, p)
end

@testset "PassthroughDecorator" begin
    M = PassthroughDecorator(ManifoldsBase.DefaultManifold(2))

    q = [0.0, 0.0]
    p = [0.0, 0.0]
    X = [1.0, 2.0]
    Y = [0.0, 0.0]
    @test inverse_retract!(M, q, p, X) == [1.0, 2.0]
    @test retract(M, p, X) == [1.0, 2.0]
    @test ManifoldsBase.retract_fused(M, p, X, 1.0) == [1.0, 2.0]
    @test retract!(M, Y, p, X) == [1.0, 2.0]
    @test ManifoldsBase.retract_fused!(M, Y, p, X, 1.0) == [1.0, 2.0]
    @test rand(M) isa Vector{Float64}
    p2 = similar(p)
    rand!(M, p2)
    @test p2 isa Vector{Float64}
    @test rand(MersenneTwister(123), M) isa Vector{Float64}
    p2 = similar(p)
    rand!(MersenneTwister(123), M, p2)
    @test p2 == rand(MersenneTwister(123), M)
end
