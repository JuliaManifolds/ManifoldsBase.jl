
using ManifoldsBase
using Test

using ManifoldsBase: TraitList, merge_traits

struct PassthoughTrait <: AbstractTrait end

struct PassthroughDecorator{𝔽,MT<:AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}
    manifold::MT
end

function ManifoldsBase.active_traits(f, ::PassthroughDecorator, ::Any...)
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
    return log!(M.manifold, q, p, X)
end

@testset "PassthroughDecorator" begin
    M = PassthroughDecorator(ManifoldsBase.DefaultManifold(2))

    q = [0.0, 0.0]
    p = [0.0, 0.0]
    X = [1.0, 2.0]
    Y = [0.0, 0.0]
    @test inverse_retract!(M, q, p, X) == [1.0, 2.0]
    @test retract(M, p, q) == [1.0, 2.0]
    @test retract!(M, Y, p, q) == [1.0, 2.0]
end