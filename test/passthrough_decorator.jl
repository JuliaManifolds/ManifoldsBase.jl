
using ManifoldsBase
using Test

using ManifoldsBase: NestedTrait, merge_traits

struct PassthoughTrait <: AbstractTrait end

struct PassthroughDecorator{𝔽,MT<:AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}
    manifold::MT
end

ManifoldsBase.active_traits(::PassthroughDecorator, ::Any...) = merge_traits(PassthoughTrait())

function ManifoldsBase.log!(
    ::NestedTrait{PassthoughTrait},
    M::AbstractDecoratorManifold,
    X,
    p,
    q,
)
    return log!(M.manifold, X, p, q)
end
function ManifoldsBase.exp!(
    ::NestedTrait{PassthoughTrait},
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


