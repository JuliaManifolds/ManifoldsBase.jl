using ManifoldsBase
using Test
using ManifoldsBase

struct TestDecorator{M<:AbstractManifold{ℝ}} <: AbstractDecoratorManifold{ℝ}
    manifold::M
end
