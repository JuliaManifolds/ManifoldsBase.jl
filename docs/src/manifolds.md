# Manifolds

While the interface `ManifoldsBase.jl` does not cover concrete manifolds, it provides a few
helpers to build or create manifolds based on existing manifolds

## (Abstract) Power Manifold

```@autodocs
Modules = [ManifoldsBase]
Pages = ["src/PowerManifold.jl"]
Order = [:macro, :type, :function]
```

## ValidationManifold

[`ValidationManifold`](@ref) is a simple decorator using the [`AbstractDecoratorManifold`](@ref) that “decorates” a manifold with tests that all involved points and vectors are valid for the wrapped manifold.
For example involved input and output paratemers are checked before and after running a function, repectively.
This is done by calling [`is_point`](@ref) or [`is_vector`](@ref) whenever applicable.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["ValidationManifold.jl"]
Order = [:macro, :type, :function]
```

## DefaultManifold

[`DefaultManifold`](@ref ManifoldsBase.DefaultManifold) is a simplified version of [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) and demonstrates a basic interface implementation.
It can be used to perform simple tests.
Since when using `Manifolds.jl` the [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) is available, the `DefaultManifold` itself is not exported.

```@docs
ManifoldsBase.DefaultManifold
```
