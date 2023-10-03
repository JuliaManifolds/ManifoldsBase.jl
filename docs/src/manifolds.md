# Manifolds

While the interface `ManifoldsBase.jl` does not cover concrete manifolds, it provides a few
helpers to build or create manifolds based on existing manifolds

## A default manifold

[`DefaultManifold`](@ref ManifoldsBase.DefaultManifold) is a simplified version of [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) and demonstrates a basic interface implementation.
It can be used to perform simple tests.
Since when using `Manifolds.jl` the [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) is available, the `DefaultManifold` itself is not exported.

```@docs
ManifoldsBase.DefaultManifold
```

## [Embedded manifold](@id sec-embedded-manifold)

The embedded manifold is a manifold ``\mathcal M`` which is modelled _explicitly_ specifying its embedding ``\mathcal N`` in which the points and tangent vectors are represented.
Most prominently [`is_point`](@ref) and [`is_vector`](@ref) of an embedded manifold are implemented to check whether the point is a valid point in the embedding. This can of course still be extended by further tests.
`ManifoldsBase.jl` provides two possibilities of easily introducing this in order to dispatch some functions to the embedding.

### [Implicit case: the `IsEmbeddedManifold` Trait](@id subsec-implicit-embedded)

For the implicit case, your manifold has to be a subtype of the [`AbstractDecoratorManifold`](@ref).
Adding a method to the [`active_traits`](@ref ManifoldsBase.active_traits) function for a manifold that returns an [`AbstractTrait`](@ref)
[`IsEmbeddedManifold`](@ref), makes that manifold an embedded manifold. You just have to also define [`get_embedding`](@ref) so that appropriate functions are passed on to that embedding.
This is the implicit case, since the manifold type itself does not carry any information about the embedding, just the trait and the function definition do.

### [Explicit case: the `EmbeddedManifold`](@id subsec-explicit-embedded)

The [`EmbeddedManifold`](@ref) itself is an [`AbstractDecoratorManifold`](@ref) so it is a case of the implicit embedding itself, but internally stores both the original manifold and the embedding.
They are also parameters of the type.
This way, an additional embedding of one manifold in another can be modelled. That is, if the manifold is implemented using the implicit embedding approach from before but can also be implemented using a _different_ embedding, then this method should be chosen, since you can dispatch functions that you want to implement in this embedding then on the type which explicitly has the manifold and its embedding as parameters.

Hence this case should be used for any further embedding after the first or if the default implementation works without an embedding and the alternative needs one.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["EmbeddedManifold.jl"]
Order = [:type, :macro, :function]
```

## Metrics

Most metric-related functionality is currently defined in [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/latest/) but a few basic types are defined here.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["metric.jl"]
Order = [:type, :function]
```


## A manifold for validation

[`ValidationManifold`](@ref) is a simple decorator using the [`AbstractDecoratorManifold`](@ref) that “decorates” a manifold with tests that all involved points and vectors are valid for the wrapped manifold.
For example involved input and output paratemers are checked before and after running a function, repectively.
This is done by calling [`is_point`](@ref) or [`is_vector`](@ref) whenever applicable.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["ValidationManifold.jl"]
Order = [:macro, :type, :function]
```
