# Manifolds

While the interface `ManifoldsBase.jl` does not cover concrete manifolds, it provides a few
helpers to build or create manifolds based on existing manifolds


## Abstract Power Manifold

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

## [EmbeddedManifold](@id EmbeddedmanifoldSec)

Some manifolds can easily be defined by using a certain embedding.
For example the [`Sphere`](@ref)`(n)` is embedded in [`Euclidean`](@ref)`(n+1)`.
Similar to the metric and [`MetricManifold`](@ref), an embedding is often implicitly assumed.
We introduce the embedded manifolds hence as an [`AbstractDecoratorManifold`](@ref).

This decorator enables to use such an embedding in an transparent way.
Different types of embeddings can be distinguished using the [`AbstractEmbeddingType`](@ref),
which is an [`AbstractDecoratorType`](@ref).

### Isometric Embeddings

For isometric embeddings the type [`AbstractIsometricEmbeddingType`](@ref) can be used to avoid reimplementing the metric.
See [`Sphere`](@ref) or [`Hyperbolic`](@ref) for example.
Here, the exponential map, the logarithmic map, the retraction and its inverse
are set to `:intransparent`, i.e. they have to be implemented.

Furthermore, the [`TransparentIsometricEmbedding`](@ref) type even states that the exponential
and logarithmic maps as well as retractions and vector transports of the embedding can be
used for the embedded manifold as well.
See [`SymmetricMatrices`](@ref) for an example.

In both cases of course [`check_point`](@ref) and [`check_vector`](@ref) have to be implemented.

### Further Embeddings

A first embedding can also just be given implementing [`embed!`](@ref) ann [`project!`](@ref)
for a manifold. This is considered to be the most usual or default embedding.

If you have two different embeddings for your manifold, a second one can be specified using
the [`EmbeddedManifold`](@ref), a type that “couples” two manifolds, more precisely a
manifold and its embedding, to define embedding and projection functions between these
two manifolds.

### Types

```@autodocs
Modules = [ManifoldsBase]
Pages = ["EmbeddedManifold.jl"]
Order = [:type]
```

### Functions

```@autodocs
Modules = [ManifoldsBase]
Pages = ["EmbeddedManifold.jl"]
Order = [:function]
```

## DefaultManifold

[`DefaultManifold`](@ref ManifoldsBase.DefaultManifold) is a simplified version of [`Euclidean`](@ref) and demonstrates a basic interface implementation.
It can be used to perform simple tests.
Since when using `Manifolds.jl` the [`Euclidean`](@ref) is available, the `DefaultManifold` itself is not exported.

```@docs
ManifoldsBase.DefaultManifold
```

## Error Messages

especially to collect and display errors on [`AbstractPowerManifold`](@ref)s the following
component and collection error messages are available.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["errors.jl"]
Order = [:type]
```
