# A Decorator for manifolds

A decorator for a manifold is the main scheme to explicitly add or change a property to a manifold. This can be done by defining a new type that has the [`AbstractDecoratorManifold`](@ref) as a supertype.

In general there are two scenarios where this is used.

1. Dispatching a certain default property to another manifold.
2. Introducing an alternate variant of a property.

Let‘s make this more precise. For the first case imagine that a manifold ``\mathcal M``
is represented by points in the embedding by default.
If it is isometrically embedded, the inner product can be inherited from the embedding.
This is a case for the first point.
We provide a default coupling to another manifold to dispatch certain functions over.
Which ones are passed on is regulated by traits and their dispatch functions.

Assume that a manifold has two different embeddings or a default representation of its points/vectors works without an embedding and there is a second way to represent the manifold in some embedding.
This is the situation for the second scenario, see [`EmbeddedManifold`](@ref) to provide a concrete “coupling” between a manifold and its embedding.

implicitly or explicitly add properties _

```@autodocs
Modules = [ManifoldsBase]
Pages = ["decorator_trait.jl"]
Order = [:type, :macro, :function]
```
