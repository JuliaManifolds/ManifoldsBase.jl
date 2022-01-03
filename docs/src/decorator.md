# A Decorator for manifolds

Several properties of a manifold are often implicitly assumed, for example the choice of the (Riemannian) metric, the group structure or the embedding. The latter shall serve as an example how to either implicitly or explicitly specify the embedding to avoid reimplemenations and/or distinguish different embeddings.

## The abstract decorator

When first implementing a manifold, it might be beneficial to dispatch certain computations to already existing manifolds.
For an embedded manifold that is isometrically embedded this might be the [`inner`](@ref) the manifold inherits in each tangent space from its embedding.

This means we disptach the default implementation of a function to some other manifold.
We refer to this as implicit decoration, since one can not “see” explicitly that a certain manifold inherits this property.
As a small example consider the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html), which in every tangent space inherits its metric from the embedding. Since in the default implementation in [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) points are represented by unit vectors and tangent vectors as vectors orthogonal to a point, we can just dispatch the inner product to the embedding without having to re-implement this.
The manifold using such an implicit dispatch just has to have [`AbstractDecoratorManifold`](@ref) as its super type.

## Traits with a inheritance hierarchy

The properties mentioned above might form a hierarchy.
For embedded manifolds, again, we might have just a manifold whose points are represented in some embedding.
If the manifold is even isometrically embedded, it is embedded but also inherits the Riemannian metric (by restriction). But it also inherits the functions form the plain embedding.
If it is even a submanifold, also further functions are inherited.

We use a variation of Tim Holys Traits Trick (THTT) which takes into account this nestedness of traits

```@autodocs
Modules = [ManifoldsBase]
Pages = ["nested_trait.jl"]
Order = [:type, :macro, :function]
```

Then the following functions and macros introduce the decoraator traits

```@autodocs
Modules = [ManifoldsBase]
Pages = ["decorator_trait.jl"]
Order = [:type, :macro, :function]
```

## [Concrete decorators: The EmbeddedManifold](@id subsec-embeddedmanifold)

While the implicit decorator ddescribed until now provides a way to dispatch functions of a manifold to other manifold, there is also the case, that certain properties might have more than one variant they appear in. For example there might be different metrics or different embeddings.

These are modelled in `ManifoldsBase.jl` as own manifolds that follow the `AbstractDecorator` manifold, but they explicitly “couple” a manifold with this new (second occurence of a) property.

If a manifold ``\mathcal M`` is implemented as an embedded manifold as described above, but it can also be implemented using a different embedding (into another manifold ``\mathcal N``), then these two manifolds can be _explicitly_ coupled using the [`EmbeddedManifold`](@ref),
which itself is an [`AbstractDecoratorManifold`](@ref) and dispatches functions to either ``\mathcal M`` or ``\mathcal N`` again based on the type of embedding specified.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["EmbeddedManifold.jl"]
Order = [:type, :macro, :function]
```
