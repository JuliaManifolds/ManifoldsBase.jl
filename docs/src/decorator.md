# A Decorator for manifolds

Several properties of a manifold are often implicitly assumed, for example the choice of the (Riemannian) metric, the group structure or the embedding. The latter shall serve as an example how to either implicitly or explicitly specify the embedding to avoid re-implementations and/or distinguish different embeddings.

## The abstract decorator

When first implementing a manifold, it might be beneficial to dispatch certain computations to already existing manifolds.
For an embedded manifold that is isometrically embedded this might be the [`inner`](@ref) the manifold inherits in each tangent space from its embedding.

This means we would like to dispatch the default implementation of a function to some other manifold.
We refer to this as implicit decoration, since one can not “see” explicitly that a certain manifold inherits this property.
As an example consider the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html). At each point the tangent space can be identified with a subspace of the tangent space in the embedding, the [Euclidean](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) manifold which the unit vectors of the sphere belong to. Thus every tangent space inherits its metric from the embedding.
Since in the default implementation in [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) points are represented by unit vectors and tangent vectors at a point as vectors orthogonal to that point, we can just dispatch the inner product to the embedding without having to re-implement this.
The manifold using such an implicit dispatch just has to be a subtype of [`AbstractDecoratorManifold`](@ref).

## Traits with an inheritance hierarchy

The properties mentioned above might form a hierarchy.
For embedded manifolds, again, we might have just a manifold whose points are represented in some embedding.
If the manifold is even isometrically embedded, it is embedded but also inherits the Riemannian metric by restricting the metric from the embedding to the corresponding tangent space under consideration.
But it also inherits the functions defined for the plain embedding, for example checking some conditions for the validity of points and vectors.
If it is even a submanifold, also further functions are inherited like the [`shortest_geodesic`](@ref).

We use a variation of [Tim Holy's Traits Trick](https://github.com/JuliaLang/julia/issues/2345#issuecomment-54537633) (THTT) which takes into account this nestedness of traits.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["nested_trait.jl"]
Order = [:type, :macro, :function]
```

## The Manifold decorator

The idea of the decorator for a manifold is to allow to exchange certain functionality by
a dispatch layer. For example that for an embedded manifold some functions are passed to the embedding,
or that for a metric manifold decorator, all functions unrelated to the metric are passed to the original manifold. The following types, functions, and macros introduce the decorator trait which allows to decorate an arbitrary `<: `[`AbstractDecoratorManifold`](@ref) with further features.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["decorator_trait.jl"]
Order = [:type, :macro, :function]
```

For an example see the [(implicit) embedded manifold](@ref subsec-implicit-embedded).