# A Decorator for manifolds

Several properties of a manifold are often implicitly assumed, for example the choice of the (Riemannian) metric, the group structure or the embedding. The latter shall serve as an example how to either implicitly or explicitly specify the embedding to avoid re-implementations and/or distinguish different embeddings.

## The abstract decorator

When first implementing a manifold, it might be beneficial to dispatch certain computations to already existing manifolds.
For an embedded manifold that is isometrically embedded this might be the [`inner`](@ref) the manifold inherits in each tangent space from its embedding.

This means we dispatch the default implementation of a function to some other manifold.
We refer to this as implicit decoration, since one can not “see” explicitly that a certain manifold inherits this property.
As a small example consider the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html), which in every tangent space inherits its metric from the embedding. Since in the default implementation in [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) points are represented by unit vectors and tangent vectors as vectors orthogonal to a point, we can just dispatch the inner product to the embedding without having to re-implement this.
The manifold using such an implicit dispatch just has to have [`AbstractDecoratorManifold`](@ref) as its super type.

## Traits with a inheritance hierarchy

The properties mentioned above might form a hierarchy.
For embedded manifolds, again, we might have just a manifold whose points are represented in some embedding.
If the manifold is even isometrically embedded, it is embedded but also inherits the Riemannian metric (by restriction). But it also inherits the functions form the plain embedding.
If it is even a submanifold, also further functions are inherited.

We use a variation of [Tim Holy's Traits Trick](https://github.com/JuliaLang/julia/issues/2345#issuecomment-54537633) (THTT) which takes into account this nestedness of traits

```@autodocs
Modules = [ManifoldsBase]
Pages = ["nested_trait.jl"]
Order = [:type, :macro, :function]
```

Then the following functions and macros introduce the decorator traits

```@autodocs
Modules = [ManifoldsBase]
Pages = ["decorator_trait.jl"]
Order = [:type, :macro, :function]
```

For an example see the [(implicit) embedded manifold](@ref subsec-implicit-embedded).