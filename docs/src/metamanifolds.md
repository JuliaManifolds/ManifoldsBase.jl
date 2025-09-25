# Meta Manifolds

While the interface does not provide concrete manifolds itself, it does provide several manifolds that can be build based on a given [`AbstractManifold`](@ref) instance.

## [(Abstract) power manifold](@id sec-power-manifold)

A power manifold is constructed like higher dimensional vector spaces are formed from the real line, just that for every point ``p = (p_1,\ldots,p_n) ∈ \mathcal M^n`` on the power manifold ``\mathcal M^n`` the entries of ``p`` are points ``p_1,\ldots,p_n ∈ \mathcal M`` on some manifold ``\mathcal M``. Note that ``n`` can also be replaced by multiple values, such that ``p`` is not a vector but a matrix or a multi-index array of points.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["src/PowerManifold.jl"]
Order = [:macro, :type, :function]
```

## [Product Manifold](@id ProductManifold)

```@autodocs
Modules = [ManifoldsBase]
Pages = ["src/ProductManifold.jl"]
Order = [:macro, :type, :function]
```

## Quotient manifolds

A manifold ``\mathcal M`` is then a a quotient manifold of another manifold ``\mathcal N``,
i.e. for an [equivalence relation](https://en.wikipedia.org/wiki/Equivalence_relation) ``∼``
on ``\mathcal N`` we have

```math
    \mathcal M = \mathcal N / ∼ = \bigl\{ [p] : p ∈ \mathcal N \bigr\},
```

where ``[p] ≔ \{ q ∈ \mathcal N : q ∼ p\}`` denotes the equivalence class containing ``p``.
For more details see Subsection 3.4.1 [AbsilMahonySepulchre:2008](@cite).

This manifold type models an explicit quotient structure.
This should be done if either the default implementation of ``\mathcal M``
uses another representation different from the quotient structure or if
it provides a (default) quotient structure that is different from the one introduced here.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["quotientmanifold.jl"]
Order = [:macro, :type, :function]
```

## Fiber

```@autodocs
Modules = [ManifoldsBase]
Pages = ["Fiber.jl"]
Order = [:macro, :type, :function]
```

## Tangent Space

```@autodocs
Modules = [ManifoldsBase]
Pages = ["TangentSpace.jl"]
Order = [:macro, :type, :function]
```
