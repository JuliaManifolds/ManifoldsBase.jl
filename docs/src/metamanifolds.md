# Meta Manifolds

While the interface does not provide concrete manifolds itself, it does provide several manifolds that can be build based on a given [`AbstractManifold`](@ref) instance.

## [(Abstract) power manifold](@id sec-power-manifold)

A power manifold is constructed like higher dimensional vector spaces are formed from the real line, just that for every point ``p = (p_1,\ldots,p_n) ∈ \mathcal M^n`` on the power manifold ``\mathcal M^n`` the entries of ``p`` are points ``p_1,\ldots,p_n ∈ \mathcal M`` on some manifold ``\mathcal M``. Note that ``n`` can also be replaced by multiple values, such that ``p`` is not a vector but a matrix or a multi-index array of points.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["src/PowerManifold.jl"]
Order = [:macro, :type, :function]
```

## VectorBundle and Vector Fibre

```@autodocs
Modules = [ManifoldsBase]
Pages = ["VectorBundle.jl"]
Order = [:macro, :type, :function]
Public=true
Private=false
```

```@autodocs
Modules = [ManifoldsBase]
Pages = ["VectorFiber.jl"]
Order = [:macro, :type, :function]
Public=true
Private=false
```

## Fiber and FiberBundle

```@autodocs
Modules = [ManifoldsBase]
Pages = ["Fiber.jl"]
Order = [:macro, :type, :function]
```

```@autodocs
Modules = [ManifoldsBase]
Pages = ["FiberBundle.jl"]
Order = [:macro, :type, :function]
```
