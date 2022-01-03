## [Retractions and inverse Retractions](@id sec-retractions)

The [exponential and logarithmic map](@ref exp-and-log) might be too expensive to evaluate or not be available in a very stable numerical way. Retractions provide a possibly cheap, fast and stable alternative.

The following figure compares the exponential map [`exp`](@ref)`(M, p, X)` on the [Circle](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/circle.html) `(ℂ)` (or [`Sphere`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html)`(1)` embedded in $ℝ^2$ with one possible retraction, the one based on projections. Note especially that ``\mathrm{dist}(p,q)=\lVert X\rVert_p`` while this is not the case for ``q'``.

![A comparson of the exponential map and a retraction on the Circle.](assets/images/retraction_illustration_600.png)

```@autodocs
Modules = [ManifoldsBase]
Pages = ["retractions.jl"]
Order = [:function]
```

To distinguish different types of retractions, the last argument of the (inverse) retraction
specifies a type. The following ones are available.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["retractions.jl"]
Order = [:type]
```
