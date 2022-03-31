## [Retractions and inverse Retractions](@id sec-retractions)

The [exponential and logarithmic map](@ref exp-and-log) might be too expensive to evaluate or not be available in a very stable numerical way on certain manifolds ``\mathcal M``.
Retractions provide a possibly cheap, fast and stable alternative.

A retraction ``\operatorname{retr}_p: T_p\mathcal M → \mathcal M`` is a smooth map that fulfils (for all ``p\in\mathcal M``) that

1. ``\operatorname{retr}_p(0) = p``
2. ``D\operatorname{retr}_p(0): T_p\mathcal M \to T_p\mathcal M`` is the identity map,
i.e. ``D\operatorname{retr}_p(0)[X]=X`` holds for all ``X\in T_p\mathcal M``,

where ``D\operatorname{retr}_p`` denotes the differential of the retraction.

A retraction ``\operatorname{retr}_p`` can be interpreted as a first order approximation to the exponential map ``\exp_p``.

The retraction is called of second order if for all ``X`` the curves ``c(t) = R_p(tX)``
have a zero acceleration at ``t=0``, i.e. ``c''(0) = 0``.

The following figure compares the exponential map [`exp`](@ref)`(M, p, X)` on the [Circle](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/circle.html) `(ℂ)` (or [`Sphere`](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html)`(1)` embedded in ``ℝ^2`` with one possible retraction, the one based on projections.
Note especially that ``\operatorname{dist}(p,q)=\lVert X\rVert_p`` while this is not the case for the result ``\operatorname{retr}_p(X) = q'``.

![A comparison of the exponential map and a retraction on the Circle.](assets/images/retraction_illustration_600.png)

Similar to the [`exp`](@ref)onential map the [`retract`](@ref)ion might not be globally invertible, but locally it is.
So locally one can define the inverse retraction ``\operatorname{retr}_p^{-1}\colon \mathcal M \to T_p\mathcal M``, which
can be seen as a first order approximation to the [`log`](@ref)arithmic map. Within the `ManifoldsBase.jl` interface the inverse retraction is called [`inverse_retract`](@ref).

The general interface looks as follows.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["retractions.jl"]
Order = [:function]
Private = false
Public = true
```

## Types of Retractions

To distinguish different types of retractions, the last argument of the retraction as well as its inverse
specifies a type. The following ones are available.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["retractions.jl"]
Order = [:type]
```

## The functions on layer 3

While you should always add your documentation to [`retract`](@ref) or [`retract!`](@ref) when implementing new manifolds, the actual implementation happens on the following functions on [layer III](@ref design-layer3).

```@autodocs
Modules = [ManifoldsBase]
Pages = ["retractions.jl"]
Order = [:function]
Public = false
Private = true
```
