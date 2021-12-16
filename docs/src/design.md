# Main Design Principles

The core types are the [`AbstractManifold`](@ref) and – if necessary – types for points.
Another important part

Before starting with these, a short prequel on the Number systems to represent is required.

## Interface design

The interface for a manifold is defined to be as generic as possible, such that applications can be implemented as independently as possible from an actual manifold.
This way, algorithms like those from [`Manopt.jl`](https://manoptjl.org) can be implemented on _arbitrary_ manifolds.

The main design criteria for the interface are:

* Aims to also provide _efficient_ _global state-free_, both _in-place_ and _out-of-place_ computations whenever possible.
* Provide a high level interface that is easy to use.

Therefore this interface has 3 main features, that we will explain using two (related)
concepts, the [exponential map](https://en.wikipedia.org/wiki/Exponential_map_(Riemannian_geometry)) that maps a tangent vector ``X`` at a point ``p`` to a point ``q`` or mathematically ``\exp_p:T_p\mathcal M \to \mathcal M`` and its generalisation, a [`retract`](@ref)ion ``\operatorname{retr}_p`` with same domain and range.

You do not need to know their exact definition at this point, just that there is _one_ exponential map on a Riemannian manifold, and several retractions, where one of them is the exponential map (sometime called exponential retraction for completeness). Every retraction has its own subtype of the [`AbstractRetractionMethod`](@ref) that uniquely defines it.

The following three design patterns aim to fulfill the criteria from above, while
also avoiding ambiguities in multiple dispatch using the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) approach.

### General order of parameters

Since the central element for functions on a manifold is the manifold itself, it should always be the first parameter, even for mutating functions.

### Mutating and allocating functions

Every function, where this is applicable should provide a mutating and an allocating variant.
For example for the exponential map `exp(M,p,x)` returns a _new_ point `q` where the result is computed in.
On the other hand `exp!(M, q, p, X)` computes the result in place of `q`, where the design of the implementation
should keep in mind that also `exp!(M,p,p,X)` should correctly overwrite `p`.

The interface provides a way to determine the allocation type and a result to compute/allocate
the resulting memory, such that the default implementation allocating functions, like [`exp`](@ref) is to allocate the resulting memory and call [`exp!`](@ref).

!!! note
    it might be useful to provide two distinct implementations, for example when using AD schemes.
    The default is meant for ease of use (concerning implementation), since then one has to just implement the mutating variants.

### The higher level interface and ease of use

The higher level interface should aim for a convenience layer that resolves defaults and
creates fallbacks for certain input parameters, that have these properties.
It usually should not dispatch on a certain manifold nor on certain point or (co- or tangent) vector types.

This layer should also not resolve/dispatch from the allocating to the mutating variant.

This is maybe best illustrated by two examples

1. the exponential map usually has a long form, where one can specify a fraction (of `X`) where the evaluation should be. This is generically implemented a

  ```julia
    exp(M::AbstractManifold, p, X, t::Real) = exp(M, p, t * X)
  ```

  On this level neither the manifold _nor_ the points should be too strictly typed, points and vectors should – for best of cases – never be types.

2. for the retraction, a default retraction (usually exp) is specified/defined via [`default_retraction_method`](@ref).
  This also means, that the last parameter of

  ```julia
  retract(M::AbstractManifold, p, X, ::AbstractRetractionMethod=default_retraction_method(M))
  ```

  is optional and for a concrete type the dispatch on a certain retraction is done next.
  To avoid ambiguities, this concrete type should always be the first argument we dispatch on:
  The `ExponentialRetractionMethod` calls `exp`, any other retraction calls a function of different name without this last parameter,
  for example the `PolarRetractionMethod` by default calls `retract_polar(M,p,X)`, which is actualy a function from the lower level, see next section.

### The lower level interface – performance

This lower level aims for performance, that is, any function should have as few as possible optional and keyword arguments
and be typed as concrete as possible/necessary. This means

* the function name should be similar to its high level parent (for example `retract` and `retract_polar`from above)
* The manifold type in method signature should always be as narrow as possible.
* the points/vectors should either be untyped (for the default representation of if there is only one) or provide all types concretely.

The first thing to do on this level is the aforementioned default to pass from allocating to mutating functions.

(TODO/Discuss - dispatch for decorators here instead of with the current/modular decorator pattern).
