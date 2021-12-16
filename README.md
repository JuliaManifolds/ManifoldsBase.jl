# ManifoldsBase.jl
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html)
[![Build Status](https://travis-ci.org/JuliaManifolds/ManifoldsBase.jl.svg?branch=master)](https://travis-ci.org/JuliaManifolds/ManifoldsBase.jl/)
[![codecov.io](http://codecov.io/github/JuliaManifolds/ManifoldsBase.jl/coverage.svg?branch=master)](https://codecov.io/gh/JuliaManifolds/ManifoldsBase.jl/)
[![arXiv](https://img.shields.io/badge/arXiv%20CS.MS-2106.08777-blue.svg)](https://arxiv.org/abs/2106.08777)

Basic interface for manifolds in Julia.

The project [`Manifolds.jl`](https://github.com/JuliaManifolds/Manifolds.jl)
is based on this interface and provides a variety of manifolds.

## Main types

This packages provides an  [`AbstractManifold`](https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html#ManifoldsBase.AbstractManifold) type, that is parametrized by its field, i.e., whether elements on the manifold are represented using real (`ℝ`) or complex numbers (`ℂ`) or even quaternions (`ℍ`).
A manifold is a set of points. These are usually represented by (real-, complex-, quarternion-) valued arrays, most prominently matrices.
Sometimes they are also vectors, higher-dimensional arrays or even more complex data structures.
For more advanced points, you can also use subtypes [`AbstractManifoldPoint`](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html#ManifoldsBase.AbstractManifoldPoint).
Sometimes points on a manifold can be represented in different forms.
The first (default) representation can use just arrays or any other arbitrary type and methods can be implemented without type annotations.
Any further representation _must_ use a specific type for methods to dispatch on. Then the default case can (and should) also have a default implementation that just unwraps to the internal value (see [`@default_manifold_fallbacks`]()).

## Interface design

The interface for a manifold is defined to be as generic as possible, such that applications can be implemented as independently as possible from an actual manifold.
This way, algorithms like those from [`Manopt.jl`](https://manoptjl.org) can be implemented on _arbitrary_ manifolds.

The main design criteria for the interface are:

* Aims to also provide _efficient_ _global state-free_, both _in-place_ and _out-of-place_ computations whenever possible.
* Provide a high level interface that is easy to use.

Therefore this interface has 3 main features, that we will explain using two (related)
concepts, the [exponential map](https://en.wikipedia.org/wiki/Exponential_map_(Riemannian_geometry)) that maps a tangent vector ``X`` at a point ``p`` to a point ``q`` or mathematically ``\exp_p:T_p\mathcal M \to \mathcal M`` and its generalisation, a [retraction]() ``\operatorname{retr}_p`` with same domain and range.

You do not need to know their exact definition at this point, just that there is _one_ exponential map on a Riemannian manifold, and several retractions, where one of them is the exponential map (sometime called exponential retraction for completeness). Every retraction has its own subtype of the [`AbstractRetractionMethod`]() that uniquely defines it.

The following three design patterns aim to fulfill the criteria from above, while
also avoiding ambioguities in multiple dispatch using the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) approach.

### General order of parameters

Since the central element for functions on a manifold is the manifold itself, it should always be the first parameter, even for mutating functions.

### Mutating and allocating functions

Every function, where this is applicable should provide a mutating and an allocating variant.
For example for the exponential map `exp(M,p,x)` returns a _new_ point `q` where the result is computed in.
On the other hand `exp!(M, q, p, X)` computes the result in place of `q`, where the design of the implementation
should keep in mind that also `exp!(M,p,p,X)` should correctly overwrite `p`.

The interface provides a way to determine the allocation type and a result (see [`allocate_result`]()) to compute/allocate
the resulting memory, such that the default implementation allocating functions, like `exp` is to allocate the resulting memory and call `exp!`.

!!! note
   it might be useful to provide two distinct implementations, for example when using AD schemes.
   The default is meant for ease of use (concerning implementation), since then one has to just implement the mutating variants.

### The higher level interface and ease of use

The higher level interface should aim for a convenience layer that resolves defaults and
creates fallbacks for certain input parameters, that have these properties.
It usually should not dispatch on a certain manifold nor on certain point or (co- or tangent) vector types.

This layer should also not resolve/dispatch from the allocating to the mutating variant.

This is maybe best illustrated by two examples

1. the exponential map usually has a long form, where one can specify a fraction (of `X`) where the evaluation should be. This is generically impülemented a

  ```julia
    exp(M::AbstractManifold, p, X, t::Real) = exp(M, p, t * X)
  ```

  On this level neither the manifold _nor_ the points should be too strictly typed, points and vectors should – for best of cases – never be types.

2. for the retraction, a default retraction (usually exp) is specified/defined via [`default_retraction_method`]().
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

## Bases

Several different types of bases for a tangent space at `p` on a [`AbstractManifold`](https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html#ManifoldsBase.AbstractManifold) are provided.
Methods are provided to obtain such a basis, to represent a tangent vector in a basis and to reconstruct a tangent vector from coefficients with respect to a basis.
The last two can be performed without computing the complete basis.
Further a basis can be cached and hence be reused, see [`CachedBasis`](https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html#ManifoldsBase.CachedBasis).

## (Abstract) Example Manifolds

While this package does not provide _actual_, _concrete_ manifolds (see [`Manifolds.jl`](https://github.com/JuliaManifolds/Manifolds.jl) for these),
it provides a few generic (meta) manifolds that are helpful.

### `DecoratorManifold`

The decorator manifold enhances a manifold by certain, in most cases implicitly
assumed to have a standard case, properties, see for example the `EmbeddedManifold`.
The decorator acts semi transparently, i.e. `:transparent` for all functions not affected by that
decorator and `:intransparent` otherwise. Another possibility is, that the decorator just
passes to `:parent` in order to fill default values.

### `EmbeddedManifold`

The embedded manifold models the embedding of a manifold into another manifold.
This way a manifold can benefit from existing implementations.
One example is the `TransparentIsometricEmbeddingType` where a manifold uses the metric,
`inner`, from its embedding.

### `ValidationManifold`

The `ValidationManifold` further illustrates how one can also used types to
represent points on a manifold, tangent vectors, and cotangent vectors,
where values are encapsulated in a certain type.

In general, `ValidationManifold` might be used for manifolds where these three types are represented
by more complicated data structures or when it is necessary to distinguish these
by type.

This adds a semantic layer to the interface, and the default implementation of
`ValidationManifold` adds checks to all inputs and outputs of typed data.

### `DefaultManifold`

This interface includes a simple `DefaultManifold`, which is a reduced version
of the [`Euclidean`](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/euclidean.html)
manifold from [`Manifolds.jl`](https://github.com/JuliaManifolds/Manifolds.jl).
This is indeed a concrete manifold, mainly to illustrate how a manifold is implemented
as well as for testing the interface functions.

## Citation

If you use `ManifoldsBase.jl` in your work, please cite the following

```biblatex
@online{2106.08777,
Author = {Seth D. Axen and Mateusz Baran and Ronny Bergmann and Krzysztof Rzecki},
Title = {Manifolds.jl: An Extensible Julia Framework for Data Analysis on Manifolds},
Year = {2021},
Eprint = {2106.08777},
Eprinttype = {arXiv},
}
```
