# What are Manifolds?

```@contents
Pages = ["what-are-manifolds-tutorial2.md"]
Depth = 2
```

## Introduction

In [JuliaManifolds](https://github.com/JuliaManifolds) we develop tools to define and work with manifolds. This package here, [`ManifoldsBase.jl`](@ref ManifoldsBase.jl), defines an API to define manifolds. But also when someone recommended to work with one of the manifolds in [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) a first question that might come up is

**What is a manifold?**

In this tutorial, we aim to start with the most general concepts to those that provide most structure.
We always provide references both into the [JuliaManifolds](https://github.com/JuliaManifolds) ecosystem but also into literature.
We focus in this tutorial on presenting concepts and ideas that are necessary to work with manifolds and gain a good intuition.
While a basic knowledge of multivariate calculus, linear algebra and ordinary differential equations is required to follow this tutorial, we try to also link for these concepts to corresponding definitions and further literature.

Let's start with an informal idea based on a concrete example. One of the first manifolds one might come in contact with is that of the 2-[`Sphere`](@extref Manifolds.Sphere), for example as model of our planet earth. Mathematically we could consider all unit-norm vectors in ``\mathbb R^3``, informally all points that are the same distance from the center of the earth, where we choose to define this distance to be ``1``.
Compared to introductory calculus, there is one main difference for this set of points that build the surface of the earth: we can not add two of these points and stay on the sphere: North pole “plus” north pole is a point twice as far away from the earth's center as the north pole itself, so not a point on the sphere.

Still, there is an important observation, that we make from local observations: locally the sphere “feels like” the usual ``\mathbb R^2`` – a plane. For example by considering a city map, we can convince ourselves, we always have two directions we can walk into: the north/south and the east/west axes.
While globally the sphere is not the same as ``\mathbb R^2``, due to the lack of addition, we can collect these local charts into an atlas to get an “overal picture” of the shere.

For more than two dimensions, we can formalize this as  [AbsilMahonySepulchre:2008; Section 3.1](@cite):

> A ``d``-dimensional manifold can be informally defined as a set M covered
> with a “suitable” collection of coordinate patches, or charts, that identify
> certain subsets of M with open subsets of ``\mathbb R^2``.

## A Topological manifold

The first concept that provides tools we can use is that of a topological manifold ``\mathcal M``.
Such a manifold ``M`` is described by its dimension `d=`[`manifold_dimension`](@ref)`(M)` and a set of functions called an atlas ``\{φ_i\}_{i\in I}`` where ``φ_i\colon U_i \to \mathbb{R}^d`` are charts indexed by ``i`` from some index set ``I`` [Lee:2012; p. 4](@cite) and ``U_i \subseteq \mathcal{M}``. [^1] [^2]

For each point ``p \in \mathcal{M}`` there is ``i\in I`` such that ``p \in U_i``.
There are also some regularity conditions which we skip here because they are rarely relevant. [^3]

[JuliaManifolds](https://github.com/JuliaManifolds) has a few functions for working at this level.
First, [`manifold_dimension`](@ref) returns the number `d` for a given manifold.
Next, [`get_chart_index`](@extref Manifolds :jl:method:`Manifolds.get_chart_index-Tuple{AbstractManifold, AbstractAtlas, Any, Any}`) points to one of the charts such that ``p`` is in its domain.
The value of chart on a point can be calculated using [`get_parameters`](@extref `Manifolds.get_parameters`) and its inverse using [`get_point`](@extref `Manifolds.get_point`).
When we have two charts ``φ_i, φ_j``, the composition ``φ_j \circ φ_i^{-1}`` is called the transition map from ``φ_i`` to ``φ_j``, see [`transition_map`](@extref `Manifolds.transition_map`).
More details are discussed in [this page](https://juliamanifolds.github.io/Manifolds.jl/stable/features/atlases/) and [this tutorial](https://juliamanifolds.github.io/Manifolds.jl/stable/tutorials/working-in-charts/) demonstrates a use case.
Often additional restrictions are imposed, for example only atlases with only [differentiable](https://en.wikipedia.org/wiki/Differentiable_manifold), smooth or [analytic](https://en.wikipedia.org/wiki/Analytic_manifold) transition maps are considered.

An even stronger restriction holds for complex manifolds: the charts are complex valued (``φ_i\colon U_i \to \mathbb{C}^n``) and transition maps are required to be holomorphic functions.
Note that not any manifold represented using complex numbers is a complex manifold.
In particular, no manifold of odd (real) dimension can be a complex manifold.
Complex representation is a feature of [CR manifolds](https://en.wikipedia.org/wiki/CR_manifold).
Instances of [`AbstractManifold`](@ref) with complex number system `𝔽` are not required to be complex manifolds, though they are expected to be CR manifolds.

## Literature

```@bibliography
Pages = ["tutorials/what-are-manifolds-tutorial2.md"]
Canonical=false
```

## Footnotes

[^1]: In JuliaManifolds we either have a finite number of charts (for example, spheres require only two) or one chart for each point on the manifold.

[^2]: Sometimes other number systems are considered for the codomain of charts, most notably complex numbers. This discussion is restricted to the real case because it’s general enough for practical purposes. Complex atlases can be represented as real atlases with real and imaginary parts separated. Quaternionic manifolds are most easily expressed though fiber bundles. [Other generalizations](https://math.stackexchange.com/a/581087) often lead to spaces that are no longer manifolds.

[^3]: Specifically manifolds are required to be second-countable Hausdorff spaces, see [Lee:2012; p. 3](@cite) for more details.
