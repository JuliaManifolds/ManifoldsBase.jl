# Functions on manifolds

This page collects several basic functions on manifolds.

## Validation

ince points and tangent vectors are represented usually as multidimensional arrays or for more complex cases as structs, there might be values, which invalidate a point of tangent vector. Here the interface provides two [high level functions](@ref design-layer1).

```@docs
is_point
is_vector
```

These are mapped to the [lower level functions](@ref design-layer3)

```@docs
ManifoldsBase.check_point
ManifoldsBase.check_vector
ManifoldsBase.check_size
```

## [The exponential and the logarithmic map, and geodesics](@id exp-and-log)

Geodesics are the generalizations of a straight line to manifolds, i.e. their intrinsic acceleration is zero.
Together with geodesics one also obtains the exponential map and its inverse, the logarithmic map.
Informally speaking, the exponential map takes a vector (think of a direction and a length) at one point and returns another point,
which lies towards this direction at distance of the specified length. The logarithmic map does the inverse, i.e. given two points, it tells which vector “points towards” the other point.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["exp_log_geo.jl"]
Order = [:function]
```


## Vector transport

There are three main functions for vector transport:

* [`vector_transport_along`](@ref)
* [`vector_transport_direction`](@ref)
* [`vector_transport_to`](@ref)

Different types of vector transport are implemented using subtypes of [`AbstractVectorTransportMethod`](@ref):

* [`ParallelTransport`](@ref)
* [`PoleLadderTransport`](@ref)
* [`ProjectionTransport`](@ref)
* [`SchildsLadderTransport`](@ref)

```@autodocs
Modules = [ManifoldsBase]
Pages = ["vector_transport.jl"]
Order = [:type, :function]
```

## Projections

A manifold might be embedded in some space.
Often this is implicitly assumed, for example the complex [Circle](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/circle.html) is embedded in the complex plane.
Let‘s keep the circle in mind in the following as a simple example.
For the general case see of explicitly stating an embedding and/or distinguising several, different embeddings, see [Embedded Manifolds](@ref subsec-embeddedmanifold) below.

To make this a little more concrete, let‘s assume we have a manifold ``\mathcal M`` which is embedded in some manifold ``\mathcal N`` and the image ``i(\mathcal M)`` of the embedding function ``i`` is a closed set (with respect to the topology on ``\mathcal N``). Then we can do two kinds of projections.

To make this concrete in an example for the Circle ``\mathcal M=\mathcal C := \{ p ∈ ℂ | |p| = 1\}``
the embedding can be chosen to be the manifold ``N = ℂ`` and due to our representation of ``\mathcal C`` as complex numbers already, we have ``i(p) = p`` the identity as the embedding function.

1. Given a point ``p∈\mathcal N`` we can look for the closest point on the manifold ``\mathcal M`` formally as

```math
  \operatorname*{arg\,min}_{q\in \mathcal M} d_{\mathcal N}(i(q),p)
```

And this resulting ``q`` we call the projection of ``p`` onto the manifold ``\mathcal M``.

2. Given a point ``p∈\mathcal M`` and a vector in ``X\inT_{i(p)}\mathcal N`` in the embedding we can similarly look for the closest point to ``Y∈ T_p\mathcal M`` using the push forward ``\mathrm{d}i_p`` of the embedding.

```math
  \operatorname*{arg\,min}_{Y\in T_p\mathcal M} \lVert \mathrm{d}i(p)[Y] - X \rVert_{i(p)}
```

And we call the resulting ``Y`` the projection of ``X`` onto the tangent space ``T_p\mathcal M`` at ``p``.

Let‘s look at the little more concrete example of the complex Circle again.
Here, the closest point of ``p ∈ ℂ`` is just the projection onto the circle, or in other words ``q = \frac{p}{\lvert p \rvert}``. A tangent space ``T_p\mathcal C`` in the embedding is the line orthogonal to a point ``p∈\mathcal C`` through the origin.
This can be better visualized by looking at ``p+T_p\mathcal C`` which is actually the line tangent to ``p``. Note that this shift does not change the resulting projection relative to the origin of the tangent space.

Here the projection can be computed as the classical projection onto the line, i.e.  ``Y = X - ⟨X,p⟩X``.

this is illustrated in the following figure

![An example illustrating the two kinds of projections on the Circle.](assets/images/projection_illustration_600.png)

```@autodocs
Modules = [ManifoldsBase]
Pages = ["projections.jl"]
Order = [:function]
```

## Further functions

```@autodocs
Modules = [ManifoldsBase]
Pages = ["ManifoldsBase.jl"]
Order = [:type, :function]
```

## Error Messages

especially to collect and display errors on [`AbstractPowerManifold`](@ref)s the following
component and collection error messages are available.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["errors.jl"]
Order = [:type]
```