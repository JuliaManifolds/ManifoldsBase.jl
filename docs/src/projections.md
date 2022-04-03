### Projections

A manifold might be embedded in some space.
Often this is implicitly assumed, for example the complex [Circle](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/circle.html) is embedded in the complex plane.
Let‘s keep the circle in mind in the following as a simple example.
For the general case of explicitly stating an embedding and/or to distinguish several, different embeddings, see [Embedded Manifolds](@ref sec-embedded-manifold) below.

To make this a little more concrete, let‘s assume we have a manifold ``\mathcal M`` which is embedded in some manifold ``\mathcal N`` and the image ``i(\mathcal M)`` of the embedding function ``i`` is a closed set (with respect to the topology on ``\mathcal N``).
Then we can do two kinds of projections.

To make this concrete in an example for the Circle ``\mathcal M=\mathcal C := \{ p ∈ ℂ | |p| = 1\}``
the embedding can be chosen to be the manifold ``\mathcal N = ℂ`` and due to our representation of ``\mathcal C`` as complex numbers already, we have ``i(p) = p``, that is the identity as the embedding function.

The first projection we can consider is for a given a point ``p∈\mathcal N`` in the embedding we can look for the closest point on the manifold ``\mathcal M``, i.e.

```math
  \operatorname*{arg\,min}_{q\in \mathcal M}\ d_{\mathcal N}(i(q),p)
```

And this resulting ``q`` we call the projection of ``p`` onto the manifold ``\mathcal M``.

The second projection we can look at is for a given a point ``p∈\mathcal M`` and a vector in ``X\in T_{i(p)}\mathcal N`` in the embedding,
where we can similarly look for the closest tangent vector ``Y∈ T_p\mathcal M``, which we have to embed itself before itself.
Embedding a tangent vector is usually the same as using the pushforward ``\mathrm{d}i_p`` of the embedding (at ``p``).
We obtain

```math
  \operatorname*{arg\,min}_{Y\in T_p\mathcal M}\ \bigl\lVert \mathrm{d}i(p)[Y] - X \bigr\rVert_{i(p)}
```

And we call the resulting ``Y`` the projection of ``X`` onto the tangent space ``T_p\mathcal M`` at ``p``.

Let‘s look at the little more concrete example of the complex circle again.
Here, the closest point of ``p ∈ ℂ`` is just the projection onto the circle, or in other words ``q = \frac{p}{\lvert p \rvert}``, as long as ``p\neq 0``. For ``p=0`` the projection is not defined.
A tangent space ``T_p\mathcal C`` in the embedding is the line through the origin that is orthogonal to a point ``p∈\mathcal C``.
This can be better visualized by looking at ``p+T_p\mathcal C`` which is actually the line tangent to ``p`` on the unit circle.
Note that this shift does not change the resulting projection relative to the origin of the tangent space.

Here the projection can be computed as the classical projection onto the line, i.e.  ``Y = X - ⟨X,p⟩X``.

Both projections onto ``\mathcal C`` and onto ``T_p\mathcal C`` are illustrated in the following figure.

![An example illustrating the two kinds of projections on the Circle.](assets/images/projection_illustration_600.png)

The functions provided in this interface are the following.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["projections.jl"]
Order = [:function]
```
