# [How to implement your own manifold](@id manifold-tutorial)

```@meta
CurrentModule = ManifoldsBase
```

## Introduction

This tutorial explains, how to implement a manifold using the `ManifoldsBase.jl` interface.
We assume that you are familiar with the basic terminology on Riemannian manifolds, especially
the dimension of a manifold, the exponential map, and the inner product on tangent spaces.
To read more about this you can for example check [[do Carmo, 1992](#doCarmo1992)], Chapter 3, first.

Furthermore, we will look into a manifold that is isometrically embedded into e Euclidean space.

In general you need just a datatype (`struct`) that inherits from [`AbstractManifold`](@ref) to define a manifold. No function is _per se_ required to be implemented.
However, it is a good idea to provide functions that might be useful to others, for example [`check_point`](@ref check_point) and [`check_vector`](@ref check_point), as we do in this tutorial.

We start with two technical preliminaries. If you want to start directly, you can [skip](@ref manifold-tutorial-task) this paragraph and revisit it for two of the implementation details.

After that, we will

* [model](@ref manifold-tutorial-task) the manifold
* [implement](@ref manifold-tutorial-checks) two tests, so that points and tangent vectors can be checked for validity, for example also within [`ValidationManifold`](@ref),
* [implement](@ref manifold-tutorial-fn) two functions, the exponential map and the manifold dimension.
* [decorate](@ref manifold-tutorial-emb) the manifold with an embedding to gain further features.

## [Technical preliminaries](@id manifold-tutorial-prel)

There are only two small technical things we need to explain at this point.
First of all our [`AbstractManifold`](@ref)`{𝔽}` has a parameter `𝔽`.
This parameter indicates the [`number_system`](@ref) the manifold is based on, for example `ℝ` for rel manifolds, which is short for `RealNumbers()`.
This indicates that the manifold is a real manifold.

## [Startup](@id manifold-tutorial-startup)

As a start, let's load `ManifoldsBase.jl` and import the functions we consider throughout this tutorial.

```@example manifold-tutorial
using ManifoldsBase, LinearAlgebra, Test
import ManifoldsBase: check_point, check_vector, manifold_dimension, exp!, inner, representation_size, get_embedding
import Base: show
```

We load `LinearAlgebra` for some computations. `Test` is only loaded for illustrations in the examples.

We import the mutating variant of the [`exp`](@ref)onential map, see [the design section on mutating and nonmutating functions](@ref mutating-and-nonmutating).

## [The manifold](@id manifold-tutorial-task)

The manifold we want to implement here a sphere, with a radius $r$.
Since this radius is a property inherent to the manifold, it will become a field of the manifold.
The second information, we want to store is the dimension of the sphere, for example whether it's the 1-sphere, i.e. the circle, represented by vectors $p\in\mathbb R^2$ or norm $r$ or the 2-sphere in $\mathbb R^3$ of radius $r$.
Since the latter might be something we want to [dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch) on, we model it as a parameter of the type.
In general the `struct` of a manifold should provide information about the manifold, which are inherent to the manifold or has to be available without a specific point or tangent vector present.
This is -- most prominently -- a way to determine the manifold dimension.

Note that this a slightly more general manifold than the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/sphere.html) in [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/index.html)

For our example we define the following struct.
While a first implementation might also just take [`AbstractManifold`](@ref)`{ℝ}` as supertype, we directly take
[`AbstractDecoratorManifold`](@ref)`{ℝ}, which will be useful later. For now it does not make a difference.

```@example manifold-tutorial
"""
    ScaledSphere{N} <: AbstractDecoratorManifold{ℝ}

Define an `N`-sphere of radius `r`. Construct by `ScaledSphere(radius,n)`.
"""
struct ScaledSphere{N} <: AbstractDecoratorManifold{ManifoldsBase.ℝ} where {N}
    radius::Float64
end
ScaledSphere(radius, n) = ScaledSphere{n}(radius)
Base.show(io::IO, M::ScaledSphere{n}) where {n} = print(io, "ScaledSphere($(M.radius),$n)")
nothing #hide
```

Here, the last line just provides a nicer print of a variable of that type.
Now we can already initialize our manifold that we will use later, the $2$-sphere of radius $1.5$.

```@example manifold-tutorial
S = ScaledSphere(1.5, 2)
```

## [Checking points and tangents](@id manifold-tutorial-checks)

If we have now a point, represented as an array, we would first like to check, that it is a valid point on the manifold.
For this one can use the easy interface [`is_point`](@ref is_point(M::AbstractManifold, p; kwargs...)).
This is a function on [layer 1](@ref design-layer1) which handles special cases as well cases, so it should not be implemented.
The actual functions where we dispatch per manifold are on [layer 3](@ref design-layer3). Generically, for both  [`is_point`](@ref is_point(M::AbstractManifold, p; kwargs...)) and  [`is_vector`](@ref is_vector(M::AbstractManifold, p, X; kwargs...)) the size is checked using [`check_size`](@ref ManifoldsBase.check_size)
For the test of points the function to implement is [`check_point`](@ref ManifoldsBase.check_point) which we actually will implement, analogously there exists also [`check_vector`](@ref ManifoldsBase.check_vector).
These functions return `nothing` if the point (vector, size) is a correct/valid and returns an error (but not throw it) otherwise.
This is usually a `DomainError`.

We have to check two things: that a point `p` is a vector with `N+1` entries and its norm is the desired radius.
To spare a few lines, we can use [short-circuit evaluation](https://docs.julialang.org/en/v1/manual/control-flow/#Short-Circuit-Evaluation-1) instead of `if` statements.

A first thing we have to specify is how points and tangent vectors are represented, that is we have to specify their [`representation_size`](@ref)

```@example manifold-tutorial
representation_size(::ScaledSphere{N}) where {N} = N+1
nothing #hide
```

This already finishes the size check which [`check_size`](@ref ManifoldsBase.check_size) performs (based on the representation size).

If something has to only hold up to precision, we can pass that down, too using the `kwargs...`, so all three functions we currenlty discuss have these in their signature usually.

```@example manifold-tutorial
function check_point(M::ScaledSphere{N}, p; kwargs...) where {N}
    if !isapprox(norm(p), M.radius; kwargs...)
        return DomainError(norm(p), "The norm of $p is not $(M.radius).")
    end
    return nothing
end
nothing #hide
```

Similarly, we can verify, whether a tangent vector `X` is valid.
It has to fulfill the same size requirements and it has to be orthogonal to `p`.
We can again use the `kwargs`, but also provide a way to check `p`, too.

```@example manifold-tutorial
function check_vector(M::ScaledSphere, p, X; kwargs...)
    if !isapprox(dot(p,X), 0.0; kwargs...)
        return DomainError(dot(p,X), "The tangent $X is not orthogonal to $p.")
    end
    return nothing
end
nothing #hide
```

to test points we can now use

```@example manifold-tutorial
is_point(S, [1.0,0.0,0.0]) # norm 1, so not on S, returns false
@test_throws DomainError is_point(S, [1.5,0.0], true) # only on R^2, throws an error.
p = [1.5,0.0,0.0]
X = [0.0,1.0,0.0]
# The following two tests return true
[ is_point(S, p); is_vector(S,p,X) ]
```

## [Functions on the manifold](@id manifold-tutorial-fn)

For the [`manifold_dimension`](@ref manifold_dimension(M::AbstractManifold)) we have to just return the `N` parameter

```@example manifold-tutorial
manifold_dimension(::ScaledSphere{N}) where {N} = N
manifold_dimension(S)
```

Note that we can even omit the variable name in the first line since we do not have to access any field or use the variable otherwise.

To implement the exponential map, we have to implement the formula for great arcs, given a start point `p` and a direction `X` on the $n$-sphere of radius $r$ the formula reads

````math
\exp_p X = \cos(\frac{1}{r}\lVert X \rVert)p + \sin(\frac{1}{r}\lVert X \rVert)\frac{r}{\lVert X \rVert}X.
````

Note that with this choice we for example implicitly assume a certain metric. This is completely fine. We only have to think about specifying a metric explicitly, when we have (at least) two different metrics on the same manifold.

An implementation of the mutation version, see the [technical note](@ref manifold-tutorial-prel), reads

```@example manifold-tutorial
function exp!(M::ScaledSphere{N}, q, p, X) where {N}
    nX = norm(X)
    if nX == 0
        q .= p
    else
        q .= cos(nX/M.radius)*p + M.radius*sin(nX/M.radius) .* (X./nX)
    end
    return q
end
nothing #hide
```

A first easy check can be done taking `p` from above and any vector `X` of length `1.5π` from its tangent space.
The resulting point is opposite of `p`, i.e. `-p`.

```@example manifold-tutorial
q = exp(S, p, [0.0,1.5π,0.0])
[isapprox(p, -q); is_point(S, q)]
```

## [Adding an isometric embedding](@id manifold-tutorial-emb)

Since the sphere is isometrically embedded, we do not have to implement the [`inner`](@ref)`(M,p,X,Y)`for tangent vectors, but we can “delegate” it to the embedding. The embedding is the [Euclidean](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/essentialmanifold.html) which is also available in `ManifoldsBase` as `DefaultManifold` for testing purposes.

```@example manifold-tutorial
using ManifoldsBase: DefaultManifold, IsIsometricEmbeddedManifold
import ManifoldsBase: active_traits, merge_traits, get_embedding
```

Now we can activate a decorator by specifying that the sphere has the [`IsIsometricEmbeddedManifold`](@ref) trait for the functions `f` on our scaled sphere manifold by writing

```@example manifold-tutorial
active_traits(f, ::ScaledSphere, args...) = merge_traits(IsIsometricEmbeddedManifold())
nothing #hide
```

and then specifying that said embedding is the default manifold

```@example manifold-tutorial
get_embedding(::ScaledSphere{N}) where {N} = DefaultManifold(N+1)
nothing #hide
```

Now metric related functions are passed to this embedding, so the inner product works by using the embedding

Now we can compute the inner product by calling [`inner`](@ref)

```@example manifold-tutorial
X = [0.0, 0.1, 3.0]
Y = [0.0, 4.0, 0.2]
    # returns 1.0 by calling the inner product in DefaultManifold(3)
inner(S, p, X, Y)
```

## [Conclusion](@id manifold-tutorial-outlook)

You can now just continue implementing further functions from `ManifoldsBase.jl`.
but with just [`exp!`](@ref exp!(M::AbstractManifold, q, p, X)) you for example already have

* [`geodesic`](@ref geodesic(M::AbstractManifold, p, X)) the (not necessarily shortest) geodesic emanating from `p` in direction `X`.
* the [`ExponentialRetraction`](@ref), that the [`retract`](@ref retract(M::AbstractManifold, p, X)) function uses by default.

For the [`shortest_geodesic`](@ref shortest_geodesic(M::AbstractManifold, p, q)) the implementation of a logarithm [`log`](@ref ManifoldsBase.log(M::AbstractManifold, p, q)), again better a [`log!`](@ref log!(M::AbstractManifold, X, p, q)) is necessary.

Sometimes a default implementation is provided; for example if you implemented [`inner`](@ref inner(M::AbstractManifold, p, X, Y)), the [`norm`](@ref norm(M, p, X)) is defined. You should overwrite it, if you can provide a more efficient version. For a start the default should suffice.
With [`log!`](@ref log!(M::AbstractManifold, X, p, q)) and [`inner`](@ref inner(M::AbstractManifold, p, X, Y)) you get the [`distance`](@ref distance(M::AbstractManifold, p, q)), and so.

In summary with just these few functions you can already explore the first things on your own manifold. Whenever a function from `Manifolds.jl` requires another function to be specifically implemented, you get a reasonable error message.

## Literature

```@raw html
<ul>
<li id="doCarmo1992">
    [<a>doCarmo, 1992</a>]
    M. P. do Carmo,
    <emph>Riemannian Geometry</emph>,
    Birkhäuser Boston, 1992,
    ISBN: 0-8176-3490-8.
</li>
</ul>
```
