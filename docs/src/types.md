# The Manifold interface

## The `AbstractManifold`

The main type is the [`AbstractManifold`](@ref). It represents the manifold per se.
Throughout the documentation of `ManifoldsBase.jl` we might use the [Euclidean Space](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) and the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html) (both implemented in [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl)) as easy examples to illustrate properties and features of this interface on concrete examples.

```@docs
AbstractManifold
```

which should store information about the manifold, for example parameters inherent to the manifold.

## Points on a manifold

Points do not necessarily have to be typed.
Usually one can just use any type. When a manifold has multiple representations, these should be distinguished by point and vector types.

```@docs
AbstractManifoldPoint
```

## Tangent spaces

```@autodocs
Modules = [ManifoldsBase]
Pages = ["vector_spaces.jl"]
Order = [:type, :function]
```

This interface also covers a large variety how to [model bases in tangent spaces](@ref bases)

## Macros for automatic forwards for simple points/tangent vectors

When distinguishing different representations of points or tangent vectors on one manifold,
it might happen that both a subtype of [`AbstractManifoldPoint`](@ref) and a subtype of [`TVector`](@ref)
are just encapsulating a value

This is taken into account by the following macros, that forward several actions just to this field. Most prominently vector operations for the tangent vectors.
If there is still a default case, a macro sets this type to be equivalent to calling the manifold functions just with the types field that carries the value.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["point_vector_fallbacks.jl"]
Order = [:type, :function, :macro]
```

## [Number Systems](@id number-system)

The [`AbstractManifold`](@ref) has one parameter to distinguish the number system a manifold is based on.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["numbers.jl"]
Order = [:type, :function]
```
