# The Manifold Type

## Number systems

```@autodocs
Modules = [ManifoldsBase]
Pages = ["numbers.jl"]
Order = [:type, :function]
```

## The main type: The `AbstractManifold`

The main type is the [`AbstractManifold`](@ref). It represents the manifold per se.
During the documentation we will use the [Euclidean Space](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) and the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html) (both implemented in [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl)) as easy examples to often illustrate properties and features of this interface

```@autodocs
Modules = [ManifoldsBase]
Pages = ["maintypes.jl"]
Order = [:type, :function]
```

which should store information about the manifold, for example parameters inherent to the manifold.

## Points and Tangent Vectors

(TODO)

```@autodocs
Modules = [ManifoldsBase]
Pages = ["point_vector_fallbacks.jl"]
Order = [:type, :function, :macro]
```
