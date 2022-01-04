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

Points and tangent vectors do not necessarily have to be typed. Usually
one can just use any type. When a manifold has multiple representations, these should be distinguished by point and vector types. Then it might be that types just encapsulate a vector value.
This is taken into account by the following macros, that forward several actions just to this field. Most prominently vector operations for the tangent vectors.
If there is still a default case, a macro sets this type to be equivalent to calling the manifold functions just with the types field that carries the value.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["point_vector_fallbacks.jl"]
Order = [:type, :function, :macro]
```
