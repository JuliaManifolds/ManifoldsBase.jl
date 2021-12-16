# ManifoldsBase.jl

`ManifoldsBase.jl` is a lightweight interface for manifolds.

You can easily implement your algorithms and even your own manifolds just using the interface.
All manifolds from the package here are also based on this interface, so any project based on the interface can benefit from all manifolds, as soon as a certain manifold provides implementations of the functions a project requires.

The main type is the [`AbstractManifold`](@ref). It represents the manifold per se.
During the documentation we will use the [Euclidean Space](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/euclidean.html) and the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/sphere.html) (both implemented in [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl)) as easy examples to often illustrate properties and features of this interface

```@autodocs
Modules = [ManifoldsBase]
Pages = ["maintypes.jl"]
Order = [:type, :function]
```

which should store information about the manifold, for example parameters inherent to the manifold.
