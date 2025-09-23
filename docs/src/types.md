# The Manifold interface

## The `AbstractManifold`

The main type is the [`AbstractManifold`](@ref). It represents the manifold per se.
Throughout the documentation of `ManifoldsBase.jl` we might use the [Euclidean Space](@extref `Manifolds.Euclidean`) and the [Sphere](@extref `Manifolds.Sphere`) (both implemented in [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl)) as easy examples to illustrate properties and features of this interface on concrete examples.

```@docs
AbstractManifold
```

which should store information about the manifold, for example parameters inherent to the manifold.

In order of parameters for any subtype of [`AbstractManifold`](@ref) it would be good if the first parameter is the number type that the abstract type has as well.
The second one should be the type of field type for size information. This might be the dimension of the manifold, like for the [`Sphere`](@extref `Manifolds.Sphere`) or any other number(s) determining the manifolds representation or dimension, like the matrix size for the [`SymmetricPositiveDefinite`](@extref `Manifolds.SymmetricPositiveDefinite`) manifold.

## Points on a manifold

Points do not necessarily have to be typed.
Usually one can just use any type. When a manifold has multiple representations, these should be distinguished by point and vector types.

```@docs
AbstractManifoldPoint
```

Converting points between different representations can be performed using the `convert` function with either two or three arguments (`convert(T, M, p)` or `convert(T, p)`). For some manifolds providing `M` may be necessary. The first variant falls back to the second variant.

## Tangent and Cotangent spaces

```@autodocs
Modules = [ManifoldsBase]
Pages = ["vector_spaces.jl"]
Order = [:type, :function]
```

This interface also covers a large variety how to [model bases in tangent spaces](@ref bases).

Converting tangent vectors between different representations can be performed using the `convert` function with either three or four arguments (`convert(T, M, p, X)` or `convert(T, p, X)`). For some manifolds providing `M` may be necessary. The first variant falls back to the second variant.

## Macros for automatic forwards for simple points/tangent vectors

When distinguishing different representations of points or tangent vectors on one manifold,
it might happen that both a subtype of [`AbstractManifoldPoint`](@ref) and a subtype of [`AbstractTangentVector`](@ref)
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

## [Type Parameter](@id type-parameter)

Concrete [`AbstractManifold`](@ref)s usually correspond to families of manifolds that are parameterized by some numbers, for example determining their [`manifold_dimension`](@ref). Those numbers can either be stored in a field or as a type parameter of the structure. The [`TypeParameter`](@ref ManifoldsBase.TypeParameter) offers the flexibility
to have this parameter either as type parameter or a field.

```@docs
ManifoldsBase.TypeParameter
ManifoldsBase.wrap_type_parameter
```
