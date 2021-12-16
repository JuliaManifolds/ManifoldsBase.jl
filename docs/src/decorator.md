# A Decorator for manifolds

A decorator manifold extends the functionality of a [`AbstractManifold`](@ref) in a semi-transparent way.
It internally stores the [`AbstractManifold`](@ref) it extends and by default for functions defined in the [`ManifoldsBase`](interface.md) it acts transparently in the sense that it passes all functions through to the base except those that it actually affects.
For example, because the [`ValidationManifold`](@ref) affects nearly all functions, it overwrites nearly all functions, except a few like [`manifold_dimension`](@ref).

By default, i.e. for a plain new decorator, all functions are transparent, i.e. passed down to the manifold the [`AbstractDecoratorManifold`](@ref) decorates.
To implement a method for a decorator that behaves differently from the method of the same function for the internal manifold, two steps are required.
Let's assume the function is called `f(M, arg1, arg2)`, and our decorator manifold `DM` of type `OurDecoratorManifold` decorates `M`.
Then

1. set `decorator_transparent_dispatch(f, M::OurDecoratorManifold, args...) = Val(:intransparent)`
2. implement `f(DM::OurDecoratorManifold, arg1, arg2)`

This makes it possible to extend a manifold or all manifolds with a feature or replace a feature of the original manifold.

A final technical note – if several manifolds have similar transparency rules concerning functions from the interface, the last parameter `T` of the [`AbstractDecoratorManifold`](@ref)`{𝔽,T<:`[`AbstractDecoratorType`](@ref)`}` can be used to dispatch on different transparency schemes.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["DecoratorManifold.jl"]
Order = [:type, :macro, :function]
```
