# Main Design Principles

The interface for a manifold is defined to be as generic as possible, such that applications can be implemented as independently as possible from an actual manifold.
This way, algorithms like those from [`Manopt.jl`](https://manoptjl.org) can be implemented on _arbitrary_ manifolds.

The main design criteria for the interface are:

* Aims to also provide _efficient_ _global state-free_, both _in-place_ and _out-of-place_ computations whenever possible.
* Provide a high level interface that is easy to use.

Therefore this interface has 3 main features, that we will explain using two (related)
concepts, the [exponential map](https://en.wikipedia.org/wiki/Exponential_map_(Riemannian_geometry)) that maps a tangent vector ``X`` at a point ``p`` to a point ``q`` or mathematically ``\exp_p:T_p\mathcal M \to \mathcal M`` and its generalization, a [`retract`](@ref)ion ``\operatorname{retr}_p`` with the same domain and range.

You do not need to know their exact definition at this point, just that there is _one_ exponential map on a Riemannian manifold, and several retractions, where one of them is the exponential map (called [`ExponentialRetraction`](@ref) for completeness). Every retraction has its own subtype of the [`AbstractRetractionMethod`](@ref) that uniquely defines it.

The following three design patterns aim to fulfil the criteria from above, while
also avoiding ambiguities in multiple dispatch using the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) approach.

## General order of parameters

Since the central element for functions on a manifold is the manifold itself, it should always be the first parameter, even for mutating functions. Then the classical parametzers of a function (for example a point and a tangent vector for the retraction) follow and the final part are parameters to further dispatch on, which usually have their defaults.

## A 3-Layer architecture for dispatch

The general architecture consists of three layers

* The high level interface for ease of use – and to dispatch on other manifolds.
* The intermediate layer to dispatch on different parameters in the last section, e.g. type of retraction or vector transport.
* The lowest layer for specific manifolds to dispatch on different types of points and tangent vectors. Usually this layer with a specific manifold and no optional parameters.

These three layers are described in more detail in the following.
The main motivation to introduce these layers is, that it reduces method ambiguities.
It also provides a good structure where to implement extensions to this interface.

### [Layer I: The high level interface and ease of use](@id design-layer1)

The highest layer for convenience of decorators.
A usual scheme is, that a manifold might assume several things implicitly, for example the default implementation of the sphere $\mathbb S^n$ using unit vectors in $\mathbb R^{n+1}$.
The embedding can be explicitly used to avoid re-implementations – the inner product can be “passed on” to its embedding.

To do so, we “decorate” the manifold by making it an [`AbstractDecoratorManifold`](@ref) and activating the right traits see [the example](@ref manifold-tutorial).

The explicit case of the [`EmbeddedManifold`](@ref) can be used to distinguish different embeddings of a manifold, but also their dispatch (onto the manifold or its embedding, depending on the type of embedding) happens here.

Note that all other parameters of a function should be as least typed as possible for all parameters besides the manifold.
With respect to the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) paradigm, this layer dispatches the _manifold first_.
We also stay as abstract as possible, for example on the [`AbstractManifold`](@ref) level if possible.


If a function has optional positional arguments, (like [`retract`](@ref)) their default values might be filled/provided on this layer.
This layer ends usually in calling the same functions like [`retract`](@ref) but prefixed with a `_` to enter [Layer II](@ref design-layer2).

!!! note
    Usually only functions from this layer are exported from the interface, since these are the ones one should use for generic implementations.
    If you implement your own manifold, `import` the necessary lower layer functions as needed.

### [Layer II: An internal dispatch interface for parameters](@id design-layer2)

This layer is an interims layer to dispatch on the (optional/default) parameters of a function.
For example the last parameter of retraction: [`retract`](@ref) determines the type (variant) to be used.
The last function in the previous layer calls `_retract`, which is an internal function.
These parameters are usually the last parameters of a function.

On this layer, e.g. for `_retract` only these last parameters should be typed, the manifold should stay at the [`AbstractManifold`](@ref) level.
The layer dispatches on different functions per existing parameter type (and might pass this one further on, if it has fields).
Function definitions on this layer should only be extended when introducing new such parameter types, for example when introducing a new type of a retraction.

The functions from this layer should never be called directly, are hence also not exported and carry the `_` prefix.
They should only be called as the final step in the previous layer.

If the default parameters are not dispatched per type, using `_` might be skipped.
The same holds for functions that do not have these parameters.
The following resolution might even be seen as a last step in layer I or the resolution here in layer II.

```julia
exp(M::AbstractManifold, p, X, t::Real) = exp(M, p, t * X)
```

When there is no dispatch for different types of the optional parameter (here `t`), the `_` might be skipped.
One could hence see the last code line as a definition on Layer I that passes directly to Layer III, since there are not parameter to dispatch on.

To close this section, let‘s look at an example.
The high level (or [Layer I](@ref design-layer1)) definition of the retraction is given by

```julia
retract(M::AbstractManifold, p, X, m::AbstractRetractionMethod=default_retraction_method(M)) = _retract(M, p, X, m)
```

This level now dispatches on different retraction types `m`.
It usually passes to specific functions implemented in [Layer III](@ref design-layer3), here for example

```julia
_retract(M::AbstractManifold, p, X, m::Exponentialretraction) = exp(M, p, X)
_retract(M::AbstractManifold, p, X, m::PolarRetraction) = retract_polar(M, p, X)
```

where the [`ExponentialRetraction`](@ref) is resolved by again calling a function on [Layer I](@ref design-layer1) (to fill futher default values if these exist). The [`PolarRetraction`](@ref) is dispatched to [`retract_polar`](@ref ManifoldsBase.retract_polar), a function on [Layer III](@ref design-layer3).

For further details and dispatches, see [retractions and inverse retractions](@ref sec-retractions) for an overview.

!!! note
    The documentation should be attached to the high level functions, since this again fosters ease of use.
    If you implement a polar retraction, you should write a method of function [`retract_polar`](@ref ManifoldsBase.retract_polar) but the doc string should be attached to `retract(::M, ::P, ::V, ::PolarRetraction)` for your types `::M, ::P, ::V` of the manifold, points and vectors, respectively.

To summarize, with respect to the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) paradigm, this layer dispatches the (optional) _parameters second_.

### [Layer III: The base layer with focus on implementations](@id design-layer3)

This lower level aims for the actual implementation of the function avoiding ambiguities.
It should have as few as possible optional parameters and as concrete as possible types for these.

This means

* the function name should be similar to its high level parent (for example [`retract`](@ref) and [`retract_polar`](@ref ManifoldsBase.retract_polar)  above)
* The manifold type in method signature should always be as narrow as possible.
* The points/vectors should either be untyped (for the default representation or if there is only one implementation) or provide all type bounds (for second representations or when using [`AbstractManifoldPoint`](@ref) and [`TVector`](@ref TVector), respectively).

The first step that often happens on this level is memory allocation and calling the mutating function. If faster, it might also implement the function at hand itself.

Usually functions from this layer are not exported, when they have an analogue on the first layer. For example the function [`retract_polar`](@ref ManifoldsBase.retract_polar)`(M, p, X)` is not exported, since when using the interface one would use the [`PolarRetraction`](@ref) or to be precise call [`retract`](@ref)`(M, p, X, PolarRetraction())`.
When implementing your own manifold, you have to import functions like these anyway.

To summarize, with respect to the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) paradigm, this layer dispatches the _concrete manifold and point/vector types last_.

## [Mutating and allocating functions](@id mutating-and-nonmutating)

Every function, where this is applicable should provide a mutating and an allocating variant.
For example for the exponential map `exp(M, p, X)` returns a _new_ point `q` where the result is computed in.
On the other hand `exp!(M, q, p, X)` computes the result in place of `q`, where the design of the implementation
should keep in mind that also `exp!(M, p, p, X)` should correctly overwrite `p`.

The interface provides a way to determine the allocation type and a result to compute/allocate
the resulting memory, such that the default implementation allocating functions, like [`exp`](@ref) is to allocate the resulting memory and call [`exp!`](@ref).

!!! note
    it might be useful to provide two distinct implementations, for example when using AD schemes.
    The default is meant for ease of use (concerning implementation), since then one has to just implement the mutating variants.

Non-mutating functions in `ManifoldsBase.jl` are typically implemented using mutating variants after a suitable allocation of memory.

Not that this allocation usually takes place only on [Layer III](@ref design-layer3) when dispatching on points.
Both [Layer I](@ref design-layer1) and [Layer II](@ref design-layer1) are usually implemented for both variants in parallel.

### Allocation of new points and vectors

The [`allocate`](@ref) function behaves like `similar` for simple representations of points and vectors (for example `Array{Float64}`).
For more complex types, such as nested representations of [`PowerManifold`](@ref) (see [`NestedPowerRepresentation`](@ref)), checked types like [`ValidationMPoint`](@ref) and more it operates differently.
While `similar` only concerns itself with the higher level of nested structures, `allocate` maps itself through all levels of nesting until a simple array of numbers is reached and then calls `similar`.
The difference can be most easily seen in the following example:

```julia
julia> x = similar([[1.0], [2.0]])
2-element Array{Array{Float64,1},1}:
 #undef
 #undef

julia> y = allocate([[1.0], [2.0]])
2-element Array{Array{Float64,1},1}:
 [6.90031725726027e-310]
 [6.9003678131654e-310]

julia> x[1]
ERROR: UndefRefError: access to undefined reference
Stacktrace:
 [1] getindex(::Array{Array{Float64,1},1}, ::Int64) at ./array.jl:744
 [2] top-level scope at REPL[12]:1

julia> y[1]
1-element Array{Float64,1}:
 6.90031725726027e-310
```

The function [`allocate_result`](@ref ManifoldsBase.allocate_result) allocates a correct return value. It takes into account the possibility that different arguments may have different numeric [`number_eltype`](@ref) types thorough the [`allocate_result_type`](@ref ManifoldsBase.allocate_result_type) function.
The most prominent example of the usage of this function is the logarithmic function [`log`](@ref) when used with typed points.
Lets assume on a manifold `M` the have points of type `P` and corresponding tangent vector types `V`.
then the logarithmic map has the signature

```julia
log(::M, ::P, ::P)
```

but the return type would be ``V``, whose internal sizes (fields/arrays) will depend on the concrete type of one of the points. This is accomplished by implementing a method `allocate_result(::M, ::typeof(log), ::P, ::P)` that returns the concrete variable for the result. This way, even with specific types, one just has to implement `log!` and the one line for the allocation.

!!! note
    This dispatch from the allocating to the mutating variant happens in Layer III, that is, functions like `exp` or [`retract_polar`](@ref ManifoldsBase.retract_polar) (but not [`retract`](@ref) itself) allocate their result (using `::typeof(retract)` for the second function)
    and call the mutating variant `exp!` and [`retract_polar!`](@ref ManifoldsBase.retract_polar!) afterwards.
