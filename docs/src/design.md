# Main Design Principles

The interface for a manifold is defined to be as generic as possible, such that applications can be implemented as independently as possible from an actual manifold.
This way, algorithms like those from [`Manopt.jl`](https://manoptjl.org) can be implemented on _arbitrary_ manifolds.

The main design criteria for the interface are:

* Aims to also provide _efficient_ _global state-free_, both _in-place_ and _out-of-place_ computations whenever possible.
* Provide a high level interface that is easy to use.

Therefore this interface has 3 main features, that we will explain using two (related)
concepts, the [exponential map](https://en.wikipedia.org/wiki/Exponential_map_(Riemannian_geometry)) that maps a tangent vector ``X`` at a point ``p`` to a point ``q`` or mathematically ``\exp_p:T_p\mathcal M \to \mathcal M`` and its generalization, a [`retract`](@ref)ion ``\operatorname{retr}_p`` with same domain and range.

You do not need to know their exact definition at this point, just that there is _one_ exponential map on a Riemannian manifold, and several retractions, where one of them is the exponential map (sometime called exponential retraction for completeness). Every retraction has its own subtype of the [`AbstractRetractionMethod`](@ref) that uniquely defines it.

The following three design patterns aim to fulfill the criteria from above, while
also avoiding ambiguities in multiple dispatch using the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) approach.

## General order of parameters

Since the central element for functions on a manifold is the manifold itself, it should always be the first parameter, even for mutating functions. Then the classical parametzers of a function (for example a point and a tangent vector for the retraction) follow and the final part are parameters to further dispatch on, which usually have their defaults.

## A 3-Layer architecture for dispatch

The general architecture consists of three layers

* The high level interface for ease of use – and to dispatch on other manifolds
* an interims layer to dispatch on different parameters in the last section
* the lowest layer with a specific manifold and no optional parameters for performance

These three layers are described in more detail in the following.

### Layer I: The high level interface and ease of use

THe highest layer for convenience of decorators.
A usual scheme is, that a manifold might assume several things implicitly, for example the default implementation of the sphere $\mathbb S^n$ using unit vectors in $\mathbb R^{n+1}$.
The embedding can be explicitly used to avoid reimplementations – the inner product can be “passed on” to its embedding.

To do so, we “decorate” the manifold by making it an [`AbstractDecoratorManifold`](@ref) and activating the right traits see [the example](@ref example).

The explicit case of the [`EmbeddedManifold`](@ref) can be used to distinguish different embeddings of a manifold, but also their dispatch (onto the manifold or its embedding, depending on the type of embedding) happens here.

Note that all other parameters of a function should be as unspecific as possible on this layer.

With respect to the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) paradigm, this layer dispatches the _manifold first_, but here we stay oon an abstract type level.

This layer ends usually in calling the same functions like `retract` but prefixed with a `_` to enter Layer II.

### Layer II: An internal dispatch interface for parameters

This layer is an interims layer to dispatch on the (optional/default) parameters of a function like the retraction:
[`retract`](@ref) has a last parameter that determines the type.
The last function in the previous layer calls `_retract`, which is an internal function.

On this layer, e.g. for `_retract` only these last parameters should be typed, the manifold should stay at the [`AbstractManifold`](@ref) level. It dispatches on different functions per existing parameter type (and might pass this one further on, if it has fields).

Note that this layer is an internal one. It is automatically called for functions with parameters to dispatch on.

It should only be extended when introducing new such parameter types, for example when introducing a new type of a retraction.

The functions from this layer should never be called directly, are hence also not exported and carry the `_` prefix.
They should only be called as the final step in the previous layer.

If the default parameters are not dispatched per type, using `_` might be skipped.
The following resolution might even be seen as a last step in layer I or the resolution here in layer II.

```julia
exp(M::AbstractManifold, p, X, t::Real) = exp(M, p, t * X)
```

When there is no dispatch for different types of the optional parameter (here `t`), the `_` might be skipped.
One could hence see the last code line as a definition on Layer I that passes directly to Layer III, since there are not parameter to dispatch on.

To close this section, let‘s look at an example. The high level (or level I) definition of the retraction is given by

```julia
retract(M::AbstractManifold, p, X, m::AbstractRetractionMethod=default_retraction_method(M)) = _retract(M,p,X,m)
```

This level now dispatches on different retraction types. It usually passes to specific functions implemented in Layer III,
here for example

```julia
_retract(M::AbstractManifold, p, X, m::Exponentialretraction) = exp(M,p,X)
_retract(M::AbstractManifold, p, X, m::PolarRetractionMethod) = retract_polar(M,p,X)
```

or the [`PolarRetractionMethod`](@ref) which dispatches to [`retract_polar`](@ref).

For further details and dispatches, see the [appendix](@ref subsec_appendix_retr) for an overview.

!!! note
    The documentation should be attached to the high level functions, since this again fosters ease of use.
    If yuo implement a polar retraction, you should write a function `polar_retract` but the doc string should be attached to `retract(::M, ::P, ::V, ::PolarRetraction)` for your types `::M, ::P, ::V` of the manifold, points and vectors, respectively.

To summarize, with respect to the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) paradigm, this layer dispatches the (optional) _parameters second_.

### Layer III: The lower level interface to gain performance

This lower level aims for performance, that is, any function should have as few as possible optional and keyword arguments
and be typed as concrete as possible/necessary. This means

* the function name should be similar to its high level parent (for example `retract` and `retract_polar`from above)
* The manifold type in method signature should always be as narrow as possible.
* the points/vectors should either be untyped (for the default representation of if there is only one) or provide all types concretely.

The first thing to do on this level is the aforementioned default to pass from allocating to mutating functions.

Note that not all of these functions are exported, but if you implement for example an existing retraction on a new manifold, you will have to import this function.

To summarize, with respect to the [dispatch on one argument at a time](https://docs.julialang.org/en/v1/manual/methods/#Dispatch-on-one-argument-at-a-time) paradigm, this layer dispatches the _concrete manifold and point/vector types last_.

## Mutating and allocating functions

Every function, where this is applicable should provide a mutating and an allocating variant.
For example for the exponential map `exp(M,p,x)` returns a _new_ point `q` where the result is computed in.
On the other hand `exp!(M, q, p, X)` computes the result in place of `q`, where the design of the implementation
should keep in mind that also `exp!(M,p,p,X)` should correctly overwrite `p`.

The interface provides a way to determine the allocation type and a result to compute/allocate
the resulting memory, such that the default implementation allocating functions, like [`exp`](@ref) is to allocate the resulting memory and call [`exp!`](@ref).

!!! note
    it might be useful to provide two distinct implementations, for example when using AD schemes.
    The default is meant for ease of use (concerning implementation), since then one has to just implement the mutating variants.

Non-mutating functions in `ManifoldsBase.jl` are typically implemented using mutating variants.
Allocation of new points is performed using a custom mechanism that relies on the following functions:

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

but the return type would be ``V``, whose internal sizes (fields/arrays) will depend on the concrete type of one of the points. This is accomplished by omplementing a `allocate_result(::M, ::typeof(log), ::P, ::P)`that returns the concrete variable for the return. This way, even with specific types, one just has to implement `log!` and the one line for the allocation.

!!! note
    This dispatch from the allocating to the mutating variant happens in Layer III, that is, functions like `exp` or `retract_polar` (but not `retract` itself) allocate their result (using `::typeof(retract)` for the second function)
    and call the mutating variant `exp!` and `retract_polar!` afterwards.

## Appendix

### Validations

The function [`is_point`](@ref) internally calls the lower level function [`check_point`](@ref). Similarly [`is_vector`](@ref) calls [`check_vector`](@ref), which assumes that the (base) point is correct.

### [Inverse Retractions](@id subsec_appendix_inv_retr)

The high level function `inverse_retract(::M, p, X, m::AbstractInverseRetractionMethod)`
as well as its mutating variant first dispatch on the lower level, before the non-mutating variant (of the name below) allocates memory and calls its mutating variant.

The following table provides an overview of the currently available types and their lower level functions.

| Name | default lower level function | comment |
| :--- | :----------------------------- | :----- |
| [`PolarInverseRetraction`](@ref) | `inverse_retract_polar` |
| [`ProjectionInverseRetraction`](@ref) | `inverse_retract_project` |
| [`QRInverseRetraction`](@ref) | `inverse_retract_qr` |
| [`NLSolveInverseRetraction`](@ref) | `inverse_retract_nlsolve` | the `m` is also passed on here. |

### [Retractions](@id subsec_appendix_retr)

The high level function `retract(::M, p, X, m::AbstractRetractionMethod)`
as well as its mutating variant first dispatch on the lower level, before the non-mutating variant (of the name below) allocates memory and calls its mutating variant.

The following table provides an overview of the currently available types and their lower level functions.

| Name | default lower level function | comment |
| :--- | :----------------------------- | :----- |
| [`PolarRetraction`](@ref) | `inverse_retract_polar` |
| [`ProjectionRetraction`](@ref) | `inverse_retract_project` |
| [`QRRetraction`](@ref) | `inverse_retract_qr` |

### Vector transport

In the follwing table the `V` in the function names stand for `along`, `to`. Nte that `vector_transport_direction` by default uses the default retraction and `to`.

| Name | default lower level function | comment |
| :--- | :----------------------------- | :----- |
| [`DifferentiatedRetractionVectorTransport`](@ref) | `vector_transport_V_diff` | the inner retraction is passed on as last argument |
| [`ParallelTransport`](@ref) | `parallel_transport_V` |
| [`ProjectionTransport`](@ref) | `vector_transport_V_project` |
