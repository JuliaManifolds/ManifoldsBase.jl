# Functions on manifolds

This page collects several basic functions on manifolds.

## [The exponential map, the logarithmic map, and geodesics](@id exp-and-log)

Geodesics are the generalizations of a straight line to manifolds, i.e. their intrinsic acceleration is zero.
Together with geodesics one also obtains the exponential map and its inverse, the logarithmic map.
Informally speaking, the exponential map takes a vector (think of a direction and a length) at one point and returns another point,
which lies towards this direction at distance of the specified length. The logarithmic map does the inverse, i.e. given two points, it tells which vector “points towards” the other point.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["exp_log_geo.jl"]
Order = [:function]
```

## [Parallel transport](@id subsec-parallel-transport)

While moving vectors from one base point to another is the identity in the Euclidean space – or in other words all tangent spaces (directions one can “walk” into) are the same. This is different on a manifold.

If we have two points ``p,q ∈ \mathcal M``, we take a ``c: [0,1] → \mathcal M`` connecting the two points, i.e. ``c(0) = p`` and ``c(1) = q``. this could be a (or the) geodesic.
If we further consider a vector field ``X: [0,1] → T\mathcal M``, i.e. where ``X(t) ∈ T_{c(t)}\mathcal M``.
Then the vector field is called _parallel_ if its covariant derivative ``\frac{\mathrm{D}}{\mathrm{d}t}X(t) = 0`` for all ``t∈ |0,1]``.

If we now impose a value for ``X=X(0) ∈ T_p\mathcal M``, we obtain an ODE with an initial condition.
The resulting value ``X(1) ∈ T_q\mathcal M`` is called the _parallel transport_ of `X` along ``c``
or in case of a geodesic the _parallel transport of `X` from `p` to `q`.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["parallel_transport.jl"]
Order = [:function]
```

## Further functions on manifolds

### General functions provided by the interface

```@autodocs
Modules = [ManifoldsBase]
Pages = ["ManifoldsBase.jl"]
Order = [:type, :function]
Public=true
Private=false
```

### Internal functions

While you should always add your documentation to functions from the last section, some of the functions dispatch onto functions on [layer III](@ref design-layer3). These are the ones
you usually implement for your manifold – unless there is no lower level function called, like for the [`manifold_dimension`](@ref).

```@autodocs
Modules = [ManifoldsBase]
Pages = ["ManifoldsBase.jl"]
Order = [:function]
Public=false
Private=true
```

## Error Messages

This interface introduces a small set of own error messages.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["errors.jl"]
Order = [:type]
```
