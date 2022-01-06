## Vector transport

Similar to the [exponential and logarithmic map](@ref exp-and-log) also the [parallel transport](@ref subsec-parallel-transport) might be costly to compute, especially when there is no closed form solution known and it has to be approximated with numerical methods.

Similar to the [retraction and its inverse](@ref sec-retractions) the generalisation of the parallel transport can be phrased as follows

A _vector transport_ is a way to transport a vector between two tangent spaces.
Let ``p,q ∈ \mathcal M`` be given, ``c`` the curve along which we want to transport (cf. [parallel transport](@ref subsec-parallel-transport), for example a geodesic or a geodesic or curve given by a retraction, and ``X ∈ T_p\mathcal M`` be a tangent vector.
Then ``T_{p→q}X`` is a vector transport if

1. Consistency. ``T_{p→p}`` is the identiy on ``T_p\mathcal M``.
2. Underlying curve ``T_{p→c(t)}X ∈ T_{c(t)}\mathcal M`` for ``t∈[0,1]``
3. Linearity. ``T_{p→q}(X+Y)=T_{p→q}X + T_{p→q}Y``

hold.

Currently the following types of vector transport are defined in `ManifoldsBase.jl`.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["vector_transport.jl"]
Order = [:function]
Public=true
Private=false
```

## Types of Retractions

To distinguish different types of retractions, the last argument of the (inverse) retraction
specifies a type. The following ones are available.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["vector_transport.jl"]
Order = [:type]
```

## The lower layer functions

While you should always add your documentation to the first layer vector transport methods above when implementing new manifolds, the actual implementation happens on the following functions on [the lower layer](@ref design-layer3).

```@autodocs
Modules = [ManifoldsBase]
Pages = ["vector_transport.jl"]
Order = [:function]
Public = false
Private = true
```
