## Vector transport

Similar to the [exponential and logarithmic map](@ref exp-and-log) also the [parallel transport](@ref subsec-parallel-transport) might be costly to compute, especially when there is no closed form solution known and it has to be approximated with numerical methods.

Similar to the [retraction and its inverse](@ref sec-retractions) the generalisation of the parallel transport can be phrased as follows

A _vector transport_ is a way to transport a vector between two tangent spaces.
Let ``p,q ∈ \mathcal M`` be given, ``c`` the curve along which we want to transport (cf. [parallel transport](@ref subsec-parallel-transport), for example a geodesic or a geodesic or curve given by a retraction. We can speficy the geodesic or curve a retraction realises for example by a direction ``d``.

More precisely using [^AbsilMahonySepulchre2008], Def. 8.1.1, a vector transport
``T_{p,d}: T_p\mathcal M \to T_q\mathcal M``, ``p∈ \mathcal M``, ``Y∈ T_p\mathcal M`` is a smooth mapping
associated to a retraction ``\operatorname{retr}_p(Y) = q`` such that

1. (associated retraction) ``\mathcal T_{p,d}X ∈ T_q\mathcal M`` if and only if ``q = \operatorname{retr}_p(d)``.
2. (consistency) ``\mathcal T_{p,0_p}X = X`` for all ``X∈T_p\mathcal M``
3. (linearity) ``\mathcal T_{p,d}(αX+βY) = \mathcal T_{p,d}αX + \mathcal T_{p,d}βY``

hold.

Currently the following methods for vector transport are defined in `ManifoldsBase.jl`.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["vector_transport.jl"]
Order = [:function]
Public=true
Private=false
```

## Types of vector transports

To distinguish different types of vector transport we introduce the [`AbstractVectorTransportMethod`](@ref). The following concrete types are available.

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
