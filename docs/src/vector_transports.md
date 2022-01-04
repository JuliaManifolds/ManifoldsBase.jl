## Vector transport

TODO

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
