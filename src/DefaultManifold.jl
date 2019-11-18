"""
    DefaultManifold <: Manifold

This default manifold illustrates the main features of the interface and provides
a good starting point for an own manifold. It is a simplified/shortened variant
of `Euclidean` from `Manifolds.jl`.

this manifold further illustrates, how to type your manifold points
and tangent vectors. Note that the interface does not require this, but it
might be handy in debugging and educative situations to verify correctness of
involved variabes.
"""
struct DefaultManifold{T<:Tuple} <: Manifold where {T} end
DefaultManifold(n::Vararg{Int,N}) where N = DefaultManifold{Tuple{n...}}()

@generated representation_size(::DefaultManifold{T}) where {T} = Tuple(T.parameters...)
@generated manifold_dimension(::DefaultManifold{T}) where {T} = *(T.parameters...)
@inline inner(::DefaultManifold, x, v, w) = dot(v, w)
distance(::DefaultManifold, x, y) = norm(x-y)
norm(::DefaultManifold, x, v) = norm(v)
exp!(M::DefaultManifold, y, x, v) = (y .= x .+ v)
log!(M::DefaultManifold, v, x, y) = (v .= y .- x)
zero_tangent_vector!(M::DefaultManifold, v, x) = fill!(v, 0)
project_point!(M::DefaultManifold, y, x) = (y .= x)
project_tangent!(M::DefaultManifold, w, x, v) = (w .= v)
vector_transport_to!(M::DefaultManifold, vto, x, v, y, ::ParallelTransport) = (vto .= v)