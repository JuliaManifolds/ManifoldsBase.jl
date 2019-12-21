"""
    DefaultManifold <: Manifold

This default manifold illustrates the main features of the interface and
provides a skeleton to build one's own manifold. It is a simplified/shortened
variant of `Euclidean` from `Manifolds.jl`.

This manifold further illustrates how to type your manifold points and tangent
vectors. Note that the interface does not require this, but it might be handy in
debugging and educative situations to verify correctness of involved variabes.
"""
struct DefaultManifold{T<:Tuple} <: Manifold where {T} end
DefaultManifold(n::Vararg{Int,N}) where {N} = DefaultManifold{Tuple{n...}}()

@generated representation_size(::DefaultManifold{T}) where {T} = Tuple(T.parameters...)
@generated manifold_dimension(::DefaultManifold{T}) where {T} = *(T.parameters...)
@inline inner(::DefaultManifold, x, v, w) = dot(v, w)
distance(::DefaultManifold, x, y) = norm(x - y)
norm(::DefaultManifold, x, v) = norm(v)
exp!(::DefaultManifold, y, x, v) = (y .= x .+ v)
log!(::DefaultManifold, v, x, y) = (v .= y .- x)
zero_tangent_vector!(::DefaultManifold, v, x) = fill!(v, 0)
project_point!(::DefaultManifold, y, x) = copyto!(y, x)
project_tangent!(::DefaultManifold, w, x, v) = copyto!(w, v)
function vector_transport_to!(::DefaultManifold, vto, x, v, y, ::ParallelTransport)
    return copyto!(vto, v)
end
injectivity_radius(::DefaultManifold) = Inf
