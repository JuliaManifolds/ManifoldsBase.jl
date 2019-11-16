struct DefaultManifold{T<:Tuple} <: Manifold where {T} end
DefaultManifold(n::Vararg{Int,N}) where N = DefaultManifold{Tuple{n...}}()
@generated representation_size(::DefaultManifold{T}) where {T} = Tuple(T.parameters...)
@generated manifold_dimension(::DefaultManifold{T}) where {T} = *(T.parameters...)
@inline inner(::DefaultManifold, x, v, w) = dot(v, w)
distance(::DefaultManifold, x, y) = norm(x-y)
norm(::DefaultManifold, x, v) = norm(v)
exp!(M::DefaultManifold, y, x, v) = (y .= x .+ v)
log!(M::DefaultManifold, v, x, y) = (v .= y .- x)
function zero_tangent_vector!(M::DefaultManifold, v, x)
    fill!(v, 0)
    return v
end
function project_point!(M::DefaultManifold, y, x)
    y .= x
    return y
end
function project_tangent!(M::DefaultManifold, w, x, v)
    w .= v
    return w
end
function vector_transport_to!(M::DefaultManifold, vto, x, v, y, ::ParallelTransport)
    vto .= v
    return vto
end