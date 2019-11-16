@doc doc"""
    Euclidean{T<:Tuple} <: Manifold

Euclidean vector space $\mathbb R^n$.

# Constructor

    Euclidean(n)

generates the $n$-dimensional vector space $\mathbb R^n$.

   Euclidean(m, n)

generates the $mn$-dimensional vector space $\mathbb R^{m \times n}$, whose
elements are interpreted as $m \times n$ matrices.
"""
struct Euclidean{T<:Tuple} <: Manifold where {T} end

Euclidean(n::Int) = Euclidean{Tuple{n}}()
Euclidean(m::Int, n::Int) = Euclidean{Tuple{m,n}}()

function representation_size(::Euclidean{Tuple{n}}) where {n}
    return (n,)
end

function representation_size(::Euclidean{Tuple{m,n}}) where {m,n}
    return (m,n)
end

@generated manifold_dimension(::Euclidean{T}) where {T} = *(T.parameters...)

@inline inner(::Euclidean, x, v, w) = dot(v, w)

distance(::Euclidean, x, y) = norm(x-y)
norm(::Euclidean, x, v) = norm(v)

exp!(M::Euclidean, y, x, v) = (y .= x .+ v)

log!(M::Euclidean, v, x, y) = (v .= y .- x)

function zero_tangent_vector!(M::Euclidean, v, x)
    fill!(v, 0)
    return v
end

function project_point!(M::Euclidean, y, x)
    copyto!(y, x)
    return y
end

function project_tangent!(M::Euclidean, w, x, v)
    w .= v
    return w
end
