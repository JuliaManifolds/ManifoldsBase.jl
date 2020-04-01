"""
    DefaultManifold <: Manifold

This default manifold illustrates the main features of the interface and provides a skeleton
to build one's own manifold. It is a simplified/shortened variant of `Euclidean` from
`Manifolds.jl`.

This manifold further illustrates how to type your manifold points and tangent vectors. Note
that the interface does not require this, but it might be handy in debugging and educative
situations to verify correctness of involved variabes.
"""
struct DefaultManifold{T<:Tuple, 𝔽} <: Manifold{𝔽} end
DefaultManifold(n::Vararg{Int,N}; field = ℝ) where {N} = DefaultManifold{Tuple{n...}, field}()

function check_manifold_point(M::DefaultManifold, p; kwargs...)
    if size(p) != representation_size(M)
        return DomainError(
            size(p),
            "The point $(p) does not lie on $M, since its size is not $(representation_size(M)).",
        )
    end
    return nothing
end

function check_tangent_vector(
    M::DefaultManifold,
    p,
    X;
    check_base_point = true,
    kwargs...,
)
    if check_base_point
        perr = check_manifold_point(M, p)
        perr === nothing || return perr
    end
    if size(X) != representation_size(M)
        return DomainError(
            size(X),
            "The vector $(X) is not a tangent to a point on $M since its size does not match $(representation_size(M)).",
        )
    end
    return nothing
end

distance(::DefaultManifold, x, y) = norm(x - y)

embed!(::DefaultManifold, y, x) = copyto!(y, x)

embed!(::DefaultManifold, w, x, v) = copyto!(w, v)


exp!(::DefaultManifold, y, x, v) = (y .= x .+ v)

function get_basis(M::DefaultManifold, p, B::DefaultOrthonormalBasis)
    return CachedBasis(B, [_euclidean_basis_vector(p, i) for i in eachindex(p)])
end
function get_basis(M::DefaultManifold, p, B::DefaultOrthogonalBasis)
    return CachedBasis(B, [_euclidean_basis_vector(p, i) for i in eachindex(p)])
end
function get_basis(M::DefaultManifold, p, B::DefaultBasis)
    return CachedBasis(B, [_euclidean_basis_vector(p, i) for i in eachindex(p)])
end
function get_basis(M::DefaultManifold, p, B::DiagonalizingOrthonormalBasis)
    vecs = get_vectors(M, p, get_basis(M, p, DefaultOrthonormalBasis()))
    eigenvalues = zeros(real(eltype(p)), manifold_dimension(M))
    return CachedBasis(B, DiagonalizingBasisData(B.frame_direction, eigenvalues, vecs))
end

function get_coordinates!(M::DefaultManifold, Y, p, X, B::DefaultOrthonormalBasis)
    copyto!(Y, reshape(X, manifold_dimension(M)))
    return Y
end

function get_vector!(M::DefaultManifold, Y, p, X, B::DefaultOrthonormalBasis)
    copyto!(Y, reshape(X, representation_size(M)))
    return Y
end

injectivity_radius(::DefaultManifold) = Inf

@inline inner(::DefaultManifold, x, v, w) = dot(v, w)

log!(::DefaultManifold, v, x, y) = (v .= y .- x)

@generated manifold_dimension(::DefaultManifold{T,𝔽}) where {T,𝔽} = *(T.parameters...)*real_dimension(𝔽)

number_system(::DefaultManifold{T,𝔽}) where {T,𝔽} = 𝔽

norm(::DefaultManifold, x, v) = norm(v)

project!(::DefaultManifold, y, x) = copyto!(y, x)
project!(::DefaultManifold, w, x, v) = copyto!(w, v)

@generated representation_size(::DefaultManifold{T}) where {T} = Tuple(T.parameters...)

function vector_transport_along!(
    ::DefaultManifold,
    vto,
    x,
    v,
    c,
    ::AbstractVectorTransportMethod,
)
    return copyto!(vto, v)
end

function vector_transport_to!(::DefaultManifold, vto, x, v, y, ::ParallelTransport)
    return copyto!(vto, v)
end

zero_tangent_vector!(::DefaultManifold, v, x) = fill!(v, 0)
