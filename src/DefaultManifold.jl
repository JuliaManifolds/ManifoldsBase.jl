"""
    DefaultManifold <: AbstractManifold

This default manifold illustrates the main features of the interface and provides a skeleton
to build one's own manifold. It is a simplified/shortened variant of `Euclidean` from
`Manifolds.jl`.

This manifold further illustrates how to type your manifold points and tangent vectors. Note
that the interface does not require this, but it might be handy in debugging and educative
situations to verify correctness of involved variabes.
"""
struct DefaultManifold{T<:Tuple,𝔽} <: AbstractManifold{𝔽} end
function DefaultManifold(n::Vararg{Int,N}; field = ℝ) where {N}
    return DefaultManifold{Tuple{n...},field}()
end



function check_approx(M::DefaultManifold, p, q; kwargs...)
    res = isapprox(p, q; kwargs...)
    res && return nothing
    v = distance(M, p, q)
    s = "The two points $p and $q on $M are not (approximately) equal."
    return ApproximatelyError(v, s)
end

function check_approx(M::DefaultManifold, p, X, Y; kwargs...)
    res = isapprox(X, Y; kwargs...)
    res && return nothing
    v = norm(M, p, X - Y)
    s = "The two tangent vectors $X and $Y in the tangent space at $p on $M are not (approximately) equal."
    return ApproximatelyError(v, s)
end

distance(::DefaultManifold, p, q) = norm(p - q)

embed!(::DefaultManifold, q, p) = copyto!(q, p)

embed!(::DefaultManifold, Y, p, X) = copyto!(Y, X)

exp!(::DefaultManifold, q, p, X) = (q .= p .+ X)

function get_basis_orthonormal(::DefaultManifold, p, N::RealNumbers)
    return CachedBasis(
        DefaultOrthonormalBasis(N),
        [_euclidean_basis_vector(p, i) for i in eachindex(p)],
    )
end
function get_basis_orthogonal(::DefaultManifold, p, N::RealNumbers)
    return CachedBasis(
        DefaultOrthogonalBasis(N),
        [_euclidean_basis_vector(p, i) for i in eachindex(p)],
    )
end
function get_basis_default(::DefaultManifold, p, N::RealNumbers)
    return CachedBasis(
        DefaultBasis(N),
        [_euclidean_basis_vector(p, i) for i in eachindex(p)],
    )
end
function get_basis_diagonalizing(M::DefaultManifold, p, B::DiagonalizingOrthonormalBasis)
    vecs = get_vectors(M, p, get_basis(M, p, DefaultOrthonormalBasis()))
    eigenvalues = zeros(real(eltype(p)), manifold_dimension(M))
    return CachedBasis(B, DiagonalizingBasisData(B.frame_direction, eigenvalues, vecs))
end

# Complex manifold, real basis -> coefficients c are complesx -> reshape
# Real manifold, real basis -> reshape
function get_coordinates_orthonormal!(M::DefaultManifold, c, p, X, ::RealNumbers)
    return copyto!(c, reshape(X, number_of_coordinates(M, ℝ)))
end
function get_coordinates_diagonalizing!(
    M::DefaultManifold,
    c,
    p,
    X,
    ::DiagonalizingOrthonormalBasis{ℝ},
)
    return copyto!(c, reshape(X, number_of_coordinates(M, ℝ)))
end
function get_coordinates_orthonormal!(::DefaultManifold, c, p, X, ::ComplexNumbers)
    m = length(X)
    return copyto!(c, [reshape(real(X), m); reshape(imag(X), m)])
end
function get_vector_orthonormal!(M::DefaultManifold, Y, p, c, ::RealNumbers)
    return copyto!(Y, reshape(c, representation_size(M)))
end
function get_vector_diagonalizing!(
    M::DefaultManifold,
    Y,
    p,
    c,
    ::DiagonalizingOrthonormalBasis{ℝ},
)
    return copyto!(Y, reshape(c, representation_size(M)))
end
function get_vector_orthonormal!(
    M::DefaultManifold{T,ℂ},
    Y,
    p,
    c,
    ::ComplexNumbers,
) where {T}
    n = div(length(c), 2)
    return copyto!(Y, reshape(c[1:n] + c[(n + 1):(2n)] * 1im, representation_size(M)))
end

injectivity_radius(::DefaultManifold) = Inf

@inline inner(::DefaultManifold, p, X, Y) = dot(X, Y)

is_flat(::DefaultManifold) = true

log!(::DefaultManifold, Y, p, q) = (Y .= q .- p)

@generated function manifold_dimension(::DefaultManifold{T,𝔽}) where {T,𝔽}
    return length(T.parameters) == 0 ? 1 : *(T.parameters...) * real_dimension(𝔽)
end

number_system(::DefaultManifold{T,𝔽}) where {T,𝔽} = 𝔽

norm(::DefaultManifold, p, X) = norm(X)

project!(::DefaultManifold, q, p) = copyto!(q, p)
project!(::DefaultManifold, Y, p, X) = copyto!(Y, X)

@generated representation_size(::DefaultManifold{T}) where {T} = tuple(T.parameters...)

function Base.show(io::IO, ::DefaultManifold{N,𝔽}) where {N,𝔽}
    return print(io, "DefaultManifold($(join(N.parameters, ", ")); field = $(𝔽))")
end

function parallel_transport_to!(::DefaultManifold, Y, p, X, q)
    return copyto!(Y, X)
end

function Random.rand!(::DefaultManifold, pX; σ = one(eltype(pX)), vector_at = nothing)
    pX .= randn(size(pX)) .* σ
    return pX
end
function Random.rand!(
    rng::AbstractRNG,
    ::DefaultManifold,
    pX;
    σ = one(eltype(pX)),
    vector_at = nothing,
)
    pX .= randn(rng, size(pX)) .* σ
    return pX
end

function riemann_tensor!(::DefaultManifold, Xresult, p, X, Y, Z)
    return fill!(Xresult, 0)
end

zero_vector(::DefaultManifold, p) = zero(p)

zero_vector!(::DefaultManifold, Y, p) = fill!(Y, 0)
