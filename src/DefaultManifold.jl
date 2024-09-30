"""
    DefaultManifold <: AbstractManifold

This default manifold illustrates the main features of the interface and provides a skeleton
to build one's own manifold. It is a simplified/shortened variant of `Euclidean` from
`Manifolds.jl`.

This manifold further illustrates how to type your manifold points and tangent vectors. Note
that the interface does not require this, but it might be handy in debugging and educative
situations to verify correctness of involved variables.

# Constructor

    DefaultManifold(n::Int...; field = ℝ, parameter::Symbol = :field)


Arguments:

- `n`: shape of array representing points on the manifold.
- `field`: field over which the manifold is defined. Either `ℝ`, `ℂ` or `ℍ`.
- `parameter`: whether a type parameter should be used to store `n`. By default size
  is stored in a field. Value can either be `:field` or `:type`.
"""
struct DefaultManifold{𝔽,T} <: AbstractManifold{𝔽}
    size::T
end
function DefaultManifold(n::Vararg{Int}; field = ℝ, parameter::Symbol = :field)
    size = wrap_type_parameter(parameter, n)
    return DefaultManifold{field,typeof(size)}(size)
end

function allocation_promotion_function(
    ::DefaultManifold{ℂ},
    ::Union{typeof(get_vector),typeof(get_coordinates)},
    ::Tuple,
)
    return complex
end

change_representer!(M::DefaultManifold, Y, ::EuclideanMetric, p, X) = copyto!(M, Y, p, X)

change_metric!(M::DefaultManifold, Y, ::EuclideanMetric, p, X) = copyto!(M, Y, p, X)

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

for fname in [:get_basis_orthonormal, :get_basis_orthogonal, :get_basis_default]
    BD = Dict(
        :get_basis_orthonormal => :DefaultOrthonormalBasis,
        :get_basis_orthogonal => :DefaultOrthogonalBasis,
        :get_basis_default => :DefaultBasis,
    )
    BT = BD[fname]
    @eval function $fname(::DefaultManifold{ℝ}, p, N::RealNumbers)
        return CachedBasis($BT(N), [_euclidean_basis_vector(p, i) for i in eachindex(p)])
    end
    @eval function $fname(::DefaultManifold{ℂ}, p, N::ComplexNumbers)
        return CachedBasis($BT(N), [_euclidean_basis_vector(p, i) for i in eachindex(p)])
    end
end
function get_basis_diagonalizing(M::DefaultManifold, p, B::DiagonalizingOrthonormalBasis)
    vecs = get_vectors(M, p, get_basis(M, p, DefaultOrthonormalBasis()))
    eigenvalues = zeros(real(eltype(p)), manifold_dimension(M))
    return CachedBasis(B, DiagonalizingBasisData(B.frame_direction, eigenvalues, vecs))
end

# Complex manifold, real basis -> coefficients c are complex -> reshape
# Real manifold, real basis -> reshape
function get_coordinates_orthonormal!(M::DefaultManifold, c, p, X, N::AbstractNumbers)
    return copyto!(c, reshape(X, number_of_coordinates(M, N)))
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
function get_coordinates_orthonormal!(::DefaultManifold{ℂ}, c, p, X, ::RealNumbers)
    m = length(X)
    return copyto!(c, [reshape(real(X), m); reshape(imag(X), m)])
end
function get_vector_orthonormal!(M::DefaultManifold, Y, p, c, ::AbstractNumbers)
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
function get_vector_orthonormal!(M::DefaultManifold{ℂ}, Y, p, c, ::RealNumbers)
    n = div(length(c), 2)
    return copyto!(Y, reshape(c[1:n] + c[(n + 1):(2n)] * 1im, representation_size(M)))
end

injectivity_radius(::DefaultManifold) = Inf

@inline inner(::DefaultManifold, p, X, Y) = dot(X, Y)

is_flat(::DefaultManifold) = true

log!(::DefaultManifold, Y, p, q) = (Y .= q .- p)

function manifold_dimension(M::DefaultManifold{𝔽}) where {𝔽}
    size = get_parameter(M.size)
    return prod(size) * real_dimension(𝔽)
end

number_system(::DefaultManifold{𝔽}) where {𝔽} = 𝔽

norm(::DefaultManifold, p, X) = norm(X)

project!(::DefaultManifold, q, p) = copyto!(q, p)
project!(::DefaultManifold, Y, p, X) = copyto!(Y, X)

representation_size(M::DefaultManifold) = get_parameter(M.size)

function Base.show(io::IO, M::DefaultManifold{𝔽,<:TypeParameter}) where {𝔽}
    return print(
        io,
        "DefaultManifold($(join(get_parameter(M.size), ", ")); field = $(𝔽), parameter = :type)",
    )
end
function Base.show(io::IO, M::DefaultManifold{𝔽}) where {𝔽}
    return print(io, "DefaultManifold($(join(get_parameter(M.size), ", ")); field = $(𝔽))")
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

sectional_curvature_max(::DefaultManifold) = 0.0
sectional_curvature_min(::DefaultManifold) = 0.0

Weingarten!(::DefaultManifold, Y, p, X, V) = fill!(Y, 0)

zero_vector(::DefaultManifold, p) = zero(p)

zero_vector!(::DefaultManifold, Y, p) = fill!(Y, 0)
