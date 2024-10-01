"""
    VectorSpaceType

Abstract type for tangent spaces, cotangent spaces, their tensor products,
exterior products, etc.

Every vector space `fiber` is supposed to provide:
* a method of constructing vectors,
* basic operations: addition, subtraction, multiplication by a scalar
  and negation (unary minus),
* `zero_vector(fiber, p)` to construct zero vectors at point `p`,
* `allocate(X)` and `allocate(X, T)` for vector `X` and type `T`,
* `copyto!(X, Y)` for vectors `X` and `Y`,
* `number_eltype(v)` for vector `v`,
* [`vector_space_dimension`](@ref).

Optionally:
* inner product via `inner` (used to provide Riemannian metric on vector
  bundles),
* [`flat`](https://juliamanifolds.github.io/Manifolds.jl/stable/features/atlases.html#Manifolds.flat-Tuple{AbstractManifold,%20Any,%20Any}) and [`sharp`](https://juliamanifolds.github.io/Manifolds.jl/stable/features/atlases.html#Manifolds.sharp-Tuple{AbstractManifold,%20Any,%20Any}),
* `norm` (by default uses `inner`),
* [`project`](@ref) (for embedded vector spaces),
* [`representation_size`](@ref),
* broadcasting for basic operations.
"""
abstract type VectorSpaceType <: FiberType end

"""
    struct TangentSpaceType <: VectorSpaceType end

A type that indicates that a [`Fiber`](@ref) is a [`TangentSpace`](@ref).
"""
struct TangentSpaceType <: VectorSpaceType end

"""
    struct CotangentSpaceType <: VectorSpaceType end

A type that indicates that a [`Fiber`](@ref) is a [`CotangentSpace`](@ref).
"""
struct CotangentSpaceType <: VectorSpaceType end

TCoTSpaceType = Union{TangentSpaceType,CotangentSpaceType}

"""
    AbstractBasis{𝔽,VST<:VectorSpaceType}

Abstract type that represents a basis of vector space of type `VST` on a manifold or
a subset of it.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

# See also

[`VectorSpaceType`](@ref)
"""
abstract type AbstractBasis{𝔽,VST<:VectorSpaceType} end

"""
    DefaultBasis{𝔽,VST<:VectorSpaceType}

An arbitrary basis of vector space of type `VST` on a manifold. This will usually
be the fastest basis available for a manifold.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

# See also

[`VectorSpaceType`](@ref)
"""
struct DefaultBasis{𝔽,VST<:VectorSpaceType} <: AbstractBasis{𝔽,VST}
    vector_space::VST
end
function DefaultBasis(𝔽::AbstractNumbers = ℝ, vs::VectorSpaceType = TangentSpaceType())
    return DefaultBasis{𝔽,typeof(vs)}(vs)
end
function DefaultBasis{𝔽}(vs::VectorSpaceType = TangentSpaceType()) where {𝔽}
    return DefaultBasis{𝔽,typeof(vs)}(vs)
end
function DefaultBasis{𝔽,TangentSpaceType}() where {𝔽}
    return DefaultBasis{𝔽,TangentSpaceType}(TangentSpaceType())
end

"""
    AbstractOrthogonalBasis{𝔽,VST<:VectorSpaceType}

Abstract type that represents an orthonormal basis of vector space of type `VST` on a
manifold or a subset of it.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

# See also

[`VectorSpaceType`](@ref)
"""
abstract type AbstractOrthogonalBasis{𝔽,VST<:VectorSpaceType} <: AbstractBasis{𝔽,VST} end

"""
    DefaultOrthogonalBasis{𝔽,VST<:VectorSpaceType}

An arbitrary orthogonal basis of vector space of type `VST` on a manifold. This will usually
be the fastest orthogonal basis available for a manifold.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

# See also

[`VectorSpaceType`](@ref)
"""
struct DefaultOrthogonalBasis{𝔽,VST<:VectorSpaceType} <: AbstractOrthogonalBasis{𝔽,VST}
    vector_space::VST
end
function DefaultOrthogonalBasis(
    𝔽::AbstractNumbers = ℝ,
    vs::VectorSpaceType = TangentSpaceType(),
)
    return DefaultOrthogonalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthogonalBasis{𝔽}(vs::VectorSpaceType = TangentSpaceType()) where {𝔽}
    return DefaultOrthogonalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthogonalBasis{𝔽,TangentSpaceType}() where {𝔽}
    return DefaultOrthogonalBasis{𝔽,TangentSpaceType}(TangentSpaceType())
end


struct VeeOrthogonalBasis{𝔽} <: AbstractOrthogonalBasis{𝔽,TangentSpaceType} end
VeeOrthogonalBasis(𝔽::AbstractNumbers = ℝ) = VeeOrthogonalBasis{𝔽}()

"""
    AbstractOrthonormalBasis{𝔽,VST<:VectorSpaceType}

Abstract type that represents an orthonormal basis of vector space of type `VST` on a
manifold or a subset of it.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

# See also

[`VectorSpaceType`](@ref)
"""
abstract type AbstractOrthonormalBasis{𝔽,VST<:VectorSpaceType} <:
              AbstractOrthogonalBasis{𝔽,VST} end

"""
    DefaultOrthonormalBasis(𝔽::AbstractNumbers = ℝ, vs::VectorSpaceType = TangentSpaceType())

An arbitrary orthonormal basis of vector space of type `VST` on a manifold. This will usually
be the fastest orthonormal basis available for a manifold.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

# See also

[`VectorSpaceType`](@ref)
"""
struct DefaultOrthonormalBasis{𝔽,VST<:VectorSpaceType} <: AbstractOrthonormalBasis{𝔽,VST}
    vector_space::VST
end

function DefaultOrthonormalBasis(
    𝔽::AbstractNumbers = ℝ,
    vs::VectorSpaceType = TangentSpaceType(),
)
    return DefaultOrthonormalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthonormalBasis{𝔽}(vs::VectorSpaceType = TangentSpaceType()) where {𝔽}
    return DefaultOrthonormalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthonormalBasis{𝔽,TangentSpaceType}() where {𝔽}
    return DefaultOrthonormalBasis{𝔽,TangentSpaceType}(TangentSpaceType())
end

"""
    ProjectedOrthonormalBasis(method::Symbol, 𝔽::AbstractNumbers = ℝ)

An orthonormal basis that comes from orthonormalization of basis vectors
of the ambient space projected onto the subspace representing the tangent space
at a given point.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

Available methods:
  - `:gram_schmidt` uses a modified Gram-Schmidt orthonormalization.
  - `:svd` uses SVD decomposition to orthogonalize projected vectors.
    The SVD-based method should be more numerically stable at the cost of
    an additional assumption (local metric tensor at a point where the
    basis is calculated has to be diagonal).
"""
struct ProjectedOrthonormalBasis{Method,𝔽} <: AbstractOrthonormalBasis{𝔽,TangentSpaceType} end

function ProjectedOrthonormalBasis(method::Symbol, 𝔽::AbstractNumbers = ℝ)
    return ProjectedOrthonormalBasis{method,𝔽}()
end

@doc raw"""
    GramSchmidtOrthonormalBasis{𝔽} <: AbstractOrthonormalBasis{𝔽}

An orthonormal basis obtained from a basis.

# Constructor

    GramSchmidtOrthonormalBasis(𝔽::AbstractNumbers = ℝ)
"""
struct GramSchmidtOrthonormalBasis{𝔽} <: AbstractOrthonormalBasis{𝔽,TangentSpaceType} end
GramSchmidtOrthonormalBasis(𝔽::AbstractNumbers = ℝ) = GramSchmidtOrthonormalBasis{𝔽}()

@doc raw"""
    DiagonalizingOrthonormalBasis{𝔽,TV} <: AbstractOrthonormalBasis{𝔽,TangentSpaceType}

An orthonormal basis `Ξ` as a vector of tangent vectors (of length determined by
[`manifold_dimension`](@ref)) in the tangent space that diagonalizes the curvature
tensor ``R(u,v)w`` and where the direction `frame_direction` ``v`` has curvature `0`.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
as coefficients in linear combinations of the basis vectors.

# Constructor

    DiagonalizingOrthonormalBasis(frame_direction, 𝔽::AbstractNumbers = ℝ)
"""
struct DiagonalizingOrthonormalBasis{𝔽,TV} <: AbstractOrthonormalBasis{𝔽,TangentSpaceType}
    frame_direction::TV
end
function DiagonalizingOrthonormalBasis(X, 𝔽::AbstractNumbers = ℝ)
    return DiagonalizingOrthonormalBasis{𝔽,typeof(X)}(X)
end
struct DiagonalizingBasisData{D,V,ET}
    frame_direction::D
    eigenvalues::ET
    vectors::V
end

const DefaultOrDiagonalizingBasis{𝔽} =
    Union{DefaultOrthonormalBasis{𝔽,TangentSpaceType},DiagonalizingOrthonormalBasis{𝔽}}

"""
    CachedBasis{𝔽,V,<:AbstractBasis{𝔽}} <: AbstractBasis{𝔽}

A cached version of the given `basis` with precomputed basis vectors. The basis vectors
are stored in `data`, either explicitly (like in cached variants of
[`ProjectedOrthonormalBasis`](@ref)) or implicitly.

# Constructor

    CachedBasis(basis::AbstractBasis, data)
"""
struct CachedBasis{𝔽,B,V} <:
       AbstractBasis{𝔽,TangentSpaceType} where {B<:AbstractBasis{𝔽,TangentSpaceType},V}
    data::V
end
function CachedBasis(::B, data::V) where {V,𝔽,B<:AbstractBasis{𝔽,<:TangentSpaceType}}
    return CachedBasis{𝔽,B,V}(data)
end
function CachedBasis(basis::CachedBasis) # avoid double encapsulation
    return basis
end
function CachedBasis(
    basis::DiagonalizingOrthonormalBasis,
    eigenvalues::ET,
    vectors::T,
) where {ET<:AbstractVector,T<:AbstractVector}
    data = DiagonalizingBasisData(basis.frame_direction, eigenvalues, vectors)
    return CachedBasis(basis, data)
end

# forward declarations
function get_coordinates end
function get_vector end

const all_uncached_bases{T} = Union{
    AbstractBasis{<:Any,T},
    DefaultBasis{<:Any,T},
    DefaultOrthogonalBasis{<:Any,T},
    DefaultOrthonormalBasis{<:Any,T},
}

function allocate_on(M::AbstractManifold, ::TangentSpaceType)
    return similar(Array{Float64}, representation_size(M))
end
function allocate_on(M::AbstractManifold, ::TangentSpaceType, T::Type{<:AbstractArray})
    return similar(T, representation_size(M))
end

"""
    allocate_coordinates(M::AbstractManifold, p, T, n::Int)

Allocate vector of coordinates of length `n` of type `T` of a vector at point `p`
on manifold `M`.
"""
allocate_coordinates(::AbstractManifold, p, T, n::Int) = allocate(p, T, n)
function allocate_coordinates(M::AbstractManifold, p::Int, T, n::Int)
    return (representation_size(M) == () && n == 0) ? zero(T) : zeros(T, n)
end
function allocate_result(
    M::AbstractManifold,
    f::typeof(get_coordinates),
    p,
    X,
    basis::AbstractBasis{𝔽},
) where {𝔽}
    T = coordinate_eltype(M, p, 𝔽)
    return allocate_coordinates(M, p, T, number_of_coordinates(M, basis))
end

@inline function allocate_result_type(
    M::AbstractManifold,
    f::typeof(get_vector),
    args::Tuple{Any,Vararg{Any}},
)
    apf = allocation_promotion_function(M, f, args)
    return apf(
        invoke(allocate_result_type, Tuple{AbstractManifold,Any,typeof(args)}, M, f, args),
    )
end

"""
    allocation_promotion_function(M::AbstractManifold, f, args::Tuple)

Determine the function that must be used to ensure that the allocated representation is of
the right type. This is needed for [`get_vector`](@ref) when a point on a complex manifold
is represented by a real-valued vectors with a real-coefficient basis, so that
a complex-valued vector representation is allocated.
"""
allocation_promotion_function(M::AbstractManifold, f, args::Tuple) = identity

"""
    change_basis(M::AbstractManifold, p, c, B_in::AbstractBasis, B_out::AbstractBasis)

Given a vector with coordinates `c` at point `p` from manifold `M` in basis `B_in`,
compute coordinates of the same vector in basis `B_out`.
"""
function change_basis(M::AbstractManifold, p, c, B_in::AbstractBasis, B_out::AbstractBasis)
    return get_coordinates(M, p, get_vector(M, p, c, B_in), B_out)
end

function change_basis!(
    M::AbstractManifold,
    c_out,
    p,
    c,
    B_in::AbstractBasis,
    B_out::AbstractBasis,
)
    return get_coordinates!(M, c_out, p, get_vector(M, p, c, B_in), B_out)
end

function combine_allocation_promotion_functions(f::T, ::T) where {T}
    return f
end
function combine_allocation_promotion_functions(::typeof(complex), ::typeof(identity))
    return complex
end
function combine_allocation_promotion_functions(::typeof(identity), ::typeof(complex))
    return complex
end

"""
    coordinate_eltype(M::AbstractManifold, p, 𝔽::AbstractNumbers)

Get the element type for 𝔽-field coordinates of the tangent space at a point `p` from
manifold `M`. This default assumes that usually complex bases of complex manifolds have
real coordinates but it can be overridden by a more specific method.
"""
@inline function coordinate_eltype(::AbstractManifold, p, 𝔽::ComplexNumbers)
    return complex(number_eltype(p))
end
@inline function coordinate_eltype(::AbstractManifold, p, ::RealNumbers)
    return real(number_eltype(p))
end

@doc raw"""
    dual_basis(M::AbstractManifold, p, B::AbstractBasis)

Get the dual basis to `B`, a basis of a vector space at point `p` from manifold `M`.

The dual to the ``i``th vector ``v_i`` from basis `B` is a vector ``v^i`` from the dual space
such that ``v^i(v_j) = δ^i_j``, where ``δ^i_j`` is the Kronecker delta symbol:
````math
δ^i_j = \begin{cases}
1 & \text{ if } i=j, \\
0 & \text{ otherwise.}
\end{cases}
````
"""
dual_basis(M::AbstractManifold, p, B::AbstractBasis) = _dual_basis(M, p, B)

function _dual_basis(
    ::AbstractManifold,
    p,
    ::DefaultOrthonormalBasis{𝔽,TangentSpaceType},
) where {𝔽}
    return DefaultOrthonormalBasis{𝔽}(CotangentSpaceType())
end
function _dual_basis(
    ::AbstractManifold,
    p,
    ::DefaultOrthonormalBasis{𝔽,CotangentSpaceType},
) where {𝔽}
    return DefaultOrthonormalBasis{𝔽}(TangentSpaceType())
end

function _euclidean_basis_vector(p::StridedArray, i)
    X = zero(p)
    X[i] = 1
    return X
end
function _euclidean_basis_vector(p, i)
    # when p is for example a SArray
    X = similar(p)
    copyto!(X, zero(p))
    X[i] = 1
    return X
end

"""
    get_basis(M::AbstractManifold, p, B::AbstractBasis; kwargs...) -> CachedBasis

Compute the basis vectors of the tangent space at a point on manifold `M`
represented by `p`.

Returned object derives from [`AbstractBasis`](@ref) and may have a field `.vectors`
that stores tangent vectors or it may store them implicitly, in which case
the function [`get_vectors`](@ref) needs to be used to retrieve the basis vectors.

See also: [`get_coordinates`](@ref), [`get_vector`](@ref)
"""
function get_basis(M::AbstractManifold, p, B::AbstractBasis; kwargs...)
    return _get_basis(M, p, B; kwargs...)
end

function _get_basis(::AbstractManifold, ::Any, B::CachedBasis)
    return B
end
function _get_basis(M::AbstractManifold, p, B::ProjectedOrthonormalBasis{:svd,ℝ})
    S = representation_size(M)
    PS = prod(S)
    dim = manifold_dimension(M)
    # projection
    # TODO: find a better way to obtain a basis of the ambient space
    Xs = [
        convert(Vector, reshape(project(M, p, _euclidean_basis_vector(p, i)), PS)) for
        i in eachindex(p)
    ]
    O = reduce(hcat, Xs)
    # orthogonalization
    # TODO: try using rank-revealing QR here
    decomp = svd(O)
    rotated = Diagonal(decomp.S) * decomp.Vt
    vecs = [collect(reshape(rotated[i, :], S)) for i in 1:dim]
    # normalization
    for i in 1:dim
        i_norm = norm(M, p, vecs[i])
        vecs[i] /= i_norm
    end
    return CachedBasis(B, vecs)
end
function _get_basis(
    M::AbstractManifold,
    p,
    B::ProjectedOrthonormalBasis{:gram_schmidt,ℝ};
    kwargs...,
)
    E = [project(M, p, _euclidean_basis_vector(p, i)) for i in eachindex(p)]
    V = gram_schmidt(M, p, E; kwargs...)
    return CachedBasis(B, V)
end
function _get_basis(M::AbstractManifold, p, B::VeeOrthogonalBasis)
    return get_basis_vee(M, p, number_system(B))
end
function get_basis_vee(M::AbstractManifold, p, N)
    return get_basis(M, p, DefaultOrthogonalBasis(N))
end

function _get_basis(M::AbstractManifold, p, B::DefaultBasis)
    return get_basis_default(M, p, number_system(B))
end
function get_basis_default(M::AbstractManifold, p, N)
    return get_basis(M, p, DefaultOrthogonalBasis(N))
end

function _get_basis(M::AbstractManifold, p, B::DefaultOrthogonalBasis)
    return get_basis_orthogonal(M, p, number_system(B))
end
function get_basis_orthogonal(M::AbstractManifold, p, N)
    return get_basis(M, p, DefaultOrthonormalBasis(N))
end

function _get_basis(M::AbstractManifold, p, B::DiagonalizingOrthonormalBasis)
    return get_basis_diagonalizing(M, p, B)
end
function get_basis_diagonalizing end

function _get_basis(M::AbstractManifold, p, B::DefaultOrthonormalBasis)
    return get_basis_orthonormal(M, p, number_system(B))
end

function get_basis_orthonormal(M::AbstractManifold, p, N::AbstractNumbers; kwargs...)
    B = DefaultOrthonormalBasis(N)
    dim = number_of_coordinates(M, B)
    Eltp = coordinate_eltype(M, p, N)
    p0 = zero(Eltp)
    p1 = one(Eltp)
    return CachedBasis(
        B,
        [get_vector(M, p, [ifelse(i == j, p1, p0) for j in 1:dim], B) for i in 1:dim],
    )
end

@doc raw"""
    get_coordinates(M::AbstractManifold, p, X, B::AbstractBasis)
    get_coordinates(M::AbstractManifold, p, X, B::CachedBasis)

Compute a one-dimensional vector of coefficients of the tangent vector `X`
at point denoted by `p` on manifold `M` in basis `B`.

Depending on the basis, `p` may not directly represent a point on the manifold.
For example if a basis transported along a curve is used, `p` may be the coordinate
along the curve. If a [`CachedBasis`](@ref) is provided, their stored vectors are used,
otherwise the user has to provide a method to compute the coordinates.

For the [`CachedBasis`](@ref) keep in mind that the reconstruction with [`get_vector`](@ref)
requires either a dual basis or the cached basis to be selfdual, for example orthonormal

See also: [`get_vector`](@ref), [`get_basis`](@ref)
"""
function get_coordinates(M::AbstractManifold, p, X, B::AbstractBasis)
    return _get_coordinates(M, p, X, B)
end

function _get_coordinates(M::AbstractManifold, p, X, B::VeeOrthogonalBasis)
    return get_coordinates_vee(M, p, X, number_system(B))
end
function get_coordinates_vee(M::AbstractManifold, p, X, N)
    return get_coordinates(M, p, X, DefaultOrthogonalBasis(N))
end

function _get_coordinates(M::AbstractManifold, p, X, B::DefaultBasis)
    return get_coordinates_default(M, p, X, number_system(B))
end
function get_coordinates_default(M::AbstractManifold, p, X, N::AbstractNumbers)
    return get_coordinates(M, p, X, DefaultOrthogonalBasis(N))
end

function _get_coordinates(M::AbstractManifold, p, X, B::DefaultOrthogonalBasis)
    return get_coordinates_orthogonal(M, p, X, number_system(B))
end
function get_coordinates_orthogonal(M::AbstractManifold, p, X, N)
    return get_coordinates_orthonormal(M, p, X, N)
end

function _get_coordinates(M::AbstractManifold, p, X, B::DefaultOrthonormalBasis)
    return get_coordinates_orthonormal(M, p, X, number_system(B))
end
function get_coordinates_orthonormal(M::AbstractManifold, p, X, N)
    c = allocate_result(M, get_coordinates, p, X, DefaultOrthonormalBasis(N))
    return get_coordinates_orthonormal!(M, c, p, X, N)
end

function _get_coordinates(M::AbstractManifold, p, X, B::DiagonalizingOrthonormalBasis)
    return get_coordinates_diagonalizing(M, p, X, B)
end
function get_coordinates_diagonalizing(
    M::AbstractManifold,
    p,
    X,
    B::DiagonalizingOrthonormalBasis,
)
    c = allocate_result(M, get_coordinates, p, X, B)
    return get_coordinates_diagonalizing!(M, c, p, X, B)
end

function _get_coordinates(M::AbstractManifold, p, X, B::CachedBasis)
    return get_coordinates_cached(M, number_system(M), p, X, B, number_system(B))
end
function get_coordinates_cached(
    M::AbstractManifold,
    ::ComplexNumbers,
    p,
    X,
    B::CachedBasis,
    ::ComplexNumbers,
)
    return map(vb -> conj(inner(M, p, X, vb)), get_vectors(M, p, B))
end
function get_coordinates_cached(
    M::AbstractManifold,
    ::𝔽,
    p,
    X,
    C::CachedBasis,
    ::RealNumbers,
) where {𝔽}
    return map(vb -> real(inner(M, p, X, vb)), get_vectors(M, p, C))
end

function get_coordinates!(M::AbstractManifold, Y, p, X, B::AbstractBasis)
    return _get_coordinates!(M, Y, p, X, B)
end
function _get_coordinates!(M::AbstractManifold, Y, p, X, B::VeeOrthogonalBasis)
    return get_coordinates_vee!(M, Y, p, X, number_system(B))
end
function get_coordinates_vee!(M::AbstractManifold, Y, p, X, N)
    return get_coordinates!(M, Y, p, X, DefaultOrthogonalBasis(N))
end

function _get_coordinates!(M::AbstractManifold, Y, p, X, B::DefaultBasis)
    return get_coordinates_default!(M, Y, p, X, number_system(B))
end
function get_coordinates_default!(M::AbstractManifold, Y, p, X, N)
    return get_coordinates!(M, Y, p, X, DefaultOrthogonalBasis(N))
end

function _get_coordinates!(M::AbstractManifold, Y, p, X, B::DefaultOrthogonalBasis)
    return get_coordinates_orthogonal!(M, Y, p, X, number_system(B))
end
function get_coordinates_orthogonal!(M::AbstractManifold, Y, p, X, N)
    return get_coordinates!(M, Y, p, X, DefaultOrthonormalBasis(N))
end

function _get_coordinates!(M::AbstractManifold, Y, p, X, B::DefaultOrthonormalBasis)
    return get_coordinates_orthonormal!(M, Y, p, X, number_system(B))
end
function get_coordinates_orthonormal! end

function _get_coordinates!(M::AbstractManifold, Y, p, X, B::DiagonalizingOrthonormalBasis)
    return get_coordinates_diagonalizing!(M, Y, p, X, B)
end
function get_coordinates_diagonalizing! end

function _get_coordinates!(M::AbstractManifold, Y, p, X, B::CachedBasis)
    return get_coordinates_cached!(M, number_system(M), Y, p, X, B, number_system(B))
end
function get_coordinates_cached!(
    M::AbstractManifold,
    ::ComplexNumbers,
    Y,
    p,
    X,
    B::CachedBasis,
    ::ComplexNumbers,
)
    map!(vb -> conj(inner(M, p, X, vb)), Y, get_vectors(M, p, B))
    return Y
end
function get_coordinates_cached!(
    M::AbstractManifold,
    ::𝔽,
    Y,
    p,
    X,
    C::CachedBasis,
    ::RealNumbers,
) where {𝔽}
    map!(vb -> real(inner(M, p, X, vb)), Y, get_vectors(M, p, C))
    return Y
end

"""
    X = get_vector(M::AbstractManifold, p, c, B::AbstractBasis)

Convert a one-dimensional vector of coefficients in a basis `B` of
the tangent space at `p` on manifold `M` to a tangent vector `X` at `p`.

Depending on the basis, `p` may not directly represent a point on the manifold.
For example if a basis transported along a curve is used, `p` may be the coordinate
along the curve.

For the [`CachedBasis`](@ref) keep in mind that the reconstruction from [`get_coordinates`](@ref)
requires either a dual basis or the cached basis to be selfdual, for example orthonormal

See also: [`get_coordinates`](@ref), [`get_basis`](@ref)
"""
@inline function get_vector(M::AbstractManifold, p, c, B::AbstractBasis)
    return _get_vector(M, p, c, B)
end

@inline function _get_vector(M::AbstractManifold, p, c, B::VeeOrthogonalBasis)
    return get_vector_vee(M, p, c, number_system(B))
end
@inline function get_vector_vee(M::AbstractManifold, p, c, N)
    return get_vector(M, p, c, DefaultOrthogonalBasis(N))
end

@inline function _get_vector(M::AbstractManifold, p, c, B::DefaultBasis)
    return get_vector_default(M, p, c, number_system(B))
end
@inline function get_vector_default(M::AbstractManifold, p, c, N)
    return get_vector(M, p, c, DefaultOrthogonalBasis(N))
end

@inline function _get_vector(M::AbstractManifold, p, c, B::DefaultOrthogonalBasis)
    return get_vector_orthogonal(M, p, c, number_system(B))
end
@inline function get_vector_orthogonal(M::AbstractManifold, p, c, N)
    return get_vector_orthonormal(M, p, c, N)
end

function _get_vector(M::AbstractManifold, p, c, B::DefaultOrthonormalBasis)
    return get_vector_orthonormal(M, p, c, number_system(B))
end
function get_vector_orthonormal(M::AbstractManifold, p, c, N)
    B = DefaultOrthonormalBasis(N)
    Y = allocate_result(M, get_vector, p, c)
    return get_vector!(M, Y, p, c, B)
end

@inline function _get_vector(M::AbstractManifold, p, c, B::DiagonalizingOrthonormalBasis)
    return get_vector_diagonalizing(M, p, c, B)
end
function get_vector_diagonalizing(
    M::AbstractManifold,
    p,
    c,
    B::DiagonalizingOrthonormalBasis,
)
    Y = allocate_result(M, get_vector, p, c)
    return get_vector!(M, Y, p, c, B)
end

@inline function _get_vector(M::AbstractManifold, p, c, B::CachedBasis)
    return get_vector_cached(M, p, c, B)
end
_get_vector_cache_broadcast(::Any) = Val(true)
function get_vector_cached(M::AbstractManifold, p, X, B::CachedBasis)
    # quite convoluted but:
    #  1) preserves the correct `eltype`
    #  2) guarantees a reasonable array type `Y`
    #     (for example scalar * `SizedValidation` is an `SArray`)
    bvectors = get_vectors(M, p, B)
    if _get_vector_cache_broadcast(bvectors[1]) === Val(false)
        Xt = X[1] * bvectors[1]
        for i in 2:length(X)
            copyto!(Xt, Xt + X[i] * bvectors[i])
        end
    else
        Xt = X[1] .* bvectors[1]
        for i in 2:length(X)
            Xt .+= X[i] .* bvectors[i]
        end
    end
    return Xt
end
@inline function get_vector!(M::AbstractManifold, Y, p, c, B::AbstractBasis)
    return _get_vector!(M, Y, p, c, B)
end

@inline function _get_vector!(M::AbstractManifold, Y, p, c, B::VeeOrthogonalBasis)
    return get_vector_vee!(M, Y, p, c, number_system(B))
end
@inline get_vector_vee!(M, Y, p, c, N) = get_vector!(M, Y, p, c, DefaultOrthogonalBasis(N))

@inline function _get_vector!(M::AbstractManifold, Y, p, c, B::DefaultBasis)
    return get_vector_default!(M, Y, p, c, number_system(B))
end
@inline function get_vector_default!(M::AbstractManifold, Y, p, c, N)
    return get_vector!(M, Y, p, c, DefaultOrthogonalBasis(N))
end

@inline function _get_vector!(M::AbstractManifold, Y, p, c, B::DefaultOrthogonalBasis)
    return get_vector_orthogonal!(M, Y, p, c, number_system(B))
end
@inline function get_vector_orthogonal!(M::AbstractManifold, Y, p, c, N)
    return get_vector!(M, Y, p, c, DefaultOrthonormalBasis(N))
end

@inline function _get_vector!(M::AbstractManifold, Y, p, c, B::DefaultOrthonormalBasis)
    return get_vector_orthonormal!(M, Y, p, c, number_system(B))
end
function get_vector_orthonormal! end

@inline function _get_vector!(
    M::AbstractManifold,
    Y,
    p,
    c,
    B::DiagonalizingOrthonormalBasis,
)
    return get_vector_diagonalizing!(M, Y, p, c, B)
end
function get_vector_diagonalizing! end

@inline function _get_vector!(M::AbstractManifold, Y, p, c, B::CachedBasis)
    return get_vector_cached!(M, Y, p, c, B)
end
function get_vector_cached!(M::AbstractManifold, Y, p, X, B::CachedBasis)
    # quite convoluted but:
    #  1) preserves the correct `eltype`
    #  2) guarantees a reasonable array type `Y`
    #     (for example scalar * `SizedValidation` is an `SArray`)
    bvectors = get_vectors(M, p, B)
    return if _get_vector_cache_broadcast(bvectors[1]) === Val(false)
        Xt = X[1] * bvectors[1]
        copyto!(Y, Xt)
        for i in 2:length(X)
            copyto!(Y, Y + X[i] * bvectors[i])
        end
        return Y
    else
        Xt = X[1] .* bvectors[1]
        copyto!(Y, Xt)
        for i in 2:length(X)
            Y .+= X[i] .* bvectors[i]
        end
        return Y
    end
end

"""
    get_vectors(M::AbstractManifold, p, B::AbstractBasis)

Get the basis vectors of basis `B` of the tangent space at point `p`.
"""
function get_vectors(M::AbstractManifold, p, B::AbstractBasis)
    return _get_vectors(M, p, B)
end
function _get_vectors(::AbstractManifold, ::Any, B::CachedBasis)
    return _get_vectors(B)
end
_get_vectors(B::CachedBasis{𝔽,<:AbstractBasis,<:AbstractArray}) where {𝔽} = B.data
function _get_vectors(B::CachedBasis{𝔽,<:AbstractBasis,<:DiagonalizingBasisData}) where {𝔽}
    return B.data.vectors
end

@doc raw"""
    gram_schmidt(M::AbstractManifold{𝔽}, p, B::AbstractBasis{𝔽}) where {𝔽}
    gram_schmidt(M::AbstractManifold, p, V::AbstractVector)

Compute an ONB in the tangent space at `p` on the [`AbstractManifold`](@ref} `M` from either an
[`AbstractBasis`](@ref) basis ´B´ or a set of (at most) [`manifold_dimension`](@ref)`(M)`
many vectors.
Note that this method requires the manifold and basis to work on the same
[`AbstractNumbers`](@ref) `𝔽`, i.e. with real coefficients.

The method always returns a basis, i.e. linearly dependent vectors are removed.

# Keyword arguments

* `warn_linearly_dependent` (`false`) – warn if the basis vectors are not linearly
  independent
* `skip_linearly_dependent` (`false`) – whether to just skip (`true`) a vector that
  is linearly dependent to the previous ones or to stop (`false`, default) at that point
* `return_incomplete_set` (`false`) – throw an error if the resulting set of vectors is not
  a basis but contains less vectors

further keyword arguments can be passed to set the accuracy of the independence test.
Especially `atol` is raised slightly by default to `atol = 5*1e-16`.

# Return value

When a set of vectors is orthonormalized a set of vectors is returned.
When an [`AbstractBasis`](@ref) is orthonormalized, a [`CachedBasis`](@ref) is returned.
"""
function gram_schmidt(
    M::AbstractManifold{𝔽},
    p,
    B::AbstractBasis{𝔽};
    warn_linearly_dependent = false,
    return_incomplete_set = false,
    skip_linearly_dependent = false,
    kwargs...,
) where {𝔽}
    V = gram_schmidt(
        M,
        p,
        get_vectors(M, p, B);
        warn_linearly_dependent = warn_linearly_dependent,
        skip_linearly_dependent = skip_linearly_dependent,
        return_incomplete_set = return_incomplete_set,
        kwargs...,
    )
    return CachedBasis(GramSchmidtOrthonormalBasis(𝔽), V)
end
function gram_schmidt(
    M::AbstractManifold,
    p,
    V::AbstractVector;
    atol = eps(number_eltype(first(V))),
    warn_linearly_dependent = false,
    return_incomplete_set = false,
    skip_linearly_dependent = false,
    kwargs...,
)
    N = length(V)
    Ξ = empty(V)
    dim = manifold_dimension(M)
    N < dim && @warn "Input only has $(N) vectors, but manifold dimension is $(dim)."
    @inbounds for n in 1:N
        Ξₙ = copy(M, p, V[n])
        for k in 1:length(Ξ)
            Ξₙ .-= real(inner(M, p, Ξ[k], Ξₙ)) .* Ξ[k]
        end
        nrmΞₙ = norm(M, p, Ξₙ)
        if isapprox(nrmΞₙ / dim, 0; atol = atol, kwargs...)
            warn_linearly_dependent &&
                @warn "Input vector $(n) lies in the span of the previous ones."
            !skip_linearly_dependent && throw(
                ErrorException("Input vector $(n) lies in the span of the previous ones."),
            )
        else
            push!(Ξ, Ξₙ ./ nrmΞₙ)
        end
        if length(Ξ) == dim
            (n < N) &&
                @warn "More vectors ($(N)) entered than the dimension of the manifold ($dim). All vectors after the $n th ignored."
            return Ξ
        end
    end
    if return_incomplete_set # if we rech this point - length(Ξ) < dim
        return Ξ
    else
        throw(
            ErrorException(
                "gram_schmidt found only $(length(Ξ)) orthonormal basis vectors, but manifold dimension is $(dim).",
            ),
        )
    end
end


@doc raw"""
    hat(M::AbstractManifold, p, Xⁱ)

Given a basis ``e_i`` on the tangent space at a point `p` and tangent
component vector ``X^i ∈ ℝ``, compute the equivalent vector representation
``X=X^i e_i``, where Einstein summation notation is used:

````math
∧ : X^i ↦ X^i e_i
````

For array manifolds, this converts a vector representation of the tangent
vector to an array representation. The [`vee`](@ref) map is the `hat` map's
inverse.
"""
@inline hat(M::AbstractManifold, p, X) = get_vector(M, p, X, VeeOrthogonalBasis(ℝ))
@inline hat!(M::AbstractManifold, Y, p, X) = get_vector!(M, Y, p, X, VeeOrthogonalBasis(ℝ))

"""
    number_of_coordinates(M::AbstractManifold, B::AbstractBasis)
    number_of_coordinates(M::AbstractManifold, ::𝔾)

Compute the number of coordinates in basis of field type `𝔾` on a manifold `M`.
This also corresponds to the number of vectors represented by `B`,
or stored within `B` in case of a [`CachedBasis`](@ref).
"""
function number_of_coordinates(M::AbstractManifold, ::AbstractBasis{𝔾}) where {𝔾}
    return number_of_coordinates(M, 𝔾)
end
function number_of_coordinates(M::AbstractManifold, f::𝔾) where {𝔾}
    return div(manifold_dimension(M), real_dimension(f))
end

"""
    number_system(::AbstractBasis)

The number system for the vectors of the given basis.
"""
number_system(::AbstractBasis{𝔽}) where {𝔽} = 𝔽

"""
    requires_caching(B::AbstractBasis)

Return whether basis `B` can be used in [`get_vector`](@ref) and [`get_coordinates`](@ref)
without calling [`get_basis`](@ref) first.
"""
requires_caching(::AbstractBasis) = true
requires_caching(::CachedBasis) = false
requires_caching(::DefaultBasis) = false
requires_caching(::DefaultOrthogonalBasis) = false
requires_caching(::DefaultOrthonormalBasis) = false

function _show_basis_vector(io::IO, X; pre = "", head = "")
    sX = sprint(show, "text/plain", X, context = io, sizehint = 0)
    sX = replace(sX, '\n' => "\n$(pre)")
    return print(io, head, pre, sX)
end
function _show_basis_vector_range(io::IO, Ξ, range; pre = "", sym = "E")
    for i in range
        _show_basis_vector(io, Ξ[i]; pre = pre, head = "\n$(sym)$(i) =\n")
    end
    return nothing
end
function _show_basis_vector_range_noheader(io::IO, Ξ; max_vectors = 4, pre = "", sym = "E")
    nv = length(Ξ)
    return if nv ≤ max_vectors
        _show_basis_vector_range(io, Ξ, 1:nv; pre = "  ", sym = " E")
    else
        halfn = div(max_vectors, 2)
        _show_basis_vector_range(io, Ξ, 1:halfn; pre = "  ", sym = " E")
        print(io, "\n ⋮")
        _show_basis_vector_range(io, Ξ, (nv - halfn + 1):nv; pre = "  ", sym = " E")
    end
end

function show(io::IO, ::DefaultBasis{𝔽}) where {𝔽}
    return print(io, "DefaultBasis($(𝔽))")
end
function show(io::IO, ::DefaultOrthogonalBasis{𝔽}) where {𝔽}
    return print(io, "DefaultOrthogonalBasis($(𝔽))")
end
function show(io::IO, ::DefaultOrthonormalBasis{𝔽}) where {𝔽}
    return print(io, "DefaultOrthonormalBasis($(𝔽))")
end
function show(io::IO, ::GramSchmidtOrthonormalBasis{𝔽}) where {𝔽}
    return print(io, "GramSchmidtOrthonormalBasis($(𝔽))")
end
function show(io::IO, ::ProjectedOrthonormalBasis{method,𝔽}) where {method,𝔽}
    return print(io, "ProjectedOrthonormalBasis($(repr(method)), $(𝔽))")
end
function show(io::IO, ::MIME"text/plain", onb::DiagonalizingOrthonormalBasis)
    println(
        io,
        "DiagonalizingOrthonormalBasis($(number_system(onb))) with eigenvalue 0 in direction:",
    )
    sk = sprint(show, "text/plain", onb.frame_direction, context = io, sizehint = 0)
    sk = replace(sk, '\n' => "\n ")
    return print(io, sk)
end
function show(
    io::IO,
    ::MIME"text/plain",
    B::CachedBasis{𝔽,T,D},
) where {𝔽,T<:AbstractBasis,D}
    try
        vectors = _get_vectors(B)
        print(
            io,
            "Cached basis of type $T with $(length(vectors)) basis vector$(length(vectors) == 1 ? "" : "s"):",
        )
        return _show_basis_vector_range_noheader(
            io,
            vectors;
            max_vectors = 4,
            pre = "  ",
            sym = " E",
        )
    catch e
        # in case _get_vectors(B) is not defined
        print(io, "Cached basis of type $T")
    end
end
function show(
    io::IO,
    ::MIME"text/plain",
    B::CachedBasis{𝔽,T,D},
) where {𝔽,T<:DiagonalizingOrthonormalBasis,D<:DiagonalizingBasisData}
    vectors = _get_vectors(B)
    nv = length(vectors)
    sk = sprint(show, "text/plain", T(B.data.frame_direction), context = io, sizehint = 0)
    sk = replace(sk, '\n' => "\n ")
    print(io, sk)
    println(io, "\nand $(nv) basis vector$(nv == 1 ? "" : "s").")
    print(io, "Basis vectors:")
    _show_basis_vector_range_noheader(io, vectors; max_vectors = 4, pre = "  ", sym = " E")
    println(io, "\nEigenvalues:")
    sk = sprint(show, "text/plain", B.data.eigenvalues, context = io, sizehint = 0)
    sk = replace(sk, '\n' => "\n ")
    return print(io, ' ', sk)
end

@doc raw"""
    vee(M::AbstractManifold, p, X)

Given a basis ``e_i`` on the tangent space at a point `p` and tangent
vector `X`, compute the vector components ``X^i ∈ ℝ``, such that ``X = X^i e_i``, where
Einstein summation notation is used:

````math
\vee : X^i e_i ↦ X^i
````

For array manifolds, this converts an array representation of the tangent
vector to a vector representation. The [`hat`](@ref) map is the `vee` map's
inverse.
"""
vee(M::AbstractManifold, p, X) = get_coordinates(M, p, X, VeeOrthogonalBasis(ℝ))
function vee!(M::AbstractManifold, Y, p, X)
    return get_coordinates!(M, Y, p, X, VeeOrthogonalBasis(ℝ))
end
