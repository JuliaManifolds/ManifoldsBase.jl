"""
    VectorSpaceType

Abstract type for tangent spaces, cotangent spaces, their tensor products,
exterior products, etc.

Every vector space `fiber` is supposed to provide:
* a method of constructing vectors,
* basic operations: addition, subtraction, multiplication by a scalar
  and negation (unary minus),
* [`zero_vector(fiber, p)`](@ref Main.Manifolds.zero_vector) to construct zero vectors at point `p`,
* `allocate(X)` and `allocate(X, T)` for vector `X` and type `T`,
* `copyto!(X, Y)` for vectors `X` and `Y`,
* `number_eltype(v)` for vector `v`,
* [`vector_space_dimension`](@ref).

Optionally:
* inner product via `inner` (used to provide Riemannian metric on vector
  bundles),
* [`flat`](@ref Main.Manifolds.flat) and [`sharp`](@ref Main.Manifolds.sharp),
* `norm` (by default uses `inner`),
* [`project`](@ref) (for embedded vector spaces),
* [`representation_size`](@ref) (if support for [`ProductArray`](@ref Main.Manifolds.ProductArray) is desired),
* broadcasting for basic operations.
"""
abstract type VectorSpaceType end

struct TangentSpaceType <: VectorSpaceType end

struct CotangentSpaceType <: VectorSpaceType end

TCoTSpaceType = Union{TangentSpaceType,CotangentSpaceType}

const TangentSpace = TangentSpaceType()
const CotangentSpace = CotangentSpaceType()

"""
    AbstractBasis{𝔽,VST<:VectorSpaceType}

Abstract type that represents a basis of vector space of type `VST` on a manifold or
a subset of it.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

# See also

[`VectorSpaceType`](@ref)
"""
abstract type AbstractBasis{𝔽,VST<:VectorSpaceType} end

"""
    DefaultBasis{𝔽,VST<:VectorSpaceType}

An arbitrary basis of vector space of type `VST` on a manifold. This will usually
be the fastest basis available for a manifold.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

# See also

[`VectorSpaceType`](@ref)
"""
struct DefaultBasis{𝔽,VST<:VectorSpaceType} <: AbstractBasis{𝔽,VST}
    vector_space::VST
end
function DefaultBasis(𝔽::AbstractNumbers = ℝ, vs::VectorSpaceType = TangentSpace)
    return DefaultBasis{𝔽,typeof(vs)}(vs)
end
function DefaultBasis{𝔽}(vs::VectorSpaceType = TangentSpace) where {𝔽}
    return DefaultBasis{𝔽,typeof(vs)}(vs)
end
function DefaultBasis{𝔽,TangentSpaceType}() where {𝔽}
    return DefaultBasis{𝔽,TangentSpaceType}(TangentSpace)
end

"""
    AbstractOrthogonalBasis{𝔽,VST<:VectorSpaceType}

Abstract type that represents an orthonormal basis of vector space of type `VST` on a
manifold or a subset of it.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

# See also

[`VectorSpaceType`](@ref)
"""
abstract type AbstractOrthogonalBasis{𝔽,VST<:VectorSpaceType} <: AbstractBasis{𝔽,VST} end

"""
    DefaultOrthogonalBasis{𝔽,VST<:VectorSpaceType}

An arbitrary orthogonal basis of vector space of type `VST` on a manifold. This will usually
be the fastest orthogonal basis available for a manifold.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

# See also

[`VectorSpaceType`](@ref)
"""
struct DefaultOrthogonalBasis{𝔽,VST<:VectorSpaceType} <: AbstractOrthogonalBasis{𝔽,VST}
    vector_space::VST
end
function DefaultOrthogonalBasis(𝔽::AbstractNumbers = ℝ, vs::VectorSpaceType = TangentSpace)
    return DefaultOrthogonalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthogonalBasis{𝔽}(vs::VectorSpaceType = TangentSpace) where {𝔽}
    return DefaultOrthogonalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthogonalBasis{𝔽,TangentSpaceType}() where {𝔽}
    return DefaultOrthogonalBasis{𝔽,TangentSpaceType}(TangentSpace)
end


struct VeeOrthogonalBasis{𝔽} <: AbstractOrthogonalBasis{𝔽,TangentSpaceType} end
VeeOrthogonalBasis(𝔽::AbstractNumbers = ℝ) = VeeOrthogonalBasis{𝔽}()

"""
    AbstractOrthonormalBasis{𝔽,VST<:VectorSpaceType}

Abstract type that represents an orthonormal basis of vector space of type `VST` on a
manifold or a subset of it.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

# See also

[`VectorSpaceType`](@ref)
"""
abstract type AbstractOrthonormalBasis{𝔽,VST<:VectorSpaceType} <:
              AbstractOrthogonalBasis{𝔽,VST} end

"""
    DefaultOrthonormalBasis(𝔽::AbstractNumbers = ℝ, vs::VectorSpaceType = TangentSpace)

An arbitrary orthonormal basis of vector space of type `VST` on a manifold. This will usually
be the fastest orthonormal basis available for a manifold.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

# See also

[`VectorSpaceType`](@ref)
"""
struct DefaultOrthonormalBasis{𝔽,VST<:VectorSpaceType} <: AbstractOrthonormalBasis{𝔽,VST}
    vector_space::VST
end

function DefaultOrthonormalBasis(𝔽::AbstractNumbers = ℝ, vs::VectorSpaceType = TangentSpace)
    return DefaultOrthonormalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthonormalBasis{𝔽}(vs::VectorSpaceType = TangentSpace) where {𝔽}
    return DefaultOrthonormalBasis{𝔽,typeof(vs)}(vs)
end
function DefaultOrthonormalBasis{𝔽,TangentSpaceType}() where {𝔽}
    return DefaultOrthonormalBasis{𝔽,TangentSpaceType}(TangentSpace)
end

"""
    ProjectedOrthonormalBasis(method::Symbol, 𝔽::AbstractNumbers = ℝ)

An orthonormal basis that comes from orthonormalization of basis vectors
of the ambient space projected onto the subspace representing the tangent space
at a given point.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

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
tensor $R(u,v)w$ and where the direction `frame_direction` $v$ has curvature `0`.

The type parameter `𝔽` denotes the [`AbstractNumbers`](@ref) that will be used
for the vectors elements.

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
function CachedBasis(::B, data::V) where {V,𝔽,B<:AbstractBasis{𝔽,TangentSpaceType}}
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
const DISAMBIGUATION_BASIS_TYPES = [
    CachedBasis,
    DefaultBasis,
    DefaultBasis{<:Any,TangentSpaceType},
    DefaultOrthonormalBasis,
    DefaultOrthonormalBasis{<:Any,TangentSpaceType},
    DefaultOrthogonalBasis,
    DefaultOrthogonalBasis{<:Any,TangentSpaceType},
    DiagonalizingOrthonormalBasis,
    ProjectedOrthonormalBasis{:svd,ℝ},
    ProjectedOrthonormalBasis{:gram_schmidt,ℝ},
    VeeOrthogonalBasis,
]
const DISAMBIGUATION_COTANGENT_BASIS_TYPES = [
    DefaultBasis{<:Any,CotangentSpaceType},
    DefaultOrthonormalBasis{<:Any,CotangentSpaceType},
    DefaultOrthogonalBasis{<:Any,CotangentSpaceType},
]

@inline function allocate_result_type(
    M::AbstractManifold,
    f::Union{typeof(get_coordinates),typeof(get_vector)},
    args::Tuple,
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

function combine_allocation_promotion_functions(f::T, ::T) where {T}
    return f
end
function combine_allocation_promotion_functions(::typeof(complex), ::typeof(identity))
    return complex
end
function combine_allocation_promotion_functions(::typeof(identity), ::typeof(complex))
    return complex
end

@doc raw"""
    dual_basis(M::AbstractManifold, p, B::AbstractBasis)

Get the dual basis to `B`, a basis of a vector space at point `p` from manifold `M`.

The dual to the $i$th vector $v_i$ from basis `B` is a vector $v^i$ from the dual space
such that $v^i(v_j) = δ^i_j$, where $δ^i_j$ is the Kronecker delta symbol:
````math
δ^i_j = \begin{cases}
1 & \text{ if } i=j, \\
0 & \text{ otherwise.}
\end{cases}
````
"""
dual_basis(M::AbstractManifold, p, B::AbstractBasis)

function dual_basis(
    ::AbstractManifold,
    p,
    ::DefaultOrthonormalBasis{𝔽,TangentSpaceType},
) where {𝔽}
    return DefaultOrthonormalBasis{𝔽}(CotangentSpace)
end
function dual_basis(
    ::AbstractManifold,
    p,
    ::DefaultOrthonormalBasis{𝔽,CotangentSpaceType},
) where {𝔽}
    return DefaultOrthonormalBasis{𝔽}(TangentSpace)
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
    get_basis(M::AbstractManifold, p, B::AbstractBasis) -> CachedBasis

Compute the basis vectors of the tangent space at a point on manifold `M`
represented by `p`.

Returned object derives from [`AbstractBasis`](@ref) and may have a field `.vectors`
that stores tangent vectors or it may store them implicitly, in which case
the function [`get_vectors`](@ref) needs to be used to retrieve the basis vectors.

See also: [`get_coordinates`](@ref), [`get_vector`](@ref)
"""
get_basis(M::AbstractManifold, p, B::AbstractBasis)
@decorator_transparent_signature get_basis(
    M::AbstractDecoratorManifold,
    p,
    B::AbstractBasis,
)
function decorator_transparent_dispatch(::typeof(get_basis), ::AbstractManifold, args...)
    return Val(:parent)
end

function get_basis(
    M::AbstractManifold,
    p,
    B::DefaultOrthonormalBasis{<:Any,TangentSpaceType},
)
    dim = manifold_dimension(M)
    return CachedBasis(
        B,
        [get_vector(M, p, [ifelse(i == j, 1, 0) for j in 1:dim], B) for i in 1:dim],
    )
end
function get_basis(::AbstractManifold, ::Any, B::CachedBasis)
    return B
end
function get_basis(M::AbstractManifold, p, B::ProjectedOrthonormalBasis{:svd,ℝ})
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
function get_basis(
    M::AbstractManifold,
    p,
    B::ProjectedOrthonormalBasis{:gram_schmidt,ℝ};
    warn_linearly_dependent = false,
    return_incomplete_set = false,
    kwargs...,
)
    E = [project(M, p, _euclidean_basis_vector(p, i)) for i in eachindex(p)]
    V = gram_schmidt(
        M,
        p,
        E;
        warn_linearly_dependent = warn_linearly_dependent,
        return_incomplete_set = return_incomplete_set,
        kwargs...,
    )
    return CachedBasis(B, V)
end
for BT in DISAMBIGUATION_BASIS_TYPES
    eval(
        quote
            @decorator_transparent_signature get_basis(
                M::AbstractDecoratorManifold,
                p,
                B::$BT,
            )
        end,
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
    Y = allocate_result_coords_like(M, get_coordinates, p, X)
    return get_coordinates!(M, Y, p, X, B)
end
@decorator_transparent_signature get_coordinates(
    M::AbstractDecoratorManifold,
    p,
    X,
    B::AbstractBasis,
)
function decorator_transparent_dispatch(
    ::typeof(get_coordinates),
    ::AbstractManifold,
    args...,
)
    return Val(:parent)
end

function get_coordinates!(M::AbstractManifold, Y, p, X, B::AbstractBasis)
    return error(
        "get_coordinates! not implemented for manifold of type $(typeof(M)) coordinates of type $(typeof(Y)), a point of type $(typeof(p)), tangent vector of type $(typeof(X)) and basis of type $(typeof(B)).",
    )
end
@decorator_transparent_signature get_coordinates!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    B::AbstractBasis,
)
for BT in [DISAMBIGUATION_BASIS_TYPES..., DISAMBIGUATION_COTANGENT_BASIS_TYPES...]
    eval(
        quote
            @decorator_transparent_signature get_coordinates!(
                M::AbstractDecoratorManifold,
                Y,
                p,
                X,
                B::$BT,
            )
        end,
    )
end

function get_coordinates!(M::AbstractManifold, Y, p, X, B::VeeOrthogonalBasis)
    return get_coordinates!(M, Y, p, X, DefaultOrthogonalBasis(number_system(B)))
end
function get_coordinates!(M::AbstractManifold, Y, p, X, B::DefaultBasis)
    return get_coordinates!(M, Y, p, X, DefaultOrthogonalBasis(number_system(B)))
end
function get_coordinates!(M::AbstractManifold, Y, p, X, B::DefaultOrthogonalBasis)
    return get_coordinates!(M, Y, p, X, DefaultOrthonormalBasis(number_system(B)))
end
function get_coordinates!(M::AbstractManifold, Y, p, X, B::CachedBasis)
    return _get_coordinates!(M, number_system(M), Y, p, X, B, number_system(B))
end
function _get_coordinates!(
    M::AbstractManifold,
    ::ComplexNumbers,
    Y,
    p,
    X,
    B::CachedBasis,
    ::RealNumbers,
)
    map!(vb -> conj(inner(M, p, X, vb)), Y, get_vectors(M, p, B))
    return Y
end
function _get_coordinates!(
    M::AbstractManifold,
    a::𝔽,
    Y,
    p,
    X,
    C::CachedBasis,
    b::𝔽,
) where {𝔽}
    map!(vb -> real(inner(M, p, X, vb)), Y, get_vectors(M, p, C))
    return Y
end

"""
    get_vector(M::AbstractManifold, p, X, B::AbstractBasis)

Convert a one-dimensional vector of coefficients in a basis `B` of
the tangent space at `p` on manifold `M` to a tangent vector `X` at `p`.

Depending on the basis, `p` may not directly represent a point on the manifold.
For example if a basis transported along a curve is used, `p` may be the coordinate
along the curve.

For the [`CachedBasis`](@ref) keep in mind that the reconstruction from [`get_coordinates`](@ref)
requires either a dual basis or the cached basis to be selfdual, for example orthonormal

See also: [`get_coordinates`](@ref), [`get_basis`](@ref)
"""
function get_vector(M::AbstractManifold, p, X, B::AbstractBasis)
    Y = allocate_result_vector(M, get_vector, p)
    return get_vector!(M, Y, p, X, B)
end
@decorator_transparent_signature get_vector(
    M::AbstractDecoratorManifold,
    p,
    X,
    B::AbstractBasis,
)
function decorator_transparent_dispatch(::typeof(get_vector), ::AbstractManifold, args...)
    return Val(:parent)
end

function get_vector!(M::AbstractManifold, Y, p, X, B::AbstractBasis)
    return error(
        "get_vector! not implemented for manifold of type $(typeof(M)) vector of type $(typeof(Y)), a point of type $(typeof(p)), coordinates of type $(typeof(X)) and basis of type $(typeof(B)).",
    )
end
@decorator_transparent_signature get_vector!(
    M::AbstractDecoratorManifold,
    Y,
    p,
    X,
    B::AbstractBasis,
)
for BT in [DISAMBIGUATION_BASIS_TYPES..., DISAMBIGUATION_COTANGENT_BASIS_TYPES...]
    eval(
        quote
            @decorator_transparent_signature get_vector!(
                M::AbstractDecoratorManifold,
                Y,
                p,
                X,
                B::$BT,
            )
        end,
    )
end

_get_vector_cache_broadcast(::Any) = Val(true)

function get_vector!(M::AbstractManifold, Y, p, X, B::VeeOrthogonalBasis)
    return get_vector!(M, Y, p, X, DefaultOrthogonalBasis(number_system(B)))
end
function get_vector!(M::AbstractManifold, Y, p, X, B::DefaultBasis)
    return get_vector!(M, Y, p, X, DefaultOrthogonalBasis(number_system(B)))
end
function get_vector!(M::AbstractManifold, Y, p, X, B::DefaultOrthogonalBasis)
    return get_vector!(M, Y, p, X, DefaultOrthonormalBasis(number_system(B)))
end
function get_vector!(M::AbstractManifold, Y, p, X, B::CachedBasis)
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
    return error(
        "get_vectors not implemented for manifold of type $(typeof(M)) a point of type $(typeof(p)) and basis of type $(typeof(B)).",
    )
end
function get_vectors(::AbstractManifold, ::Any, B::CachedBasis)
    return _get_vectors(B)
end
#internal for directly cached basis i.e. those that are just arrays – used in show
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
* `return_incomplete_set` (`false`) – throw an error if the resulting set of vectors is not
  a basis but contains less vectors

further keyword arguments can be passed to set the accuracy of the independence test.

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
    kwargs...,
) where {𝔽}
    V = gram_schmidt(
        M,
        p,
        get_vectors(M, p, B);
        warn_linearly_dependent = warn_linearly_dependent,
        return_incomplete_set = return_incomplete_set,
        kwargs...,
    )
    return CachedBasis(GramSchmidtOrthonormalBasis(𝔽), V)
end
function gram_schmidt(
    M::AbstractManifold,
    p,
    V::AbstractVector;
    warn_linearly_dependent = false,
    return_incomplete_set = false,
    kwargs...,
)
    N = length(V)
    Ξ = empty(V)
    dim = manifold_dimension(M)
    N < dim && @warn "Input only has $(N) vectors, but manifold dimension is $(dim)."
    @inbounds for n in 1:N
        Ξₙ = V[n]
        for k in 1:length(Ξ)
            Ξₙ .-= real(inner(M, p, Ξ[k], Ξₙ)) .* Ξ[k]
        end
        nrmΞₙ = norm(M, p, Ξₙ)
        if nrmΞₙ == 0
            warn_linearly_dependent && @warn "Input vector $(n) has length 0."
            @goto skip
        end
        Ξₙ ./= nrmΞₙ
        for k in 1:length(Ξ)
            if !isapprox(real(inner(M, p, Ξ[k], Ξₙ)), 0; kwargs...)
                warn_linearly_dependent &&
                    @warn "Input vector $(n) is not linearly independent of output basis vector $(k)."
                @goto skip
            end
        end
        push!(Ξ, Ξₙ)
        length(Ξ) == dim && return Ξ
        @label skip
    end
    return if return_incomplete_set
        return Ξ
    else
        error(
            "gram_schmidt found only $(length(Ξ)) orthonormal basis vectors, but manifold dimension is $(dim).",
        )
    end
end

@doc raw"""
    hat(M::AbstractManifold, p, Xⁱ)

Given a basis $e_i$ on the tangent space at a point `p` and tangent
component vector $X^i$, compute the equivalent vector representation
$X=X^i e_i$, where Einstein summation notation is used:

````math
∧ : X^i ↦ X^i e_i
````

For array manifolds, this converts a vector representation of the tangent
vector to an array representation. The [`vee`](@ref) map is the `hat` map's
inverse.
"""
hat(M::AbstractManifold, p, X) = get_vector(M, p, X, VeeOrthogonalBasis())
hat!(M::AbstractManifold, Y, p, X) = get_vector!(M, Y, p, X, VeeOrthogonalBasis())

"""
    number_of_coordinates(M::AbstractManifold, B::AbstractBasis)

Compute the number of coordinates in basis `B` of manifold `M`.
This also corresponds to the number of vectors represented by `B`,
or stored within `B` in case of a [`CachedBasis`](@ref).
"""
function number_of_coordinates(M::AbstractManifold{𝔽}, B::AbstractBasis{𝔾}) where {𝔽,𝔾}
    return div(manifold_dimension(M), real_dimension(𝔽)) * real_dimension(𝔾)
end
function number_of_coordinates(M::AbstractManifold{𝔽}, B::AbstractBasis{𝔽}) where {𝔽}
    return manifold_dimension(M)
end

"""
    number_system(::AbstractBasis)

The number system for the vectors of the given basis.
"""
number_system(::AbstractBasis{𝔽}) where {𝔽} = 𝔽

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
function show(io::IO, mime::MIME"text/plain", onb::DiagonalizingOrthonormalBasis)
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
    mime::MIME"text/plain",
    B::CachedBasis{𝔽,T,D},
) where {𝔽,T<:AbstractBasis,D}
    print(
        io,
        "$(T()) with $(length(_get_vectors(B))) basis vector$(length(_get_vectors(B)) == 1 ? "" : "s"):",
    )
    return _show_basis_vector_range_noheader(
        io,
        _get_vectors(B);
        max_vectors = 4,
        pre = "  ",
        sym = " E",
    )
end
function show(
    io::IO,
    mime::MIME"text/plain",
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

Given a basis $e_i$ on the tangent space at a point `p` and tangent
vector `X`, compute the vector components $X^i$, such that $X = X^i e_i$, where
Einstein summation notation is used:

````math
\vee : X^i e_i ↦ X^i
````

For array manifolds, this converts an array representation of the tangent
vector to a vector representation. The [`hat`](@ref) map is the `vee` map's
inverse.
"""
vee(M::AbstractManifold, p, X) = get_coordinates(M, p, X, VeeOrthogonalBasis())
vee!(M::AbstractManifold, Y, p, X) = get_coordinates!(M, Y, p, X, VeeOrthogonalBasis())

macro invoke_maker(argnum, type, sig)
    parts = ManifoldsBase._split_signature(sig)
    kwargs_list = parts[:kwargs_list]
    callargs = parts[:callargs]
    fname = parts[:fname]
    where_exprs = parts[:where_exprs]
    argnames = parts[:argnames]
    argtypes = parts[:argtypes]
    kwargs_call = parts[:kwargs_call]

    return esc(
        quote
            function ($fname)($(callargs...); $(kwargs_list...)) where {$(where_exprs...)}
                return invoke(
                    $fname,
                    Tuple{
                        $(argtypes[1:(argnum - 1)]...),
                        $type,
                        $(argtypes[(argnum + 1):end]...),
                    },
                    $(argnames...);
                    $(kwargs_call...),
                )
            end
        end,
    )
end
