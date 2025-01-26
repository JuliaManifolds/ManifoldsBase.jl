"""
    ManifoldsBaseTestUtils

A small module to collect common definitions and functions used in (several) tests of
`ManifoldsBase.jl`.

* A `TestSphere`
*
"""
module ManifoldsBaseTestUtils

using ManifoldsBase, LinearAlgebra, Random
using ManifoldsBase: ℝ, ℂ, AbstractManifold, DefaultManifold, EuclideanMetric
using ManifoldsBase: AbstractNumbers, RealNumbers, ComplexNumbers
using ManifoldsBase:
    @manifold_element_forwards, @manifold_vector_forwards, @default_manifold_fallbacks
import Base: +, *, -
#
#
# minimal implementation of the sphere – to test a few more involved Riemannian functions
struct TestSphere{N,𝔽} <: AbstractManifold{𝔽} end
TestSphere(N::Int, 𝔽 = ℝ) = TestSphere{N,𝔽}()

function ManifoldsBase.change_metric!(
    M::TestSphere,
    Y,
    ::ManifoldsBase.EuclideanMetric,
    p,
    X,
)
    return copyto!(M, Y, p, X)
end
function ManifoldsBase.change_representer!(
    M::TestSphere,
    Y,
    ::ManifoldsBase.EuclideanMetric,
    p,
    X,
)
    return copyto!(M, Y, p, X)
end
function ManifoldsBase.check_point(M::TestSphere, p; kwargs...)
    if !isapprox(norm(p), 1.0; kwargs...)
        return DomainError(
            norm(p),
            "The point $(p) does not lie on the $(M) since its norm is not 1.",
        )
    end
    return nothing
end
function ManifoldsBase.check_vector(M::TestSphere, p, X; kwargs...)
    if !isapprox(abs(real(dot(p, X))), 0.0; kwargs...)
        return DomainError(
            abs(dot(p, X)),
            "The vector $(X) is not a tangent vector to $(p) on $(M), since it is not orthogonal in the embedding.",
        )
    end
    return nothing
end
function ManifoldsBase.exp!(M::TestSphere, q, p, X)
    return ManifoldsBase.expt!(M, q, p, X, one(number_eltype(X)))
end
function ManifoldsBase.expt!(::TestSphere, q, p, X, t::Number)
    θ = abs(t) * norm(X)
    if θ == 0
        copyto!(q, p)
    else
        X_scale = t * sin(θ) / θ
        q .= p .* cos(θ) .+ X .* X_scale
    end
    return q
end
function ManifoldsBase.get_basis_diagonalizing(
    M::TestSphere{n},
    p,
    B::DiagonalizingOrthonormalBasis{ℝ},
) where {n}
    A = zeros(n + 1, n + 1)
    A[1, :] = transpose(p)
    A[2, :] = transpose(B.frame_direction)
    V = nullspace(A)
    κ = ones(n)
    if !iszero(B.frame_direction)
        # if we have a nonzero direction for the geodesic, add it and it gets curvature zero from the tensor
        V = hcat(B.frame_direction / norm(M, p, B.frame_direction), V)
        κ[1] = 0 # no curvature along the geodesic direction, if x!=y
    end
    T = typeof(similar(B.frame_direction))
    Ξ = [convert(T, V[:, i]) for i in 1:n]
    return CachedBasis(B, κ, Ξ)
end
@inline function nzsign(z, absz = abs(z))
    psignz = z / absz
    return ifelse(iszero(absz), one(psignz), psignz)
end
function ManifoldsBase.get_coordinates_orthonormal!(M::TestSphere, Y, p, X, ::RealNumbers)
    n = manifold_dimension(M)
    p1 = p[1]
    cosθ = abs(p1)
    λ = nzsign(p1, cosθ)
    pend, Xend = view(p, 2:(n + 1)), view(X, 2:(n + 1))
    factor = λ * X[1] / (1 + cosθ)
    Y .= Xend .- pend .* factor
    return Y
end
function ManifoldsBase.get_vector_orthonormal!(M::TestSphere, Y, p, X, ::RealNumbers)
    n = manifold_dimension(M)
    p1 = p[1]
    cosθ = abs(p1)
    λ = nzsign(p1, cosθ)
    pend = view(p, 2:(n + 1))
    pX = dot(pend, X)
    factor = pX / (1 + cosθ)
    Y[1] = -λ * pX
    Y[2:(n + 1)] .= X .- pend .* factor
    return Y
end
ManifoldsBase.injectivity_radius(::TestSphere) = π
ManifoldsBase.inner(::TestSphere, p, X, Y) = dot(X, Y)
function ManifoldsBase.inverse_retract_project!(M::TestSphere, X, p, q)
    X .= q .- p
    project!(M, X, p, X)
    return X
end
ManifoldsBase.is_flat(M::TestSphere) = manifold_dimension(M) == 1
function ManifoldsBase.log!(::TestSphere, X, p, q)
    cosθ = clamp(real(dot(p, q)), -1, 1)
    X .= q .- p .* cosθ
    nrmX = norm(X)
    if nrmX > 0
        X .*= acos(cosθ) / nrmX
    end
    return X
end
ManifoldsBase.manifold_dimension(::TestSphere{N}) where {N} = N
function ManifoldsBase.parallel_transport_to!(::TestSphere, Y, p, X, q)
    m = p .+ q
    mnorm2 = real(dot(m, m))
    factor = 2 * real(dot(X, q)) / mnorm2
    Y .= X .- m .* factor
    return Y
end
function ManifoldsBase.project!(::TestSphere, q, p)
    q .= p ./ norm(p)
    return q
end
function ManifoldsBase.project!(::TestSphere, Y, p, X)
    Y .= X .- p .* real(dot(p, X))
    return Y
end
function Random.rand!(M::TestSphere, pX; vector_at = nothing, σ = one(eltype(pX)))
    return rand!(Random.default_rng(), M, pX; vector_at = vector_at, σ = σ)
end
function Random.rand!(
    rng::AbstractRNG,
    M::TestSphere,
    pX;
    vector_at = nothing,
    σ = one(eltype(pX)),
)
    if vector_at === nothing
        project!(M, pX, randn(rng, eltype(pX), representation_size(M)))
    else
        n = σ * randn(rng, eltype(pX), size(pX)) # Gaussian in embedding
        project!(M, pX, vector_at, n) #project to TpM (keeps Gaussianness)
    end
    return pX
end
ManifoldsBase.representation_size(::TestSphere{N}) where {N} = (N + 1,)
function ManifoldsBase.retract_project_t!(M::TestSphere, q, p, X, t::Number)
    q .= p .+ t .* X
    project!(M, q, q)
    return q
end
function ManifoldsBase.riemann_tensor!(M::TestSphere, Xresult, p, X, Y, Z)
    innerZX = inner(M, p, Z, X)
    innerZY = inner(M, p, Z, Y)
    Xresult .= innerZY .* X .- innerZX .* Y
    return Xresult
end
function ManifoldsBase.sectional_curvature_max(M::TestSphere)
    return ifelse(manifold_dimension(M) == 1, 0.0, 1.0)
end
function ManifoldsBase.sectional_curvature_min(M::TestSphere)
    return ifelse(manifold_dimension(M) == 1, 0.0, 1.0)
end

# from Absil, Mahony, Trumpf, 2013 https://sites.uclouvain.be/absil/2013-01/Weingarten_07PA_techrep.pdf
function ManifoldsBase.Weingarten!(::TestSphere, Y, p, X, V)
    return Y .= -X * p' * V
end

#
#
# minimal implementation of SPD (for its nontrivial metric conversion)
# ---
struct TestSPD <: AbstractManifold{ℝ}
    n::Int
end
function ManifoldsBase.change_metric!(::TestSPD, Y, ::EuclideanMetric, p, X)
    Y .= p * X
    return Y
end
function ManifoldsBase.change_representer!(::TestSPD, Y, ::EuclideanMetric, p, X)
    Y .= p * X * p
    return Y
end
function ManifoldsBase.inner(::TestSPD, p, X, Y)
    F = cholesky(Symmetric(convert(AbstractMatrix, p)))
    return dot((F \ Symmetric(X)), (Symmetric(Y) / F))
end
function spd_sqrt_and_sqrt_inv(p::AbstractMatrix)
    e = eigen(Symmetric(p))
    U = e.vectors
    S = max.(e.values, floatmin(eltype(e.values)))
    Ssqrt = Diagonal(sqrt.(S))
    SsqrtInv = Diagonal(1 ./ sqrt.(S))
    return (Symmetric(U * Ssqrt * transpose(U)), Symmetric(U * SsqrtInv * transpose(U)))
end
function ManifoldsBase.log!(::TestSPD, X, p, q)
    (p_sqrt, p_sqrt_inv) = spd_sqrt_and_sqrt_inv(p)
    T = Symmetric(p_sqrt_inv * convert(AbstractMatrix, q) * p_sqrt_inv)
    e2 = eigen(T)
    Se = Diagonal(log.(max.(e2.values, eps())))
    pUe = p_sqrt * e2.vectors
    return mul!(X, pUe, Se * transpose(pUe))
end
ManifoldsBase.representation_size(M::TestSPD) = (M.n, M.n)

#
#
# A simple Manifold based on a projection onto a subspace
struct ProjManifold <: AbstractManifold{ℝ} end
ManifoldsBase.inner(::ProjManifold, p, X, Y) = dot(X, Y)
ManifoldsBase.project!(::ProjManifold, Y, p, X) = (Y .= X .- dot(p, X) .* p)
ManifoldsBase.representation_size(::ProjManifold) = (2, 3)
ManifoldsBase.manifold_dimension(::ProjManifold) = 5
ManifoldsBase.get_vector_orthonormal(::ProjManifold, p, X, N) = reverse(X)
#
#
# A second simple Manifold based on a projection onto a subspace
struct ProjectionTestManifold <: AbstractManifold{ℝ} end
ManifoldsBase.inner(::ProjectionTestManifold, ::Any, X, Y) = dot(X, Y)
function ManifoldsBase.project!(::ProjectionTestManifold, Y, p, X)
    Y .= X .- dot(p, X) .* p
    Y[end] = 0
    return Y
end
ManifoldsBase.manifold_dimension(::ProjectionTestManifold) = 100

#
# Thre Non-Things to check for the correct errors in case functions are not implemented
struct NonManifold <: AbstractManifold{ℝ} end
struct NonBasis <: ManifoldsBase.AbstractBasis{ℝ,TangentSpaceType} end
struct NonMPoint <: AbstractManifoldPoint end
struct NonTangentVector <: AbstractTangentVector end
struct NonCotangentVector <: AbstractCotangentVector end
struct NotImplementedRetraction <: AbstractRetractionMethod end
struct NotImplementedInverseRetraction <: AbstractInverseRetractionMethod end
*(t::Float64, X::NonTangentVector) = X

struct NonBroadcastBasisThing{T}
    v::T
end
+(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing) = NonBroadcastBasisThing(a.v + b.v)
*(α, a::NonBroadcastBasisThing) = NonBroadcastBasisThing(α * a.v)
-(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing) = NonBroadcastBasisThing(a.v - b.v)

function ManifoldsBase.isapprox(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing)
    return isapprox(a.v, b.v)
end

function ManifoldsBase.number_eltype(a::NonBroadcastBasisThing)
    return typeof(reduce(+, one(number_eltype(eti)) for eti in a.v))
end

ManifoldsBase.allocate(a::NonBroadcastBasisThing) = NonBroadcastBasisThing(allocate(a.v))
function ManifoldsBase.allocate(a::NonBroadcastBasisThing, ::Type{T}) where {T}
    return NonBroadcastBasisThing(allocate(a.v, T))
end
function ManifoldsBase.allocate(::NonBroadcastBasisThing, ::Type{T}, s::Integer) where {T}
    return Vector{T}(undef, s)
end

function Base.copyto!(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing)
    copyto!(a.v, b.v)
    return a
end

function ManifoldsBase.log!(
    ::DefaultManifold,
    v::NonBroadcastBasisThing,
    x::NonBroadcastBasisThing,
    y::NonBroadcastBasisThing,
)
    return copyto!(v, y - x)
end

function ManifoldsBase.exp!(
    ::DefaultManifold,
    y::NonBroadcastBasisThing,
    x::NonBroadcastBasisThing,
    v::NonBroadcastBasisThing,
)
    return copyto!(y, x + v)
end

function ManifoldsBase.get_basis_orthonormal(
    ::DefaultManifold{ℝ},
    p::NonBroadcastBasisThing,
    𝔽::RealNumbers,
)
    return CachedBasis(
        DefaultOrthonormalBasis(𝔽),
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i)) for
            i in eachindex(p.v)
        ],
    )
end
function ManifoldsBase.get_basis_orthogonal(
    ::DefaultManifold{ℝ},
    p::NonBroadcastBasisThing,
    𝔽::RealNumbers,
)
    return CachedBasis(
        DefaultOrthogonalBasis(𝔽),
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i)) for
            i in eachindex(p.v)
        ],
    )
end
function ManifoldsBase.get_basis_default(
    ::DefaultManifold{ℝ},
    p::NonBroadcastBasisThing,
    N::ManifoldsBase.RealNumbers,
)
    return CachedBasis(
        DefaultBasis(N),
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i)) for
            i in eachindex(p.v)
        ],
    )
end

function ManifoldsBase.get_coordinates_orthonormal!(
    M::DefaultManifold,
    Y,
    ::NonBroadcastBasisThing,
    X::NonBroadcastBasisThing,
    ::RealNumbers,
)
    copyto!(Y, reshape(X.v, manifold_dimension(M)))
    return Y
end

function ManifoldsBase.get_vector_orthonormal!(
    M::DefaultManifold,
    Y::NonBroadcastBasisThing,
    ::NonBroadcastBasisThing,
    X,
    ::RealNumbers,
)
    copyto!(Y.v, reshape(X, representation_size(M)))
    return Y
end

function ManifoldsBase.inner(
    ::DefaultManifold,
    ::NonBroadcastBasisThing,
    X::NonBroadcastBasisThing,
    Y::NonBroadcastBasisThing,
)
    return dot(X.v, Y.v)
end

ManifoldsBase._get_vector_cache_broadcast(::NonBroadcastBasisThing) = Val(false)

DiagonalizingBasisProxy() = DiagonalizingOrthonormalBasis([1.0, 0.0, 0.0])

#
#
# Vector Space types
struct TestVectorSpaceType <: ManifoldsBase.VectorSpaceType end

struct TestFiberType <: ManifoldsBase.FiberType end

function ManifoldsBase.fiber_dimension(M::AbstractManifold, ::TestFiberType)
    return 2 * manifold_dimension(M)
end

function ManifoldsBase.vector_space_dimension(M::AbstractManifold, ::TestVectorSpaceType)
    return 2 * manifold_dimension(M)
end

#
#
# DefaultManifold with a few artificiual retrations

struct CustomDefinedRetraction <: ManifoldsBase.AbstractRetractionMethod end
struct CustomUndefinedRetraction <: ManifoldsBase.AbstractRetractionMethod end
struct CustomDefinedKeywordRetraction <: ManifoldsBase.AbstractRetractionMethod end
struct CustomDefinedKeywordInverseRetraction <:
       ManifoldsBase.AbstractInverseRetractionMethod end
struct CustomDefinedInverseRetraction <: ManifoldsBase.AbstractInverseRetractionMethod end

struct DefaultPoint{T} <: AbstractManifoldPoint
    value::T
end
Base.convert(::Type{DefaultPoint{T}}, v::T) where {T} = DefaultPoint(v)
Base.eltype(p::DefaultPoint) = eltype(p.value)
Base.length(p::DefaultPoint) = length(p.value)
struct DefaultTangentVector{T} <: AbstractTangentVector
    value::T
end
function Base.convert(::Type{DefaultTangentVector{T}}, ::DefaultPoint, v::T) where {T}
    return DefaultTangentVector(v)
end
Base.eltype(X::DefaultTangentVector) = eltype(X.value)
function Base.fill!(X::DefaultTangentVector, x)
    fill!(X.value, x)
    return X
end
function ManifoldsBase.allocate_result_type(
    ::DefaultManifold,
    ::typeof(log),
    ::Tuple{DefaultPoint,DefaultPoint},
)
    return DefaultTangentVector
end
function ManifoldsBase.allocate_result_type(
    ::DefaultManifold,
    ::typeof(inverse_retract),
    ::Tuple{DefaultPoint,DefaultPoint},
)
    return DefaultTangentVector
end

ManifoldsBase.@manifold_element_forwards DefaultPoint value
ManifoldsBase.@manifold_vector_forwards DefaultTangentVector value
ManifoldsBase.@default_manifold_fallbacks ManifoldsBase.DefaultManifold DefaultPoint DefaultTangentVector value value

function ManifoldsBase._injectivity_radius(::DefaultManifold, ::CustomDefinedRetraction)
    return 10.0
end
function ManifoldsBase._retract_t!(
    M::DefaultManifold,
    q,
    p,
    X,
    t::Number,
    ::CustomDefinedKeywordRetraction;
    kwargs...,
)
    return retract_custom_kw_t!(M, q, p, X, t; kwargs...)
end
function retract_custom_kw_t!(
    ::DefaultManifold,
    q::DefaultPoint,
    p::DefaultPoint,
    X::DefaultTangentVector,
    t::Number;
    scale = 2.0,
)
    q.value .= scale .* p.value .+ t .* X.value
    return q
end
function ManifoldsBase._inverse_retract!(
    M::DefaultManifold,
    X,
    p,
    q,
    ::CustomDefinedKeywordInverseRetraction;
    kwargs...,
)
    return inverse_retract_custom_kw!(M, X, p, q; kwargs...)
end
function inverse_retract_custom_kw!(
    ::DefaultManifold,
    X::DefaultTangentVector,
    p::DefaultPoint,
    q::DefaultPoint;
    scale = 2.0,
)
    X.value .= q.value - scale * p.value
    return X
end

function ManifoldsBase.retract_t!(
    ::DefaultManifold,
    q::DefaultPoint,
    p::DefaultPoint,
    X::DefaultTangentVector,
    t::Number,
    ::CustomDefinedRetraction,
)
    q.value .= 2 .* p.value .+ t * X.value
    return q
end

function ManifoldsBase.inverse_retract!(
    ::DefaultManifold,
    X::DefaultTangentVector,
    p::DefaultPoint,
    q::DefaultPoint,
    ::CustomDefinedInverseRetraction,
)
    X.value .= q.value .- 2 .* p.value
    return X
end

struct MatrixVectorTransport{T} <: AbstractVector{T}
    m::Matrix{T}
end
# dummy retractions, inverse retracions for fallback tests - mutating should be enough
ManifoldsBase.retract_polar_t!(::DefaultManifold, q, p, X, t::Number) = (q .= p .+ t .* X)
ManifoldsBase.retract_project_t!(::DefaultManifold, q, p, X, t::Number) = (q .= p .+ t .* X)
ManifoldsBase.retract_qr_t!(::DefaultManifold, q, p, X, t::Number) = (q .= p .+ t .* X)
function ManifoldsBase.retract_exp_ode_t!(
    ::DefaultManifold,
    q,
    p,
    X,
    t::Number,
    m::AbstractRetractionMethod,
    B::ManifoldsBase.AbstractBasis,
)
    return (q .= p .+ t .* X)
end

function ManifoldsBase.retract_pade_t!(
    ::DefaultManifold,
    q,
    p,
    X,
    t::Number,
    m::PadeRetraction,
)
    return (q .= p .+ t .* X)
end
function ManifoldsBase.retract_sasaki_t!(
    ::DefaultManifold,
    q,
    p,
    X,
    t::Number,
    ::SasakiRetraction,
)
    return (q .= p .+ t .* X)
end
ManifoldsBase.retract_softmax_t!(::DefaultManifold, q, p, X, t::Number) = (q .= p .+ t .* X)
ManifoldsBase.get_embedding(M::DefaultManifold) = M # dummy embedding
ManifoldsBase.inverse_retract_polar!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
ManifoldsBase.inverse_retract_project!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
ManifoldsBase.inverse_retract_qr!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
ManifoldsBase.inverse_retract_softmax!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
function ManifoldsBase.inverse_retract_nlsolve!(
    ::DefaultManifold,
    Y,
    p,
    q,
    m::NLSolveInverseRetraction,
)
    return (Y .= q .- p)
end
Base.getindex(x::MatrixVectorTransport, i) = x.m[:, i]
Base.size(x::MatrixVectorTransport) = (size(x.m, 2),)

struct TestArrayRepresentation <: AbstractPowerRepresentation end

const TestPowerManifoldMultidimensional =
    AbstractPowerManifold{𝔽,<:AbstractManifold{𝔽},TestArrayRepresentation} where {𝔽}

export CustomDefinedInverseRetraction, CustomDefinedKeywordInverseRetraction
export CustomDefinedKeywordRetraction, CustomDefinedRetraction, CustomUndefinedRetraction
export DefaultPoint, DefaultTangentVector
export DiagonalizingBasisProxy
export MatrixVectorTransport
export NonManifold, NonBasis, NonBroadcastBasisThing
export NonMPoint, NonTangentVector, NonCotangentVector
export NotImplementedRetraction, NotImplementedInverseRetraction
export ProjManifold, ProjectionTestManifold
export TestSphere, TestSPD
export TestVectorSpaceType, TestFiberType
end
