using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: DefaultManifold, ℝ, ℂ
using ManifoldsBase: CotangentSpace, CotangentSpaceType, TangentSpace, TangentSpaceType
using ManifoldsBase: FVector
using Test
import Base: +, -, *, copyto!, isapprox
import ManifoldsBase:
    allocate,
    get_vector_orthonormal!,
    get_coordinates_orthonormal!,
    get_basis_orthogonal,
    get_basis_orthonormal


struct ProjManifold <: AbstractManifold{ℝ} end

ManifoldsBase.inner(::ProjManifold, p, X, Y) = dot(X, Y)
ManifoldsBase.project!(::ProjManifold, Y, p, X) = (Y .= X .- dot(p, X) .* p)
ManifoldsBase.representation_size(::ProjManifold) = (2, 3)
ManifoldsBase.manifold_dimension(::ProjManifold) = 5
ManifoldsBase.get_vector_orthonormal(::ProjManifold, p, X, N) = reverse(X)

struct ProjectionTestManifold <: AbstractManifold{ℝ} end

ManifoldsBase.inner(::ProjectionTestManifold, ::Any, X, Y) = dot(X, Y)
function ManifoldsBase.project!(::ProjectionTestManifold, Y, p, X)
    Y .= X .- dot(p, X) .* p
    Y[end] = 0
    return Y
end
ManifoldsBase.manifold_dimension(::ProjectionTestManifold) = 100

@testset "Projected and arbitrary orthonormal basis" begin
    M = ProjManifold()
    x = [
        sqrt(2)/2 0.0 0.0
        0.0 sqrt(2)/2 0.0
    ]

    for pB in (ProjectedOrthonormalBasis(:svd), ProjectedOrthonormalBasis(:gram_schmidt))
        if pB isa ProjectedOrthonormalBasis{:gram_schmidt,ℝ}
            pb = get_basis(
                M,
                x,
                pB;
                warn_linearly_dependent = true,
                skip_linearly_dependent = true,
            ) # skip V4, wich is -V1 after proj.
            @test_throws ErrorException get_basis(M, x, pB) # error
        else
            pb = get_basis(M, x, pB) # skips V4 automatically
        end
        @test number_system(pb) == ℝ
        N = manifold_dimension(M)
        @test isa(pb, CachedBasis)
        @test CachedBasis(pb) === pb
        @test length(get_vectors(M, x, pb)) == N
        # test orthonormality
        for i in 1:N
            @test norm(M, x, get_vectors(M, x, pb)[i]) ≈ 1
            for j in (i + 1):N
                @test inner(M, x, get_vectors(M, x, pb)[i], get_vectors(M, x, pb)[j]) ≈ 0 atol =
                    1e-15
            end
        end
    end
    aonb = get_basis(M, x, DefaultOrthonormalBasis())
    @test size(get_vectors(M, x, aonb)) == (5,)
    @test get_vectors(M, x, aonb)[1] ≈ [0, 0, 0, 0, 1]

    @testset "Gram-Schmidt" begin
        # for a basis
        M = ManifoldsBase.DefaultManifold(3)
        p = zeros(3)
        V = [[2.0, 0.0, 0.0], [1.1, 2.2, 0.0], [0.0, 3.3, 4.4]]
        b1 = ManifoldsBase.gram_schmidt(M, p, V)
        b2 = ManifoldsBase.gram_schmidt(M, zeros(3), CachedBasis(DefaultBasis(), V))
        @test b1 == get_vectors(M, p, b2)
        # projected gram schmidt
        tm = ProjectionTestManifold()
        bt = ProjectedOrthonormalBasis(:gram_schmidt)
        p = [sqrt(2) / 2, 0.0, sqrt(2) / 2, 0.0, 0.0]
        @test_throws ErrorException get_basis(tm, p, bt)
        b = get_basis(
            tm,
            p,
            bt;
            return_incomplete_set = true,
            skip_linearly_dependent = true, #skips 3 and 5
        )
        @test length(get_vectors(tm, p, b)) == 3
        @test_throws ErrorException ManifoldsBase.gram_schmidt(M, p, [V[1]])
        @test_throws ErrorException ManifoldsBase.gram_schmidt(
            M,
            p,
            [V[1], V[1], V[1]];
            skip_linearly_dependent = true,
        )
    end
end

struct NonManifold <: AbstractManifold{ℝ} end
struct NonBasis <: ManifoldsBase.AbstractBasis{ℝ,TangentSpaceType} end

struct NonBroadcastBasisThing{T}
    v::T
end

+(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing) = NonBroadcastBasisThing(a.v + b.v)
*(α, a::NonBroadcastBasisThing) = NonBroadcastBasisThing(α * a.v)
-(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing) = NonBroadcastBasisThing(a.v - b.v)

isapprox(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing) = isapprox(a.v, b.v)

function ManifoldsBase.number_eltype(a::NonBroadcastBasisThing)
    return typeof(reduce(+, one(number_eltype(eti)) for eti in a.v))
end

allocate(a::NonBroadcastBasisThing) = NonBroadcastBasisThing(allocate(a.v))
function allocate(a::NonBroadcastBasisThing, ::Type{T}) where {T}
    return NonBroadcastBasisThing(allocate(a.v, T))
end
allocate(::NonBroadcastBasisThing, ::Type{T}, s::Integer) where {S,T} = Vector{T}(undef, s)

function copyto!(a::NonBroadcastBasisThing, b::NonBroadcastBasisThing)
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
    ::DefaultManifold,
    p::NonBroadcastBasisThing,
    𝔽,
)
    return CachedBasis(
        B,
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i)) for
            i in eachindex(p.v)
        ],
    )
end
function ManifoldsBase.get_basis_orthogonal(::DefaultManifold, p::NonBroadcastBasisThing, 𝔽)
    return CachedBasis(
        DefaultOrthogonalBasis(𝔽),
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i)) for
            i in eachindex(p.v)
        ],
    )
end
function ManifoldsBase.get_basis_default(
    ::DefaultManifold,
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
    𝔽,
)
    copyto!(Y, reshape(X.v, manifold_dimension(M)))
    return Y
end

function ManifoldsBase.get_vector_orthonormal!(
    M::DefaultManifold,
    Y::NonBroadcastBasisThing,
    ::NonBroadcastBasisThing,
    X,
    𝔽,
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

@testset "ManifoldsBase.jl stuff" begin

    @testset "Errors" begin
        m = NonManifold()
        onb = DefaultOrthonormalBasis()

        @test_throws MethodError get_basis(m, [0], onb)
        @test_throws MethodError get_basis(m, [0], NonBasis())
        @test_throws MethodError get_coordinates(m, [0], [0], onb)
        @test_throws MethodError get_coordinates!(m, [0], [0], [0], onb)
        @test_throws MethodError get_vector(m, [0], [0], onb)
        @test_throws MethodError get_vector!(m, [0], [0], [0], onb)
        @test_throws MethodError get_vectors(m, [0], NonBasis())
    end

    M = DefaultManifold(3)

    @testset "Constructors" begin
        @test DefaultBasis{ℂ,TangentSpaceType}() === DefaultBasis(ℂ)
        @test DefaultOrthogonalBasis{ℂ,TangentSpaceType}() === DefaultOrthogonalBasis(ℂ)
        @test DefaultOrthonormalBasis{ℂ,TangentSpaceType}() === DefaultOrthonormalBasis(ℂ)

        @test DefaultBasis{ℂ}(CotangentSpace) === DefaultBasis(ℂ, CotangentSpace)
        @test DefaultOrthogonalBasis{ℂ}(CotangentSpace) ===
              DefaultOrthogonalBasis(ℂ, CotangentSpace)
        @test DefaultOrthonormalBasis{ℂ}(CotangentSpace) ===
              DefaultOrthonormalBasis(ℂ, CotangentSpace)
    end

    _pts = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    @testset "basis representation" for BT in (
            DefaultBasis,
            DefaultOrthonormalBasis,
            DefaultOrthogonalBasis,
            DiagonalizingBasisProxy,
        ),
        pts in (_pts, map(NonBroadcastBasisThing, _pts))

        if BT == DiagonalizingBasisProxy && pts !== _pts
            continue
        end
        v1 = log(M, pts[1], pts[2])
        @test ManifoldsBase.number_of_coordinates(M, BT()) == 3

        if BT != DiagonalizingBasisProxy
            vb = get_coordinates(M, pts[1], v1, BT())
            @test isa(vb, AbstractVector)
            vbi = get_vector(M, pts[1], vb, BT())
            @test isapprox(M, pts[1], v1, vbi)
        end

        b = get_basis(M, pts[1], BT())
        if BT != DiagonalizingBasisProxy
            if pts[1] isa Array
                @test isa(b, CachedBasis{ℝ,BT{ℝ,TangentSpaceType},Vector{Vector{Float64}}})
            else
                @test isa(
                    b,
                    CachedBasis{
                        ℝ,
                        BT{ℝ,TangentSpaceType},
                        Vector{NonBroadcastBasisThing{Vector{Float64}}},
                    },
                )
            end
        end
        @test get_basis(M, pts[1], b) === b
        N = manifold_dimension(M)
        @test length(get_vectors(M, pts[1], b)) == N
        # check orthonormality
        if BT isa DefaultOrthonormalBasis && pts[1] isa Vector
            for i in 1:N
                @test norm(M, pts[1], get_vectors(M, pts[1], b)[i]) ≈ 1
                for j in (i + 1):N
                    @test inner(
                        M,
                        pts[1],
                        get_vectors(M, pts[1], b)[i],
                        get_vectors(M, pts[1], b)[j],
                    ) ≈ 0
                end
            end
            # check that the coefficients correspond to the basis
            for i in 1:N
                @test inner(M, pts[1], v1, get_vectors(M, pts[1], b)[i]) ≈ vb[i]
            end
        end

        if BT != DiagonalizingBasisProxy
            @test get_coordinates(M, pts[1], v1, b) ≈ get_coordinates(M, pts[1], v1, BT())
            @test get_vector(M, pts[1], vb, b) ≈ get_vector(M, pts[1], vb, BT())
        end

        v1c = Vector{Float64}(undef, 3)
        get_coordinates!(M, v1c, pts[1], v1, b)
        @test v1c ≈ get_coordinates(M, pts[1], v1, b)

        v1cv = allocate(v1)
        get_vector!(M, v1cv, pts[1], v1c, b)
        @test isapprox(M, pts[1], v1, v1cv)
    end
end

@testset "Complex DeaultManifold with real and complex Cached Bases" begin
    M = ManifoldsBase.DefaultManifold(3; field = ℂ)
    p = [1.0, 2.0im, 3.0]
    X = [1.2, 2.2im, 2.3im]
    b = [Matrix{Float64}(I, 3, 3)[:, i] for i in 1:3]
    Bℝ = CachedBasis(DefaultOrthonormalBasis{ℝ}(), b)
    aℝ = get_coordinates(M, p, X, Bℝ)
    Yℝ = get_vector(M, p, aℝ, Bℝ)
    @test Yℝ ≈ X
    @test ManifoldsBase.number_of_coordinates(M, Bℝ) == 3

    bℂ = [b..., (b .* 1im)...]
    Bℂ = CachedBasis(DefaultOrthonormalBasis{ℂ}(), bℂ)
    aℂ = get_coordinates(M, p, X, Bℂ)
    Yℂ = get_vector(M, p, aℂ, Bℂ)
    @test Yℂ ≈ X
    @test ManifoldsBase.number_of_coordinates(M, Bℂ) == 6
end

@testset "Basis show methods" begin
    @test sprint(show, DefaultBasis()) == "DefaultBasis(ℝ)"
    @test sprint(show, DefaultOrthogonalBasis()) == "DefaultOrthogonalBasis(ℝ)"
    @test sprint(show, DefaultOrthonormalBasis()) == "DefaultOrthonormalBasis(ℝ)"
    @test sprint(show, DefaultOrthonormalBasis(ℂ)) == "DefaultOrthonormalBasis(ℂ)"
    @test sprint(show, GramSchmidtOrthonormalBasis(ℂ)) == "GramSchmidtOrthonormalBasis(ℂ)"
    @test sprint(show, ProjectedOrthonormalBasis(:svd)) ==
          "ProjectedOrthonormalBasis(:svd, ℝ)"
    @test sprint(show, ProjectedOrthonormalBasis(:gram_schmidt, ℂ)) ==
          "ProjectedOrthonormalBasis(:gram_schmidt, ℂ)"

    diag_onb = DiagonalizingOrthonormalBasis(Float64[1, 2, 3])
    @test sprint(show, "text/plain", diag_onb) == """
    DiagonalizingOrthonormalBasis(ℝ) with eigenvalue 0 in direction:
    3-element $(sprint(show, Vector{Float64})):
      1.0
      2.0
      3.0"""

    M = DefaultManifold(2, 3)
    x = collect(reshape(1.0:6.0, (2, 3)))
    pb = get_basis(M, x, DefaultOrthonormalBasis())

    @test sprint(show, "text/plain", pb) == """
    DefaultOrthonormalBasis(ℝ) with 6 basis vectors:
     E1 =
      2×3 $(sprint(show, Matrix{Float64})):
       1.0  0.0  0.0
       0.0  0.0  0.0
     E2 =
      2×3 $(sprint(show, Matrix{Float64})):
       0.0  0.0  0.0
       1.0  0.0  0.0
     ⋮
     E5 =
      2×3 $(sprint(show, Matrix{Float64})):
       0.0  0.0  1.0
       0.0  0.0  0.0
     E6 =
      2×3 $(sprint(show, Matrix{Float64})):
       0.0  0.0  0.0
       0.0  0.0  1.0"""

    b = DiagonalizingOrthonormalBasis(get_vectors(M, x, pb)[1])
    dpb = CachedBasis(b, Float64[1, 2, 3, 4, 5, 6], get_vectors(M, x, pb))
    @test sprint(show, "text/plain", dpb) == """
    DiagonalizingOrthonormalBasis(ℝ) with eigenvalue 0 in direction:
     2×3 $(sprint(show, Matrix{Float64})):
       1.0  0.0  0.0
       0.0  0.0  0.0
    and 6 basis vectors.
    Basis vectors:
     E1 =
      2×3 $(sprint(show, Matrix{Float64})):
       1.0  0.0  0.0
       0.0  0.0  0.0
     E2 =
      2×3 $(sprint(show, Matrix{Float64})):
       0.0  0.0  0.0
       1.0  0.0  0.0
     ⋮
     E5 =
      2×3 $(sprint(show, Matrix{Float64})):
       0.0  0.0  1.0
       0.0  0.0  0.0
     E6 =
      2×3 $(sprint(show, Matrix{Float64})):
       0.0  0.0  0.0
       0.0  0.0  1.0
    Eigenvalues:
     6-element $(sprint(show, Vector{Float64})):
      1.0
      2.0
      3.0
      4.0
      5.0
      6.0"""

    M = DefaultManifold(1, 1, 1)
    x = reshape(Float64[1], (1, 1, 1))
    pb = get_basis(M, x, DefaultOrthonormalBasis())
    @test sprint(show, "text/plain", pb) == """
    DefaultOrthonormalBasis(ℝ) with 1 basis vector:
     E1 =
      1×1×1 $(sprint(show, Array{Float64,3})):
      [:, :, 1] =
       1.0"""

    dpb = CachedBasis(
        DiagonalizingOrthonormalBasis(get_vectors(M, x, pb)),
        Float64[1],
        get_vectors(M, x, pb),
    )

    @test sprint(show, "text/plain", dpb) == """
    DiagonalizingOrthonormalBasis(ℝ) with eigenvalue 0 in direction:
     1-element $(sprint(show, Vector{Array{Float64,3}})):
       $(sprint(show, dpb.data.frame_direction[1]))
    and 1 basis vector.
    Basis vectors:
     E1 =
      1×1×1 $(sprint(show, Array{Float64,3})):
      [:, :, 1] =
       1.0
    Eigenvalues:
     1-element $(sprint(show, Vector{Float64})):
      1.0"""
end

@testset "Bases of cotangent spaces" begin
    b1 = DefaultOrthonormalBasis(ℝ, CotangentSpace)
    @test b1.vector_space == CotangentSpace

    b2 = DefaultOrthogonalBasis(ℝ, CotangentSpace)
    @test b2.vector_space == CotangentSpace

    b3 = DefaultBasis(ℝ, CotangentSpace)
    @test b3.vector_space == CotangentSpace

    M = DefaultManifold(2; field = ℂ)
    p = [1.0, 2.0im]
    b1_d = ManifoldsBase.dual_basis(M, p, b1)
    @test b1_d isa DefaultOrthonormalBasis
    @test b1_d.vector_space == TangentSpace

    b1_d_d = ManifoldsBase.dual_basis(M, p, b1_d)
    @test b1_d_d isa DefaultOrthonormalBasis
    @test b1_d_d.vector_space == CotangentSpace
end

@testset "FVector" begin
    @test sprint(show, TangentSpace) == "TangentSpace"
    @test sprint(show, CotangentSpace) == "CotangentSpace"
    tvs = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    fv_tvs = map(v -> TFVector(v, DefaultOrthonormalBasis()), tvs)
    fv1 = fv_tvs[1]
    tv1s = allocate(fv_tvs[1])
    @test isa(tv1s, FVector)
    @test tv1s.type == TangentSpace
    @test size(tv1s.data) == size(tvs[1])
    @test number_eltype(tv1s) == number_eltype(tvs[1])
    @test number_eltype(tv1s) == number_eltype(typeof(tv1s))
    @test isa(fv1 + fv1, FVector)
    @test (fv1 + fv1).type == TangentSpace
    @test isa(fv1 - fv1, FVector)
    @test (fv1 - fv1).type == TangentSpace
    @test isa(-fv1, FVector)
    @test (-fv1).type == TangentSpace
    @test isa(2 * fv1, FVector)
    @test (2 * fv1).type == TangentSpace
    tv1s_32 = allocate(fv_tvs[1], Float32)
    @test isa(tv1s, FVector)
    @test eltype(tv1s_32.data) === Float32
    copyto!(tv1s, fv_tvs[2])
    @test isapprox(tv1s.data, fv_tvs[2].data)

    @test sprint(show, fv1) == "TFVector([1.0, 0.0, 0.0], $(fv1.basis))"

    cofv1 = CoTFVector(tvs[1], DefaultOrthonormalBasis(ℝ, CotangentSpace))
    @test cofv1 isa CoTFVector
    @test sprint(show, cofv1) == "CoTFVector([1.0, 0.0, 0.0], $(fv1.basis))"
end

@testset "vector_space_dimension" begin
    M = ManifoldsBase.DefaultManifold(3)
    MC = ManifoldsBase.DefaultManifold(3; field = ℂ)
    @test ManifoldsBase.vector_space_dimension(M, TangentSpace) == 3
    @test ManifoldsBase.vector_space_dimension(M, CotangentSpace) == 3
    @test ManifoldsBase.vector_space_dimension(MC, TangentSpace) == 6
    @test ManifoldsBase.vector_space_dimension(MC, CotangentSpace) == 6
end
