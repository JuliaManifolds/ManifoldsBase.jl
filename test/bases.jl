using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: DefaultManifold, ℝ, ℂ
using Test
import Base: +, -, *, copyto!, isapprox
import ManifoldsBase: allocate

struct ProjManifold <: Manifold{ℝ} end

ManifoldsBase.inner(::ProjManifold, x, w, v) = dot(w, v)
ManifoldsBase.project!(::ProjManifold, w, x, v) = (w .= v .- dot(x, v) .* x)
ManifoldsBase.representation_size(::ProjManifold) = (2, 3)
ManifoldsBase.manifold_dimension(::ProjManifold) = 5
ManifoldsBase.get_vector(::ProjManifold, x, v, ::DefaultOrthonormalBasis) = reverse(v)

@testset "Dispatch" begin
    @test ManifoldsBase.decorator_transparent_dispatch(
        get_basis,
        DefaultManifold(3),
        [0.0, 0.0, 0.0],
        DefaultBasis(),
    ) === Val(:parent)
    @test ManifoldsBase.decorator_transparent_dispatch(
        get_coordinates,
        DefaultManifold(3),
        [0.0, 0.0, 0.0],
    ) === Val(:parent)
    @test ManifoldsBase.decorator_transparent_dispatch(
        get_coordinates!,
        DefaultManifold(3),
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ) === Val(:transparent)
    @test ManifoldsBase.decorator_transparent_dispatch(
        get_vector,
        DefaultManifold(3),
        [0.0, 0.0, 0.0],
    ) === Val(:parent)
    @test ManifoldsBase.decorator_transparent_dispatch(
        get_vector!,
        DefaultManifold(3),
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ) === Val(:transparent)
end

struct ProjectionTestManifold <: Manifold{ℝ} end

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

    for pb in (ProjectedOrthonormalBasis(:svd), ProjectedOrthonormalBasis(:gram_schmidt))
        pb = get_basis(M, x, pb)
        @test number_system(pb) == ℝ
        @test get_basis(M, x, pb) == pb
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
        # check projection idempotency
        for i in 1:N
            @test project(M, x, get_vectors(M, x, pb)[i]) ≈ get_vectors(M, x, pb)[i]
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
            warn_linearly_dependent = true,
        )
        @test length(get_vectors(tm, p, b)) == 3
    end
end

struct NonManifold <: Manifold{ℝ} end
struct NonBasis <: ManifoldsBase.AbstractBasis{ℝ} end

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

function ManifoldsBase.get_basis(
    M::DefaultManifold,
    p::NonBroadcastBasisThing,
    B::DefaultOrthonormalBasis,
)
    return CachedBasis(
        B,
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i))
            for i in eachindex(p.v)
        ],
    )
end
function ManifoldsBase.get_basis(
    M::DefaultManifold,
    p::NonBroadcastBasisThing,
    B::DefaultOrthogonalBasis,
)
    return CachedBasis(
        B,
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i))
            for i in eachindex(p.v)
        ],
    )
end
function ManifoldsBase.get_basis(
    M::DefaultManifold,
    p::NonBroadcastBasisThing,
    B::DefaultBasis,
)
    return CachedBasis(
        B,
        [
            NonBroadcastBasisThing(ManifoldsBase._euclidean_basis_vector(p.v, i))
            for i in eachindex(p.v)
        ],
    )
end

function ManifoldsBase.get_coordinates!(
    M::DefaultManifold,
    Y,
    p::NonBroadcastBasisThing,
    X::NonBroadcastBasisThing,
    B::DefaultOrthonormalBasis,
)
    copyto!(Y, reshape(X.v, manifold_dimension(M)))
    return Y
end

function ManifoldsBase.get_vector!(
    M::DefaultManifold,
    Y::NonBroadcastBasisThing,
    p::NonBroadcastBasisThing,
    X,
    B::DefaultOrthonormalBasis,
)
    copyto!(Y.v, reshape(X, representation_size(M)))
    return Y
end

function ManifoldsBase.inner(
    ::DefaultManifold,
    x::NonBroadcastBasisThing,
    v::NonBroadcastBasisThing,
    w::NonBroadcastBasisThing,
)
    return dot(v.v, w.v)
end

ManifoldsBase._get_vector_cache_broadcast(::NonBroadcastBasisThing) = Val(false)

DiagonalizingBasisProxy() = DiagonalizingOrthonormalBasis([1.0, 0.0, 0.0])

@testset "ManifoldsBase.jl stuff" begin

    @testset "Errors" begin
        m = NonManifold()
        onb = DefaultOrthonormalBasis()

        @test_throws ErrorException get_basis(m, [0], onb)
        @test_throws ErrorException get_basis(m, [0], NonBasis())
        @test_throws ErrorException get_coordinates(m, [0], [0], onb)
        @test_throws ErrorException get_coordinates!(m, [0], [0], [0], onb)
        @test_throws ErrorException get_vector(m, [0], [0], onb)
        @test_throws ErrorException get_vector!(m, [0], [0], [0], onb)
        @test_throws ErrorException get_vectors(m, [0], NonBasis())
    end

    M = DefaultManifold(3)

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
                @test isa(b, CachedBasis{ℝ,BT{ℝ},Vector{Vector{Float64}}})
            else
                @test isa(
                    b,
                    CachedBasis{ℝ,BT{ℝ},Vector{NonBroadcastBasisThing{Vector{Float64}}}},
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

    @test sprint(show, "text/plain", DiagonalizingOrthonormalBasis(Float64[1, 2, 3])) == """
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
       [1.0]
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
