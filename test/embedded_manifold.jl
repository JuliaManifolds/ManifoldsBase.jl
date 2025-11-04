using LinearAlgebra, ManifoldsBase, Test

using ManifoldsBase: DefaultManifold, ℝ
#
# A first artificial (not real) manifold that is modelled as a submanifold
# half plane with euclidean metric is not a manifold but should test all things correctly here
#
struct HalfPlaneManifold <: AbstractDecoratorManifold{ℝ} end
struct PosQuadrantManifold <: AbstractDecoratorManifold{ℝ} end

function ManifoldsBase.get_embedding_type(::HalfPlaneManifold)
    return ManifoldsBase.EmbeddedSubmanifoldType()
end
function ManifoldsBase.get_embedding_type(::PosQuadrantManifold)
    return ManifoldsBase.EmbeddedSubmanifoldType()
end

ManifoldsBase.get_embedding(::HalfPlaneManifold) = ManifoldsBase.DefaultManifold(1, 3)
ManifoldsBase.decorated_manifold(::HalfPlaneManifold) = ManifoldsBase.DefaultManifold(2)
ManifoldsBase.representation_size(::HalfPlaneManifold) = (1, 3)

ManifoldsBase.get_embedding(::PosQuadrantManifold) = HalfPlaneManifold()
ManifoldsBase.representation_size(::PosQuadrantManifold) = (3,)

function ManifoldsBase.check_point(::HalfPlaneManifold, p)
    return p[1] > 0 ? nothing : DomainError(p[1], "p[1] ≤ 0")
end
function ManifoldsBase.check_vector(::HalfPlaneManifold, p, X)
    return X[1] > 0 ? nothing : DomainError(X[1], "X[1] ≤ 0")
end
function ManifoldsBase.check_point(::PosQuadrantManifold, p)
    return p[2] > 0 ? nothing : DomainError(p[1], "p[2] ≤ 0")
end
function ManifoldsBase.check_vector(::PosQuadrantManifold, p, X)
    return X[2] > 0 ? nothing : DomainError(X[1], "X[2] ≤ 0")
end

ManifoldsBase.embed(::HalfPlaneManifold, p) = reshape(p, 1, :)
ManifoldsBase.embed(::HalfPlaneManifold, p, X) = reshape(X, 1, :)

ManifoldsBase.project!(::HalfPlaneManifold, q, p) = (q .= [p[1] p[2] 0.0])
ManifoldsBase.project!(::HalfPlaneManifold, Y, p, X) = (Y .= [X[1] X[2] 0.0])

function ManifoldsBase.get_coordinates_orthonormal!(
        ::HalfPlaneManifold,
        Y,
        p,
        X,
        ::ManifoldsBase.RealNumbers,
    )
    return (Y .= [X[1], X[2]])
end
function ManifoldsBase.get_vector_orthonormal!(
        ::HalfPlaneManifold,
        Y,
        p,
        c,
        ::ManifoldsBase.RealNumbers,
    )
    return (Y .= [c[1] c[2] 0.0])
end

#
# A second manifold that is modelled as just isometrically embedded but not a submanifold
#
struct AnotherHalfPlaneManifold <: AbstractDecoratorManifold{ℝ} end

ManifoldsBase.get_embedding(::AnotherHalfPlaneManifold) = ManifoldsBase.DefaultManifold(3)
function ManifoldsBase.decorated_manifold(::AnotherHalfPlaneManifold)
    return ManifoldsBase.DefaultManifold(2)
end
ManifoldsBase.representation_size(::AnotherHalfPlaneManifold) = (2,)

function ManifoldsBase.get_embedding_type(::AnotherHalfPlaneManifold)
    return ManifoldsBase.IsometricallyEmbeddedManifoldType(ManifoldsBase.DirectEmbedding())
end

function ManifoldsBase.embed!(::AnotherHalfPlaneManifold, q, p)
    q[1:2] .= p
    q[3] = 0
    return q
end
function ManifoldsBase.embed!(::AnotherHalfPlaneManifold, Y, p, X)
    Y[1:2] .= X
    Y[3] = 0
    return Y
end
function ManifoldsBase.project!(::AnotherHalfPlaneManifold, q, p)
    return q .= [p[1], p[2]]
end
function ManifoldsBase.project!(::AnotherHalfPlaneManifold, Y, p, X)
    return Y .= [X[1], X[2]]
end
function ManifoldsBase.exp!(::AnotherHalfPlaneManifold, q, p, X)
    return q .= p .+ X
end
#
# Third example - explicitly mention an embedding.
#
function ManifoldsBase.embed!(
        ::EmbeddedManifold{𝔽, DefaultManifold{𝔽, nL}, DefaultManifold{𝔽2, mL}},
        q,
        p,
    ) where {nL, mL, 𝔽, 𝔽2}
    n = size(p)
    ln = length(n)
    m = size(q)
    lm = length(m)
    (length(n) > length(m)) && throw(
        DomainError(
            "Invalid embedding, since Euclidean dimension ($(n)) is longer than embedding dimension $(m).",
        ),
    )
    any(n .> m[1:ln]) && throw(
        DomainError(
            "Invalid embedding, since Euclidean dimension ($(n)) has entry larger than embedding dimensions ($(m)).",
        ),
    )
    fill!(q, 0)
    q[map(ind_n -> Base.OneTo(ind_n), n)..., ntuple(_ -> 1, lm - ln)...] .= p
    return q
end

function ManifoldsBase.project!(
        ::EmbeddedManifold{𝔽, DefaultManifold{𝔽, nL}, DefaultManifold{𝔽2, mL}},
        q,
        p,
    ) where {nL, mL, 𝔽, 𝔽2}
    n = size(p)
    ln = length(n)
    m = size(q)
    lm = length(m)
    (length(n) < length(m)) && throw(
        DomainError(
            "Invalid embedding, since Euclidean dimension ($(n)) is longer than embedding dimension $(m).",
        ),
    )
    any(n .< m[1:ln]) && throw(
        DomainError(
            "Invalid embedding, since Euclidean dimension ($(n)) has entry larger than embedding dimensions ($(m)).",
        ),
    )
    #  fill q with the „top left edge“ of p.
    q .= p[map(i -> Base.OneTo(i), m)..., ntuple(_ -> 1, lm - ln)...]
    return q
end

#
# A manifold that is a submanifold but otherwise has not implementations
#
struct NotImplementedEmbeddedSubManifold <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.get_embedding(::NotImplementedEmbeddedSubManifold)
    return ManifoldsBase.DefaultManifold(3)
end
function ManifoldsBase.decorated_manifold(::NotImplementedEmbeddedSubManifold)
    return ManifoldsBase.DefaultManifold(2)
end

function ManifoldsBase.get_embedding_type(::NotImplementedEmbeddedSubManifold)
    return ManifoldsBase.EmbeddedSubmanifoldType()
end

#
# A manifold that is isometrically embedded and `embed` is indicated as not an identity but has no implementations
#
struct NotImplementedIsometricEmbeddedManifoldNE <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.get_embedding_type(::NotImplementedIsometricEmbeddedManifoldNE)
    return ManifoldsBase.IsometricallyEmbeddedManifoldType(ManifoldsBase.DirectEmbedding())
end

#
# A manifold that is isometrically embedded and `embed` is indicated as an identity but has no implementations
#
struct NotImplementedIsometricEmbeddedManifoldDNE <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.get_embedding_type(::NotImplementedIsometricEmbeddedManifoldDNE)
    return ManifoldsBase.IsometricallyEmbeddedManifoldType(ManifoldsBase.IndirectEmbedding())
end

#
# A manifold that is an embedded manifold and `embed` is indicated as not an identity but not isometric and has no other implementation
#
struct NotImplementedEmbeddedManifoldNE <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.get_embedding_type(::NotImplementedEmbeddedManifoldNE)
    return ManifoldsBase.EmbeddedManifoldType(ManifoldsBase.DirectEmbedding())
end

#
# A manifold that is an embedded manifold and `embed` is indicated as an identity but not isometric and has no other implementation
#
struct NotImplementedEmbeddedManifoldDNE <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.get_embedding_type(::NotImplementedEmbeddedManifoldDNE)
    return ManifoldsBase.EmbeddedManifoldType(ManifoldsBase.IndirectEmbedding())
end

struct SimpleEmbeddedTestManifold <: AbstractDecoratorManifold{ℝ} end
ManifoldsBase.get_embedding(::SimpleEmbeddedTestManifold) = DefaultManifold(3)
function ManifoldsBase.get_forwarding_type(::SimpleEmbeddedTestManifold, ::Any, P::Type = Nothing)
    return ManifoldsBase.EmbeddedForwardingType(ManifoldsBase.DirectEmbedding())
end

struct EmbeddedTestManifold <: AbstractDecoratorManifold{ℝ} end
ManifoldsBase.get_embedding(::EmbeddedTestManifold) = DefaultManifold(3)
function ManifoldsBase.get_forwarding_type(::EmbeddedTestManifold, ::Any, P::Type = Nothing)
    return ManifoldsBase.EmbeddedForwardingType(ManifoldsBase.IndirectEmbedding())
end
ManifoldsBase.embed(::EmbeddedTestManifold, p) = p
ManifoldsBase.embed!(::EmbeddedTestManifold, q, p) = (q .= p)
ManifoldsBase.embed(::EmbeddedTestManifold, p, X) = X
ManifoldsBase.embed!(::EmbeddedTestManifold, Y, p, X) = (Y .= X)
ManifoldsBase.project(::EmbeddedTestManifold, p) = p
ManifoldsBase.project!(::EmbeddedTestManifold, q, p) = (q .= p)
ManifoldsBase.project(::EmbeddedTestManifold, p, X) = X
ManifoldsBase.project!(::EmbeddedTestManifold, Y, p, X) = (Y .= X)

#
# A Manifold with a fallback
#
struct FallbackManifold <: AbstractDecoratorManifold{ℝ} end
ManifoldsBase.decorated_manifold(::FallbackManifold) = DefaultManifold(3)

function ManifoldsBase.get_forwarding_type(::FallbackManifold, ::Any, P::Type = Nothing)
    return ManifoldsBase.SimpleForwardingType()
end

@testset "Embedded Manifolds" begin
    @testset "EmbeddedManifold basic tests" begin
        M = EmbeddedManifold(
            ManifoldsBase.DefaultManifold(2),
            ManifoldsBase.DefaultManifold(3),
        )
        @test repr(M) ==
            "EmbeddedManifold($(sprint(show, M.manifold)), $(sprint(show, M.embedding)))"
        @test base_manifold(M) == ManifoldsBase.DefaultManifold(2)
        @test base_manifold(M, Val(0)) == M
        @test base_manifold(M, Val(1)) == ManifoldsBase.DefaultManifold(2)
        @test base_manifold(M, Val(2)) == ManifoldsBase.DefaultManifold(2)
        @test get_embedding(M) == ManifoldsBase.DefaultManifold(3)
        @test get_embedding(M, Vector{Int}) == ManifoldsBase.DefaultManifold(3)
    end

    @testset "HalfPlaneManifold" begin
        M = HalfPlaneManifold()
        N = PosQuadrantManifold()
        @test repr(M) == "HalfPlaneManifold()"
        @test get_embedding(M) == ManifoldsBase.DefaultManifold(1, 3)
        @test representation_size(M) == (1, 3)
        # Check point checks using embedding
        @test is_point(M, [1 0.1 0.1], true)
        @test !is_point(M, [-1, 0, 0]) #wrong dim (3,1)
        @test_throws ManifoldDomainError is_point(M, [-1, 0, 0], true)
        @test !is_point(M, [-1, 0, 0])
        @test !is_point(M, [1, 0.1]) # size (from embedding)
        @test_throws ManifoldDomainError is_point(M, [1, 0.1], true)
        @test !is_point(M, [1, 0.1])
        @test is_point(M, [1 0 0], true)
        @test !is_point(M, [-1 0 0]) # right size but <0 1st
        @test_throws DomainError is_point(M, [-1 0 0], true) # right size but <0 1st
        @test !is_vector(M, [1 0 0], [1]) # right point, wrong size vector
        @test_throws ManifoldDomainError is_vector(M, [1 0 0], [1]; error = :error)
        @test !is_vector(M, [1 0 0], [1])
        @test_throws DomainError is_vector(M, [1 0 0], [-1 0 0]; error = :error) # right point, vec 1st <0
        @test !is_vector(M, [1 0 0], [-1 0 0])
        @test is_vector(M, [1 0 0], [1 0 1], true)
        @test !is_vector(M, [-1, 0, 0], [0, 0, 0])
        @test_throws ManifoldDomainError is_vector(M, [-1, 0, 0], [0, 0, 0]; error = :error)
        @test_throws ManifoldDomainError is_vector(M, [1, 0, 0], [-1, 0, 0]; error = :error)
        @test !is_vector(M, [-1, 0, 0], [0, 0, 0])
        @test !is_vector(M, [1, 0, 0], [-1, 0, 0])
        # check manifold domain error from embedding to obtain ManifoldDomainErrors
        @test !is_point(N, [0, 0, 0])
        @test_throws ManifoldDomainError is_point(N, [0, 0, 0]; error = :error)
        @test !is_vector(N, [1, 1, 0], [0, 0, 0])
        @test_throws ManifoldDomainError is_vector(N, [1, 1, 0], [0, 0, 0]; error = :error)
        p = [1.0 1.0 0.0]
        q = [1.0 0.0 0.0]
        X = q - p
        @test ManifoldsBase.check_size(M, p) === nothing
        @test ManifoldsBase.check_size(M, p, X) === nothing
        @test ManifoldsBase.check_size(M, [1, 2]) isa ManifoldDomainError
        @test ManifoldsBase.check_size(M, [1 2 3 4]) isa ManifoldDomainError
        @test ManifoldsBase.check_size(M, p, [1, 2]) isa ManifoldDomainError
        @test ManifoldsBase.check_size(M, p, [1 2 3 4]) isa ManifoldDomainError
        @test embed(M, p) == p
        pE = similar(p)
        embed!(M, pE, p)
        @test pE == p
        P = [1.0 1.0 2.0]
        Q = similar(P)
        @test project!(M, Q, P) == project!(M, Q, P)
        @test project!(M, Q, P) == embed_project!(M, Q, P)
        @test project!(M, Q, P) == [1.0 1.0 0.0]
        @test isapprox(M, p, zero_vector(M, p), [0 0 0])
        XZ = similar(X)
        zero_vector!(M, XZ, p)
        @test isapprox(M, p, XZ, [0 0 0])

        XE = similar(X)
        embed!(M, XE, p, X)
        XE2 = embed(M, p, X)
        @test X == XE
        @test XE == XE2

        @test log(M, p, q) == q - p
        Y = similar(p)
        log!(M, Y, p, q)
        @test Y == q - p
        @test exp(M, p, X) == q
        @test ManifoldsBase.exp_fused(M, p, X, 1.0) == q
        r = similar(p)
        exp!(M, r, p, X)
        @test r == q
        ManifoldsBase.exp_fused!(M, r, p, X, 1.0)
        @test r == q
        @test distance(M, p, r) == norm(r - p)

        @test retract(M, p, X) == q
        @test ManifoldsBase.retract_fused(M, p, X, 1.0) == q
        q2 = similar(q)
        @test retract!(M, q2, p, X) == q
        @test ManifoldsBase.retract_fused!(M, q2, p, X, 1.0) == q
        @test q2 == q
        @test inverse_retract(M, p, q) == X
        Y = similar(X)
        @test inverse_retract!(M, Y, p, q) == X
        @test Y == X

        @test parallel_transport_direction(M, p, X, X) == X
        @test parallel_transport_direction!(M, Y, p, X, X) == X
        @test parallel_transport_to(M, p, X, q) == X
        @test parallel_transport_to!(M, Y, p, X, q) == X

        @test get_basis(M, p, DefaultOrthonormalBasis()) isa CachedBasis
        Xc = [X[1], X[2]]
        Yc = similar(Xc)
        @test get_coordinates(M, p, X, DefaultOrthonormalBasis()) == Xc
        @test get_coordinates(M, p, X) == Xc
        @test get_coordinates!(M, Yc, p, X, DefaultOrthonormalBasis()) == Xc
        @test get_coordinates!(M, Yc, p, X) == Xc
        @test get_vector(M, p, Xc, DefaultOrthonormalBasis()) == X
        @test get_vector(M, p, Xc) == X
        @test get_vector!(M, Y, p, Xc, DefaultOrthonormalBasis()) == X
        @test get_vector!(M, Y, p, Xc) == X
    end

    @testset "AnotherHalfPlaneManifold" begin
        M = AnotherHalfPlaneManifold()
        p = [1.0, 2.0]
        pe = embed(M, p)
        @test pe == [1.0, 2.0, 0.0]
        X = [2.0, 3.0]
        Xe = embed(M, pe, X)
        @test Xe == [2.0, 3.0, 0.0]
        @test project(M, pe) == p
        @test embed_project(M, p) == p
        @test embed_project(M, p, X) == X
        Xs = similar(X)
        @test embed_project!(M, Xs, p, X) == X
        @test project(M, pe, Xe) == X
        # injectivity_radius shouldn't pass through
        @test_throws MethodError injectivity_radius(M)
        @test_throws MethodError injectivity_radius(M, p)
        @test_throws MethodError injectivity_radius(M, p, ExponentialRetraction())
        @test_throws MethodError injectivity_radius(M, ExponentialRetraction())

        @test inner(M, p, X, X) ≈ dot(X, X)
        @test norm(M, p, X) ≈ norm(X)

        q = similar(p)
        copyto!(M, q, p)
        @test p == q
        Y = similar(X)
        copyto!(M, Y, p, X)
        @test Y == X
        @test isapprox(M, p, X, Y)
        @test exp(M, p, X) ≈ p + X
        exp!(M, q, p, X)
        @test isapprox(M, q, p + X)

        @test has_components(M)

        # test vector transports in the embedding
        m = EmbeddedVectorTransport(ParallelTransport())
        q = [2.0, 2.0]
        @test vector_transport_to(M, p, X, q, m) == X #since its PT on R^3
        @test vector_transport_direction(M, p, X, q, m) == X
    end

    @testset "Test nonimplemented fallbacks" begin
        @testset "Submanifold Embedding Fallbacks & Error Tests" begin
            M = NotImplementedEmbeddedSubManifold()
            A = zeros(2)
            # for a submanifold quite a lot of functions are passed on
            @test ManifoldsBase.check_point(M, [1, 2]) === nothing
            @test ManifoldsBase.check_vector(M, [1, 2], [3, 4]) === nothing
            @test norm(M, [1, 2], [2, 3]) ≈ sqrt(13)
            @test distance(M, [1, 2], [3, 4]) ≈ sqrt(8)
            @test inner(M, [1, 2], [2, 3], [2, 3]) == 13
            @test manifold_dimension(M) == 2 # since base is defined is defined
            @test_throws MethodError project(M, [1, 2])
            @test_throws MethodError project(M, [1, 2], [2, 3]) == [2, 3]
            @test_throws MethodError project!(M, A, [1, 2], [2, 3])
            @test vector_transport_direction(M, [1, 2], [2, 3], [3, 4]) == [2, 3]
            vector_transport_direction!(M, A, [1, 2], [2, 3], [3, 4])
            @test A == [2, 3]
            @test vector_transport_to(M, [1, 2], [2, 3], [3, 4]) == [2, 3]
            vector_transport_to!(M, A, [1, 2], [2, 3], [3, 4])
            @test A == [2, 3]
            @test @inferred !isapprox(M, [1, 2], [2, 3])
            @test @inferred !isapprox(M, [1, 2], [2, 3], [4, 5])

            @test ManifoldsBase.get_forwarding_type_embedding(ManifoldsBase.EmbeddedSubmanifoldType{ManifoldsBase.DirectEmbedding}(), M, exp) === EmbeddedForwardingType()
        end
        @testset "Isometric Embedding Fallbacks & Error Tests" begin
            for M2 in [NotImplementedIsometricEmbeddedManifoldNE(), NotImplementedIsometricEmbeddedManifoldDNE()]
                @test base_manifold(M2) == M2
                A = zeros(2)
                if M2 isa NotImplementedIsometricEmbeddedManifoldDNE
                    @test_throws MethodError ManifoldsBase.allocate_result(M2, zero_vector, A)
                else
                    @test size(ManifoldsBase.allocate_result(M2, zero_vector, A)) == size(A)
                end
                # Check that all of these report not to be implemented, i.e.
                @test_throws MethodError exp(M2, [1, 2], [2, 3])
                @test_throws MethodError ManifoldsBase.exp_fused(M2, [1, 2], [2, 3], 1.0)
                @test_throws MethodError exp!(M2, A, [1, 2], [2, 3])
                @test_throws MethodError ManifoldsBase.exp_fused!(M2, A, [1, 2], [2, 3], 1.0)
                @test_throws MethodError retract(M2, [1, 2], [2, 3])
                @test_throws MethodError ManifoldsBase.retract_fused(M2, [1, 2], [2, 3], 1.0)
                @test_throws MethodError retract!(M2, A, [1, 2], [2, 3])
                @test_throws MethodError ManifoldsBase.retract_fused!(
                    M2,
                    A,
                    [1, 2],
                    [2, 3],
                    1.0,
                )
                @test_throws MethodError log(M2, [1, 2], [2, 3])
                @test_throws MethodError log!(M2, A, [1, 2], [2, 3])
                @test_throws MethodError inverse_retract(M2, [1, 2], [2, 3])
                @test_throws MethodError inverse_retract!(M2, A, [1, 2], [2, 3])
                @test_throws MethodError distance(M2, [1, 2], [2, 3])
                @test_throws MethodError project(M2, [1, 2])
                @test_throws MethodError project!(M2, A, [1, 2])
                @test_throws MethodError project(M2, [1, 2], [2, 3])
                @test_throws MethodError project!(M2, A, [1, 2], [2, 3])
                @test_throws MethodError vector_transport_direction(M2, [1, 2], [2, 3], [3, 4])
                @test_throws MethodError vector_transport_direction!(
                    M2,
                    A,
                    [1, 2],
                    [2, 3],
                    [3, 4],
                )
                @test_throws MethodError vector_transport_to(M2, [1, 2], [2, 3], [3, 4])
                @test_throws MethodError vector_transport_to!(M2, A, [1, 2], [2, 3], [3, 4])
            end
        end
        @testset "Nonisometric Embedding Fallback Error Tests" begin
            for M3 in [NotImplementedEmbeddedManifoldNE(), NotImplementedEmbeddedManifoldDNE()]
                @test_throws MethodError inner(M3, [1, 2], [2, 3], [2, 3])
                @test_throws MethodError distance(M3, [1, 2], [2, 3])
                @test_throws MethodError norm(M3, [1, 2], [2, 3])
                @test_throws MethodError zero_vector(M3, [1, 2])
            end
        end
    end
    @testset "Explicit Embeddings using EmbeddedManifold" begin
        M = DefaultManifold(3, 3)
        N = DefaultManifold(4, 4)
        O = EmbeddedManifold(M, N)
        # first test with same length of sizes
        p = ones(3, 3)
        q = zeros(4, 4)
        qT = zeros(4, 4)
        qT[1:3, 1:3] .= 1.0
        embed!(O, q, p)
        @test norm(qT - q) == 0
        qM = embed(O, p)
        @test norm(project(O, qM) - p) == 0
        @test norm(qT - qM) == 0
        # test with different sizes, check that it only fills first element
        q2 = zeros(4, 4, 3)
        q2T = zeros(4, 4, 3)
        q2T[1:3, 1:3, 1] .= 1.0
        embed!(O, q2, p)
        @test norm(q2T - q2) == 0
        O2 = EmbeddedManifold(M, DefaultManifold(4, 4, 3))
        q2M = embed(O2, p)
        @test norm(q2T - q2M) == 0
        # wrong size error checks
        @test_throws DomainError embed!(O, zeros(3, 3), zeros(3, 3, 5))
        @test_throws DomainError embed!(O, zeros(3, 3), zeros(4, 4))
        @test_throws DomainError project!(O, zeros(3, 3, 5), zeros(3, 3))
        @test_throws DomainError project!(O, zeros(4, 4), zeros(3, 3))
    end
    @testset "Explicit Fallback" begin
        M = FallbackManifold()
        # test the explicit fallback to DefaultManifold(3)
        @test inner(M, [1, 0, 0], [1, 2, 3], [0, 1, 0]) == 2
        @test is_point(M, [1, 0, 0])
        @test ManifoldsBase.allocate_result(M, exp, [1.0, 0.0, 2.0]) isa Vector
    end
    @testset "SimpleEmbeddedTestManifold" begin
        M = SimpleEmbeddedTestManifold()
        p = [1.0, 2.0, 3.0]
        X = [2.0, 3.0, 4.0]
        vf = [1.0, 2.0, 3.0, 4.0]
        is_point(M, p; error = :error)
        is_vector(M, p, X; error = :error)
        @test_throws ManifoldDomainError is_point(M, vf; error = :error)
        @test_throws ManifoldDomainError is_vector(M, p, vf; error = :error)
        @test_throws ManifoldDomainError is_vector(M, vf, X; error = :error)
        # Triggert actual test in embedding to fail
        @test_throws ManifoldDomainError is_point(M, [1.0, 2.0im, 3.0]; error = :error)
        @test_throws ManifoldDomainError is_vector(M, [1.0, 2.0im, 3.0], X; error = :error)
        @test_throws ManifoldDomainError is_vector(M, p, [1.0, 2.0im, 3.0]; error = :error)
    end
    @testset "EmbeddedTestManifold" begin
        M = EmbeddedTestManifold()
        p = [1.0, 2.0, 3.0]
        X = [2.0, 3.0, 4.0]
        vf = [1.0, 2.0, 3.0, 4.0]
        is_point(M, p; error = :error)
        is_vector(M, p, X; error = :error)
        @test_throws ManifoldDomainError is_point(M, vf; error = :error)
        @test_throws ManifoldDomainError is_vector(M, p, vf; error = :error)
        @test_throws ManifoldDomainError is_vector(M, vf, X; error = :error)
        # Triggert actual test in embedding to fail
        @test_throws ManifoldDomainError is_point(M, [1.0, 2.0im, 3.0]; error = :error)
        @test_throws ManifoldDomainError is_vector(M, [1.0, 2.0im, 3.0], X; error = :error)
        @test_throws ManifoldDomainError is_vector(M, p, [1.0, 2.0im, 3.0]; error = :error)
    end
end
