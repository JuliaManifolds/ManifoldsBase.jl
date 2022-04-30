using ManifoldsBase, Test

using ManifoldsBase: DefaultManifold, ℝ
#
# A first artificial (not real) manifold that is modelled as a submanifold
# half plane with euclidean metric is not a manifold but should test all things correctly here
#
struct HalfPlanemanifold <: AbstractDecoratorManifold{ℝ} end

ManifoldsBase.get_embedding(::HalfPlanemanifold) = ManifoldsBase.DefaultManifold(1, 3)
ManifoldsBase.decorated_manifold(::HalfPlanemanifold) = ManifoldsBase.DefaultManifold(2)
ManifoldsBase.representation_size(::HalfPlanemanifold) = (2,)

function ManifoldsBase.check_point(::HalfPlanemanifold, p)
    return p[1] > 0 ? nothing : DomainError(p[1], "p[1] ≤ 0")
end
function ManifoldsBase.check_vector(::HalfPlanemanifold, p, X)
    return X[1] > 0 ? nothing : DomainError(X[1], "X[1] ≤ 0")
end
ManifoldsBase.embed(::HalfPlanemanifold, p) = reshape(p, 1, :)
ManifoldsBase.embed(::HalfPlanemanifold, p, X) = reshape(X, 1, :)

ManifoldsBase.project!(::HalfPlanemanifold, q, p) = (q .= [p[1] p[2] 0.0])
ManifoldsBase.project!(::HalfPlanemanifold, Y, p, X) = (Y .= [X[1] X[2] 0.0])

function ManifoldsBase.get_coordinates_orthonormal!(
    ::HalfPlanemanifold,
    Y,
    p,
    X,
    ::ManifoldsBase.RealNumbers,
)
    return (Y .= [X[1], X[2]])
end
function ManifoldsBase.get_vector_orthonormal!(
    ::HalfPlanemanifold,
    Y,
    p,
    c,
    ::ManifoldsBase.RealNumbers,
)
    return (Y .= [c[1] c[2] 0.0])
end

function ManifoldsBase.active_traits(f, ::HalfPlanemanifold, args...)
    return ManifoldsBase.merge_traits(ManifoldsBase.IsEmbeddedSubmanifold())
end

#
# A second manifold that is modelled as just isometrically embedded but not a submanifold
#
struct AnotherHalfPlanemanifold <: AbstractDecoratorManifold{ℝ} end

ManifoldsBase.get_embedding(::AnotherHalfPlanemanifold) = ManifoldsBase.DefaultManifold(3)
function ManifoldsBase.decorated_manifold(::AnotherHalfPlanemanifold)
    return ManifoldsBase.DefaultManifold(2)
end
ManifoldsBase.representation_size(::AnotherHalfPlanemanifold) = (2,)

function ManifoldsBase.active_traits(f, ::AnotherHalfPlanemanifold, args...)
    return ManifoldsBase.merge_traits(ManifoldsBase.IsIsometricEmbeddedManifold())
end

function ManifoldsBase.embed!(::AnotherHalfPlanemanifold, q, p)
    q[1:2] .= p
    q[3] = 0
    return q
end
function ManifoldsBase.embed!(::AnotherHalfPlanemanifold, Y, p, X)
    Y[1:2] .= X
    Y[3] = 0
    return Y
end
function ManifoldsBase.project!(::AnotherHalfPlanemanifold, q, p)
    return q .= [p[1], p[2]]
end
function ManifoldsBase.project!(::AnotherHalfPlanemanifold, Y, p, X)
    return Y .= [X[1], X[2]]
end

#
# Third example - explicitly mention an embedding.
#
function ManifoldsBase.embed!(
    ::EmbeddedManifold{𝔽,DefaultManifold{nL,𝔽},DefaultManifold{mL,𝔽2}},
    q,
    p,
) where {nL,mL,𝔽,𝔽2}
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
    ::EmbeddedManifold{𝔽,DefaultManifold{nL,𝔽},DefaultManifold{mL,𝔽2}},
    q,
    p,
) where {nL,mL,𝔽,𝔽2}
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
function ManifoldsBase.active_traits(f, ::NotImplementedEmbeddedSubManifold, args...)
    return ManifoldsBase.merge_traits(ManifoldsBase.IsEmbeddedSubmanifold())
end

#
# A manifold that is isometrically embedded but has no implementations
#
struct NotImplementedIsometricEmbeddedManifold <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.active_traits(f, ::NotImplementedIsometricEmbeddedManifold, args...)
    return ManifoldsBase.merge_traits(ManifoldsBase.IsIsometricEmbeddedManifold())
end

#
# A manifold that is an embedded manifold but not isometric and has no other implementation
#
struct NotImplementedEmbeddedManifold <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.active_traits(f, ::NotImplementedEmbeddedManifold, args...)
    return ManifoldsBase.merge_traits(ManifoldsBase.IsEmbeddedManifold())
end

#
# A Manifold with a fallback
#
struct FallbackManifold <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.active_traits(f, ::FallbackManifold, args...)
    return ManifoldsBase.merge_traits(ManifoldsBase.IsExplicitDecorator())
end
ManifoldsBase.decorated_manifold(::FallbackManifold) = DefaultManifold(3)

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
        @test get_embedding(M, [1, 2, 3]) == ManifoldsBase.DefaultManifold(3)
    end

    @testset "HalfPlanemanifold" begin
        M = HalfPlanemanifold()
        @test repr(M) == "HalfPlanemanifold()"
        @test get_embedding(M) == ManifoldsBase.DefaultManifold(1, 3)
        @test representation_size(M) == (2,)
        # Check point checks using embedding
        @test is_point(M, [1 0.1 0.1], true)
        @test_throws DomainError is_point(M, [-1, 0, 0], true) #wrong dim (3,1)
        @test !is_point(M, [-1, 0, 0])
        @test_throws ManifoldDomainError is_point(M, [1, 0.1], true) # size
        @test is_point(M, [1 0 0], true)
        @test !is_point(M, [-1 0 0]) # right size but <0 1st
        @test_throws DomainError is_point(M, [-1 0 0], true) # right size but <0 1st
        @test_throws ManifoldDomainError is_vector(M, [1 0 0], [1], true) # right point, wrong size vector
        @test !is_vector(M, [1 0 0], [1])
        @test_throws DomainError is_vector(M, [1 0 0], [-1 0 0], true) # right point, vec 1st <0
        @test !is_vector(M, [1 0 0], [-1 0 0])
        @test is_vector(M, [1 0 0], [1 0 1], true)
        @test_throws DomainError is_vector(M, [-1, 0, 0], [0, 0, 0], true)
        @test_throws DomainError is_vector(M, [1, 0, 0], [-1, 0, 0], true)
        @test !is_vector(M, [-1, 0, 0], [0, 0, 0])
        @test !is_vector(M, [1, 0, 0], [-1, 0, 0])
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
        r = similar(p)
        exp!(M, r, p, X)
        @test r == q
        @test distance(M, p, r) == norm(r - p)

        @test retract(M, p, X) == q
        q2 = similar(q)
        @test retract!(M, q2, p, X) == q
        @test q2 == q
        @test inverse_retract(M, p, q) == X
        Y = similar(X)
        @test inverse_retract!(M, Y, p, q) == X
        @test Y == X

        @test vector_transport_along(M, p, X, []) == X
        @test vector_transport_along!(M, Y, p, X, []) == X
        @test parallel_transport_along(M, p, X, []) == X
        @test parallel_transport_along!(M, Y, p, X, []) == X
        @test parallel_transport_direction(M, p, X, X) == X
        @test parallel_transport_direction!(M, Y, p, X, X) == X
        @test parallel_transport_to(M, p, X, q) == X
        @test parallel_transport_to!(M, Y, p, X, q) == X

        @test get_basis(M, p, DefaultOrthonormalBasis()) isa CachedBasis
        Xc = [X[1], X[2]]
        Yc = similar(Xc)
        @test get_coordinates(M, p, X, DefaultOrthonormalBasis()) == Xc
        @test get_coordinates!(M, Yc, p, X, DefaultOrthonormalBasis()) == Xc
        @test get_vector(M, p, Xc, DefaultOrthonormalBasis()) == X
        @test get_vector!(M, Y, p, Xc, DefaultOrthonormalBasis()) == X
    end

    @testset "AnotherHalfPlanemanifold" begin
        M = AnotherHalfPlanemanifold()
        p = [1.0, 2.0]
        pe = embed(M, p)
        @test pe == [1.0, 2.0, 0.0]
        X = [2.0, 3.0]
        Xe = embed(M, pe, X)
        @test Xe == [2.0, 3.0, 0.0]
        @test project(M, pe) == p
        @test project(M, pe, Xe) == X
        # isometric passtthrough
        @test injectivity_radius(M) == Inf
        @test injectivity_radius(M, p) == Inf
        @test injectivity_radius(M, p, ExponentialRetraction()) == Inf
        @test injectivity_radius(M, ExponentialRetraction()) == Inf
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
        end
        @testset "Isometric Embedding Fallbacks & Error Tests" begin
            M2 = NotImplementedIsometricEmbeddedManifold()
            @test base_manifold(M2) == M2
            A = zeros(2)
            # Check that all of these report not to be implemented, i.e.
            @test_throws MethodError exp(M2, [1, 2], [2, 3])
            @test_throws MethodError exp!(M2, A, [1, 2], [2, 3])
            @test_throws MethodError retract(M2, [1, 2], [2, 3])
            @test_throws MethodError retract!(M2, A, [1, 2], [2, 3])
            @test_throws MethodError log(M2, [1, 2], [2, 3])
            @test_throws MethodError log!(M2, A, [1, 2], [2, 3])
            @test_throws MethodError inverse_retract(M2, [1, 2], [2, 3])
            @test_throws MethodError inverse_retract!(M2, A, [1, 2], [2, 3])
            @test_throws MethodError distance(M2, [1, 2], [2, 3])
            @test_throws StackOverflowError manifold_dimension(M2)
            @test_throws MethodError project(M2, [1, 2])
            @test_throws MethodError project!(M2, A, [1, 2])
            @test_throws MethodError project(M2, [1, 2], [2, 3])
            @test_throws MethodError project!(M2, A, [1, 2], [2, 3])
            @test_throws MethodError vector_transport_along(M2, [1, 2], [2, 3], [[1, 2]])
            @test_throws MethodError vector_transport_along(
                M2,
                [1, 2],
                [2, 3],
                [[1, 2]],
                ParallelTransport(),
            )
            @test vector_transport_along!(M2, A, [1, 2], [2, 3], []) == [2, 3]
            @test A == [2, 3]
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
        @testset "Nonisometric Embedding Fallback Error Rests" begin
            M3 = NotImplementedEmbeddedManifold()
            @test_throws MethodError inner(M3, [1, 2], [2, 3], [2, 3])
            @test_throws StackOverflowError manifold_dimension(M3)
            @test_throws MethodError distance(M3, [1, 2], [2, 3])
            @test_throws MethodError norm(M3, [1, 2], [2, 3])
            @test_throws MethodError embed(M3, [1, 2], [2, 3])
            @test_throws MethodError embed(M3, [1, 2])
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
    end
end
