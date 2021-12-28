using ManifoldsBase, Test

using ManifoldsBase: DefaultManifold, ℝ
import ManifoldsBase:
    is_embedded_manifold, is_embedded_submanifold, is_isometric_embedded_manifold
#
# A first manifold that is modelled as a submanifold
#
struct PlaneManifold <: AbstractDecoratorManifold{ℝ} end

ManifoldsBase.get_embedding(::PlaneManifold) = ManifoldsBase.DefaultManifold(1, 3)
ManifoldsBase.base_manifold(::PlaneManifold) = ManifoldsBase.DefaultManifold(2)

ManifoldsBase.project!(::PlaneManifold, q, p) = (q .= [p[1] p[2] 0.0])
ManifoldsBase.project!(::PlaneManifold, Y, p, X) = (Y .= [X[1] X[2] 0.0])

ManifoldsBase.is_embedded_submanifold(::Type{<:PlaneManifold}) = true

function ManifoldsBase.activate_traits(::PlaneManifold, args...)
    return merge_traits(IsEmbeddedSubManifold())
end

#
# A second manifold that is modelled as just isometrically embedded but not a submanifold
#
struct AnotherPlaneManifold <: AbstractDecoratorManifold{ℝ} end

ManifoldsBase.get_embedding(::AnotherPlaneManifold) = ManifoldsBase.DefaultManifold(3)
ManifoldsBase.base_manifold(::AnotherPlaneManifold) = ManifoldsBase.DefaultManifold(2)

is_isometric_embedded_manifold(::Type{<:AnotherPlaneManifold}) = true
function ManifoldsBase.activate_traits(::AnotherPlaneManifold, args...)
    return merge_traits(IsIsometricEmbeddedManifold())
end

function ManifoldsBase.embed!(::AnotherPlaneManifold, q, p)
    q[1:2] .= p
    q[3] = 0
    return q
end
function ManifoldsBase.embed!(::AnotherPlaneManifold, Y, p, X)
    Y[1:2] .= X
    Y[3] = 0
    return Y
end
function ManifoldsBase.project!(::AnotherPlaneManifold, q, p)
    return q .= [p[1], p[2]]
end
function ManifoldsBase.project!(::AnotherPlaneManifold, Y, p, X)
    return Y .= [X[1], X[2]]
end

#
# This approach - explicitly mention an embedding.
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
struct NotImplementedEmbeddedManifold <: AbstractDecoratorManifold{ℝ} end
function ManifoldsBase.get_embedding(::NotImplementedEmbeddedManifold)
    return ManifoldsBase.DefaultManifold(3)
end
function ManifoldsBase.base_manifold(::NotImplementedEmbeddedManifold)
    return ManifoldsBase.DefaultManifold(2)
end
is_embedded_submanifold(::Type{<:NotImplementedEmbeddedManifold}) = true

#
# A manifold that is isometrically embedded but has no implementations
#
struct NotImplementedEmbeddedManifold2 <: AbstractDecoratorManifold{ℝ} end
is_isometric_embedded_manifold(::Type{<:NotImplementedEmbeddedManifold2}) = true

#
# A manifold that is an embedded manifold but not isometric and has no other implementation
#
struct NotImplementedEmbeddedManifold3 <: AbstractDecoratorManifold{ℝ} end
is_embedded_manifold(::Type{<:NotImplementedEmbeddedManifold3}) = true

@testset "Embedded Manifolds" begin
    @testset "EmbeddedManifold basic tests" begin
        M = EmbeddedManifold(
            ManifoldsBase.DefaultManifold(2),
            ManifoldsBase.DefaultManifold(3),
        )
        @test repr(M) ==
              "EmbeddedManifold($(sprint(show, M.manifold)), $(sprint(show, M.embedding)))"
        @test base_manifold(M) == ManifoldsBase.DefaultManifold(2)
        @test get_embedding(M) == ManifoldsBase.DefaultManifold(3)
    end

    @testset "PlaneManifold" begin
        M = PlaneManifold()
        @test repr(M) == "PlaneManifold()"
        @test get_embedding(M) == ManifoldsBase.DefaultManifold(1, 3)
        # Check fallbacks to check embed->check_manifoldpoint Defaults
        @test_throws DomainError is_point(M, [1, 0, 0], true)
        @test_throws DomainError is_point(M, [1 0], true)
        @test is_point(M, [1 0 0], true)
        @test_throws DomainError is_vector(M, [1 0 0], [1], true)
        @test_throws DomainError is_vector(M, [1 0 0], [0 0 0 0], true)
        @test is_vector(M, [1 0 0], [1 0 1], true)
        p = [1.0 1.0 0.0]
        q = [1.0 0.0 0.0]
        X = q - p
        @test check_size(M, p) === nothing
        @test check_size(M, p, X) === nothing
        @test check_size(M, [1, 2]) isa DomainError
        @test check_size(M, [1 2 3 4]) isa DomainError
        @test check_size(M, p, [1, 2]) isa DomainError
        @test check_size(M, p, [1 2 3 4]) isa DomainError
        @test embed(M, p) == p
        pE = similar(p)
        embed!(M, pE, p)
        @test pE == p
        P = [1.0 1.0 2.0]
        Q = similar(P)
        @test project!(M, Q, P) == project!(M, Q, P)
        @test project!(M, Q, P) == [1.0 1.0 0.0]

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
    end

    @testset "AnotherPlaneManifold" begin
        M = AnotherPlaneManifold()
        p = [1.0, 2.0]
        pe = embed(M, p)
        @test pe == [1.0, 2.0, 0.0]
        X = [2.0, 3.0]
        Xe = embed(M, pe, X)
        @test Xe == [2.0, 3.0, 0.0]
        @test project(M, pe) == p
        @test project(M, pe, Xe) == X
    end

    @testset "Test nonimplemented fallbacks" begin
        @testset "Default Isometric Embedding Fallback Error Tests" begin
            M = NotImplementedEmbeddedManifold()
            A = zeros(2)
            # without any extra tests just the embedding is asked
            @test check_point(M, [1, 2]) === nothing
            @test check_vector(M, [1, 2], [3, 4]) === nothing
            @test norm(M, [1, 2], [2, 3]) ≈ sqrt(13)
            @test distance(M, [1, 2], [3, 4]) ≈ sqrt(8)
            @test inner(M, [1, 2], [2, 3], [2, 3]) ≈ 13
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
        @testset "General Isometric Embedding Fallback Error Tests" begin
            M2 = NotImplementedEmbeddedManifold2()
            @test base_manifold(M2) == M2
            A = zeros(2)
            # Check that all of these report not to be implemented, i.e.
            @test_throws MethodError exp(M2, [1, 2], [2, 3])
            @test_throws MethodError exp!(M2, A, [1, 2], [2, 3])
            @test_throws MethodError log(M2, [1, 2], [2, 3])
            @test_throws MethodError log!(M2, A, [1, 2], [2, 3])
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
            @test_throws MethodError vector_transport_along!(M2, A, [1, 2], [2, 3], [])
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
            M3 = NotImplementedEmbeddedManifold3()
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
end
