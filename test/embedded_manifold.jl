struct PlaneManifold <: AbstractEmbeddedManifold{ℝ,TransparentIsometricEmbedding} end

ManifoldsBase.decorated_manifold(::PlaneManifold) = ManifoldsBase.DefaultManifold(3)
ManifoldsBase.base_manifold(::PlaneManifold) = ManifoldsBase.DefaultManifold(2)

ManifoldsBase.embed!(::PlaneManifold, q, p) = copyto!(q, p)
ManifoldsBase.embed!(::PlaneManifold, Y, p, X) = copyto!(Y, X)
ManifoldsBase.project!(::PlaneManifold, q, p) = (q .= [p[1], p[2], 0.0])
ManifoldsBase.project!(::PlaneManifold, Y, p, X) = (Y .= [X[1], X[2], 0.0])

struct AnotherPlaneManifold <: AbstractEmbeddedManifold{ℝ,DefaultIsometricEmbeddingType} end

ManifoldsBase.decorated_manifold(::AnotherPlaneManifold) = ManifoldsBase.DefaultManifold(3)
ManifoldsBase.base_manifold(::AnotherPlaneManifold) = ManifoldsBase.DefaultManifold(2)

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
    q .= [p[1], p[2]]
end
function ManifoldsBase.project!(::AnotherPlaneManifold, Y, p, X)
    Y .= [X[1], X[2]]
end

struct NotImplementedEmbeddedManifold <:
       AbstractEmbeddedManifold{ℝ,TransparentIsometricEmbedding} end
function ManifoldsBase.decorated_manifold(::NotImplementedEmbeddedManifold)
    return ManifoldsBase.DefaultManifold(2)
end
function ManifoldsBase.base_manifold(::NotImplementedEmbeddedManifold)
    return ManifoldsBase.DefaultManifold(2)
end

struct NotImplementedEmbeddedManifold2 <:
       AbstractEmbeddedManifold{ℝ,DefaultIsometricEmbeddingType} end

struct NotImplementedEmbeddedManifold3 <: AbstractEmbeddedManifold{ℝ,DefaultEmbeddingType} end

@testset "Embedded Manifolds" begin
    @testset "EmbeddedManifold basic tests" begin
        M = EmbeddedManifold(
            ManifoldsBase.DefaultManifold(2),
            ManifoldsBase.DefaultManifold(3),
        )
        @test repr(M) ==
              "EmbeddedManifold(DefaultManifold{Tuple{2},ℝ}(), DefaultManifold{Tuple{3},ℝ}(), TransparentIsometricEmbedding())"
        @test base_manifold(M) == ManifoldsBase.DefaultManifold(2)
        @test ManifoldsBase.decorated_manifold(M) == ManifoldsBase.DefaultManifold(3)
        @test ManifoldsBase.default_embedding_dispatch(M) === Val{false}()
        @test ManifoldsBase.default_decorator_dispatch(M) ===
              ManifoldsBase.default_embedding_dispatch(M)
    end
    @testset "PlaneManifold" begin
        M = PlaneManifold()
        @test repr(M) == "PlaneManifold()"
        @test ManifoldsBase.default_decorator_dispatch(M) === Val{false}()
        @test get_embedding(M) == ManifoldsBase.DefaultManifold(3)
        # Check fallbacks to check embed->check_manifoldpoint Defaults
        @test_throws DomainError is_manifold_point(M, [1, 0], true)
        @test_throws DomainError is_tangent_vector(M, [1, 0], [1, 0, 0], true)
        @test_throws DomainError is_tangent_vector(M, [1, 0, 0], [0, 0], true)
        p = [1.0, 1.0, 0.0]
        q = [1.0, 0.0, 0.0]
        X = q - p
        @test embed(M, p) == p
        pE = similar(p)
        embed!(M, pE, p)
        @test pE == p
        P = [1.0, 1.0, 2.0]
        Q = similar(P)
        @test project!(M, Q, P) == project!(M, Q, P)
        @test project!(M, Q, P) == [1.0, 1.0, 0.0]

        @test log(M, p, q) == q - p
        Y = similar(p)
        log!(M, Y, p, q)
        @test Y == q - p
        @test exp(M, p, X) == q
        r = similar(p)
        exp!(M, r, p, X)
        @test r == q
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
            @test norm(M, [1, 2], [2, 3]) ≈ sqrt(13)
            @test inner(M, [1, 2], [2, 3], [2, 3]) ≈ 13
            @test_throws ErrorException manifold_dimension(M)
            # without any implementation the projections are the identity
            @test project(M, [1, 2]) == [1, 2]
            @test project(M, [1, 2], [2, 3]) == [2, 3]
            project!(M, A, [1, 2], [2, 3])
            @test A == [2, 3]
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
            @test_throws ErrorException exp(M2, [1, 2], [2, 3])
            @test_throws ErrorException exp!(M2, A, [1, 2], [2, 3])
            @test_throws ErrorException log(M2, [1, 2], [2, 3])
            @test_throws ErrorException log!(M2, A, [1, 2], [2, 3])
            @test_throws ErrorException manifold_dimension(M2)
            @test_throws ErrorException project(M2, [1, 2])
            @test_throws ErrorException project!(M2, A, [1, 2])
            @test_throws ErrorException project(M2, [1, 2], [2, 3])
            @test_throws ErrorException project!(M2, A, [1, 2], [2, 3])
            @test_throws ErrorException vector_transport_along(M2, [1, 2], [2, 3], [[1, 2]])
            @test_throws ErrorException vector_transport_along(
                M2,
                [1, 2],
                [2, 3],
                [[1, 2]],
                ParallelTransport(),
            )
            @test_throws ErrorException vector_transport_along!(M2, A, [1, 2], [2, 3], [])
            @test_throws ErrorException vector_transport_direction(
                M2,
                [1, 2],
                [2, 3],
                [3, 4],
            )
            @test_throws ErrorException vector_transport_direction!(
                M2,
                A,
                [1, 2],
                [2, 3],
                [3, 4],
            )
            @test_throws ErrorException vector_transport_to(M2, [1, 2], [2, 3], [3, 4])
            @test_throws ErrorException vector_transport_to!(M2, A, [1, 2], [2, 3], [3, 4])
        end
        @testset "Nonisometric Embedding Fallback Error Rests" begin
            M3 = NotImplementedEmbeddedManifold3()
            @test_throws ErrorException inner(M3, [1, 2], [2, 3], [2, 3])
            @test_throws ErrorException manifold_dimension(M3)
            @test_throws ErrorException norm(M3, [1, 2], [2, 3])
            @test_throws ErrorException embed(M3, [1, 2], [2, 3])
            @test_throws ErrorException embed(M3, [1, 2])
        end
    end
    @testset "EmbeddedManifold decorator dispatch" begin
        TM = NotImplementedEmbeddedManifold() # transparently iso
        IM = NotImplementedEmbeddedManifold2() # iso
        AM = NotImplementedEmbeddedManifold3() # general
        for f in [embed, exp, get_basis, get_coordinates, get_vector, inverse_retract, log]
            @test ManifoldsBase.decorator_transparent_dispatch(f, AM) === Val(:parent)
        end
        for f in
            [project, retract, inverse_retract!, retract!, get_coordinates!, get_vector!]
            @test ManifoldsBase.decorator_transparent_dispatch(f, AM) === Val(:parent)
        end
        for f in [vector_transport_along, vector_transport_direction, vector_transport_to]
            @test ManifoldsBase.decorator_transparent_dispatch(f, AM) === Val(:parent)
        end
        for f in [mid_point, mid_point!]
            @test ManifoldsBase.decorator_transparent_dispatch(f, AM) === Val(:parent)
        end
        for f in [check_manifold_point, check_tangent_vector, exp!, inner, embed!]
            @test ManifoldsBase.decorator_transparent_dispatch(f, AM) ===
                  Val(:intransparent)
        end
        for f in [log!, norm, manifold_dimension, project!]
            @test ManifoldsBase.decorator_transparent_dispatch(f, AM) ===
                  Val(:intransparent)
        end
        @test ManifoldsBase.decorator_transparent_dispatch(vector_transport_along!, AM) ===
              Val(:intransparent)
        @test ManifoldsBase.decorator_transparent_dispatch(
            vector_transport_to!,
            AM,
            1,
            1,
            1,
            1,
            1,
        ) === Val(:intransparent)
        @test ManifoldsBase.decorator_transparent_dispatch(
            vector_transport_direction!,
            AM,
        ) === Val(:parent)

        for f in [inner, norm]
            @test ManifoldsBase.decorator_transparent_dispatch(f, IM) === Val(:transparent)
        end
        for f in [inverse_retract!, retract!, mid_point!]
            @test ManifoldsBase.decorator_transparent_dispatch(f, IM) === Val(:parent)
        end

        for f in [exp, inverse_retract, log, project, retract, mid_point]
            @test ManifoldsBase.decorator_transparent_dispatch(f, TM) === Val(:transparent)
        end
        for f in [exp!, inverse_retract!, log!, project!, retract!, mid_point!]
            @test ManifoldsBase.decorator_transparent_dispatch(f, TM) === Val(:transparent)
        end
        for f in [vector_transport_along, vector_transport_direction, vector_transport_to]
            @test ManifoldsBase.decorator_transparent_dispatch(f, TM) === Val(:transparent)
        end
        for f in
            [vector_transport_along!, vector_transport_direction!, vector_transport_to!]
            @test ManifoldsBase.decorator_transparent_dispatch(f, TM) === Val(:transparent)
        end
        for t in [PoleLadderTransport(), SchildsLadderTransport()], M in [AM, TM, IM]
            @test ManifoldsBase.decorator_transparent_dispatch(
                vector_transport_to!,
                M,
                1,
                1,
                1,
                1,
                PoleLadderTransport(),
            ) === Val(:parent)
        end
        @test ManifoldsBase.decorator_transparent_dispatch(embed, TM) === Val{:parent}()
        @test ManifoldsBase.decorator_transparent_dispatch(embed!, TM) ===
              Val(:intransparent)
    end
end
