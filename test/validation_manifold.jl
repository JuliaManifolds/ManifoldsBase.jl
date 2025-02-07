using ManifoldsBase, LinearAlgebra, Random, Test

s = joinpath(@__DIR__, "ManifoldsBaseTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using ManifoldsBaseTestSuite

@testset "Validation manifold" begin
    M = ManifoldsBase.DefaultManifold(3)
    A = ValidationManifold(M)
    @test A == ValidationManifold(M, A) # A equals the VM with same defaults
    x = [1.0, 0.0, 0.0]
    y = 1 / sqrt(2) * [1.0, 1.0, 0.0]
    z = [0.0, 1.0, 0.0]
    v = log(M, x, y)
    x2 = ValidationMPoint(x)
    y2 = ValidationMPoint(y)
    v2 = log(A, x, y) # auto convert
    y2 = exp(A, x, v2)
    w = log(M, x, z)
    w2 = log(A, x, z; atol = 10^(-15))
    @testset "Checks forward" begin
        @test ManifoldsBase.check_size(A, x2, v2) === ManifoldsBase.check_size(M, x, v)
        @test ManifoldsBase.check_size(A, x2) === ManifoldsBase.check_size(M, x)
    end
    @testset "is_point / is_vector error." begin
        @test is_point(A, x; error = :error)
        @test_throws DomainError is_point(A, [1, 2, 3, 4]; error = :error)
        @test is_vector(A, x, v; error = :error)
        @test_throws DomainError is_vector(A, x, [1, 2, 3, 4]; error = :error)
        # Test that we can ignore point contexts
        A2a = ValidationManifold(M; ignore_contexts = [:Point])
        @test is_point(A2a, [1, 2, 3, 4])
        @test_throws DomainError !is_vector(A2a, x, [1, 2, 3, 4])
        A2b = ValidationManifold(M; ignore_contexts = [:Vector])
        @test_throws DomainError !is_point(A2b, [1, 2, 3, 4])
        @test is_vector(A2b, x, [1, 2, 3, 4])
        A3a = ValidationManifold(M; ignore_functions = Dict(exp => :All))
        @test is_point(A3a, [1, 2, 3, 4]; within = exp)
        @test_throws DomainError is_point(A3a, [1, 2, 3, 4]; within = log)
        A3b = ValidationManifold(M; ignore_functions = Dict(exp => [:Point, :Vector]))
        @test is_point(A3b, [1, 2, 3, 4]; within = exp)
        @test_throws DomainError is_point(A3b, [1, 2, 3, 4]; within = log)
        @test is_vector(A3b, x, [1, 2, 3, 4]; within = exp)
        @test_throws DomainError is_vector(A3b, x, [1, 2, 3, 4]; within = log)
    end
    @testset "Types and Conversion" begin
        @test convert(typeof(M), A) == M
        # equality does not work, since we have mutable fields, but the type agrees
        @test typeof(convert(typeof(A), M)) == typeof(A)
        @test base_manifold(A) == M
        @test base_manifold(base_manifold(A)) == base_manifold(A)
        @test ManifoldsBase.representation_size(A) == ManifoldsBase.representation_size(M)
        @test ManifoldsBase.representation_size(A) == (3,)
        @test manifold_dimension(A) == manifold_dimension(M)
        @test manifold_dimension(A) == 3
        for T in [ValidationMPoint, ValidationTangentVector, ValidationCotangentVector]
            p = T(x)
            @test convert(typeof(x), p) == x
            @test convert(typeof(p), y) == T(y)
            @test number_eltype(typeof(p)) == eltype(x)
            @test number_eltype(p) == eltype(x)
            @test typeof(allocate(p)) == typeof(p)
            @test typeof(allocate(p, eltype(x))) == typeof(p)
            if T === ValidationMPoint
                @test typeof(allocate(p, eltype(x), (3, 1))) == T{Matrix{Float64}}
            else
                @test typeof(allocate(p, eltype(x), (3, 1))) == T{Matrix{Float64},Nothing}
            end
            @test allocate(p) isa T
            @test allocate(p, Float32) isa T
            @test number_eltype(allocate(p, Float32)) == Float32
            @test similar(p) isa T
            @test similar(p, Float32) isa T
            @test number_eltype(similar(p, Float32)) == Float32
            q = allocate(p)
            if T == ValidationMPoint
                copyto!(A, q, p)
            else
                copyto!(A, q, ValidationMPoint(x), p) # generate base point “on the fly”
            end
            @test isapprox(A, q, p)
            @test ManifoldsBase.internal_value(p) == x
            @test ManifoldsBase.internal_value(x) == x
            @test copy(p) == p
            q = allocate(p)
            copyto!(q, p)
            @test q == p
            @test isapprox(A, q, p)
        end
    end
    @testset "Vector functions" begin
        for T in [ValidationTangentVector, ValidationCotangentVector]
            a = T(v)
            b = T(w)
            @test isapprox(A, a + b, T(v + w))
            @test isapprox(A, (a - b), T(v - w))
            @test isapprox(A, -b, T(-w))
            @test isapprox(A, +b, T(+w))
            @test isapprox(A, 2 * a, T(2 .* v))
            @test isapprox(A, a * 2, T(v .* 2))
            @test isapprox(A, 2 \ a, T(2 .\ v))
            @test isapprox(A, a / 2, T(v ./ 2))
            @test zero(a) == T(zero(v))
            @test isapprox(A, 2 .* a .+ b, T(2 .* v .+ w))
            c = similar(a)
            c .= a .+ b
            @test isapprox(A, c, a .+ b)
        end
    end
    @testset "AbstractManifold functions" begin
        @test manifold_dimension(A) == manifold_dimension(M)
        @test isapprox(y2.value, y)
        @test distance(A, x, y) == distance(M, x, y)
        @test norm(A, x, v) == norm(M, x, v)
        @test inner(A, x, v2, w2; atol = 10^(-15)) == inner(M, x, v, w)
        @test isapprox(A, x2, y2) == isapprox(M, x, y)
        @test isapprox(A, x, y) == isapprox(A, x2, y2)
        @test isapprox(A, x, v2, v2) == isapprox(M, x, v, v)
        v2s = similar(v2)
        project!(A, v2s, x2, v2)
        @test isapprox(A, v2, v2s)
        y2s = similar(y2)
        exp!(A, y2s, x2, v2)
        @test isapprox(A, y2s, y2)
        log!(A, v2s, x, y)
        @test mid_point(A, x, y) == 0.5 * (x .+ y)
        mp = similar(x)
        mid_point!(A, mp, x, y)
        @test mp == 0.5 * (x .+ y)
        @test isapprox(A, x, v2s, v2)
        @test isapprox(A, exp(A, x, v), y2)
        @test isapprox(A, zero_vector(A, x), zero_vector(M, x))
        vector_transport_to!(A, v2s, x2, v2, y2)
        @test isapprox(A, x2, v2, v2s)
        @test isapprox(A, x2, v2, vector_transport_to(A, x2, v2, y2))
        slt = ManifoldsBase.SchildsLadderTransport()
        vector_transport_to!(A, v2s, x2, v2, y2, slt)
        @test isapprox(A, x2, v2, v2s)
        @test isapprox(A, x2, v2, vector_transport_to(A, x2, v2, y2, slt))
        plt = ManifoldsBase.PoleLadderTransport()
        vector_transport_to!(A, v2s, x2, v2, y, plt)
        @test isapprox(A, x2, v2, vector_transport_to(A, x2, v2, y2, plt))
        @test isapprox(A, x2, v2, v2s)
        pt = ManifoldsBase.ProjectionTransport()
        vector_transport_to!(A, v2s, x2, v2, y2, pt)
        @test isapprox(A, x2, v2, v2s)
        @test isapprox(A, x2, v2, vector_transport_to(A, x2, v2, y2, pt))
        zero_vector!(A, v2s, x)
        @test isapprox(A, x, v2s, zero_vector(M, x))
        c2 = [x2]
        v3 = similar(v2)
        @test injectivity_radius(A) == Inf
        @test injectivity_radius(A, x) == Inf
        @test injectivity_radius(A, ManifoldsBase.ExponentialRetraction()) == Inf
        @test injectivity_radius(A, x, ManifoldsBase.ExponentialRetraction()) == Inf
        @test injectivity_radius(
            A,
            ManifoldsBaseTestSuite.CustomValidationManifoldRetraction(),
        ) == 10
        @test injectivity_radius(
            A,
            x,
            ManifoldsBaseTestSuite.CustomValidationManifoldRetraction(),
        ) == 11
    end
    @testset "ValidationManifold basis" begin
        b = [Matrix(I, 3, 3)[:, i] for i in 1:3]
        for BT in (DefaultBasis, DefaultOrthonormalBasis, DefaultOrthogonalBasis)
            @testset "Basis $(BT)" begin
                cb = BT()
                @test b == get_vectors(M, x, get_basis(A, x, cb))
                v = similar(x)
                v2 = similar(x)
                @test_throws ArgumentError get_vector(A, x, [1.0], cb)
                @test_throws DomainError get_vector(A, [1.0], [1.0, 0.0, 0.0], cb)
                @test_throws ArgumentError get_vector!(A, v, x, [], cb)
                @test_throws DomainError get_vector!(A, v, [1.0], [1.0, 0.0, 0.0], cb)
                @test_throws DomainError get_coordinates(A, x, [1.0], cb)
                @test_throws DomainError get_coordinates!(A, v, x, [], cb)
                @test_throws DomainError get_coordinates!(A, v, [1.0], [1.0, 0.0, 0.0], cb)
                @test get_vector(A, x, [1, 2, 3], cb) ≈ get_vector(M, x, [1, 2, 3], cb)
                @test get_vector!(A, v2, x, [1, 2, 3], cb) ≈
                      get_vector!(M, v, x, [1, 2, 3], cb)
                @test get_coordinates(A, x, [1, 2, 3], cb) ≈
                      get_coordinates(M, x, [1, 2, 3], cb)
                @test get_coordinates!(A, v2, x, [1, 2, 3], cb) ≈
                      get_coordinates!(M, v, x, [1, 2, 3], cb)

                @test_throws ErrorException get_basis(A, x, CachedBasis(cb, [x]))
                @test_throws ErrorException get_basis(A, x, CachedBasis(cb, [x, x, x]))
                @test_throws ErrorException get_basis(A, x, CachedBasis(cb, [2 * x, x, x]))
                if BT <: ManifoldsBase.AbstractOrthonormalBasis
                    @test_throws ArgumentError get_basis(
                        A,
                        x,
                        CachedBasis(
                            cb,
                            [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        ),
                    )
                elseif BT <: ManifoldsBase.AbstractOrthogonalBasis
                    @test_throws ArgumentError get_basis(
                        A,
                        x,
                        CachedBasis(
                            cb,
                            [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        ),
                    )
                end
            end
        end
    end
    @testset "Output distance an norm – on different message types" begin
        Dm = ManifoldsBaseTestSuite.ValidationDummyManifold()
        Ad = ValidationManifold(Dm)
        @test_throws DomainError distance(Ad, [], [])
        @test_throws DomainError norm(Ad, [], [])
        AdI = ValidationManifold(Dm; error = :info)
        @test_logs (:info,) distance(AdI, [], [])
        @test_logs (:info,) norm(AdI, [], [])
        AdW = ValidationManifold(Dm; error = :warn)
        @test_logs (:warn,) distance(AdW, [], [])
        @test_logs (:warn,) norm(AdW, [], [])
        AdO = ValidationManifold(Dm; ignore_contexts = [:Output])
        @test distance(AdO, [], []) == -1.0
        @test norm(AdO, [], []) == -1.0
        AdN = ValidationManifold(Dm; error = :None)
        @test distance(AdN, [], []) == -1.0
        @test norm(AdN, [], []) == -1.0
    end
    @testset "rand" begin
        Random.seed!(42)
        p = rand(A)
        @test is_point(A, p)
        X = rand(A; vector_at = p)
        @test is_vector(A, p, X)
    end
    @testset "embed and project" begin
        Dm = ManifoldsBaseTestSuite.ValidationDummyManifold()
        Ad = ValidationManifold(Dm)
        p = [0.0, 0.0, 1.0]
        X = [1.0, 0.0, 0.0]
        q = embed(Ad, p)
        @test q isa ValidationMPoint
        q2 = copy(Ad, q)
        embed!(Ad, q2, p)
        @test q2 == q
        Y = embed(Ad, p, X)
        @test Y isa ValidationTangentVector
        Y2 = copy(Ad, q, Y)
        embed!(Ad, Y2, p, X)
        @test Y2 == Y
        #
        q3 = embed_project(Ad, p)
        @test q3 isa ValidationMPoint
        q4 = copy(Ad, q3)
        embed_project!(Ad, q4, p)
        @test q3 == q4
        Y3 = embed_project(Ad, p, X)
        @test Y3 isa ValidationTangentVector
        Y4 = copy(Ad, p, Y3)
        embed_project!(Ad, Y4, p, X)
        @test Y3 == Y4
    end
    @testset "riemannian_tensor" begin
        r = riemann_tensor(M, x, x, y, z)
        r2 = riemann_tensor(A, x, x, y, z)
        @test r2 isa ValidationTangentVector
        r3 = copy(A, r2)
        r3 = riemann_tensor!(A, r3, x, x, y, z)
        @test r2 == r3
    end
    @testset "_msg defaults w/strings" begin
        @test_logs (:warn, "msg") ManifoldsBase._msg(A, "msg"; error = :warn)
        @test_logs (:info, "msg") ManifoldsBase._msg(A, "msg"; error = :info)
    end
    @testset "_update_basepoint" begin
        v = ValidationTangentVector([1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        ManifoldsBase._update_basepoint!(A, v, [0.0, 0.0, 1.0])
        @test v.point == [0.0, 0.0, 1.0]
    end
    @testset "show" begin
        As = ValidationManifold(
            M;
            ignore_contexts = [:Point],
            ignore_functions = Dict(exp => :All),
        )
        sAs = repr(As)
        @test contains(sAs, "ValidationManifold")
        @test contains(sAs, repr(Dict(exp => :All)))
        @test contains(sAs, ":Point")
    end
end
