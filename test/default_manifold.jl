using ManifoldsBase

using DoubleFloats
using ForwardDiff
using LinearAlgebra
using Random
using ReverseDiff
using StaticArrays
using Test

s = @__DIR__
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using ManifoldsBaseTestUtils

@testset "Testing Default (Euclidean)" begin
    M = ManifoldsBase.DefaultManifold(3)
    types = [
        Vector{Float64},
        SizedVector{3, Float64, Vector{Float64}},
        MVector{3, Float64},
        Vector{Float32},
        SizedVector{3, Float32, Vector{Float32}},
        MVector{3, Float32},
        Vector{Double64},
        MVector{3, Double64},
        SizedVector{3, Double64, Vector{Double64}},
        DefaultPoint{Vector{Float64}},
    ]

    @test repr(M) == "DefaultManifold(3; field = ℝ)"
    @test isa(manifold_dimension(M), Integer)
    @test manifold_dimension(M) ≥ 0
    @test base_manifold(M) == M
    @test has_components(M)
    @test number_system(M) == ManifoldsBase.ℝ
    @test ManifoldsBase.representation_size(M) == (3,)

    p = zeros(3)
    m = PolarRetraction()
    @test injectivity_radius(M) == Inf
    @test injectivity_radius(M, p) == Inf
    @test injectivity_radius(M, m) == Inf
    @test injectivity_radius(M, p, m) == Inf
    @test default_retraction_method(M) == ExponentialRetraction()
    @test default_inverse_retraction_method(M) == LogarithmicInverseRetraction()

    @test is_flat(M)

    rm = ManifoldsBase.ExponentialRetraction()
    irm = ManifoldsBase.LogarithmicInverseRetraction()

    # Representation sizes not equal
    @test ManifoldsBase.check_size(M, zeros(3, 3)) isa DomainError
    @test ManifoldsBase.check_size(M, zeros(3), zeros(3, 3)) isa DomainError
    rm2 = CustomDefinedRetraction()
    rm3 = CustomUndefinedRetraction()

    for T in types
        @testset "Type $T" begin
            pts = convert.(Ref(T), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

            @test injectivity_radius(M, pts[1]) == Inf
            @test injectivity_radius(M, pts[1], rm) == Inf
            @test injectivity_radius(M, rm) == Inf
            @test injectivity_radius(M, rm2) == 10
            @test injectivity_radius(M, pts[1], rm2) == 10
            # ManifoldsBase._injectivity_radius always requires the method to be passed
            # to reduce the number of ambiguities
            @test_throws MethodError ManifoldsBase._injectivity_radius(M, pts[1])

            tv1 = log(M, pts[1], pts[2])

            for pt in pts
                @test is_point(M, pt)
            end
            @test is_vector(M, pts[1], tv1; atol = eps(eltype(pts[1])))

            tv2 = log(M, pts[2], pts[1])
            tv3 = log(M, pts[2], pts[3])
            @test isapprox(M, pts[2], exp(M, pts[1], tv1))
            @test_logs (:info,) !isapprox(M, pts[1], pts[2]; error = :info)
            @test isapprox(M, pts[1], pts[1]; error = :info)
            @test_logs (:info,) !isapprox(
                M,
                pts[1],
                convert(T, [NaN, NaN, NaN]);
                error = :info,
            )
            @test isapprox(M, pts[1], ManifoldsBase.exp_fused(M, pts[1], tv1, 0))
            @test isapprox(M, pts[2], ManifoldsBase.exp_fused(M, pts[1], tv1, 1))
            @test isapprox(M, pts[1], exp(M, pts[2], tv2))
            if T <: Array
                @test_throws ApproximatelyError isapprox(M, pts[1], pts[2]; error = :error)
                @test_throws ApproximatelyError isapprox(
                    M,
                    pts[2],
                    tv2,
                    tv3;
                    error = :error,
                )
            end
            # test lower level fallbacks
            @test ManifoldsBase.check_approx(M, pts[1], pts[2]) isa ApproximatelyError
            @test ManifoldsBase.check_approx(M, pts[1], pts[1]) === nothing
            @test ManifoldsBase.check_approx(M, pts[2], tv2, tv2) === nothing
            @test is_point(M, retract(M, pts[1], tv1))
            @test isapprox(M, pts[1], ManifoldsBase.retract_fused(M, pts[1], tv1, 0))

            @test is_point(M, retract(M, pts[1], tv1, rm))
            @test isapprox(M, pts[1], ManifoldsBase.retract_fused(M, pts[1], tv1, 0, rm))

            new_pt = exp(M, pts[1], tv1)
            retract!(M, new_pt, pts[1], tv1)
            @test is_point(M, new_pt)
            @test !isapprox(M, pts[1], [1, 2, 3], [3, 2, 4]; error = :other)

            @test isapprox(
                M,
                pts[1],
                ManifoldsBase.retract_fused!(
                    M,
                    new_pt,
                    pts[1],
                    tv1,
                    0,
                    ExponentialRetraction(),
                ),
            )
            for p in pts
                X_p_zero = zero_vector(M, p)
                X_p_nan = NaN * X_p_zero
                @test isapprox(M, p, X_p_zero, log(M, p, p); atol = eps(eltype(p)))
                if T <: Array
                    @test_logs (:info,) !isapprox(
                        M, p, X_p_zero, X_p_nan; atol = eps(eltype(p)), error = :info,
                    )
                    @test isapprox(
                        M, p, X_p_zero, log(M, p, p); atol = eps(eltype(p)), error = :info,
                    )
                end
                @test isapprox(
                    M, p, X_p_zero, inverse_retract(M, p, p); atol = eps(eltype(p)),
                )
                @test isapprox(
                    M, p, X_p_zero, inverse_retract(M, p, p, irm); atol = eps(eltype(p)),
                )
            end
            zero_vector!(M, tv1, pts[1])
            @test isapprox(M, pts[1], tv1, zero_vector(M, pts[1]))
            log!(M, tv1, pts[1], pts[2])
            @test norm(M, pts[1], tv1) ≈ sqrt(inner(M, pts[1], tv1, tv1))

            @test isapprox(M, ManifoldsBase.exp_fused(M, pts[1], tv1, 1), pts[2])
            @test isapprox(M, ManifoldsBase.exp_fused(M, pts[1], tv1, 0), pts[1])

            @test distance(M, pts[1], pts[2]) ≈ norm(M, pts[1], tv1)
            @test distance(M, pts[1], pts[2], LogarithmicInverseRetraction()) ≈ norm(M, pts[1], tv1)

            @test mid_point(M, pts[1], pts[2]) == convert(T, [0.5, 0.5, 0.0])
            midp = allocate(pts[1])
            @test mid_point!(M, midp, pts[1], pts[2]) === midp
            @test midp == convert(T, [0.5, 0.5, 0.0])

            @test riemann_tensor(M, pts[1], tv1, tv2, tv1) == zero(tv1)
            tv_rt = allocate(tv1)
            @test riemann_tensor!(M, tv_rt, pts[1], tv1, tv2, tv1) === tv_rt
            @test tv_rt == zero(tv1)

            @test sectional_curvature(M, pts[1], tv1, log(M, pts[1], pts[3])) == 0.0
            @test sectional_curvature_max(M) == 0.0
            @test sectional_curvature_min(M) == 0.0

            q = copy(M, pts[1])
            Ts = [0.0, 1.0 / 2, 1.0]
            @testset "Geodesic interface test" begin
                @test isapprox(M, geodesic(M, pts[1], tv1)(0.0), pts[1])
                @test isapprox(M, geodesic(M, pts[1], tv1)(1.0), pts[2])
                g! = geodesic!(M, pts[1], tv1)
                g!(q, 0.0)
                isapprox(M, q, pts[1])
                g!(q, 0.0)
                isapprox(M, q, pts[2])
                geodesic!(M, q, pts[1], tv1, 1.0 / 2)
                isapprox(M, q, midp)
                @test isapprox(M, geodesic(M, pts[1], tv1, 1.0), pts[2])
                @test isapprox(M, geodesic(M, pts[1], tv1, 1.0 / 2), midp)
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2])(0.0), pts[1])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2])(1.0), pts[2])
                sg! = shortest_geodesic!(M, pts[1], pts[2])
                sg!(q, 0.0)
                isapprox(M, q, pts[1])
                sg!(q, 1.0)
                isapprox(M, q, pts[2])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2], 0.0), pts[1])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2], 1.0), pts[2])
                shortest_geodesic!(M, q, pts[1], pts[2], 0.5)
                isapprox(M, q, midp)
                @test all(
                    isapprox.(Ref(M), geodesic(M, pts[1], tv1, Ts), [pts[1], midp, pts[2]]),
                )
                @test all(
                    isapprox.(
                        Ref(M), shortest_geodesic(M, pts[1], pts[2], Ts), [pts[1], midp, pts[2]],
                    ),
                )
                Q = [copy(M, q), copy(M, q), copy(M, q)]
                geodesic!(M, Q, pts[1], tv1, Ts)
                @test all(isapprox.(Ref(M), Q, [pts[1], midp, pts[2]]))
                shortest_geodesic!(M, Q, pts[1], pts[2], Ts)
                @test all(isapprox.(Ref(M), Q, [pts[1], midp, pts[2]]))
            end

            @testset "basic linear algebra in tangent space" begin
                @test isapprox(
                    M, pts[1], 0 * tv1, zero_vector(M, pts[1]); atol = eps(eltype(pts[1])),
                )
                @test isapprox(M, pts[1], 2 * tv1, tv1 + tv1)
                @test isapprox(M, pts[1], 0 * tv1, tv1 - tv1)
                @test isapprox(M, pts[1], (-1) * tv1, -tv1)
            end

            @testset "Change Representer and Metric" begin
                G = ManifoldsBase.EuclideanMetric()
                p = pts[1]
                for X in [tv1, tv2, tv3]
                    @test change_representer(M, G, p, X) == X
                    Y = similar(X)
                    change_representer!(M, Y, G, p, X)
                    @test isapprox(M, p, Y, X)
                    @test change_metric(M, G, p, X) == X
                    Z = similar(X)
                    change_metric!(M, Z, G, p, X)
                    @test isapprox(M, p, Z, X)
                end
            end

            @testset "Hat and vee in the tangent space" begin
                X = log(M, pts[1], pts[2])
                a = vee(M, pts[1], X)
                b = similar(a)
                vee!(M, b, pts[1], X)
                Y = hat(M, pts[1], a)
                Z = similar(Y)
                hat!(M, Z, pts[1], a)
                @test a == b
                @test X == Y
                @test Z == X
                @test a == ((T <: DefaultPoint) ? vec(X.value) : vec(X))
            end

            @testset "broadcasted linear algebra in tangent space" begin
                @test isapprox(M, pts[1], 3 * tv1, 2 .* tv1 .+ tv1)
                @test isapprox(M, pts[1], -tv1, tv1 .- 2 .* tv1)
                @test isapprox(M, pts[1], -tv1, .-tv1)
                v = similar(tv1)
                v .= 2 .* tv1 .+ tv1
                @test isapprox(M, pts[1], v, 3 * tv1)
            end

            @testset "project test" begin
                # point
                @test isapprox(M, pts[1], project(M, pts[1]))
                pt = similar(pts[1])
                project!(M, pt, pts[1])
                @test isapprox(M, pt, pts[1])

                @test isapprox(M, pts[1], embed(M, pts[1]))
                pt = similar(pts[1])
                embed!(M, pt, pts[1])
                @test isapprox(M, pt, pts[1])

                # tangents
                @test isapprox(M, pts[1], tv1, project(M, pts[1], tv1))
                tv = similar(tv1)
                project!(M, tv, pts[1], tv1)
                @test isapprox(M, pts[1], tv, tv1)

                @test isapprox(M, pts[1], tv1, embed(M, pts[1], tv1))
                tv = similar(tv1)
                embed!(M, tv, pts[1], tv1)
                @test isapprox(M, pts[1], tv, tv1)
                @test Weingarten(M, pts[1], tv, tv) == zero_vector(M, pts[1])
            end

            @testset "random tests" begin
                Random.seed!(23)
                pr = rand(M)
                @test is_point(M, pr, true)
                Xr = rand(M; vector_at = pr)
                @test is_vector(M, pr, Xr, true)
                rng = MersenneTwister(42)
                P = rand(M, 3)
                @test length(P) == 3
                @test all([is_point(M, pj) for pj in P])
                Xv = rand(M, 3; vector_at = pr)
                @test length(Xv) == 3
                @test all([is_point(M, Xj) for Xj in Xv])
                # and the same again with rng upfront
                rng = MersenneTwister(42)
                pr = rand(rng, M)
                @test is_point(M, pr, true)
                Xr = rand(rng, M; vector_at = pr)
                @test is_vector(M, pr, Xr, true)
                rng = MersenneTwister(42)
                P = rand(rng, M, 3)
                @test length(P) == 3
                @test all([is_point(M, pj) for pj in P])
                Xv = rand(rng, M, 3; vector_at = pr)
                @test length(Xv) == 3
                @test all([is_point(M, Xj) for Xj in Xv])
            end

            @testset "vector transport" begin
                # test constructor
                @test default_vector_transport_method(M) == ParallelTransport()
                X1 = log(M, pts[1], pts[2])
                X2 = log(M, pts[1], pts[3])
                X1t1 = vector_transport_to(M, pts[1], X1, pts[3])
                X1t2 = zero(X1t1)
                vector_transport_to!(M, X1t2, pts[1], X1, X2, ProjectionTransport())
                X1t3 = vector_transport_direction(M, pts[1], X1, X2)
                @test ManifoldsBase.is_vector(M, pts[3], X1t1)
                @test ManifoldsBase.is_vector(M, pts[3], X1t3)
                @test isapprox(M, pts[3], X1t1, X1t3)

                # On Euclidean Space Schild & Pole are identity
                @test vector_transport_to(
                    M, pts[1], X2, pts[2], SchildsLadderTransport(),
                ) == X2
                @test vector_transport_to(M, pts[1], X2, pts[2], PoleLadderTransport()) ==
                    X2
                @test vector_transport_to(
                    M, pts[1], X2, pts[2], ScaledVectorTransport(ParallelTransport()),
                ) == X2

                # along is also the identity
                c = [
                    mid_point(M, pts[1], pts[2]), pts[2], mid_point(M, pts[2], pts[3]), pts[3],
                ]
                # check mutating ones with defaults
                p = allocate(pts[1])
                ManifoldsBase.pole_ladder!(M, p, pts[1], pts[2], pts[3])
                # -log_p3 p == log_p1 p2
                @test isapprox(M, pts[3], -log(M, pts[3], p), log(M, pts[1], pts[2]))
                ManifoldsBase.schilds_ladder!(M, p, pts[1], pts[2], pts[3])
                @test isapprox(M, pts[3], log(M, pts[3], p), log(M, pts[1], pts[2]))

                @test repr(ParallelTransport()) == "ParallelTransport()"
                @test repr(ScaledVectorTransport(ParallelTransport())) ==
                    "ScaledVectorTransport(ParallelTransport())"
            end

            @testset "ForwardDiff support" begin
                exp_f(t) = distance(M, pts[1], exp(M, pts[1], t * tv1))
                d12 = distance(M, pts[1], pts[2])
                for t in 0.1:0.1:0.9
                    @test d12 ≈ ForwardDiff.derivative(exp_f, t)
                end
                retract_f(t) = distance(M, pts[1], retract(M, pts[1], t * tv1))
                for t in 0.1:0.1:0.9
                    @test ForwardDiff.derivative(retract_f, t) ≥ 0
                end
            end

            isa(pts[1], Union{Vector, SizedVector}) && @testset "ReverseDiff support" begin
                exp_f(t) = distance(M, pts[1], exp(M, pts[1], t[1] * tv1))
                d12 = distance(M, pts[1], pts[2])
                for t in 0.1:0.1:0.9
                    @test d12 ≈ ReverseDiff.gradient(exp_f, [t])[1]
                end

                retract_f(t) = distance(M, pts[1], retract(M, pts[1], t[1] * tv1))
                for t in 0.1:0.1:0.9
                    @test ReverseDiff.gradient(retract_f, [t])[1] ≥ 0
                end
            end
        end
    end

    @testset "mid_point on 0-index arrays" begin
        M = ManifoldsBase.DefaultManifold(1)
        p1 = fill(0.0)
        p2 = fill(1.0)
        @test isapprox(M, fill(0.5), mid_point(M, p1, p2))
    end

    @testset "Retraction" begin
        a = NLSolveInverseRetraction(ExponentialRetraction())
        @test a.retraction isa ExponentialRetraction
    end

    @testset "copy of points and vectors" begin
        M = ManifoldsBase.DefaultManifold(2)
        for (p, X) in (
                ([2.0, 3.0], [4.0, 5.0]),
                (DefaultPoint([2.0, 3.0]), DefaultTangentVector([4.0, 5.0])),
            )
            q = similar(p)
            copyto!(M, q, p)
            @test p == q
            r = copy(M, p)
            @test r == p
            Y = similar(X)
            copyto!(M, Y, p, X)
            @test Y == X
            Z = copy(M, p, X)
            @test Z == X
        end

        p1 = DefaultPoint([2.0, 3.0])
        p2 = copy(p1)
        @test (p1 == p2) && (p1 !== p2)
    end

    @testset "further vector and point automatic forwards" begin
        M = ManifoldsBase.DefaultManifold(3)
        p = DefaultPoint([1.0, 0.0, 0.0])
        q = DefaultPoint([0.0, 0.0, 0.0])
        X = DefaultTangentVector([0.0, 1.0, 0.0])
        Y = DefaultTangentVector([1.0, 0.0, 0.0])
        @test angle(M, p, X, Y) ≈ π / 2
        @test inverse_retract(M, p, q, LogarithmicInverseRetraction()) == -Y
        @test retract(M, q, Y, CustomDefinedRetraction()) == p
        @test ManifoldsBase.retract_fused(M, q, Y, 1.0, CustomDefinedRetraction()) == p
        @test retract(M, q, Y, ExponentialRetraction()) == p
        @test ManifoldsBase.retract_fused(M, q, Y, 1.0, ExponentialRetraction()) == p
        # rest not implemented - so they also fall back even onto mutating
        Z = similar(Y)
        r = similar(p)
        # test passthrough using the dummy implementations
        for retr in [
                PolarRetraction(),
                ProjectionRetraction(),
                QRRetraction(),
                SoftmaxRetraction(),
                PadeRetraction(2),
                EmbeddedRetraction(ExponentialRetraction()),
                SasakiRetraction(5),
                ApproximateExponentialRetraction((;)),
            ]
            @test retract(M, q, Y, retr) == DefaultPoint(q.value + Y.value)
            @test ManifoldsBase.retract_fused(M, q, Y, 0.5, retr) ==
                DefaultPoint(q.value + 0.5 * Y.value)
            @test retract!(M, r, q, Y, retr) == DefaultPoint(q.value + Y.value)
            @test ManifoldsBase.retract_fused!(M, r, q, Y, 0.5, retr) ==
                DefaultPoint(q.value + 0.5 * Y.value)
        end

        mRK = RetractionWithKeywords(CustomDefinedKeywordRetraction(); scale = 3.0)
        pRK = allocate(p, eltype(p.value), size(p.value))
        @test retract(M, p, X, mRK) == DefaultPoint(3 * p.value + X.value)
        @test ManifoldsBase.retract_fused(M, p, X, 0.5, mRK) ==
            DefaultPoint(3 * p.value + 0.5 * X.value)
        @test retract!(M, pRK, p, X, mRK) == DefaultPoint(3 * p.value + X.value)
        @test ManifoldsBase.retract_fused!(M, pRK, p, X, 0.5, mRK) ==
            DefaultPoint(3 * p.value + 0.5 * X.value)
        mIRK = InverseRetractionWithKeywords(
            CustomDefinedKeywordInverseRetraction();
            scale = 3.0,
        )
        XIRK = allocate(X, eltype(X.value), size(X.value))
        @test inverse_retract(M, p, pRK, mIRK) ==
            DefaultTangentVector(pRK.value - 3 * p.value)
        @test inverse_retract!(M, XIRK, p, pRK, mIRK) ==
            DefaultTangentVector(pRK.value - 3 * p.value)
        p2 = allocate(p, eltype(p.value), size(p.value))
        @test size(p2.value) == size(p.value)
        X2 = allocate(X, eltype(X.value), size(X.value))
        @test size(X2.value) == size(X.value)
        X3 = ManifoldsBase.allocate_result(M, log, p, q)
        @test log!(M, X3, p, q) == log(M, p, q)
        @test X3 == log(M, p, q)
        @test log!(M, X3, p, q) == log(M, p, q)
        @test X3 == log(M, p, q)
        @test inverse_retract(M, p, q, CustomDefinedInverseRetraction()) == -2 * Y
        @test distance(M, p, q, CustomDefinedInverseRetraction()) == 2.0
        X4 = ManifoldsBase.allocate_result(M, inverse_retract, p, q)
        @test inverse_retract!(M, X4, p, q) == inverse_retract(M, p, q)
        @test X4 == inverse_retract(M, p, q)
        # rest not implemented but check passthrough
        for r in [
                PolarInverseRetraction(),
                ProjectionInverseRetraction(),
                QRInverseRetraction(),
                SoftmaxInverseRetraction(),
                ApproximateLogarithmicInverseRetraction((;)),
            ]
            @test inverse_retract(M, q, p, r) == DefaultTangentVector(p.value - q.value)
            @test inverse_retract!(M, Z, q, p, r) == DefaultTangentVector(p.value - q.value)
        end
        @test inverse_retract(
            M, q, p, EmbeddedInverseRetraction(LogarithmicInverseRetraction()),
        ) == DefaultTangentVector(p.value - q.value)
        @test inverse_retract(M, q, p, NLSolveInverseRetraction(ExponentialRetraction())) ==
            DefaultTangentVector(p.value - q.value)
        @test inverse_retract!(
            M, Z, q, p, EmbeddedInverseRetraction(LogarithmicInverseRetraction()),
        ) == DefaultTangentVector(p.value - q.value)
        @test inverse_retract!(
            M, Z, q, p, NLSolveInverseRetraction(ExponentialRetraction()),
        ) == DefaultTangentVector(p.value - q.value)
        c = ManifoldsBase.allocate_coordinates(M, p, Float64, manifold_dimension(M))
        @test c isa Vector
        @test length(c) == 3
        @test 2.0 \ X == DefaultTangentVector(2.0 \ X.value)
        @test X + Y == DefaultTangentVector(X.value + Y.value)
        @test +X == X
        @test (Y .= X) === Y
        # vector transport pass through
        @test vector_transport_to(M, p, X, q, ProjectionTransport()) == X
        @test vector_transport_direction(M, p, X, X, ProjectionTransport()) == X
        @test vector_transport_to!(M, Y, p, X, q, ProjectionTransport()) == X
        @test vector_transport_direction!(M, Y, p, X, X, ProjectionTransport()) == X
        @test vector_transport_to(M, p, X, :q, ProjectionTransport()) == X
        @test parallel_transport_to(M, p, X, q) == X
        @test parallel_transport_direction(M, p, X, X) == X
        @test parallel_transport_to!(M, Y, p, X, q) == X
        @test parallel_transport_direction!(M, Y, p, X, X) == X

        # convert with manifold
        @test convert(typeof(p), M, p.value) == p
        @test convert(typeof(X), M, p, X.value) == X
    end

    @testset "DefaultManifold and ONB" begin
        M = ManifoldsBase.DefaultManifold(3)
        p = [1.0f0, 0.0f0, 0.0f0]
        CB = get_basis(M, p, DefaultOrthonormalBasis())
        # make sure the right type is propagated
        @test CB.data isa Vector{Vector{Float32}}
        @test CB.data ==
            [[1.0f0, 0.0f0, 0.0f0], [0.0f0, 1.0f0, 0.0f0], [0.0f0, 0.0f0, 1.0f0]]

        # test complex point -> real coordinates
        MC = ManifoldsBase.DefaultManifold(3; field = ManifoldsBase.ℂ)
        p = [1.0im, 2.0im, -1.0im]
        CB = get_basis(MC, p, DefaultOrthonormalBasis(ManifoldsBase.ℂ))
        @test CB.data == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        @test CB.data isa Vector{Vector{Float64}}
        @test ManifoldsBase.coordinate_eltype(MC, p, ManifoldsBase.ℂ) === ComplexF64
        @test ManifoldsBase.coordinate_eltype(MC, p, ManifoldsBase.ℝ) === Float64
        CBR = get_basis(MC, p, DefaultOrthonormalBasis())
        @test CBR.data == [
            [1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im],
            [0.0 + 0.0im, 1.0 + 0.0im, 0.0 + 0.0im],
            [0.0 + 0.0im, 0.0 + 0.0im, 1.0 + 0.0im],
            [0.0 + 1.0im, 0.0 + 0.0im, 0.0 + 0.0im],
            [0.0 + 0.0im, 0.0 + 1.0im, 0.0 + 0.0im],
            [0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 1.0im],
        ]
    end

    @testset "Show methods" begin
        @test repr(CayleyRetraction()) == "CayleyRetraction()"
        @test repr(PadeRetraction(2)) == "PadeRetraction(2)"
    end

    @testset "Further TestArrayRepresentation" begin
        M = ManifoldsBase.DefaultManifold(3)
        p = [1.0, 0.0, 0.0]
        X = [1.0, 0.0, 0.0]
        @test is_point(M, p, true)
        @test is_vector(M, p, X, true)
        pF = [1.0, 0.0]
        XF = [0.0, 0.0]
        m = ExponentialRetraction()
        @test_throws DomainError is_point(M, pF, true)
        @test_throws DomainError is_vector(M, p, XF; error = :error)
        @test_throws DomainError is_vector(M, pF, XF, true; error = :error)
        @test injectivity_radius(M) == Inf
        @test injectivity_radius(M, p) == Inf
        @test injectivity_radius(M, p, m) == Inf
        @test injectivity_radius(M, m) == Inf
    end

    @testset "performance" begin
        @allocated isapprox(M, SA[1, 2], SA[3, 4]) == 0
    end

    @testset "scalars" begin
        M = ManifoldsBase.DefaultManifold()
        p = 1.0
        X = 2.0
        @test copy(M, p) === p
        @test copy(M, p, X) === X
    end

    @testset "static (size in type parameter)" begin
        MS = ManifoldsBase.DefaultManifold(3; parameter = :type)
        @test (@inferred representation_size(MS)) == (3,)
        @test repr(MS) == "DefaultManifold(3; field = ℝ, parameter = :type)"
        @test_throws ArgumentError ManifoldsBase.DefaultManifold(3; parameter = :foo)
    end

    @testset "complex vee and hat" begin
        MC = ManifoldsBase.DefaultManifold(3; field = ManifoldsBase.ℂ)
        p = [1im, 2 + 2im, 3.0]
        @test isapprox(vee(MC, p, [1 + 2im, 3 + 4im, 5 + 6im]), [1, 3, 5, 2, 4, 6])
        @test isapprox(hat(MC, p, [1, 3, 5, 2, 4, 6]), [1 + 2im, 3 + 4im, 5 + 6im])
    end

    ManifoldsBase.default_approximation_method(
        ::ManifoldsBase.DefaultManifold,
        ::typeof(exp),
    ) = GradientDescentEstimation()
    @testset "Estimation Method defaults" begin
        M = ManifoldsBase.DefaultManifold(3)
        # Point generic type fallback
        @test default_approximation_method(M, exp, Float64) ==
            default_approximation_method(M, exp)
        # Retraction
        @test default_approximation_method(M, retract) == default_retraction_method(M)
        @test default_approximation_method(M, retract, DefaultPoint) ==
            default_retraction_method(M)
        # Inverse Retraction
        @test default_approximation_method(M, inverse_retract) ==
            default_inverse_retraction_method(M)
        @test default_approximation_method(M, inverse_retract, DefaultPoint) ==
            default_inverse_retraction_method(M)
        # Vector Transports – all 3: to
        @test default_approximation_method(M, vector_transport_to) ==
            default_vector_transport_method(M)
        @test default_approximation_method(M, vector_transport_to, DefaultPoint) ==
            default_vector_transport_method(M)
        @test default_approximation_method(M, vector_transport_direction) ==
            default_vector_transport_method(M)
        @test default_approximation_method(M, vector_transport_direction, DefaultPoint) ==
            default_vector_transport_method(M)
    end

    @test ManifoldsBase.tangent_vector_type(M, DefaultPoint) === DefaultTangentVector
    @test ManifoldsBase.tangent_vector_type(M, Array) === Array

    @testset "Type promotion in allocation" begin
        @test ManifoldsBase.exp_fused(M, [1, 2], [2, 3], 1.0) isa Vector{Float64}
    end
    @testset "Trait forwarding" begin
        @test ManifoldsBase.get_forwarding_type(ManifoldsBase.DefaultManifold(2), Vector{Int}) ==
            ManifoldsBase.StopForwardingType()
    end
    @testset "Error on nonnumeric types on Complex" begin
        Mc = ManifoldsBase.DefaultManifold(3; field = ManifoldsBase.ℂ)
        # Error on nonnumber points and vectors
        @test_throws DomainError is_point(Mc, ["a", "b", "c"]; error = :error)
        @test_throws DomainError is_vector(Mc, zeros(3), ["a", "b", "c"]; error = :error)
        @test_throws DomainError is_vector(Mc, ["a", "b", "c"], zeros(3); error = :error)
    end
end
