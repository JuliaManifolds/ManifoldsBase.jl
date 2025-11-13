using ManifoldsBase, Plots, Statistics, Test
# don't show plots actually
default(; show = false, reuse = true)

s = @__DIR__
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using ManifoldsBaseTestUtils

@testset "Numerical Check functions" begin
    @testset "Test retract checks" begin
        M = TestSphere(10)
        q = zeros(11)
        q[1] = 1.0
        p = zeros(11)
        p[1:4] .= 1 / sqrt(4)
        X = log(M, p, q)

        @test check_retraction(M, ProjectionRetraction(), p, X; limits = (-2.5, 0.0))
        # One call with generating a plot
        check_retraction(M, ProjectionRetraction(), p, X; limits = (-2.5, 0.0), plot = true)

        # ProjectionRetraction only works <= 1 well in stepsize
        @test_throws ErrorException check_retraction(
            M, ProjectionRetraction(), p, X; limits = (-2.5, 2.0), error = :error,
        )
        @test !check_retraction(M, ProjectionRetraction(), p, X; limits = (-2.5, 2.0))

        #test window size error
        @test_throws ErrorException ManifoldsBase.find_best_slope_window(
            zeros(2), zeros(2), 20,
        )
        @test_throws ErrorException ManifoldsBase.find_best_slope_window(
            zeros(2), zeros(2), [2, 20],
        )
        @test check_retraction(M, ExponentialRetraction(), p, X; exactness_tol = 1.0e-7)
        check_retraction(
            M, ExponentialRetraction(), p, X; plot = true, exactness_tol = 1.0e-7,
        )
    end
    @testset "Test inverse_retract checks" begin
        M = TestSphere(10)
        q = zeros(11)
        q[1] = 1.0
        p = zeros(11)
        p[1:4] .= 1 / sqrt(4)
        X = log(M, p, q)

        @test check_inverse_retraction(
            M,
            ProjectionInverseRetraction(),
            p,
            X;
            limits = (-2.5, 0.0),
        )
        # One call with generating a plot
        check_inverse_retraction(
            M,
            ProjectionInverseRetraction(),
            p,
            X;
            limits = (-2.5, 0.0),
            plot = true,
        )

        # ProjectionRetraction only works <= 1 well in stepsize
        @test_throws ErrorException check_inverse_retraction(
            M,
            ProjectionInverseRetraction(),
            p,
            X;
            limits = (-2.5, 2.0), # yields a bit too long tangents
            error = :error,
        )
        @test !check_inverse_retraction(
            M,
            ProjectionInverseRetraction(),
            p,
            X;
            limits = (-2.5, 2.0),
        )
        # Check exatness case
        @test check_inverse_retraction(
            M,
            LogarithmicInverseRetraction(),
            p,
            X;
            exactness_tol = 1.0e-7,
        )
        check_inverse_retraction(
            M,
            LogarithmicInverseRetraction(),
            p,
            X;
            plot = true,
            exactness_tol = 1.0e-7,
        )
    end
    @testset "Test vector_transport_to checks" begin
        M = TestSphere(10)
        q = zeros(11)
        q[1] = 1.0
        p = zeros(11)
        p[1:4] .= 1 / sqrt(4)
        X = log(M, p, q)
        r = zeros(11)
        r[3] = 1.0
        Y = log(M, p, r)

        @test check_vector_transport(
            M,
            ProjectionTransport(),
            p,
            X,
            Y;
            second_order = false,
        )
        # One call with generating a plot
        check_vector_transport(
            M,
            ProjectionTransport(),
            p,
            X,
            Y;
            second_order = false,
            plot = true,
        )

        # ProjectionRetraction only works <= 1 well in stepsize
        @test_throws ErrorException check_vector_transport(
            M,
            ProjectionTransport(),
            p,
            X,
            Y;
            limits = (-2.5, 2.0), # yields a bit too long tangents
            error = :error,
        )
        @test !check_vector_transport(
            M,
            ProjectionTransport(),
            p,
            X,
            Y;
            limits = (-2.5, 2.0),
        )
        # Check exactness case
        @test check_vector_transport(M, ParallelTransport(), p, X, Y; exactness_tol = 1.0e-7)
        check_vector_transport(
            M,
            ParallelTransport(),
            p,
            X,
            Y;
            plot = true,
            exactness_tol = 1.0e-7,
        )
    end
    @testset "StabilizedRetraction and its inverse" begin
        M = TestSphere(10)
        p1 = [
            -0.3676232664793671,
            0.6626844596366377,
            -0.061901832129538106,
            -0.025815090234793503,
            -0.07866889208988073,
            0.3529796828577404,
            -0.12233959915005299,
            0.004087934754632594,
            -0.37359690909547544,
            -0.2105625166208829,
            0.30253234782275296,
        ]
        p2 = copy(p1)
        p3 = copy(p1)
        p4 = copy(p1)
        p5 = copy(p1)
        X1 = [
            -1.3114742420546033,
            0.2735101165024386,
            -0.30976071299605357,
            -0.11884353912873682,
            0.5654677284268061,
            0.3525791677021648,
            0.4338068524082538,
            1.2746832300177995,
            0.7740842662015813,
            -0.06325213712780621,
            -1.460513999668164,
        ]
        X2 = copy(X1)
        X3 = copy(X1)
        X4 = copy(X1)
        X5 = copy(X1)
        # inverse
        p6 = copy(p1) # keep for checks
        p7 = exp(M, p1, X1)
        p8 = exp(M, p1, X1)
        p9 = exp(M, p1, X1)
        Y1 = copy(X1)
        Y2 = copy(X1)
        Y3 = copy(X1)
        q = similar(p1)
        SR = StabilizedRetraction()
        @test repr(SR) == "StabilizedRetraction()"
        SIR = StabilizedInverseRetraction()
        @test repr(SIR) == "StabilizedInverseRetraction()"
        err_eps = 1.0e-5
        for _ in 1:1000
            X1 .+= err_eps
            exp!(M, q, p1, X1)
            parallel_transport_to!(M, X1, p1, X1, q)
            p1 .= q
            # mutating non-fused
            X2 .+= err_eps
            retract!(M, q, p2, X2, SR)
            parallel_transport_to!(M, X2, p2, X2, q)
            p2 .= q
            # mutating fused
            X3 .+= err_eps
            ManifoldsBase.retract_fused!(M, q, p3, X3, 1.0, SR)
            parallel_transport_to!(M, X3, p3, X3, q)
            p3 .= q
            # non-mutating non-fused
            X4 .+= err_eps
            q4 = retract(M, p4, X4, SR)
            parallel_transport_to!(M, X4, p4, X4, q4)
            p4 .= q4
            # non-mutating fused
            X5 .+= err_eps
            q5 = ManifoldsBase.retract_fused(M, p5, X5, 1.0, SR)
            parallel_transport_to!(M, X5, p5, X5, q5)
            p5 .= q5

            p7 .+= err_eps
            Y1 .= log(M, p6, p7)
            p8 .+= err_eps
            Y2 .= inverse_retract(M, p6, p8, SIR)
            p9 .+= err_eps
            inverse_retract!(M, Y3, p6, p9, SIR)
        end
        @test !is_point(M, p1)
        @test is_point(M, p2; error = :error)
        @test is_point(M, p3; error = :error)
        @test is_point(M, p4; error = :error)
        @test is_point(M, p5; error = :error)
        # test the inverse as well
        @test !is_vector(M, p6, Y1)
        @test is_vector(M, p6, Y2; error = :error, atol = 1.0e-16)
        @test is_vector(M, p6, Y3; error = :error, atol = 1.0e-16)
    end
    @testset "Test check_geodesic" begin
        M = TestSphere(10)
        q = zeros(11)
        q[1] = 1.0
        p = zeros(11)
        p[1:4] .= 1 / sqrt(4)
        X = log(M, p, q)
        @test check_geodesic(M, p, X)
        # One call with generating a plot
        check_geodesic(M, p, X; plot = true)

        # Check with too large steps
        X2 = zeros(11) # but make it non-tangent to q
        X2[1] = 1.0
        X2[2:4] .= 2.0
        @test_throws ErrorException check_geodesic(M, q, X2; error = :error)
    end
end
