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
            M,
            ProjectionRetraction(),
            p,
            X;
            limits = (-2.5, 2.0),
            error = :error,
        )
        @test !check_retraction(M, ProjectionRetraction(), p, X; limits = (-2.5, 2.0))

        #test window size error
        @test_throws ErrorException ManifoldsBase.find_best_slope_window(
            zeros(2),
            zeros(2),
            20,
        )
        @test_throws ErrorException ManifoldsBase.find_best_slope_window(
            zeros(2),
            zeros(2),
            [2, 20],
        )
        @test check_retraction(M, ExponentialRetraction(), p, X; exactness_tol = 1e-7)
        check_retraction(
            M,
            ExponentialRetraction(),
            p,
            X;
            plot = true,
            exactness_tol = 1e-7,
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
            exactness_tol = 1e-7,
        )
        check_inverse_retraction(
            M,
            LogarithmicInverseRetraction(),
            p,
            X;
            plot = true,
            exactness_tol = 1e-7,
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
        @test check_vector_transport(M, ParallelTransport(), p, X, Y; exactness_tol = 1e-7)
        check_vector_transport(
            M,
            ParallelTransport(),
            p,
            X,
            Y;
            plot = true,
            exactness_tol = 1e-7,
        )
    end
    @testset "StabilizedRetraction" begin
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
        q = similar(p1)
        SR = StabilizedRetraction()
        err_eps = 1e-5
        for _ in 1:1000
            X1 .+= err_eps
            retract!(M, q, p1, X1, SR)
            parallel_transport_to!(M, X1, p1, X1, q)
            p1 .= q

            X2 .+= err_eps
            exp!(M, q, p2, X2)
            parallel_transport_to!(M, X2, p2, X2, q)
            p2 .= q

            X3 .+= err_eps
            ManifoldsBase.retract_fused!(M, q, p3, X3, 1.0, SR)
            parallel_transport_to!(M, X3, p3, X3, q)
            p3 .= q
        end
        @test is_point(M, p1)
        @test !is_point(M, p2)
        @test is_point(M, p3)
    end
end
