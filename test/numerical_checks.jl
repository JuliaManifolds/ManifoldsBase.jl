using ManifoldsBase, Plots, Statistics, Test
# don't show plots actually
default(; show = false, reuse = true)

!(pwd() in LOAD_PATH) && (push!(LOAD_PATH, pwd()))
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
end
