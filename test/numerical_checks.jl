using ManifoldsBase, Plots, Test
# don't show plots actually
default(; show = false, reuse = true)

!(pwd() in LOAD_PATH) && (push!(LOAD_PATH, pwd()))
using ManifoldsBaseTestUtils

@testset "Numerical Check functions" begin
    @testset "Test Retraction checks" begin
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
        @test !check_retraction(M, ProjectionRetraction(), p, X)

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
end
