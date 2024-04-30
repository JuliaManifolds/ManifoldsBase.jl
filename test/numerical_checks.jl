using ManifoldsBase, Plots, Test
# don't show plots actually
default(; show = false, reuse = true)

push!(LOAD_PATH, pwd())
using ManifoldsBaseTestUtils

@testset "Numerical Check functions" begin
    @testset "Test Retraction checks" begin
        M = TestSphere(10)
        q = zeros(11)
        q[1] = 1.0
        p = zeros(11)
        p[1:4] .= 1 / sqrt(4)
        X = log(M, p, q)

        # TODO: Implement projectionretraction in the test file to check it here.
        @test check_retraction(M, ExponentialRetraction(), p, X)
        check_retraction(M, ExponentialRetraction(), p, X; plot = true)

        # TODO: Implement a non-retraction to get an error here
        # @test_throws ErrorException check_gradient(M, f, grad_fb, p, X; throw_error=true)
        # @test !check_gradient(M, f, grad_fb, p, X)

        #test window size error
        @test_throws ErrorException Manopt.find_best_slope_window(zeros(2), zeros(2), 20)
        @test_throws ErrorException Manopt.find_best_slope_window(
            zeros(2),
            zeros(2),
            [2, 20],
        )
        # Exponential Map -> exact
        @test check_retraction(M, ExponentialRetraction(), p, X)
        check_retraction(M, ExponentialRetraction(), p, X; plot = true)
    end
end
