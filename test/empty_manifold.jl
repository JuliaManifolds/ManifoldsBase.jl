using ManifoldsBase

using Test
import Base: *
struct NonManifold <: Manifold end
struct NonMPoint <: MPoint end
struct NonTVector <: TVector end
struct NonCoTVector <: CoTVector end
*(t::Float64, X::NonTVector) = X
@testset "Manifold with empty implementation" begin
    M = NonManifold()
    p = NonMPoint()
    v = NonTVector()
    @test base_manifold(M) == M
    @test_throws ErrorException ManifoldsBase.representation_size(M)

    @test_throws ErrorException manifold_dimension(M)

    # by default isapprox compares given points or vectors
    @test isapprox(M, [0], [0])
    @test isapprox(M, [0], [0]; atol = 1e-6)
    @test !isapprox(M, [0], [1])
    @test !isapprox(M, [0], [1]; atol = 1e-6)
    @test isapprox(M, [0], [0], [0])
    @test isapprox(M, [0], [0], [0]; atol = 1e-6)
    @test !isapprox(M, [0], [0], [1])
    @test !isapprox(M, [0], [0], [1]; atol = 1e-6)

    exp_retr = ManifoldsBase.ExponentialRetraction()

    @test_throws ErrorException retract!(M, p, p, v)
    @test_throws ErrorException retract!(M, p, p, v, exp_retr)
    @test_throws ErrorException retract!(M, p, p, [0.0], 0.0)
    @test_throws ErrorException retract!(M, p, p, [0.0], 0.0, exp_retr)
    @test_throws ErrorException retract!(M, [0], [0], [0])
    @test_throws ErrorException retract!(M, [0], [0], [0], exp_retr)
    @test_throws ErrorException retract!(M, [0], [0], [0], 0.0)
    @test_throws ErrorException retract!(M, [0], [0], [0], 0.0, exp_retr)
    @test_throws ErrorException retract(M, [0], [0])
    @test_throws ErrorException retract(M, [0], [0], exp_retr)
    @test_throws ErrorException retract(M, [0], [0], 0.0)
    @test_throws ErrorException retract(M, [0], [0], 0.0, exp_retr)
    @test_throws ErrorException retract(M, [0.0], [0.0])
    @test_throws ErrorException retract(M, [0.0], [0.0], exp_retr)
    @test_throws ErrorException retract(M, [0.0], [0.0], 0.0)
    @test_throws ErrorException retract(M, [0.0], [0.0], 0.0, exp_retr)

    log_invretr = ManifoldsBase.LogarithmicInverseRetraction()

    @test_throws ErrorException inverse_retract!(M, p, p, p)
    @test_throws ErrorException inverse_retract!(M, p, p, p, log_invretr)
    @test_throws ErrorException inverse_retract!(M, [0], [0], [0])
    @test_throws ErrorException inverse_retract!(M, [0], [0], [0], log_invretr)
    @test_throws ErrorException inverse_retract(M, [0], [0])
    @test_throws ErrorException inverse_retract(M, [0], [0], log_invretr)

    @test_throws ErrorException inverse_retract(M, [0.0], [0.0])
    @test_throws ErrorException inverse_retract(M, [0.0], [0.0], log_invretr)

    @test_throws ErrorException project_point!(M, p, [0])
    @test_throws ErrorException project_point!(M, [0], [0])
    @test_throws ErrorException project_point(M, [0])

    @test_throws ErrorException project_tangent!(M, v, p, [0.0])
    @test_throws ErrorException project_tangent!(M, [0], [0], [0])
    @test_throws ErrorException project_tangent(M, [0], [0])
    @test_throws ErrorException project_tangent(M, [0.0], [0.0])

    @test_throws ErrorException inner(M, p, v, v)
    @test_throws ErrorException inner(M, [0], [0], [0])
    @test_throws ErrorException norm(M, p, v)
    @test_throws ErrorException norm(M, [0], [0])
    @test_throws ErrorException angle(M, p, v, v)
    @test_throws ErrorException angle(M, [0], [0], [0])

    @test_throws ErrorException distance(M, [0.0], [0.0])

    @test_throws ErrorException exp!(M, p, p, v)
    @test_throws ErrorException exp!(M, p, p, v, 0.0)
    @test_throws ErrorException exp!(M, [0], [0], [0])
    @test_throws ErrorException exp!(M, [0], [0], [0], 0.0)
    @test_throws ErrorException exp(M, [0], [0])
    @test_throws ErrorException exp(M, [0], [0], 0.0)
    @test_throws ErrorException exp(M, [0.0], [0.0])
    @test_throws ErrorException exp(M, [0.0], [0.0], 0.0)

    @test_throws ErrorException log!(M, v, p, p)
    @test_throws ErrorException log!(M, [0], [0], [0])
    @test_throws ErrorException log(M, [0.0], [0.0])

    @test_throws ErrorException vector_transport_to!(M, [0], [0], [0], [0])
    @test_throws ErrorException vector_transport_to(M, [0], [0], [0])
    @test_throws ErrorException vector_transport_to!(
        M,
        [0],
        [0],
        [0],
        ProjectionTransport(),
    )

    @test_throws ErrorException vector_transport_direction!(M, [0], [0], [0], [0])
    @test_throws ErrorException vector_transport_direction(M, [0], [0], [0])

    @test_throws ErrorException ManifoldsBase.vector_transport_along!(
        M,
        [0],
        [0],
        [0],
        x -> x,
    )
    @test_throws ErrorException ManifoldsBase.vector_transport_along(M, [0], [0], x -> x)

    @test_throws ErrorException injectivity_radius(M)
    @test_throws ErrorException injectivity_radius(M, [0])
    @test_throws ErrorException injectivity_radius(M, [0], exp_retr)

    @test_throws ErrorException zero_tangent_vector!(M, [0], [0])
    @test_throws ErrorException zero_tangent_vector(M, [0])

    @test check_manifold_point(M, [0]) === nothing
    @test_throws ErrorException check_manifold_point(M, p)
    @test is_manifold_point(M, [0])
    @test check_manifold_point(M, [0]) == nothing

    @test check_tangent_vector(M, [0], [0]) === nothing
    @test_throws ErrorException check_tangent_vector(M, p, v)
    @test is_tangent_vector(M, [0], [0])
    @test check_tangent_vector(M, [0], [0]) == nothing

    @test_throws ErrorException hat!(M,[0],[0],[0])
    @test_throws ErrorException vee!(M,[0],[0],[0])
end
