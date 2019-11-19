using ManifoldsBase

using Test
import Base: *
struct NonManifold <: Manifold end
struct NonMPoint <: MPoint end
struct NonTVector <: TVector end
struct NonCoTVector <: CoTVector end
*(t::Float64,v::NonTVector) = v
@testset "Manifold with empty implementation" begin
    m = NonManifold()
    p = NonMPoint()
    v = NonTVector()
    @test is_decorator_manifold(m) == Val(false)
    @test base_manifold(m) == m
    @test_throws ErrorException ManifoldsBase.representation_size(m)

    @test_throws ErrorException manifold_dimension(m)

    # by default isapprox compares given points or vectors
    @test isapprox(m, [0], [0])
    @test isapprox(m, [0], [0]; atol=1e-6)
    @test !isapprox(m, [0], [1])
    @test !isapprox(m, [0], [1]; atol=1e-6)
    @test isapprox(m, [0], [0], [0])
    @test isapprox(m, [0], [0], [0]; atol=1e-6)
    @test !isapprox(m, [0], [0], [1])
    @test !isapprox(m, [0], [0], [1]; atol=1e-6)

    exp_retr = ManifoldsBase.ExponentialRetraction()

    @test_throws ErrorException retract!(m, p, p, v)
    @test_throws ErrorException retract!(m, p, p, v, exp_retr)
    @test_throws ErrorException retract!(m, p, p, [0.0], 0.)
    @test_throws ErrorException retract!(m, p, p, [0.0], 0., exp_retr)
    @test_throws ErrorException retract!(m, [0], [0], [0])
    @test_throws ErrorException retract!(m, [0], [0], [0], exp_retr)
    @test_throws ErrorException retract!(m, [0], [0], [0], 0.0)
    @test_throws ErrorException retract!(m, [0], [0], [0], 0.0, exp_retr)
    @test_throws ErrorException retract(m, [0], [0])
    @test_throws ErrorException retract(m, [0], [0], exp_retr)
    @test_throws ErrorException retract(m, [0], [0], 0.0)
    @test_throws ErrorException retract(m, [0], [0], 0.0, exp_retr)
    @test_throws ErrorException retract(m, [0.0], [0.0])
    @test_throws ErrorException retract(m, [0.0], [0.0], exp_retr)
    @test_throws ErrorException retract(m, [0.0], [0.0], 0.0)
    @test_throws ErrorException retract(m, [0.0], [0.0], 0.0, exp_retr)

    log_invretr = ManifoldsBase.LogarithmicInverseRetraction()

    @test_throws ErrorException inverse_retract!(m, p, p, p)
    @test_throws ErrorException inverse_retract!(m, p, p, p, log_invretr)
    @test_throws ErrorException inverse_retract!(m, [0], [0], [0])
    @test_throws ErrorException inverse_retract!(m, [0], [0], [0], log_invretr)
    @test_throws ErrorException inverse_retract(m, [0], [0])
    @test_throws ErrorException inverse_retract(m, [0], [0], log_invretr)

    @test_throws ErrorException inverse_retract(m, [0.0], [0.0])
    @test_throws ErrorException inverse_retract(m, [0.0], [0.0], log_invretr)

    @test_throws ErrorException project_point!(m, p, [0])
    @test_throws ErrorException project_point!(m, [0], [0])
    @test_throws ErrorException project_point(m, [0])

    @test_throws ErrorException project_tangent!(m, v, p, [0.0])
    @test_throws ErrorException project_tangent!(m, [0], [0], [0])
    @test_throws ErrorException project_tangent(m, [0], [0])
    @test_throws ErrorException project_tangent(m, [0.0], [0.0])

    @test_throws ErrorException inner(m, p, v, v)
    @test_throws ErrorException inner(m, [0], [0], [0])
    @test_throws ErrorException norm(m, p, v)
    @test_throws ErrorException norm(m, [0], [0])
    @test_throws ErrorException angle(m, p, v, v)
    @test_throws ErrorException angle(m, [0], [0], [0])

    @test_throws ErrorException distance(m, [0.0], [0.0])

    @test_throws ErrorException exp!(m, p, p, v)
    @test_throws ErrorException exp!(m, p, p, v, 0.0)
    @test_throws ErrorException exp!(m, [0], [0], [0])
    @test_throws ErrorException exp!(m, [0], [0], [0], 0.0)
    @test_throws ErrorException exp(m, [0], [0])
    @test_throws ErrorException exp(m, [0], [0], 0.0)
    @test_throws ErrorException exp(m, [0], [0], [0])
    @test_throws ErrorException exp(m, [0.0], [0.0])
    @test_throws ErrorException exp(m, [0.0], [0.0], 0.0)
    @test_throws ErrorException exp(m, [0.0], [0.0], [0])

    @test_throws ErrorException log!(m, v, p, p)
    @test_throws ErrorException log!(m, [0], [0], [0])
    @test_throws ErrorException log(m, [0.0], [0.0])

    @test_throws ErrorException vector_transport_to!(m, [0], [0], [0], [0])
    @test_throws ErrorException vector_transport_to(m, [0], [0], [0])
    @test_throws ErrorException vector_transport_to!(m, [0], [0], [0], ProjectionTransport())

    @test_throws ErrorException vector_transport_direction!(m, [0], [0], [0], [0])
    @test_throws ErrorException vector_transport_direction(m, [0], [0], [0])

    @test_throws ErrorException ManifoldsBase.vector_transport_along!(m, [0], [0], [0], x -> x)
    @test_throws ErrorException ManifoldsBase.vector_transport_along(m, [0], [0], x -> x)

    @test injectivity_radius(m) === Inf
    @test injectivity_radius(m, [0]) === Inf
    @test injectivity_radius(m, [0], exp_retr) === Inf

    @test_throws ErrorException zero_tangent_vector!(m, [0], [0])
    @test_throws ErrorException zero_tangent_vector(m, [0])
    
    @test manifold_point_error(m, [0]) === nothing
    @test_throws ErrorException manifold_point_error(m,p)
    @test is_manifold_point(m, [0])
    @test check_manifold_point(m, [0]) == nothing

    @test tangent_vector_error(m, [0], [0]) === nothing
    @test_throws ErrorException tangent_vector_error(m,p,v)
    @test is_tangent_vector(m, [0], [0])
    @test check_tangent_vector(m, [0], [0]) == nothing

end
