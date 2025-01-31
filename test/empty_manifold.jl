using ManifoldsBase, Test



using Test

@testset "AbstractManifold with empty implementation" begin
    M = NonManifold()
    p = NonMPoint()
    v = NonTangentVector()
    @test base_manifold(M) === M
    @test number_system(M) === ℝ
    @test representation_size(M) === nothing
    @test !has_components(M)

    @test_throws MethodError manifold_dimension(M)

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

    @test_throws MethodError retract!(M, p, p, v)
    @test_throws MethodError retract!(M, p, p, v, exp_retr)
    @test_throws MethodError retract_fused!(M, p, p, [0.0], 0.0)
    @test_throws MethodError retract_fused!(M, p, p, [0.0], 0.0, exp_retr)
    @test_throws MethodError retract!(M, [0], [0], [0])
    @test_throws MethodError retract!(M, [0], [0], [0], exp_retr)
    @test_throws MethodError retract_fused!(M, [0], [0], [0], 0.0)
    @test_throws MethodError retract_fused!(M, [0], [0], [0], 0.0, exp_retr)
    @test_throws MethodError retract(M, [0], [0])
    @test_throws MethodError retract(M, [0], [0], exp_retr)
    @test_throws MethodError retract_fused(M, [0], [0], 0.0)
    @test_throws MethodError retract_fused(M, [0], [0], 0.0, exp_retr)
    @test_throws MethodError retract(M, [0.0], [0.0])
    @test_throws MethodError retract(M, [0.0], [0.0], exp_retr)
    @test_throws MethodError retract_fused(M, [0.0], [0.0], 0.0)
    @test_throws MethodError retract_fused(M, [0.0], [0.0], 0.0, exp_retr)
    @test_throws MethodError retract(M, [0.0], [0.0], NotImplementedRetraction())

    log_invretr = ManifoldsBase.LogarithmicInverseRetraction()

    @test_throws MethodError inverse_retract!(M, p, p, p)
    @test_throws MethodError inverse_retract!(M, p, p, p, log_invretr)
    @test_throws MethodError inverse_retract!(M, [0], [0], [0])
    @test_throws MethodError inverse_retract!(M, [0], [0], [0], log_invretr)
    @test_throws MethodError inverse_retract(M, [0], [0])
    @test_throws MethodError inverse_retract(M, [0], [0], log_invretr)
    @test_throws MethodError inverse_retract(M, [0.0], [0.0])
    @test_throws MethodError inverse_retract(M, [0.0], [0.0], log_invretr)
    @test_throws MethodError inverse_retract(
        M,
        [0.0],
        [0.0],
        NotImplementedInverseRetraction(),
    )

    @test_throws MethodError project!(M, p, [0])
    @test_throws MethodError project!(M, [0], [0])
    @test_throws MethodError project(M, [0])

    @test_throws MethodError project!(M, v, p, [0.0])
    @test_throws MethodError project!(M, [0], [0], [0])
    @test_throws MethodError project(M, [0], [0])
    @test_throws MethodError project(M, [0.0], [0.0])

    @test_throws MethodError inner(M, p, v, v)
    @test_throws MethodError inner(M, [0], [0], [0])
    @test_throws MethodError norm(M, p, v)
    @test_throws MethodError norm(M, [0], [0])
    @test_throws MethodError angle(M, p, v, v)
    @test_throws MethodError angle(M, [0], [0], [0])

    @test_throws MethodError distance(M, [0.0], [0.0])

    @test_throws MethodError exp!(M, p, p, v)
    @test_throws MethodError ManifoldsBase.exp_fused!(M, p, p, v, 0.0)
    @test_throws MethodError exp!(M, [0], [0], [0])
    @test_throws MethodError ManifoldsBase.exp_fused!(M, [0], [0], [0], 0.0)
    @test_throws MethodError exp(M, [0], [0])
    @test_throws MethodError ManifoldsBase.exp_fused(M, [0], [0], 0.0)
    @test_throws MethodError exp(M, [0.0], [0.0])
    @test_throws MethodError ManifoldsBase.exp_fused(M, [0.0], [0.0], 0.0)

    @test_throws MethodError embed!(M, p, [0]) # no copy for NoPoint p
    @test embed!(M, [0], [0]) == [0]
    @test embed(M, [0]) == [0]

    # Identity
    @test_throws MethodError embed!(M, v, p, [0.0]) # no copyto
    @test embed!(M, [0], [0], [0]) == [0]
    @test_throws MethodError embed(M, [0], v) # no copyto
    @test embed(M, [0.0], [0.0]) == [0.0]


    @test_throws MethodError log!(M, v, p, p)
    @test_throws MethodError log!(M, [0], [0], [0])
    @test_throws MethodError log(M, [0.0], [0.0])

    @test_throws MethodError vector_transport_to!(M, [0], [0], [0], [0])
    @test_throws MethodError vector_transport_to(M, [0], [0], [0])
    @test_throws MethodError vector_transport_to!(M, [0], [0], [0], ProjectionTransport())

    @test_throws MethodError vector_transport_direction!(M, [0], [0], [0], [0])
    @test_throws MethodError vector_transport_direction(M, [0], [0], [0])

    @test_throws MethodError injectivity_radius(M)
    @test_throws MethodError injectivity_radius(M, [0])
    @test_throws MethodError injectivity_radius(M, [0], exp_retr)
    @test_throws MethodError injectivity_radius(M, exp_retr)

    @test_throws MethodError zero_vector!(M, [0], [0])
    @test_throws MethodError zero_vector(M, [0])

    @test ManifoldsBase.check_point(M, [0]) === nothing
    @test ManifoldsBase.check_point(M, p) === nothing
    @test is_point(M, [0])
    @test ManifoldsBase.check_point(M, [0]) === nothing

    @test ManifoldsBase.check_vector(M, [0], [0]) === nothing
    @test ManifoldsBase.check_vector(M, p, v) === nothing
    @test is_vector(M, [0], [0])
    @test ManifoldsBase.check_vector(M, [0], [0]) === nothing

    @test_throws MethodError hat!(M, [0], [0], [0])
    @test_throws MethodError vee!(M, [0], [0], [0])
end
