using RecursiveArrayTools, ManifoldsBase, Test
using Random
using ManifoldsBase: DefaultManifold, VectorSpaceType, ℝ, Fiber

s = @__DIR__
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using ManifoldsBaseTestUtils

@testset "vector space fibers" begin
    M = DefaultManifold(3)

    p = [1.0, 0.0, 0.0]

    TpM = TangentSpace(M, p)
    @test is_flat(TpM)

    @test ManifoldsBase.fiber_dimension(M, CotangentSpaceType()) == 3

    @test ManifoldsBase.fiber_dimension(M, ManifoldsBase.TangentSpaceType()) == 3

    @testset "spaces at point" begin
        p = [1.0, 0.0, 0.0]
        q = [0.0, 2.0, 0.0]
        t_p = TangentSpace(M, p)
        t_ps = sprint(show, "text/plain", t_p)
        sp = sprint(show, "text/plain", p)
        sp = replace(sp, '\n' => "\n ")
        t_ps_test = "Tangent space to the manifold $(M) at point:\n $(sp)"
        @test t_ps == t_ps_test
        ct_p = CotangentSpace(M, p)
        ct_ps = sprint(show, "text/plain", ct_p)
        ct_ps_test = "Cotangent space to the manifold $(M) at point:\n $(sp)"
        @test ct_ps == ct_ps_test
        @test base_manifold(t_p) == M
        @test manifold_dimension(t_p) == 3
        @test t_p.manifold == M
        @test t_p.fiber_type == TangentSpaceType()
        @test t_p.point == p
        @test injectivity_radius(t_p) == Inf
        @test representation_size(t_p) == representation_size(M)
        X = [0.0, 0.0, 1.0]
        Y = [0.0, 2.0, -1.0]
        @test embed(t_p, X) == X
        @test embed(t_p, X, X) == X
        @test distance(t_p, p, q) ≈ sqrt(5)
        @test isapprox(t_p, exp(t_p, X, X), 2 * X)
        @test isapprox(t_p, log(t_p, X, Y), [0.0, 2.0, -2.0])
        @test isapprox(t_p, X, log(t_p, X, Y), [0.0, 2.0, -2.0])
        @test inner(t_p, X, X, X) ≈ 1.0
        @test norm(t_p, X) ≈ 1.0
        @test norm(t_p, X) ≈ 1.0
        @test parallel_transport_to(t_p, X, Y, X) ≈ Y
        @test vector_transport_to(t_p, X, Y, X) ≈ Y
        @test vector_transport_to(t_p, X, Y, X, ProjectionTransport()) ≈ Y
        Z = similar(X)
        @test vector_transport_to!(t_p, Z, X, Y, X, ProjectionTransport()) ≈ Y
        @test project(t_p, X, Y) ≈ Y
        @test project(t_p, Y) ≈ Y
        @test rand(t_p) isa Vector{Float64}
        @test rand(t_p; vector_at = X) isa Vector{Float64}
        @test rand(Random.default_rng(), t_p) isa Vector{Float64}
        @test rand(Random.default_rng(), t_p; vector_at = X) isa Vector{Float64}
        # generic vector space at
        X_p = Fiber(M, p, TestVectorSpaceType())
        X_ps = sprint(show, "text/plain", X_p)
        X_ps_test = "VectorSpaceFiber{ℝ, DefaultManifold{ℝ, Tuple{Int64}}, TestVectorSpaceType, Vector{Float64}}\nFiber:\n TestVectorSpaceType()DefaultManifold(3; field = ℝ)\nBase point:\n $(sp)"
        @test X_ps == X_ps_test

        for basis in
            [DefaultOrthonormalBasis(), get_basis(t_p, p, DefaultOrthonormalBasis())]
            @test length(get_vectors(t_p, p, get_basis(t_p, p, basis))) == 3
            X1c = get_coordinates(t_p, p, X, basis)
            @test isapprox(X1c, [0.0, 0.0, 1.0])
            Y1c = similar(X1c)
            get_coordinates!(t_p, Y1c, p, X, basis)
            @test isapprox(X1c, Y1c)
            @test isapprox(get_vector(t_p, p, X1c, basis), X)
            Z1 = similar(X)
            get_vector!(t_p, Z1, p, X1c, basis)
            @test isapprox(Z1, X)
        end

        @testset "scalar fiber" begin
            Ms = DefaultManifold()
            p = 1.0
            t_p = TangentSpace(Ms, p)
            @test norm(t_p, 2.0) == 2.0
        end
    end

    @testset "Weingarten Map" begin
        p0 = [1.0, 0.0, 0.0]
        M = TangentSpace(M, p0)
        p = [0.0, 1.0, 1.0]
        X = [0.0, 1.0, 0.0]
        V = [1.0, 0.0, 0.0]
        @test Weingarten(M, p, X, V) == zero_vector(M, p)
    end
end
