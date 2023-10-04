using RecursiveArrayTools, ManifoldsBase, Test
using Random
using ManifoldsBase: DefaultManifold, VectorSpaceType, VectorSpaceFiberType, ℝ, FiberAtPoint
struct TestVectorSpaceType <: VectorSpaceType end

struct TestFiberType <: ManifoldsBase.FiberType end

function ManifoldsBase.fiber_dimension(M::AbstractManifold, ::TestFiberType)
    return 2 * manifold_dimension(M)
end

function ManifoldsBase.inner(::VectorBundleFibers{TestVectorSpaceType}, p, X, Y)
    return 2 * dot(X, Y)
end

include("test_sphere.jl")

@testset "Tangent bundle" begin
    M = DefaultManifold(3)
    m_prod_retr = ManifoldsBase.FiberBundleProductRetraction()
    m_prod_invretr = ManifoldsBase.FiberBundleInverseProductRetraction()
    m_sasaki = SasakiRetraction(5)

    TB = TangentBundle(M)

    p = [1.0, 0.0, 0.0]
    @test injectivity_radius(TB) == Inf
    TpM = TangentSpaceAtPoint(M, p)
    @test sprint(show, TB) == "TangentBundle($(M))"
    @test base_manifold(TB) == M
    @test manifold_dimension(TB) == 2 * manifold_dimension(M)
    @test is_flat(TB) == is_flat(M)
    @test is_flat(TpM)
    @test representation_size(TB) === nothing
    @test default_inverse_retraction_method(TB) === m_prod_invretr
    @test default_retraction_method(TB) == m_prod_retr
    @test default_vector_transport_method(TB) isa
          ManifoldsBase.FiberBundleProductVectorTransport
    @test TB ===
          VectorBundle(TangentSpace, M, ManifoldsBase.FiberBundleProductVectorTransport())
    @test TB === TangentBundle(M, ManifoldsBase.FiberBundleProductVectorTransport())

    @test ManifoldsBase.FiberBundleProductVectorTransport(M) ===
          ManifoldsBase.FiberBundleProductVectorTransport(
        ParallelTransport(),
        ParallelTransport(),
    )

    @test sprint(show, VectorBundle(TestVectorSpaceType(), M)) ==
          "VectorBundle(TestVectorSpaceType(), $(M))"
    TVBF = ManifoldsBase.VectorBundleFibers(TestVectorSpaceType(), M)
    @test sprint(show, TVBF) == "VectorBundleFibers(TestVectorSpaceType(), $(M))"
    @test sprint(show, VectorSpaceFiberType(TestVectorSpaceType())) ==
          "VectorSpaceFiberType(TestVectorSpaceType())"
    @test ManifoldsBase.VectorBundleFibers(
        VectorSpaceFiberType(TestVectorSpaceType()),
        M,
    ) === TVBF
    TBF = ManifoldsBase.BundleFibers(TestFiberType(), M)
    @test sprint(show, TBF) == "BundleFibers(TestFiberType(), $(M))"
    @test norm(TVBF, p, [2.0, 2.0, 0.0]) ≈ 4.0
    @test ManifoldsBase.fiber_dimension(TBF) == 6
    @test ManifoldsBase.fiber_dimension(M, CotangentSpace) == 3

    @test ManifoldsBase.TangentBundleFibers(M) ===
          ManifoldsBase.BundleFibers(ManifoldsBase.TangentFiber, M)

    @test vector_space_dimension(TB.fiber) == 3
    @test ManifoldsBase.fiber_dimension(M, ManifoldsBase.TangentFiber) == 3
    @test ManifoldsBase.fiber_bundle_transport(TangentSpace, M) === ParallelTransport()

    @test representation_size(TB.fiber) == (3,)

    @testset "spaces at point" begin
        p = [1.0, 0.0, 0.0]
        q = [0.0, 2.0, 0.0]
        t_p = TangentSpaceAtPoint(M, p)
        t_p2 = TangentSpace(M, p)
        @test t_p == t_p2
        t_ps = sprint(show, "text/plain", t_p)
        sp = sprint(show, "text/plain", p)
        sp = replace(sp, '\n' => "\n ")
        t_ps_test = "Tangent space to the manifold $(M) at point:\n $(sp)"
        @test t_ps == t_ps_test
        @test base_manifold(t_p) == M
        @test manifold_dimension(t_p) == 3
        @test t_p.fiber.manifold == M
        @test t_p.fiber.fiber == VectorSpaceFiberType(TangentSpace)
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
        @test inner(t_p, X, X, X) ≈ 1.0
        @test norm(t_p, X, X) ≈ 1.0
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
        fiber = VectorBundleFibers(TestVectorSpaceType(), M)
        X_p = FiberAtPoint(fiber, p)
        X_ps = sprint(show, "text/plain", X_p)
        fiber_s = sprint(show, "text/plain", fiber)
        X_ps_test = "$(typeof(X_p))\nFiber:\n $(fiber_s)\nBase point:\n $(sp)"
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
    end

    @testset "Weingarten Map" begin
        p0 = [1.0, 0.0, 0.0]
        M = TangentSpaceAtPoint(M, p0)
        p = [0.0, 1.0, 1.0]
        X = [0.0, 1.0, 0.0]
        V = [1.0, 0.0, 0.0]
        @test Weingarten(M, p, X, V) == zero_vector(M, p)
    end

    @testset "Tangent Bundle and Array Partition" begin
        p = ArrayPartition([1.0, 0.0, 0.0], [0.0, 2.0, 4.0])
        @test p[TB, :point] === p.x[1]
        @test p[TB, :vector] === p.x[2]
        p[TB, :vector] = [0.0, 3.0, 1.0]
        @test p.x[2] == [0.0, 3.0, 1.0]
        p[TB, :point] = [0.0, 1.0, 0.0]
        @test p.x[1] == [0.0, 1.0, 0.0]
        @test_throws DomainError p[TB, :error]
        @test_throws DomainError p[TB, :error] = [1, 2, 3]

        @test view(p, TB, :point) === p.x[1]
        @test view(p, TB, :vector) === p.x[2]
        view(p, TB, :point) .= [2.0, 3.0, 5.0]
        @test p.x[1] == [2.0, 3.0, 5.0]
        view(p, TB, :vector) .= [-2.0, -3.0, -5.0]
        @test p.x[2] == [-2.0, -3.0, -5.0]
        @test_throws DomainError view(p, TB, :error)
    end
end

@testset "Tangent bundle of a sphere" begin
    M = TestSphere(2)
    TM = TangentBundle(M)
    pm = [1.0, 0.0, 0.0]
    Xm = [0.0, 1.0, -2.0]
    Ym = [0.0, -1.0, 2.0]
    p1 = ArrayPartition(pm, Xm)
    X1 = ArrayPartition(Xm, Ym)
    X2 = ArrayPartition(Ym, Xm)
    @test distance(TM.fiber, pm, Xm, Ym) ≈ sqrt(20)
    @test injectivity_radius(TM) == 0.0
    @test ManifoldsBase.fiber_dimension(TM.fiber) == 2
    @test bundle_projection(TM, p1) === pm

    qm = exp(M, pm, Xm)
    Xt = vector_transport_direction(TM, p1, X1, X2)
    Xref = ArrayPartition(
        [1.7592245393784949, -0.6172728764571667, 1.2345457529143333],
        [-1.7592245393784949, 0.6172728764571667, -1.2345457529143333],
    )
    Yt = similar(Xt)
    vector_transport_direction!(TM, Yt, p1, X1, X2)
    @test isapprox(Xt, Xref)
    @test isapprox(Yt, Xref)

    p3 = ArrayPartition(exp(M, pm, Xm), parallel_transport_direction(M, pm, Ym, Xm))
    X3 = ArrayPartition(
        [-1.7592245393784949, -0.6172728764571667, 1.2345457529143333],
        [1.7592245393784949, 0.6172728764571667, -1.2345457529143333],
    )
    Yt = similar(Xt)
    vector_transport_to!(TM, Yt, p1, X1, p3)
    @test isapprox(vector_transport_to(TM, p1, X1, p3), X3)
    @test isapprox(Yt, X3)
    @test zero_vector(TM, p1) == zero(p1)

    @test inner(TM, p1, X1, X2) ≈ -10.0
    @test inner(TM.fiber, pm, Xm, Ym) ≈ -5.0
    @test norm(TM, p1, X1) ≈ sqrt(10.0)

    for basis in [DefaultOrthonormalBasis(), get_basis(TM, p1, DefaultOrthonormalBasis())]
        @test length(get_vectors(TM, p1, get_basis(TM, p1, basis))) == 4
        X1c = get_coordinates(TM, p1, X1, basis)
        @test isapprox(X1c, [1.0, -2.0, -1.0, 2.0])
        Y1c = similar(X1c)
        get_coordinates!(TM, Y1c, p1, X1, basis)
        @test isapprox(X1c, Y1c)
        @test isapprox(get_vector(TM, p1, X1c, basis), X1)
        Z1 = similar(X1)
        get_vector!(TM, Z1, p1, X1c, basis)
        @test isapprox(TM, p1, Z1, X1)
    end

    Bd = get_basis(TM, p1, DiagonalizingOrthonormalBasis(X1))
    @test Bd.data isa ManifoldsBase.FiberBundleBasisData

    m_prod_retr = ManifoldsBase.FiberBundleProductRetraction()
    m_prod_invretr = ManifoldsBase.FiberBundleInverseProductRetraction()
    m_sasaki = SasakiRetraction(5)

    @test inverse_retract(TM, p1, p3, m_prod_invretr) ≈
          ArrayPartition([0.0, 1.0, -2.0], [0.0, -2.0, 4.0])
    Xn = similar(Xt)
    @test inverse_retract!(TM, Xn, p1, p3, m_prod_invretr) ≈
          ArrayPartition([0.0, 1.0, -2.0], [0.0, -2.0, 4.0])
    @test isapprox(
        TM,
        retract(TM, p1, X1, m_prod_retr),
        ArrayPartition(
            [-0.6172728764571667, 0.35184490787569894, -0.7036898157513979],
            [0.0, 0.0, 0.0],
        ),
    )
    qn = similar(p1)
    @test isapprox(
        TM,
        retract!(TM, qn, p1, X1, m_prod_retr),
        ArrayPartition(
            [-0.6172728764571667, 0.35184490787569894, -0.7036898157513979],
            [0.0, 0.0, 0.0],
        ),
    )

    @test retract(TM, p1, X1, m_sasaki) ≈ ArrayPartition(
        [-0.6172728764571667, 0.35184490787569894, -0.7036898157513979],
        [0.0, 0.0, 0.0],
    )
    @test retract!(TM, qn, p1, X1, m_sasaki) ≈ ArrayPartition(
        [-0.6172728764571667, 0.35184490787569894, -0.7036898157513979],
        [0.0, 0.0, 0.0],
    )

    @test rand(TM) isa ArrayPartition{Float64}
    @test rand(TM; vector_at = p1) isa ArrayPartition{Float64}
    @test rand(Random.default_rng(), TM) isa ArrayPartition{Float64}
    @test rand(Random.default_rng(), TM; vector_at = p1) isa ArrayPartition{Float64}

    @test project(TM, p1) ≈ p1
    @test project(TM, p1, X1) ≈ X1
    @test project(TM.fiber, pm, Xm) ≈ Xm

    @test ManifoldsBase.allocate_result(TM, rand) isa ArrayPartition
end
