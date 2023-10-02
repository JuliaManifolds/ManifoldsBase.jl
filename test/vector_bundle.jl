using RecursiveArrayTools, ManifoldsBase, Test
using ManifoldsBase: DefaultManifold, VectorSpaceType, VectorSpaceFiberType, ℝ, FiberAtPoint
struct TestVectorSpaceType <: VectorSpaceType end

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
    CTB = CotangentBundle(M)
    @test sprint(show, CTB) == "CotangentBundle($(M))"
    @test sprint(show, VectorBundle(TestVectorSpaceType(), M)) ==
          "VectorBundle(TestVectorSpaceType(), $(M))"

    @test vector_space_dimension(TB.fiber) == 3
    @test vector_space_dimension(CTB.fiber) == 3

    @testset "spaces at point" begin
        p = [1.0, 0.0, 0.0]
        t_p = TangentSpaceAtPoint(M, p)
        t_p2 = TangentSpace(M, p)
        @test t_p == t_p2
        ct_p = CotangentSpaceAtPoint(M, p)
        t_ps = sprint(show, "text/plain", t_p)
        sp = sprint(show, "text/plain", p)
        sp = replace(sp, '\n' => "\n ")
        t_ps_test = "Tangent space to the manifold $(M) at point:\n $(sp)"
        @test t_ps == t_ps_test
        @test base_manifold(t_p) == M
        @test base_manifold(ct_p) == M
        @test t_p.fiber.manifold == M
        @test ct_p.fiber.manifold == M
        @test t_p.fiber.fiber == VectorSpaceFiberType(TangentSpace)
        @test ct_p.fiber.fiber == VectorSpaceFiberType(CotangentSpace)
        @test t_p.point == p
        @test ct_p.point == p
        @test injectivity_radius(t_p) == Inf
        @test representation_size(t_p) == representation_size(M)
        X = [0.0, 0.0, 1.0]
        @test embed(t_p, X) == X
        @test embed(t_p, X, X) == X
        # generic vector space at
        fiber = VectorBundleFibers(TestVectorSpaceType(), M)
        X_p = FiberAtPoint(fiber, p)
        X_ps = sprint(show, "text/plain", X_p)
        fiber_s = sprint(show, "text/plain", fiber)
        X_ps_test = "$(typeof(X_p))\nFiber:\n $(fiber_s)\nBase point:\n $(sp)"
        @test X_ps == X_ps_test
    end

    @testset "tensor product" begin
        TT = ManifoldsBase.TensorProductType(TangentSpace, TangentSpace)
        TTs = "TensorProductType(TangentSpace, TangentSpace)"
        VBF = VectorBundleFibers(TT, M)
        @test sprint(show, TT) == TTs
        @test vector_space_dimension(VBF) == 9
        @test base_manifold(VBF) == M
        @test sprint(show, VBF) == "VectorBundleFibers($(TTs), $(M))"
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

    qm = exp(M, pm, Xm)
    Xt = vector_transport_direction(TM, p1, X1, X2)
    p2 = ArrayPartition(
        [1.7592245393784949, -0.6172728764571667, 1.2345457529143333],
        [-1.7592245393784949, 0.6172728764571667, -1.2345457529143333],
    )
    Yt = similar(Xt)
    vector_transport_direction!(TM, Yt, p1, X1, X2)
    @test isapprox(Xt, p2)
    @test isapprox(Yt, p2)

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

    X1c = get_coordinates(TM, p1, X1, DefaultOrthonormalBasis())
    @test isapprox(X1c, [1.0, -2.0, -1.0, 2.0])
    @test isapprox(get_vector(TM, p1, X1c, DefaultOrthonormalBasis()), X1)

end
