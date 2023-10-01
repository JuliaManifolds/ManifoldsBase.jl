using RecursiveArrayTools, ManifoldsBase, Test
using ManifoldsBase: DefaultManifold, VectorSpaceType, ℝ
struct TestVectorSpaceType <: VectorSpaceType end

@testset "Tangent bundle" begin
    M = DefaultManifold(3)
    m_prod_retr = ManifoldsBase.VectorBundleProductRetraction()
    m_prod_invretr = ManifoldsBase.VectorBundleInverseProductRetraction()
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
          ManifoldsBase.VectorBundleProductVectorTransport
    CTB = CotangentBundle(M)
    @test sprint(show, CTB) == "CotangentBundle($(M))"
    @test sprint(show, VectorBundle(TestVectorSpaceType(), M)) ==
          "VectorBundle(TestVectorSpaceType(), $(M))"

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
        @test t_p.fiber.fiber == TangentSpace
        @test ct_p.fiber.fiber == CotangentSpace
        @test t_p.point == p
        @test ct_p.point == p
        @test injectivity_radius(t_p) == Inf
        @test representation_size(t_p) == representation_size(M)
        X = [0.0, 0.0, 1.0]
        @test embed(t_p, X) == X
        @test embed(t_p, X, X) == X
        # generic vector space at
        fiber = VectorBundleFibers(TestVectorSpaceType(), M)
        X_p = VectorSpaceAtPoint(fiber, p)
        X_ps = sprint(show, "text/plain", X_p)
        fiber_s = sprint(show, "text/plain", fiber)
        X_ps_test = "$(typeof(X_p))\nFiber:\n $(fiber_s)\nBase point:\n $(sp)"
        @test X_ps == X_ps_test
        @test_throws ErrorException project(fiber, p, X)
        @test_throws ErrorException norm(fiber, p, X)
        @test_throws ErrorException distance(fiber, p, X, X)
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

    @testset "Error messages" begin
        vbf = VectorBundleFibers(TestVectorSpaceType(), M)
        @test_throws ErrorException inner(vbf, [1, 2, 3], [1, 2, 3], [1, 2, 3])
        @test_throws ErrorException ManifoldsBase.project!(
            vbf,
            [1, 2, 3],
            [1, 2, 3],
            [1, 2, 3],
        )
        @test_throws MethodError zero_vector!(vbf, [1, 2, 3], [1, 2, 3])
        @test_throws MethodError vector_space_dimension(vbf)
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
