using Test
using ManifoldsBase
using ManifoldsBase: submanifold_component, submanifold_components
using ManifoldsBase:
    AbstractNumbers, ℝ, ℂ, NestedReplacingPowerRepresentation, ProductBasisData
using LinearAlgebra
using Random
using RecursiveArrayTools

include("test_sphere.jl")

@testset "Product manifold" begin
    M1 = TestSphere(2)
    M2 = ManifoldsBase.DefaultManifold(2, 2)
    @test (@inferred ProductManifold(M1, M2)) isa ProductManifold

    M = ProductManifold(M1, M2)

    @test_throws MethodError ProductManifold()

    p1 = ArrayPartition([1, 0.0, 0.0], [4.0 5.0; 6.0 7.0])
    p2 = ArrayPartition([0.0, 1.0, 0.0], [4.0 8.0; 3.0 7.5])
    X1 = ArrayPartition([0.0, 1.0, 0.2], [4.0 0.0; 2.0 7.0])

    @test !is_flat(M)
    @test M[1] == M1
    @test M[2] == M2
    @test injectivity_radius(M) ≈ π
    @test injectivity_radius(
        M,
        ProductRetraction(ExponentialRetraction(), ExponentialRetraction()),
    ) ≈ π
    @test injectivity_radius(M, ExponentialRetraction()) ≈ π
    @test injectivity_radius(
        M,
        p1,
        ProductRetraction(ExponentialRetraction(), ExponentialRetraction()),
    ) ≈ π
    @test injectivity_radius(M, p1) ≈ π
    @test injectivity_radius(M, p1, ExponentialRetraction()) ≈ π

    @test ManifoldsBase.number_of_components(M) == 2
    @test ManifoldsBase.submanifold(M, 1) === M1
    @test ManifoldsBase.submanifold(M, [2, 1]) === ProductManifold(M2, M1)
    @test ManifoldsBase.check_point(M, p1) === nothing
    @test ManifoldsBase.check_vector(M, p1, X1) === nothing
    # test that arrays are not points
    @test_throws DomainError is_point(M, [1, 2], true)
    @test ManifoldsBase.check_point(M, [1, 2]) isa DomainError
    @test_throws DomainError is_vector(M, 1, [1, 2], true; check_base_point = false)
    @test ManifoldsBase.check_vector(M, 1, [1, 2]; check_base_point = false) isa DomainError
    #default fallbacks for check_size, Product not working with Arrays
    @test ManifoldsBase.check_size(M, zeros(2)) isa DomainError
    @test ManifoldsBase.check_size(M, zeros(2), zeros(3)) isa DomainError

    retraction_methods = [
        ManifoldsBase.ProductRetraction(
            ManifoldsBase.ExponentialRetraction(),
            ManifoldsBase.ExponentialRetraction(),
        ),
    ]
    inverse_retraction_methods = [
        ManifoldsBase.InverseProductRetraction(
            ManifoldsBase.LogarithmicInverseRetraction(),
            ManifoldsBase.LogarithmicInverseRetraction(),
        ),
    ]

    @test ManifoldsBase.default_retraction_method(M) === retraction_methods[1]
    @test ManifoldsBase.default_inverse_retraction_method(M) ===
          inverse_retraction_methods[1]
    @test ManifoldsBase.default_inverse_retraction_method(M, typeof(X1)) ===
          inverse_retraction_methods[1]

    @testset "get_component, set_component!, getindex and setindex!" begin
        @test get_component(M, p1, 1) == p1.x[1]
        @test get_component(M, p1, Val(1)) == p1.x[1]
        @test p1[M, 1] == p1.x[1]
        @test p1[M, Val(1)] == p1.x[1]
        @test p1[M, 1] isa Vector
        @test p1[M, Val(1)] isa Vector
        p2c = [5 6; 4 0]
        set_component!(M, p1, p2c, 2)
        @test get_component(M, p1, 2) == p2c
        p1[M, 2] = 2 * p2c
        @test p1[M, 2] == 2 * p2c
        p3 = [11.0 15.0; 3 3]
        set_component!(M, p1, p3, Val(2))
        @test get_component(M, p1, Val(2)) == p3
        p1[M, Val(2)] = 2 * p3
        @test p1[M, Val(2)] == 2 * p3

        p1c = copy(p1)
        p1c.x[1][1] = -123.0
        @test p1c.x[1][1] == -123.0
        @test p1.x[1][1] == 1.0
        copyto!(p1c, p1)
        @test p1c.x[1][1] == 1.0

        p1c.x[1][1] = -123.0
        copyto!(p1, p1c)
        @test p1.x[1][1] == -123.0
    end

    @testset "some ArrayPartition functions" begin
        q = allocate(p1)
        @test q.x[1] isa Vector
        p = ArrayPartition([[0.0, 1.0, 0.0]], [0.0, 0.0])
        q = allocate(p, Int)
        @test q.x[1] isa Vector{Vector{Int}}

        q = allocate(p1)
        copyto!(M, q, p1)
        @test isapprox(q, p1)
        @test ManifoldsBase.allocate_result(M, zero_vector) isa ArrayPartition
    end

    @testset "allocate on PowerManifold of ProductManifold" begin
        q = allocate([p1])
        @test q[1] isa ArrayPartition
        @test q[1].x[1] isa Vector
    end

    p1 = ArrayPartition([1.0, 0.0, 0.0], [4.0 5.0; 6.0 7.0])
    p2 = ArrayPartition([0.0, 1.0, 0.0], [4.0 8.0; 3.0 7.5])
    X1 = ArrayPartition([0.0, 1.0, 0.2], [4.0 0.0; 2.0 7.0])
    X2 = ArrayPartition([0.0, -1.0, 0.4], [3.0 1.0; -2.0 2.0])

    @testset "Basic operations" begin
        @test manifold_dimension(M) == 6
        @test representation_size(M) == (7,)
        @test distance(M, p1, p2) ≈ 4.551637188998299
        qr = similar(p1)
        exp!(M, qr, p1, X1)
        @test exp(M, p1, X1) ≈ ArrayPartition(
            [0.5235330372543839, 0.8354600062374664, 0.1670920012474933],
            [8.0 5.0; 8.0 14.0],
        )
        @test exp(M, p1, X1) ≈ qr
        @test exp(M, p1, X1, 2.0) ≈ exp(M, p1, 2 * X1)
        exp!(M, qr, p1, X1, 2.0)
        @test qr ≈ exp(M, p1, 2 * X1)
        @test qr ≈ retract(
            M,
            p1,
            2 * X1,
            ProductRetraction(ExponentialRetraction(), ExponentialRetraction()),
        )
        # single retraction gets “broadcasted”
        @test qr ≈ retract(M, p1, 2 * X1, ExponentialRetraction())
        qr2 = similar(p1)
        retract!(
            M,
            qr2,
            p1,
            2 * X1,
            ProductRetraction(ExponentialRetraction(), ExponentialRetraction()),
        )
        @test qr2 ≈ qr
        Xr = similar(X1)
        log!(M, Xr, p1, p2)
        @test log(M, p1, p2) ≈
              ArrayPartition([0.0, 1.5707963267948966, 0.0], [0.0 3.0; -3.0 0.5])
        @test log(M, p1, p2) ≈ Xr
        @test inverse_retract(
            M,
            p1,
            p2,
            InverseProductRetraction(
                LogarithmicInverseRetraction(),
                LogarithmicInverseRetraction(),
            ),
        ) ≈ Xr
        # single inverse retraction gets “broadcasted”
        @test inverse_retract(M, p1, p2, LogarithmicInverseRetraction()) ≈ Xr

        Zr = similar(X1)
        inverse_retract!(
            M,
            Zr,
            p1,
            p2,
            InverseProductRetraction(
                LogarithmicInverseRetraction(),
                LogarithmicInverseRetraction(),
            ),
        )
        @test Zr ≈ Xr

        @test inner(M, p1, X1, X2) ≈ 21.08
        @test norm(M, p1, X1) ≈ sqrt(70.04)
        @test project(M, p1) ≈ p1
        q1 = similar(p1)
        project!(M, q1, p1)
        @test q1 ≈ p1
        @test project(M, p1, X1) ≈ X1
        Y1 = similar(X1)
        project!(M, Y1, p1, X1)
        @test X1 ≈ Y1
        @test zero_vector(M, p1) == zero(X1)

        @test rand(M) isa ArrayPartition
        @test rand(M; vector_at = p1) isa ArrayPartition
        @test rand(Random.default_rng(), M) isa ArrayPartition
        @test rand(Random.default_rng(), M; vector_at = p1) isa ArrayPartition
        @test is_point(M, rand!(M, q1))
        @test is_vector(M, p1, rand!(M, Y1; vector_at = p1))
    end

    @testset "Broadcasting" begin
        br_result = p1 .+ 2.0 .* p2
        @test br_result isa ArrayPartition
        @test br_result.x[1] ≈ [1.0, 2.0, 0.0]
        @test br_result.x[2] ≈ [12.0 21.0; 12.0 22.0]

        br_result .= 2.0 .* p1 .+ p2
        @test br_result.x[1] ≈ [2.0, 1.0, 0.0]
        @test br_result.x[2] ≈ [12.0 18.0; 15.0 21.5]

        br_result .= p1
        @test br_result.x[1] ≈ p1.x[1]
        @test br_result.x[2] ≈ p1.x[2]

        @test axes(p1) == (Base.OneTo(7),)

        @test ManifoldsBase._get_vector_cache_broadcast(p1) === Val(false)

        # errors
        p3 = ArrayPartition([3.0, 4.0, 5.0], [2.0, 5.0], [3.0, 2.0])
        @test_throws DimensionMismatch p1 .+ p3
        @test_throws DimensionMismatch p1 .= p3
    end

    @testset "CompositeManifoldError" begin
        Mpr = ProductManifold(M1, M1)
        let p1 = [1.0, 0.0, 0.0],
            p2 = [0.0, 1.0, 0.0],
            X1 = [0.0, 1.0, 0.2],
            X2 = [1.0, 0.0, 0.2]

            p = ArrayPartition(p1, p2)
            X = ArrayPartition(X1, X2)
            pf = ArrayPartition(p1, X1)
            Xf = ArrayPartition(X1, p2)
            @test is_point(Mpr, p, true)
            @test_throws CompositeManifoldError is_point(Mpr, X, true)
            @test_throws ComponentManifoldError is_vector(Mpr, pf, X, true)
            @test_throws ComponentManifoldError is_vector(Mpr, p, Xf, true)
        end
    end

    @testset "Weingarten" begin
        Mpr = ProductManifold(M1, M1)
        p = [1.0, 0.0, 0.0]
        X = [0.0, 0.2, 0.0]
        V = [0.1, 0.0, 0.0] #orthogonal to TpM -> parallel to p
        @test isapprox(
            M,
            Weingarten(
                Mpr,
                ArrayPartition(p, p),
                ArrayPartition(X, X),
                ArrayPartition(V, V),
            ),
            ArrayPartition(-0.1 * X, -0.1 * X),
        )
    end

    @testset "arithmetic" begin
        @test isapprox(M, p1 + p2, ArrayPartition([1.0, 1.0, 0.0], [8.0 13.0; 9.0 14.5]))
        @test isapprox(M, p1 - p2, ArrayPartition([1.0, -1.0, 0.0], [0.0 -3.0; 3.0 -0.5]))
        @test isapprox(M, -p1, ArrayPartition([-1.0, -0.0, -0.0], [-4.0 -5.0; -6.0 -7.0]))
        @test isapprox(M, p1 * 2, ArrayPartition([2.0, 0.0, 0.0], [8.0 10.0; 12.0 14.0]))
        @test isapprox(M, 2 * p1, ArrayPartition([2.0, 0.0, 0.0], [8.0 10.0; 12.0 14.0]))
        @test isapprox(M, p1 / 2, ArrayPartition([0.5, 0.0, 0.0], [2.0 2.5; 3.0 3.5]))
    end

    @testset "Show methods" begin
        M2p = ProductManifold(M1, M1, M2, M2)
        @test sprint(show, M2p) == "ProductManifold($(M1), $(M1), $(M2), $(M2))"
        withenv("LINES" => 10, "COLUMNS" => 100) do
            @test sprint(show, "text/plain", ProductManifold(M1)) ==
                  "ProductManifold with 1 submanifold:\n $(M1)"
            @test sprint(show, "text/plain", M2p) ==
                  "ProductManifold with 4 submanifolds:\n $(M1)\n $(M1)\n $(M2)\n $(M2)"
            return nothing
        end
        withenv("LINES" => 7, "COLUMNS" => 100) do
            @test sprint(show, "text/plain", M2p) ==
                  "ProductManifold with 4 submanifolds:\n $(M1)\n ⋮\n $(M2)"
            return nothing
        end

        @test sprint(show, "text/plain", ProductManifold(M, M)) == """
        ProductManifold with 2 submanifolds:
         ProductManifold($(M1), $(M2))
         ProductManifold($(M1), $(M2))"""
    end

    @testset "Change Representer and Metric" begin
        G = ManifoldsBase.EuclideanMetric()
        @test change_representer(M, G, p1, X1) == X1
        Y1 = similar(X1)
        change_representer!(M, Y1, G, p1, X1)
        @test isapprox(M, p1, Y1, X1)
        @test change_metric(M, G, p1, X1) == X1
        Z1 = similar(X1)
        change_metric!(M, Z1, G, p1, X1)
        @test isapprox(M, p1, Z1, X1)
    end

    @testset "product vector transport" begin
        X = log(M, p1, p2)
        m = ProductVectorTransport(ParallelTransport(), ParallelTransport())
        @test default_vector_transport_method(M) === m
        @test default_vector_transport_method(M, typeof(X)) === m
        Y = vector_transport_to(M, p1, X, p2, m)
        Y2 = similar(Y)
        vector_transport_to!(M, Y2, p1, X, p2, m)
        Z = -log(M, p2, p1)
        @test isapprox(M, p2, Y, Z)
        @test isapprox(M, p2, Y2, Z)

        Y3 = vector_transport_direction(M, p1, X, X, m)
        Y4 = similar(Y)
        vector_transport_direction!(M, Y4, p1, X, X, m)
        @test isapprox(M, p2, Y3, Z)
        @test isapprox(M, p2, Y4, Z)
    end

    @testset "Implicit product vector transport" begin
        p = ArrayPartition([1, 0.0, 0.0], [4 5.0; 6 7])
        q = ArrayPartition([0.0, 1.0, 0.0], [4 8.0; 3 7.5])

        X = log(M, p, q)
        for m in [ParallelTransport(), SchildsLadderTransport(), PoleLadderTransport()]
            Y = vector_transport_to(M, p, X, q, m)
            Z1 = vector_transport_to(
                M.manifolds[1],
                submanifold_component.([p, X, q], Ref(1))...,
                m,
            )
            Z2 = vector_transport_to(
                M.manifolds[2],
                submanifold_component.([p, X, q], Ref(2))...,
                m,
            )
            Z = ArrayPartition(Z1, Z2)
            @test isapprox(M, q, Y, Z)
            Y2 = allocate(M, Y)
            vector_transport_to!(M, Y2, p, X, q, m)
            @test isapprox(M, q, Y2, Z)
        end
        for m in [ParallelTransport(), SchildsLadderTransport(), PoleLadderTransport()]
            Y = vector_transport_direction(M, p, X, X, m)
            Z1 = vector_transport_direction(
                M.manifolds[1],
                submanifold_component.([p, X, X], Ref(1))...,
                m,
            )
            Z2 = vector_transport_direction(
                M.manifolds[2],
                submanifold_component.([p, X, X], Ref(2))...,
                m,
            )
            Z = ArrayPartition(Z1, Z2)
            @test isapprox(M, q, Y, Z)
        end
    end

    @testset "Parallel transport" begin
        p = ArrayPartition([1, 0.0, 0.0], [4 5.0; 6 7])
        q = ArrayPartition([0.0, 1.0, 0.0], [4 8.0; 3 7.5])

        X = log(M, p, q)
        # to
        Y = parallel_transport_to(M, p, X, q)
        Z1 = parallel_transport_to(
            M.manifolds[1],
            submanifold_component.([p, X, q], Ref(1))...,
        )
        Z2 = parallel_transport_to(
            M.manifolds[2],
            submanifold_component.([p, X, q], Ref(2))...,
        )
        Z = ArrayPartition(Z1, Z2)
        @test isapprox(M, q, Y, Z)
        Ym = allocate(Y)
        parallel_transport_to!(M, Ym, p, X, q)
        @test isapprox(M, q, Y, Z)

        # direction
        Y = parallel_transport_direction(M, p, X, X)
        Z1 = parallel_transport_direction(
            M.manifolds[1],
            submanifold_component.([p, X, X], Ref(1))...,
        )
        Z2 = parallel_transport_direction(
            M.manifolds[2],
            submanifold_component.([p, X, X], Ref(2))...,
        )
        Z = ArrayPartition(Z1, Z2)
        @test isapprox(M, q, Y, Z)
        Ym = allocate(Y)
        parallel_transport_direction!(M, Ym, p, X, X)
        @test isapprox(M, q, Ym, Z)
    end

    @testset "ArrayPartition" begin
        @test submanifold_component(M, p1, 1) === p1.x[1]
        @test submanifold_component(M, p1, Val(1)) === p1.x[1]
        @test submanifold_component(p1, 1) === p1.x[1]
        @test submanifold_component(p1, Val(1)) === p1.x[1]
        @test submanifold_components(M, p1) === p1.x
        @test submanifold_components(p1) === p1.x
    end

    @testset "vee/hat" begin
        X = [0.1, 0.2, 0.3, -1.0, 2.0, -3.0]

        Xc = hat(M, p1, X)
        X2 = vee(M, p1, Xc)
        @test isapprox(X, X2)
    end

    @testset "get_coordinates" begin
        # make sure `get_coordinates` does not return an `ArrayPartition`
        p1 = ArrayPartition([0.0, 1.0, 0.0], [0.0 0.0; 1.0 2.0])
        X1 = ArrayPartition([1.0, 0.0, -1.0], [1.0 0.0; 2.0 1.0])
        Tp1M = TangentSpace(M, p1)
        c = get_coordinates(Tp1M, p1, X1, DefaultOrthonormalBasis())
        @test c isa Vector

        p1 = ArrayPartition([0.0, 1.0, 0.0], [0.0 0.0; 3.0 20])
        X1ap = ArrayPartition([1.0, 0.0, -1.0], [1.0 0.0; 0.0 3.0])
        Tp1M = TangentSpace(M, p1)
        cap = get_coordinates(Tp1M, p1, X1ap, DefaultOrthonormalBasis())
        @test cap isa Vector
    end

    @testset "Basis printing" begin
        p = ArrayPartition([1.0, 0.0, 0.0], [1.0 0.0; 1.0 2.0])
        B = DefaultOrthonormalBasis()
        Bc = get_basis(M, p, B)
        Bc_components_s = sprint.(show, "text/plain", Bc.data.parts)
        @test sprint(show, "text/plain", Bc) == """
        $(typeof(B)) for a product manifold
        Basis for component 1:
        $(Bc_components_s[1])
        Basis for component 2:
        $(Bc_components_s[2])
        """
    end

    @testset "Basis-related errors" begin
        a = ArrayPartition([1.0, 0.0, 0.0], [0.0 0.0; 0.0 0.0])
        B = CachedBasis(DefaultOrthonormalBasis(), ProductBasisData(([],)))
        @test_throws AssertionError get_vector!(
            M,
            a,
            ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0]),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], # this is one element too long, hence assertion error
            B,
        )
        @test_throws MethodError get_vector!(
            M,
            a,
            ArrayPartition([1.0, 0.0, 0.0], [0.0 0.0; 0.0 0.0]),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            B, # empty elements yield a submanifold MethodError
        )
    end

    @testset "Basis -- other" begin
        p1 = ArrayPartition([1.0, 0.0, 0.0], [1.0 0.0; 1.0 2.0])
        X1 = ArrayPartition([0.0, 2.0, -1.0], [1.0 0.0; 2.0 1.0])
        B = DefaultOrthonormalBasis()
        Bc = get_basis(M, p1, B)

        BcO = get_basis(M, p1, DiagonalizingOrthonormalBasis(X1))
        @test BcO.data isa ProductBasisData

        for basis in
            [DefaultOrthonormalBasis(), get_basis(M, p1, DefaultOrthonormalBasis())]
            @test length(get_vectors(M, p1, get_basis(M, p1, basis))) == 6
            X1c = get_coordinates(M, p1, X1, basis)
            @test isapprox(X1c, [2.0, -1.0, 1.0, 2.0, 0.0, 1.0])
            Y1c = allocate(M, X1c)
            get_coordinates!(M, Y1c, p1, X1, basis)
            @test isapprox(X1c, Y1c)
            @test isapprox(get_vector(M, p1, X1c, basis), X1)
            Z1 = allocate(M, X1)
            get_vector!(M, Z1, p1, X1c, basis)
            @test isapprox(X1, X1)
        end
    end

    @testset "allocation promotion" begin
        M2c = DefaultManifold(2; field = ℂ)
        Mc = ProductManifold(M1, M2c)
        @test ManifoldsBase.allocation_promotion_function(Mc, get_vector, ()) === complex
        @test ManifoldsBase.allocation_promotion_function(M, get_vector, ()) === identity
    end

    @testset "Riemann tensor" begin
        p = ArrayPartition([0.0, 1.0, 0.0], [2.0 3.0; 0.0 0.0])
        X = ArrayPartition([1.0, 0.0, 0.0], [2.0 3.0; 0.0 0.0])
        Y = ArrayPartition([0.0, 0.0, 3.0], [-2.0 3.0; 0.0 0.0])
        Z = ArrayPartition([-1.0, 0.0, 2.0], [2.0 -3.0; 0.0 0.0])
        Xresult = ArrayPartition([6.0, 0.0, 3.0], [0.0 0.0; 0.0 0.0])
        @test isapprox(riemann_tensor(M, p, X, Y, Z), Xresult)
        Xresult2 = allocate(Xresult)
        riemann_tensor!(M, Xresult2, p, X, Y, Z)
        @test isapprox(Xresult2, Xresult)
    end

    @testset "× constructors" begin
        r1 = ExponentialRetraction()
        r2 = ProjectionRetraction()
        s1 = r1 × r2
        @test s1 == ProductRetraction(r1, r2)
        @test "$(s1)" == "ProductRetraction($(r1), $(r2))"
        @test s1 × r2 == ProductRetraction(r1, r2, r2)
        @test r2 × s1 == ProductRetraction(r2, r1, r2)
        @test r1 × r1 × r1 == ProductRetraction(r1, r1, r1)

        ir1 = LogarithmicInverseRetraction()
        ir2 = ProjectionInverseRetraction()
        is1 = ir1 × ir2
        @test is1 == InverseProductRetraction(ir1, ir2)
        @test "$(is1)" == "InverseProductRetraction($(ir1), $(ir2))"
        @test is1 × ir2 == InverseProductRetraction(ir1, ir2, ir2)
        @test ir2 × is1 == InverseProductRetraction(ir2, ir1, ir2)
        @test ir1 × ir1 × ir1 == InverseProductRetraction(ir1, ir1, ir1)

        tr1 = ParallelTransport()
        tr2 = ProjectionTransport()
        ts1 = tr1 × tr2
        @test ts1 == ProductVectorTransport(tr1, tr2)
        @test "$(ts1)" == "ProductVectorTransport($(tr1), $(tr2))"
        @test ts1 × tr2 == ProductVectorTransport(tr1, tr2, tr2)
        @test tr2 × ts1 == ProductVectorTransport(tr2, tr1, tr2)
        @test tr1 × tr1 × tr1 == ProductVectorTransport(tr1, tr1, tr1)
    end


end
