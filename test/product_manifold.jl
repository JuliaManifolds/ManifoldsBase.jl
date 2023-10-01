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

    p1 = ArrayPartition([1, 0.0, 0.0], [4 5.0; 6 7])
    p2 = ArrayPartition([0.0, 1.0, 0.0], [4 8.0; 3 7.5])

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
    @test injectivity_radius(M, p1, ExponentialRetraction()) ≈ π

    @test ManifoldsBase.number_of_components(M) == 2
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
    end

    @testset "allocate on PowerManifold of ProductManifold" begin
        q = allocate([p1])
        @test q[1] isa ArrayPartition
        @test q[1].x[1] isa Vector
    end

    p1 = ArrayPartition([1, 0.0, 0.0], [4 5.0; 6 7])
    p2 = ArrayPartition([0.0, 1.0, 0.0], [4 8.0; 3 7.5])

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

        # errors
        p3 = ArrayPartition([3.0, 4.0, 5.0], [2.0, 5.0], [3.0, 2.0])
        @test_throws DimensionMismatch p1 .+ p3
        @test_throws DimensionMismatch p1 .= p3
    end

    @testset "CompositeManifoldError" begin
        Mpr = ProductManifold(TestSphere(2), TestSphere(2))
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

    @testset "arithmetic" begin
        @test isapprox(M, p1 + p2, ArrayPartition([1.0, 1.0, 0.0], [8.0 13.0; 9.0 14.5]))
        @test isapprox(M, p1 - p2, ArrayPartition([1.0, -1.0, 0.0], [0.0 -3.0; 3.0 -0.5]))
        @test isapprox(M, -p1, ArrayPartition([-1.0, -0.0, -0.0], [-4.0 -5.0; -6.0 -7.0]))
        @test isapprox(M, p1 * 2, ArrayPartition([2.0, 0.0, 0.0], [8.0 10.0; 12.0 14.0]))
        @test isapprox(M, 2 * p1, ArrayPartition([2.0, 0.0, 0.0], [8.0 10.0; 12.0 14.0]))
        @test isapprox(M, p1 / 2, ArrayPartition([0.5, 0.0, 0.0], [2.0 2.5; 3.0 3.5]))
    end

    @testset "Show methods" begin
        M2 = ProductManifold(M1, M1, M2, M2)
        @test sprint(show, M2) == "ProductManifold($(M1), $(M1), $(M2), $(M2))"
        withenv("LINES" => 10, "COLUMNS" => 100) do
            @test sprint(show, "text/plain", ProductManifold(M1)) ==
                  "ProductManifold with 1 submanifold:\n $(M1)"
            @test sprint(show, "text/plain", M2) ==
                  "ProductManifold with 4 submanifolds:\n $(M1)\n $(M1)\n $(M2)\n $(M2)"
            return nothing
        end
        withenv("LINES" => 7, "COLUMNS" => 100) do
            @test sprint(show, "text/plain", M2) ==
                  "ProductManifold with 4 submanifolds:\n $(M1)\n ⋮\n $(M2)"
            return nothing
        end

        @test sprint(show, "text/plain", ProductManifold(M, M)) == """
        ProductManifold with 2 submanifolds:
         ProductManifold(Sphere(2, ℝ), Euclidean(2; field=ℝ))
         ProductManifold(Sphere(2, ℝ), Euclidean(2; field=ℝ))"""
    end

    @testset "product vector transport" begin
        X = log(M, p1, p2)
        m = ProductVectorTransport(ParallelTransport(), ParallelTransport())
        Y = vector_transport_to(M, p1, X, p2, m)
        Z = -log(M, p2, p1)
        @test isapprox(M, p2, Y, Z)
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
        Tp1M = TangentSpaceAtPoint(M, p1)
        c = get_coordinates(Tp1M, p1, X1, DefaultOrthonormalBasis())
        @test c isa Vector

        p1 = ArrayPartition([0.0, 1.0, 0.0], [0.0 0.0; 3.0 20])
        X1ap = ArrayPartition([1.0, 0.0, -1.0], [1.0 0.0; 0.0 3.0])
        Tp1M = TangentSpaceAtPoint(M, p1)
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
            [1.0, 2.0, 3.0, 4.0, 5.0], # this is one element too long, hence assertion error
            B,
        )
        @test_throws MethodError get_vector!(
            M,
            a,
            ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0]),
            [1.0, 2.0, 3.0, 4.0],
            B, # empty elements yield a submanifold MethodError
        )
    end

    @testset "allocation promotion" begin
        M2c = DefaultManifold(2; field = ℂ)
        Mc = ProductManifold(M1, M2c)
        @test ManifoldsBase.allocation_promotion_function(Mc, get_vector, ()) === complex
        @test ManifoldsBase.allocation_promotion_function(M, get_vector, ()) === identity
    end

    @testset "default retraction, inverse retraction and VT" begin
        Mstb = ProductManifold(M1, TangentBundle(M1))
        T_p_ap = ArrayPartition{
            Float64,
            Tuple{
                Matrix{Float64},
                ArrayPartition{Float64,Tuple{Matrix{Float64},Matrix{Float64}}},
            },
        }
        @test ManifoldsBase.default_retraction_method(Mstb) === ProductRetraction(
            ExponentialRetraction(),
            ManifoldsBase.FiberBundleProductRetraction(),
        )
        @test ManifoldsBase.default_retraction_method(Mstb, T_p_ap) === ProductRetraction(
            ExponentialRetraction(),
            ManifoldsBase.FiberBundleProductRetraction(),
        )

        @test ManifoldsBase.default_inverse_retraction_method(Mstb) ===
              ManifoldsBase.InverseProductRetraction(
            LogarithmicInverseRetraction(),
            ManifoldsBase.FiberBundleInverseProductRetraction(),
        )
        @test ManifoldsBase.default_inverse_retraction_method(Mstb, T_p_ap) ===
              ManifoldsBase.InverseProductRetraction(
            LogarithmicInverseRetraction(),
            ManifoldsBase.FiberBundleInverseProductRetraction(),
        )

        @test ManifoldsBase.default_vector_transport_method(Mstb) ===
              ProductVectorTransport(
            ParallelTransport(),
            ManifoldsBase.FiberBundleProductVectorTransport(
                ParallelTransport(),
                ParallelTransport(),
            ),
        )
        @test ManifoldsBase.default_vector_transport_method(Mstb, T_p_ap) ===
              ProductVectorTransport(
            ParallelTransport(),
            ManifoldsBase.FiberBundleProductVectorTransport(
                ParallelTransport(),
                ParallelTransport(),
            ),
        )
        @test ManifoldsBase.default_vector_transport_method(Mstb, T_p_ap) ===
              ProductVectorTransport(
            ParallelTransport(),
            ManifoldsBase.FiberBundleProductVectorTransport(
                ParallelTransport(),
                ParallelTransport(),
            ),
        )
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

end
