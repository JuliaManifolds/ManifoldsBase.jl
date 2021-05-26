using Test
using ManifoldsBase
using ManifoldsBase: AbstractNumbers, ℝ, ℂ, NestedReplacingPowerRepresentation
using StaticArrays

struct DummyPowerRepresentation <: AbstractPowerRepresentation end
struct DummyDecorator{TM<:AbstractManifold{ManifoldsBase.ℝ}} <:
       AbstractDecoratorManifold{ManifoldsBase.ℝ}
    manifold::TM
end

power_array_wrapper(::Type{NestedPowerRepresentation}, ::Int) = identity
power_array_wrapper(::Type{NestedReplacingPowerRepresentation}, i::Int) = SVector{i}

@testset "Power AbstractManifold" begin
    M = ManifoldsBase.DefaultManifold(3)
    for PowerRepr in [NestedPowerRepresentation, NestedReplacingPowerRepresentation]
        N = PowerManifold(M, PowerRepr(), 2)
        wrapper_2 = power_array_wrapper(PowerRepr, 2)
        wrapper_3 = power_array_wrapper(PowerRepr, 3)
        p = wrapper_3.([zeros(3), ones(3)])
        q = wrapper_3.([ones(3), zeros(3)])

        @testset "Constructors" begin
            @test repr(N) ==
                  "PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), 2)"
            # add to type
            @test repr(N^3) ==
                  "PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), 2, 3)"
            # add to type
            @test repr(PowerManifold(N, 3)) ==
                  "PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), 2, 3)"
            # switch type
            @test repr(PowerManifold(N, DummyPowerRepresentation(), 3)) ==
                  "PowerManifold(DefaultManifold(3; field = ℝ), DummyPowerRepresentation(), 2, 3)"
            # nest
            @test repr(PowerManifold(N, PowerRepr(), 3)) ==
                  "PowerManifold(PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), 2), $(PowerRepr)(), 3)"
        end

        @test N^3 == PowerManifold(M, PowerRepr(), 2, 3)
        @test ManifoldsBase.get_iterator(N^3) == Base.product(Base.OneTo(2), Base.OneTo(3))

        @testset "point/tangent checks" begin
            pE1 = [wrapper_3(zeros(3)), wrapper_2(ones(2))] # one component wrong
            pE2 = wrapper_2.([zeros(2), ones(2)]) # both wrong
            @test is_point(N, p)
            @test is_point(N, p, true)
            @test !is_point(N, pE1)
            @test !is_point(N, pE2)
            @test_throws ComponentManifoldError is_point(N, pE1, true)
            @test_throws CompositeManifoldError is_point(N, pE2, true)
            # tangent - test base
            @test is_vector(N, p, p)
            @test is_vector(N, p, p, true)
            @test !is_vector(N, pE1, p)
            @test !is_vector(N, pE2, p)
            # tangents - with proper base
            @test is_vector(N, p, p, true)
            @test !is_vector(N, p, pE1)
            @test !is_vector(N, p, pE2)
            @test_throws ComponentManifoldError is_vector(N, p, pE1, true)
            @test_throws CompositeManifoldError is_vector(N, p, pE2, true)
        end
        @testset "specific functions" begin
            @test distance(N, p, q) == sqrt(sum(distance.(Ref(M), p, q) .^ 2))
            @test exp(N, p, q) == p .+ q
            @test retract(N, p, q) == p .+ q
            @test ManifoldsBase.get_iterator(N) == Base.OneTo(2)
            @test injectivity_radius(N) == injectivity_radius(M)
            @test injectivity_radius(N, p) == injectivity_radius(M, p)
            p2 = allocate(p)
            copyto!(N, p2, p)
            @test !(p === p2)
            @test p == p2
            q2 = allocate(q)
            copyto!(N, q2, p, q)
            @test !(q === q2)
            @test q == q2
            @test inner(N, p, q, q) == sum(inner.(Ref(M), p, q, q))
            @test isapprox(N, p, q) == (all(isapprox.(Ref(M), p, q)))
            @test isapprox(N, p, p) == (all(isapprox.(Ref(M), p, p)))
            @test isapprox(N, p, q, p) == (all(isapprox.(Ref(M), p, q, p)))
            @test isapprox(N, p, p, p) == (all(isapprox.(Ref(M), p, p, p)))
            @test log(N, p, q) == q .- p
            @test inverse_retract(N, p, q) == q .- p
            @test manifold_dimension(N) == 2 * manifold_dimension(M)
            @test mid_point(N, p, q) == mid_point.(Ref(M), p, q)
            @test sqrt(inner(N, p, q, q)) == norm(N, p, q)
            @test project(N, p) == p
            @test project(N, p, q) == q
            @test power_dimensions(N) == (2,)
            @test power_dimensions(N^3) == (2, 3)
            m = ParallelTransport()
            @test vector_transport_to(N, p, p, q) == p
            @test vector_transport_to(N, p, p, q, m) == p
            @test vector_transport_to(N, p, p, q, PowerVectorTransport(m)) == p
            q2 = [zeros(3), zeros(3)]
            vector_transport_to!(N, q2, p, q, p)
            @test q2 == q
            @test vector_transport_direction(N, p, p, q) == p
            @test vector_transport_direction(N, p, p, q, m) == p
            @test vector_transport_direction(N, p, p, q, PowerVectorTransport(m)) == p
            q2 = [zeros(3), zeros(3)]
            vector_transport_direction!(N, q2, p, q, p)
            @test q2 == q
            @test p[N, 1] == p[1]
            p[N, 1] = 2 .* ones(3)
            @test p[N, 1] == 2 .* ones(3)
            @test p[N, [2, 1]] == [p[2], p[1]]
            if PowerRepr == NestedPowerRepresentation
                @test view(p, N, 2) isa SubArray
                @test view(p, N, 2) == p[2]
            end
        end
        @testset "Basis, coordinates & vector" begin
            v = get_coordinates(N, p, q, DefaultBasis())
            @test v == [q[1]..., q[2]...]
            @test get_vector(N, p, v, DefaultBasis()) == q
            B = get_basis(N, p, DefaultBasis())
            @test get_coordinates(N, p, q, B) == v
            # the method tested below should not be used but it prevents ambiguities from occurring
            # and the test is here to make coverage happy
            @test ManifoldsBase.allocate_result(N, get_coordinates, p, q, B) isa Vector
            v2 = zeros(size(v))
            get_coordinates!(N, v2, p, q, B)
            @test v2 == v
            @test get_coordinates(N, p, q, DefaultOrthonormalBasis()) == v
            @test B.data.bases[1].data == get_basis(M, p[1], DefaultBasis()).data
            @test B.data.bases[2].data == get_basis(M, p[2], DefaultBasis()).data
            @test get_vector(N, p, v, B) == q
            @test zero_vector(N, p) == [zeros(3), zeros(3)]
            B2 = DiagonalizingOrthonormalBasis([ones(3), ones(3)])
            B3 = get_basis(N, p, B2)
            @test sprint(show, "text/plain", B) == """$(DefaultBasis()) for a power manifold
            Basis for component (1,):
            $(sprint(show, "text/plain", B.data.bases[1]))
            Basis for component (2,):
            $(sprint(show, "text/plain", B.data.bases[2]))
            """
        end

        @testset "Zero index manifold" begin
            Mzero = ManifoldsBase.DefaultManifold()
            N = PowerManifold(Mzero, PowerRepr(), 3)
            p = [1.0, 2.0, 3.0]

            @test p[N, 1] == 1.0
            @test zero_vector(N, p) == zero(p)
        end
        @testset "Decorator passthrough for getindex" begin
            Mzero = ManifoldsBase.DefaultManifold()
            N = PowerManifold(Mzero, PowerRepr(), 3)
            p = [1.0, 2.0, 3.0]
            DN = DummyDecorator(N)
            @test p[DN, 1] == p[N, 1]
        end
    end
end
