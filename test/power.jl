using Test
using ManifoldsBase
using ManifoldsBase: AbstractNumbers, ℝ, ℂ, NestedReplacingPowerRepresentation
using StaticArrays
using LinearAlgebra

power_array_wrapper(::Type{NestedPowerRepresentation}, ::Int) = identity
power_array_wrapper(::Type{NestedReplacingPowerRepresentation}, i::Int) = SVector{i}

function ManifoldsBase.allocate(
    ::ManifoldsBase.PowerManifoldNestedReplacing,
    x::AbstractArray{<:SArray},
)
    return similar(x)
end

struct TestArrayRepresentation <: AbstractPowerRepresentation end

@testset "Power Manifold" begin

    @testset "Power Manifold with a test representation" begin
        M = ManifoldsBase.DefaultManifold(3)
        N = PowerManifold(M, TestArrayRepresentation(), 2)
        O = PowerManifold(N, TestArrayRepresentation(), 3) # joins instead of nesting.
        @test repr(O) ==
              "PowerManifold(DefaultManifold(3; field = ℝ), TestArrayRepresentation(), 2, 3)"
        p = zeros(6)
        X = zeros(6)
        @test ManifoldsBase.check_power_size(N, p) === nothing
        @test ManifoldsBase.check_power_size(O, p) isa DomainError
        @test ManifoldsBase.check_power_size(N, p, X) === nothing
        @test ManifoldsBase.check_power_size(O, p, X) isa DomainError
        @test default_inverse_retraction_method(M) == default_inverse_retraction_method(N)
        @test default_retraction_method(M) == default_retraction_method(N)
        @test default_vector_transport_method(M) == default_vector_transport_method(N)
    end

    @testset "PowerManifold and allocation with empty representation size" begin
        M = ManifoldsBase.DefaultManifold()
        N = PowerManifold(M, NestedPowerRepresentation(), 2)
        p = [1, 1]
        X = [2, 2]
        # check - though only because this function exists for avoiding ambiguities.
        cm = ManifoldsBase.allocate_result(N, get_coordinates, p, X, DefaultBasis())
        @test size(X) == size(cm)
    end

    @testset "PowerManifoldNested with mutable element" begin
        M = ManifoldsBase.DefaultManifold(2, 2)
        N = PowerManifold(M, NestedPowerRepresentation(), 2)
        p = [UpperTriangular([1 2; 2 1]), UpperTriangular([1 2; 2 1])]
        q = [UpperTriangular([2 3; 3 2]), UpperTriangular([1 2; 2 1])]
        @test typeof(log(N, p, q)) === typeof(p)
    end

    @testset "PowerManifoldNestedReplacing with SArray element" begin
        M = ManifoldsBase.DefaultManifold(2, 2)
        N = PowerManifold(M, NestedReplacingPowerRepresentation(), 2)
        p = [SMatrix{2,2,Float64}([i i+1; i-1 i-2]) for i in 1:2]
        allocate(M, p) isa Vector{SMatrix{2,2,Float64,4}}
    end

    for PowerRepr in [NestedPowerRepresentation, NestedReplacingPowerRepresentation]
        @testset "PowerManifold with $(PowerRepr)" begin
            M = ManifoldsBase.DefaultManifold(3)
            for pow_size in ((2,), (2, 3))
                pow_cart = CartesianIndices(pow_size)
                N = PowerManifold(M, PowerRepr(), pow_size...)
                wrapper_2 = power_array_wrapper(PowerRepr, 2)
                wrapper_3 = power_array_wrapper(PowerRepr, 3)

                p = wrapper_3.([i[1] == 1 ? zeros(3) : ones(3) for i in pow_cart])
                q = wrapper_3.([i[1] == 1 ? ones(3) : zeros(3) for i in pow_cart])

                pow_string = join(pow_size, ", ")

                @testset "Constructors" begin
                    @test repr(N) ==
                          "PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), $pow_string)"
                    # add to type
                    @test repr(N^3) ==
                          "PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), $pow_string, 3)"
                    # add to type
                    @test repr(PowerManifold(N, 3)) ==
                          "PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), $pow_string, 3)"
                    # nest
                    @test repr(PowerManifold(N, PowerRepr(), 3)) ==
                          "PowerManifold(PowerManifold(DefaultManifold(3; field = ℝ), $(PowerRepr)(), $pow_string), $(PowerRepr)(), 3)"
                end

                @test N^3 == PowerManifold(M, PowerRepr(), pow_size..., 3)

                if pow_size == (2,)

                    @test ManifoldsBase.get_iterator(N^3) ==
                          Base.product(Base.OneTo(2), Base.OneTo(3))
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
                end
                @testset "specific functions" begin
                    @test distance(N, p, q) == sqrt(sum(distance.(Ref(M), p, q) .^ 2))
                    @test exp(N, p, q) == p .+ q

                    @test retract(N, p, q) == p .+ q
                    @test retract(N, p, q, ExponentialRetraction()) == p .+ q
                    r = allocate(p)
                    @test retract!(N, r, p, q, ExponentialRetraction()) == p .+ q
                    @test r == p .+ q
                    @test inverse_retract(N, p, r) == q
                    @test inverse_retract(N, p, r, LogarithmicInverseRetraction()) == q
                    X = allocate(p)
                    @test inverse_retract!(N, X, p, r, LogarithmicInverseRetraction()) == q
                    @test X == q

                    if pow_size == (2,)
                        @test ManifoldsBase.get_iterator(N) == Base.OneTo(2)
                    else
                        @test ManifoldsBase.get_iterator(N) ==
                              Base.product(Base.OneTo(2), Base.OneTo(3))
                    end
                    @test injectivity_radius(N) == injectivity_radius(M)
                    @test injectivity_radius(N, ExponentialRetraction()) ==
                          injectivity_radius(M)
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
                    @test manifold_dimension(N) == prod(pow_size) * manifold_dimension(M)
                    @test mid_point(N, p, q) == mid_point.(Ref(M), p, q)
                    @test sqrt(inner(N, p, q, q)) ≈ norm(N, p, q)
                    @test project(N, p) == p
                    @test project(N, p, q) == q
                    @test power_dimensions(N) == pow_size
                    @test power_dimensions(N^3) == (pow_size..., 3)
                    m = ParallelTransport()
                    @test vector_transport_to(N, p, p, q) == p
                    @test vector_transport_to(N, p, p, q, m) == p
                    @test parallel_transport_to(N, p, p, q) == p
                    q2 = [zeros(3) for _ in pow_cart]
                    vector_transport_to!(N, q2, p, q, p)
                    @test q2 == q
                    q2 = [zeros(3) for _ in pow_cart]
                    parallel_transport_to!(N, q2, p, q, p)
                    @test q2 == q
                    @test vector_transport_direction(N, p, p, q) == p
                    @test vector_transport_direction(N, p, p, q, m) == p
                    @test parallel_transport_direction(N, p, p, q) == p
                    q2 = [zeros(3) for _ in pow_cart]
                    vector_transport_direction!(N, q2, p, q, p)
                    @test q2 == q
                    q2 = [zeros(3) for _ in pow_cart]
                    parallel_transport_direction!(N, q2, p, q, p)
                    @test q2 == q
                    if pow_size == (2,)
                        @test p[N, 1] == p[1]
                        p[N, 1] = 2 .* ones(3)
                        @test p[N, 1] == 2 .* ones(3)
                        @test p[N, [2, 1]] == [p[2], p[1]]
                        if PowerRepr == NestedPowerRepresentation
                            @test view(p, N, 2) isa SubArray
                            @test view(p, N, 2) == p[2]
                        end
                    end
                end
                @testset "Basis, coordinates & vector" begin
                    v = get_coordinates(N, p, q, DefaultBasis())
                    @test v == reduce(vcat, q)
                    @test get_vector(N, p, v, DefaultBasis()) == q
                    B = get_basis(N, p, DefaultBasis())
                    @test get_coordinates(N, p, q, B) == v
                    # the method tested below should not be used but it prevents ambiguities from occurring
                    # and the test is here to make coverage happy
                    @test ManifoldsBase.allocate_result(N, get_coordinates, p, q, B) isa
                          Vector
                    v2 = similar(v)
                    get_coordinates!(N, v2, p, q, B)
                    @test v2 == v
                    @test get_coordinates(N, p, q, DefaultOrthonormalBasis()) == v
                    v3 = similar(v2)
                    @test get_coordinates!(N, v3, p, q, DefaultOrthonormalBasis()) == v
                    @test v3 == v
                    @test B.data.bases[1].data == get_basis(M, p[1], DefaultBasis()).data
                    @test B.data.bases[2].data == get_basis(M, p[2], DefaultBasis()).data
                    @test get_vector(N, p, v, B) == q
                    q2 = similar.(q)
                    @test get_vector!(N, q2, p, v, B) == q
                    @test q2 == q
                    q3 = similar.(q)
                    @test get_vector(N, p, v, DefaultOrthonormalBasis()) == q
                    @test get_vector!(N, q3, p, v, DefaultOrthonormalBasis()) == q
                    @test q3 == q
                    @test zero_vector(N, p) == [zeros(3) for _ in pow_cart]
                    B2 = DiagonalizingOrthonormalBasis([ones(3) for _ in pow_cart])
                    B3 = get_basis(N, p, B2)
                    if pow_size == (2,)
                        @test sprint(show, "text/plain", B) ==
                              """$(DefaultBasis()) for a power manifold
          Basis for component (1,):
          $(sprint(show, "text/plain", B.data.bases[1]))
          Basis for component (2,):
          $(sprint(show, "text/plain", B.data.bases[2]))
          """
                    end

                    bv = get_vectors(N, p, B)
                    @test length(bv) == number_of_coordinates(N, B)
                    for i in 1:number_of_coordinates(N, B)
                        @test bv[i] == get_vector(
                            N,
                            p,
                            ManifoldsBase._euclidean_basis_vector(v, i),
                            B,
                        )
                    end
                end
                @testset "zero argument allocation" begin
                    @test size(ManifoldsBase.allocate_result(N, rand)) == size(p)
                end
            end

            @testset "Zero index manifold" begin
                Mzero = ManifoldsBase.DefaultManifold()
                N = PowerManifold(Mzero, PowerRepr(), 3)
                p = [1.0, 2.0, 3.0]

                @test p[N, 1] == 1.0
                @test zero_vector(N, p) == zero(p)
            end
        end
    end
end
