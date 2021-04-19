using ManifoldsBase
using LinearAlgebra
using Test

struct CustomValidationManifoldRetraction <: ManifoldsBase.AbstractRetractionMethod end

function ManifoldsBase.injectivity_radius(
    ::ManifoldsBase.DefaultManifold,
    ::CustomValidationManifoldRetraction,
)
    return 10.0
end
function ManifoldsBase.injectivity_radius(
    ::ManifoldsBase.DefaultManifold,
    p,
    ::CustomValidationManifoldRetraction,
)
    return 11.0
end

@testset "Validation manifold" begin
    M = ManifoldsBase.DefaultManifold(3)
    A = ValidationManifold(M)
    x = [1.0, 0.0, 0.0]
    y = 1 / sqrt(2) * [1.0, 1.0, 0.0]
    z = [0.0, 1.0, 0.0]
    v = log(M, x, y)
    x2 = ValidationMPoint(x)
    y2 = ValidationMPoint(y)
    v2 = log(A, x, y) # auto convert
    y2 = exp(A, x, v2)
    w = log(M, x, z)
    w2 = log(A, x, z; atol = 10^(-15))
    @testset "Types and Conversion" begin
        @test convert(typeof(M), A) == M
        @test convert(typeof(A), M) == A
        @test base_manifold(A) == M
        @test base_manifold(base_manifold(A)) == base_manifold(A)
        @test ManifoldsBase.representation_size(A) == ManifoldsBase.representation_size(M)
        @test ManifoldsBase.representation_size(A) == (3,)
        @test manifold_dimension(A) == manifold_dimension(M)
        @test manifold_dimension(A) == 3
        for T in [ValidationMPoint, ValidationTVector, ValidationCoTVector]
            p = T(x)
            @test convert(typeof(x), p) == x
            @test convert(typeof(p), y) == T(y)
            @test number_eltype(typeof(p)) == eltype(x)
            @test number_eltype(p) == eltype(x)
            @test typeof(allocate(p)) == typeof(p)
            @test typeof(allocate(p, eltype(x))) == typeof(p)
            @test allocate(p) isa T
            @test allocate(p, Float32) isa T
            @test number_eltype(allocate(p, Float32)) == Float32
            @test similar(p) isa T
            @test similar(p, Float32) isa T
            @test number_eltype(similar(p, Float32)) == Float32
            q = allocate(p)
            if T == ValidationMPoint
                copyto!(A, q, p)
            else
                copyto!(A, q, ValidationMPoint(x), p) # generate base point “on the fly”
            end
            @test isapprox(A, q, p)
            @test ManifoldsBase.array_value(p) == x
            @test ManifoldsBase.array_value(x) == x
        end
    end
    @testset "Vector functions" begin
        for T in [ValidationTVector, ValidationCoTVector]
            a = T(v)
            b = T(w)
            @test isapprox(A, a + b, T(v + w))
            @test isapprox(A, (a - b), T(v - w))
            @test isapprox(A, -b, T(-w))
            @test isapprox(A, 2 * a, T(2 .* v))
        end
    end
    @testset "Manifold functions" begin
        @test manifold_dimension(A) == manifold_dimension(M)
        @test isapprox(y2.value, y)
        @test distance(A, x, y) == distance(M, x, y)
        @test norm(A, x, v) == norm(M, x, v)
        @test inner(A, x, v2, w2; atol = 10^(-15)) == inner(M, x, v, w)
        @test isapprox(A, x2, y2) == isapprox(M, x, y)
        @test isapprox(A, x, y) == isapprox(A, x2, y2)
        @test isapprox(A, x, v2, v2) == isapprox(M, x, v, v)
        v2s = similar(v2)
        project!(A, v2s, x2, v2)
        @test isapprox(A, v2, v2s)
        y2s = similar(y2)
        exp!(A, y2s, x2, v2)
        @test isapprox(A, y2s, y2)
        log!(A, v2s, x, y)
        @test mid_point(A, x, y) == 0.5 * (x .+ y)
        mp = similar(x)
        mid_point!(A, mp, x, y)
        @test mp == 0.5 * (x .+ y)
        @test isapprox(A, x, v2s, v2)
        @test isapprox(A, exp(A, x, v), y2)
        @test isapprox(A, zero_tangent_vector(A, x), zero_tangent_vector(M, x))
        vector_transport_to!(A, v2s, x2, v2, y2)
        @test isapprox(A, x2, v2, v2s)
        vector_transport_to!(A, v2s, x2, v2, y2, ManifoldsBase.SchildsLadderTransport())
        @test isapprox(A, x2, v2, v2s)
        vector_transport_to!(A, v2s, x2, v2, y, ManifoldsBase.PoleLadderTransport())
        @test isapprox(A, x2, v2, v2s)
        vector_transport_to!(A, v2s, x2, v2, y2, ManifoldsBase.ProjectionTransport())
        @test isapprox(A, x2, v2, v2s)
        zero_tangent_vector!(A, v2s, x)
        @test isapprox(A, x, v2s, zero_tangent_vector(M, x))
        c = [x2]
        v3 = similar(v2)
        @test isapprox(
            A,
            x2,
            v2,
            vector_transport_along!(A, v3, x2, v2, c, ParallelTransport()),
        )
        @test isapprox(
            A,
            x2,
            v2,
            vector_transport_along(A, x2, v2, c, ManifoldsBase.ProjectionTransport()),
        )
        @test injectivity_radius(A) == Inf
        @test injectivity_radius(A, x) == Inf
        @test injectivity_radius(A, ManifoldsBase.ExponentialRetraction()) == Inf
        @test injectivity_radius(A, x, ManifoldsBase.ExponentialRetraction()) == Inf
        @test injectivity_radius(A, CustomValidationManifoldRetraction()) == 10
        @test injectivity_radius(A, x, CustomValidationManifoldRetraction()) == 11
    end

    @testset "ValidationManifold basis" begin
        b = [Matrix(I, 3, 3)[:, i] for i in 1:3]
        for BT in (DefaultBasis, DefaultOrthonormalBasis, DefaultOrthogonalBasis)
            @testset "Basis $(BT)" begin
                cb = BT()
                @test b == get_vectors(M, x, get_basis(A, x, cb))
                v = similar(x)
                v2 = similar(x)
                @test_throws ErrorException get_vector(A, x, [1.0], cb)
                @test_throws DomainError get_vector(A, [1.0], [1.0, 0.0, 0.0], cb)
                @test_throws ErrorException get_vector!(A, v, x, [], cb)
                @test_throws DomainError get_vector!(A, v, [1.0], [1.0, 0.0, 0.0], cb)
                @test_throws DomainError get_coordinates(A, x, [1.0], cb)
                @test_throws DomainError get_coordinates!(A, v, x, [], cb)
                @test_throws DomainError get_coordinates!(A, v, [1.0], [1.0, 0.0, 0.0], cb)
                @test get_vector(A, x, [1, 2, 3], cb) ≈ get_vector(M, x, [1, 2, 3], cb)
                @test get_vector!(A, v2, x, [1, 2, 3], cb) ≈
                      get_vector!(M, v, x, [1, 2, 3], cb)
                @test get_coordinates(A, x, [1, 2, 3], cb) ≈
                      get_coordinates(M, x, [1, 2, 3], cb)
                @test get_coordinates!(A, v2, x, [1, 2, 3], cb) ≈
                      get_coordinates!(M, v, x, [1, 2, 3], cb)

                @test_throws ErrorException get_basis(A, x, CachedBasis(cb, [x]))
                @test_throws ErrorException get_basis(A, x, CachedBasis(cb, [x, x, x]))
                @test_throws ErrorException get_basis(A, x, CachedBasis(cb, [2 * x, x, x]))
                if BT <: ManifoldsBase.AbstractOrthogonalBasis
                    @test_throws ArgumentError get_basis(
                        A,
                        x,
                        CachedBasis(
                            cb,
                            [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        ),
                    )
                elseif BT <: ManifoldsBase.AbstractOrthonormalBasis
                    @test_throws ArgumentError get_basis(
                        A,
                        x,
                        CachedBasis(
                            cb,
                            [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        ),
                    )
                end
            end
        end
    end
end
