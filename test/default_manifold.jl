using ManifoldsBase

using LinearAlgebra
using DoubleFloats
using ForwardDiff
using ReverseDiff
using StaticArrays
using Test

struct CustomDefinedRetraction <: ManifoldsBase.AbstractRetractionMethod end
struct CustomUndefinedRetraction <: ManifoldsBase.AbstractRetractionMethod end

function ManifoldsBase.injectivity_radius(
    ::ManifoldsBase.DefaultManifold,
    ::CustomDefinedRetraction,
)
    return 10.0
end

struct MatrixVectorTransport{T} <: AbstractVector{T}
    m::Matrix{T}
end

Base.getindex(x::MatrixVectorTransport, i) = x.m[:, i]

Base.size(x::MatrixVectorTransport) = (size(x.m, 2),)

@testset "Testing Default (Euclidean)" begin
    M = ManifoldsBase.DefaultManifold(3)
    types = [
        Vector{Float64},
        SizedVector{3,Float64},
        MVector{3,Float64},
        Vector{Float32},
        SizedVector{3,Float32},
        MVector{3,Float32},
        Vector{Double64},
        MVector{3,Double64},
        SizedVector{3,Double64},
    ]

    @test isa(manifold_dimension(M), Integer)
    @test manifold_dimension(M) ≥ 0
    @test base_manifold(M) == M
    @test number_system(M) == ℝ
    @test ManifoldsBase.representation_size(M) == (3,)

    @test injectivity_radius(M) == Inf

    rm = ManifoldsBase.ExponentialRetraction()
    irm = ManifoldsBase.LogarithmicInverseRetraction()

    rm2 = CustomDefinedRetraction()
    rm3 = CustomUndefinedRetraction()

    for T in types
        @testset "Type $T" begin
            pts = convert.(Ref(T), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

            @test injectivity_radius(M, pts[1]) == Inf
            @test injectivity_radius(M, pts[1], rm) == Inf
            @test injectivity_radius(M, rm) == Inf
            @test injectivity_radius(M, rm2) == 10
            @test injectivity_radius(M, pts[1], rm2) == 10
            @test_throws ErrorException injectivity_radius(M, rm3)
            @test_throws ErrorException injectivity_radius(M, pts[1], rm3)

            tv1 = log(M, pts[1], pts[2])

            for pt in pts
                @test is_manifold_point(M, pt)
            end
            @test is_tangent_vector(M, pts[1], tv1; atol = eps(eltype(pts[1])))

            tv2 = log(M, pts[2], pts[1])
            @test isapprox(M, pts[2], exp(M, pts[1], tv1))
            @test isapprox(M, pts[1], exp(M, pts[1], tv1, 0))
            @test isapprox(M, pts[2], exp(M, pts[1], tv1, 1))
            @test isapprox(M, pts[1], exp(M, pts[2], tv2))
            @test is_manifold_point(M, retract(M, pts[1], tv1))
            @test isapprox(M, pts[1], retract(M, pts[1], tv1, 0))

            @test is_manifold_point(M, retract(M, pts[1], tv1, rm))
            @test isapprox(M, pts[1], retract(M, pts[1], tv1, 0, rm))

            new_pt = exp(M, pts[1], tv1)
            retract!(M, new_pt, pts[1], tv1)
            @test is_manifold_point(M, new_pt)
            for x in pts
                @test isapprox(
                    M,
                    zero_tangent_vector(M, x),
                    log(M, x, x);
                    atol = eps(eltype(x)),
                )
                @test isapprox(
                    M,
                    zero_tangent_vector(M, x),
                    inverse_retract(M, x, x);
                    atol = eps(eltype(x)),
                )
                @test isapprox(
                    M,
                    zero_tangent_vector(M, x),
                    inverse_retract(M, x, x, irm);
                    atol = eps(eltype(x)),
                )
            end
            zero_tangent_vector!(M, tv1, pts[1])
            @test isapprox(M, pts[1], tv1, zero_tangent_vector(M, pts[1]))
            log!(M, tv1, pts[1], pts[2])
            @test norm(M, pts[1], tv1) ≈ sqrt(inner(M, pts[1], tv1, tv1))

            @test isapprox(M, exp(M, pts[1], tv1, 1), pts[2])
            @test isapprox(M, exp(M, pts[1], tv1, 0), pts[1])

            @test distance(M, pts[1], pts[2]) ≈ norm(M, pts[1], tv1)

            @test mid_point(M, pts[1], pts[2]) == [0.5, 0.5, 0.0]
            midp = allocate(pts[1])
            @test mid_point!(M, midp, pts[1], pts[2]) === midp
            @test midp == [0.5, 0.5, 0.0]

            @testset "Geodesic interface test" begin
                @test isapprox(M, geodesic(M, pts[1], tv1)(0.0), pts[1])
                @test isapprox(M, geodesic(M, pts[1], tv1)(1.0), pts[2])
                @test isapprox(M, geodesic(M, pts[1], tv1, 1.0), pts[2])
                @test isapprox(M, geodesic(M, pts[1], tv1, 1.0 / 2), (pts[1] + pts[2]) / 2)
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2])(0.0), pts[1])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2])(1.0), pts[2])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2], 0.0), pts[1])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2], 1.0), pts[2])
                @test all(isapprox.(
                    Ref(M),
                    geodesic(M, pts[1], tv1, [0.0, 1.0 / 2, 1.0]),
                    [pts[1], (pts[1] + pts[2]) / 2, pts[2]],
                ))
                @test all(isapprox.(
                    Ref(M),
                    shortest_geodesic(M, pts[1], pts[2], [0.0, 1.0 / 2, 1.0]),
                    [pts[1], (pts[1] + pts[2]) / 2, pts[2]],
                ))
            end

            @testset "basic linear algebra in tangent space" begin
                @test isapprox(
                    M,
                    pts[1],
                    0 * tv1,
                    zero_tangent_vector(M, pts[1]);
                    atol = eps(eltype(pts[1])),
                )
                @test isapprox(M, pts[1], 2 * tv1, tv1 + tv1)
                @test isapprox(M, pts[1], 0 * tv1, tv1 - tv1)
                @test isapprox(M, pts[1], (-1) * tv1, -tv1)
            end

            @testset "Hat and vee in the tangent space" begin
                X = log(M, pts[1], pts[2])
                a = vee(M, pts[1], X)
                b = similar(a)
                vee!(M, b, pts[1], X)
                Y = hat(M, pts[1], a)
                Z = similar(Y)
                hat!(M, Z, pts[1], a)
                @test a == b
                @test X == Y
                @test Z == X
                @test a == vec(X)
            end

            @testset "broadcasted linear algebra in tangent space" begin
                @test isapprox(M, pts[1], 3 * tv1, 2 .* tv1 .+ tv1)
                @test isapprox(M, pts[1], -tv1, tv1 .- 2 .* tv1)
                @test isapprox(M, pts[1], -tv1, .-tv1)
                v = similar(tv1)
                v .= 2 .* tv1 .+ tv1
                @test v ≈ 3 * tv1
            end

            @testset "project test" begin
                # point
                @test isapprox(M, pts[1], project(M, pts[1]))
                pt = similar(pts[1])
                project!(M, pt, pts[1])
                @test isapprox(M, pt, pts[1])

                @test isapprox(M, pts[1], embed(M, pts[1]))
                pt = similar(pts[1])
                embed!(M, pt, pts[1])
                @test isapprox(M, pt, pts[1])

                # tangents
                @test isapprox(M, pts[1], tv1, project(M, pts[1], tv1))
                tv = similar(tv1)
                project!(M, tv, pts[1], tv1)
                @test isapprox(M, pts[1], tv, tv1)

                @test isapprox(M, pts[1], tv1, embed(M, pts[1], tv1))
                tv = similar(tv1)
                embed!(M, tv, pts[1], tv1)
                @test isapprox(M, pts[1], tv, tv1)
            end

            @testset "vector transport" begin
                v1 = log(M, pts[1], pts[2])
                v2 = log(M, pts[1], pts[3])
                v1t1 = vector_transport_to(M, pts[1], v1, pts[3])
                v1t2 = zero(v1t1)
                vector_transport_to!(M, v1t2, pts[1], v1, v2, ProjectionTransport())
                v1t3 = vector_transport_direction(M, pts[1], v1, v2)
                @test is_tangent_vector(M, pts[3], v1t1)
                @test is_tangent_vector(M, pts[3], v1t3)
                @test isapprox(M, pts[3], v1t1, v1t3)
                # along a `Vector` of points
                c = [pts[1]]
                v1t4 = vector_transport_along(M, pts[1], v1, c)
                @test isapprox(M, pts[1], v1, v1t4)
                v1t5 = allocate(v1)
                vector_transport_along!(M, v1t5, pts[1], v1, c)
                @test isapprox(M, pts[1], v1, v1t5)
                # along a custom type of points
                T = eltype(pts[1])
                c2 = MatrixVectorTransport{T}(reshape(pts[1], length(pts[1]), 1))
                v1t4c2 = vector_transport_along(M, pts[1], v1, c2)
                @test isapprox(M, pts[1], v1, v1t4c2)
                v1t5c2 = allocate(v1)
                vector_transport_along!(M, v1t5c2, pts[1], v1, c2)
                @test isapprox(M, pts[1], v1, v1t5c2)
                # On Euclidean Space Schild & Pole are identity
                @test vector_transport_to(
                    M,
                    pts[1],
                    v2,
                    pts[2],
                    SchildsLadderTransport(),
                ) == v2
                @test vector_transport_to(M, pts[1], v2, pts[2], PoleLadderTransport()) ==
                      v2
                # along is also the identity
                c = [0.5 * (pts[1] + pts[2]), pts[2], 0.5 * (pts[2] + pts[3]), pts[3]]
                @test vector_transport_along(M, pts[1], v2, c, SchildsLadderTransport()) ==
                      v2
                @test vector_transport_along(M, pts[1], v2, c, PoleLadderTransport()) == v2
                @test vector_transport_along(M, pts[1], v2, c, ParallelTransport()) == v2
                # check mutating ones with defaults
                p = allocate(pts[1])
                ManifoldsBase.pole_ladder!(M, p, pts[1], pts[2], pts[3])
                # -log_p3 p == log_p1 p2
                @test isapprox(M, -log(M, pts[3], p), log(M, pts[1], pts[2]))
                ManifoldsBase.schilds_ladder!(M, p, pts[1], pts[2], pts[3])
                @test isapprox(M, log(M, pts[3], p), log(M, pts[1], pts[2]))
            end

            @testset "ForwardDiff support" begin
                exp_f(t) = distance(M, pts[1], exp(M, pts[1], t * tv1))
                d12 = distance(M, pts[1], pts[2])
                for t in 0.1:0.1:0.9
                    @test d12 ≈ ForwardDiff.derivative(exp_f, t)
                end

                retract_f(t) = distance(M, pts[1], retract(M, pts[1], t * tv1))
                for t in 0.1:0.1:0.9
                    @test ForwardDiff.derivative(retract_f, t) ≥ 0
                end
            end

            isa(pts[1], Union{Vector,SizedVector}) && @testset "ReverseDiff support" begin
                exp_f(t) = distance(M, pts[1], exp(M, pts[1], t[1] * tv1))
                d12 = distance(M, pts[1], pts[2])
                for t in 0.1:0.1:0.9
                    @test d12 ≈ ReverseDiff.gradient(exp_f, [t])[1]
                end

                retract_f(t) = distance(M, pts[1], retract(M, pts[1], t[1] * tv1))
                for t in 0.1:0.1:0.9
                    @test ReverseDiff.gradient(retract_f, [t])[1] ≥ 0
                end
            end
        end
    end
end
