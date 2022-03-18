using ManifoldsBase
using ManifoldsBase:
    @manifold_element_forwards, @manifold_vector_forwards, @default_manifold_fallbacks
using ManifoldsBase: DefaultManifold
import ManifoldsBase:
    number_eltype,
    check_point,
    distance,
    embed!,
    exp!,
    inner,
    isapprox,
    log!,
    parallel_transport_to!,
    retract!,
    inverse_retract!
import Base: angle, convert
using LinearAlgebra
using DoubleFloats
using ForwardDiff
using ReverseDiff
using StaticArrays
using Test

struct CustomDefinedRetraction <: ManifoldsBase.AbstractRetractionMethod end
struct CustomUndefinedRetraction <: ManifoldsBase.AbstractRetractionMethod end
struct CustomDefinedInverseRetraction <: ManifoldsBase.AbstractInverseRetractionMethod end

struct DefaultPoint{T} <: AbstractManifoldPoint
    value::T
end
DefaultPoint(v::T) where {T} = DefaultPoint{T}(v)
convert(::Type{DefaultPoint{T}}, v::T) where {T} = DefaultPoint(v)
Base.size(p::DefaultPoint) = size(p.value)
Base.eltype(v::DefaultPoint) = eltype(v.value)

struct DefaultTVector{T} <: TVector
    value::T
end
DefaultTVector(v::T) where {T} = DefaultTVector{T}(v)
Base.size(X::DefaultTVector) = size(X.value)
Base.eltype(X::DefaultTVector) = eltype(X.value)

ManifoldsBase.@manifold_element_forwards DefaultPoint value
ManifoldsBase.@manifold_vector_forwards DefaultTVector value
ManifoldsBase.@default_manifold_fallbacks ManifoldsBase.DefaultManifold DefaultPoint DefaultTVector value value

function ManifoldsBase._injectivity_radius(::DefaultManifold, ::CustomDefinedRetraction)
    return 10.0
end
function ManifoldsBase._retract(M::DefaultManifold, p, X, ::CustomDefinedRetraction)
    return retract_custom(M, p, X)
end
function retract_custom(::DefaultManifold, p::DefaultPoint, X::DefaultTVector)
    return DefaultPoint(p.value + X.value)
end
function ManifoldsBase._inverse_retract(
    M::DefaultManifold,
    p,
    q,
    ::CustomDefinedInverseRetraction,
)
    return inverse_retract_custom(M, p, q)
end
function inverse_retract_custom(::DefaultManifold, p::DefaultPoint, q::DefaultPoint)
    return DefaultTVector(q.value - p.value)
end
struct MatrixVectorTransport{T} <: AbstractVector{T}
    m::Matrix{T}
end
# dummy retractions, inverse retracions for fallback tests - mutating should be enough
ManifoldsBase.retract_polar!(::DefaultManifold, q, p, X) = (q .= p .+ X)
ManifoldsBase.retract_project!(::DefaultManifold, q, p, X) = (q .= p .+ X)
ManifoldsBase.retract_qr!(::DefaultManifold, q, p, X) = (q .= p .+ X)
ManifoldsBase.retract_exp_ode!(::DefaultManifold, q, p, X, m, B) = (q .= p .+ X)
ManifoldsBase.retract_pade!(::DefaultManifold, q, p, X, i) = (q .= p .+ X)
ManifoldsBase.retract_softmax!(::DefaultManifold, q, p, X) = (q .= p .+ X)
ManifoldsBase.get_embedding(M::DefaultManifold) = M # dummy embedding
ManifoldsBase.inverse_retract_polar!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
ManifoldsBase.inverse_retract_project!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
ManifoldsBase.inverse_retract_qr!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
ManifoldsBase.inverse_retract_softmax!(::DefaultManifold, Y, p, q) = (Y .= q .- p)
ManifoldsBase.inverse_retract_nlsolve!(::DefaultManifold, Y, p, q, m) = (Y .= q .- p)
ManifoldsBase.vector_transport_along_project!(::DefaultManifold, Y, p, X, c) = (Y .= X)


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
        DefaultPoint{Vector{Float64}},
    ]

    @test repr(M) == "DefaultManifold(3; field = ℝ)"
    @test isa(manifold_dimension(M), Integer)
    @test manifold_dimension(M) ≥ 0
    @test base_manifold(M) == M
    @test number_system(M) == ManifoldsBase.ℝ
    @test ManifoldsBase.representation_size(M) == (3,)

    p = zeros(3)
    m = PolarRetraction()
    @test injectivity_radius(M) == Inf
    @test injectivity_radius(M, p) == Inf
    @test injectivity_radius(M, m) == Inf
    @test injectivity_radius(M, p, m) == Inf
    @test default_retraction_method(M) == ExponentialRetraction()
    @test default_inverse_retraction_method(M) == LogarithmicInverseRetraction()

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

            tv1 = log(M, pts[1], pts[2])

            for pt in pts
                @test is_point(M, pt)
            end
            @test is_vector(M, pts[1], tv1; atol = eps(eltype(pts[1])))

            tv2 = log(M, pts[2], pts[1])
            @test isapprox(M, pts[2], exp(M, pts[1], tv1))
            @test isapprox(M, pts[1], exp(M, pts[1], tv1, 0))
            @test isapprox(M, pts[2], exp(M, pts[1], tv1, 1))
            @test isapprox(M, pts[1], exp(M, pts[2], tv2))
            @test is_point(M, retract(M, pts[1], tv1))
            @test isapprox(M, pts[1], retract(M, pts[1], tv1, 0))

            @test is_point(M, retract(M, pts[1], tv1, rm))
            @test isapprox(M, pts[1], retract(M, pts[1], tv1, 0, rm))

            new_pt = exp(M, pts[1], tv1)
            retract!(M, new_pt, pts[1], tv1)
            @test is_point(M, new_pt)
            for x in pts
                @test isapprox(M, x, zero_vector(M, x), log(M, x, x); atol = eps(eltype(x)))
                @test isapprox(
                    M,
                    x,
                    zero_vector(M, x),
                    inverse_retract(M, x, x);
                    atol = eps(eltype(x)),
                )
                @test isapprox(
                    M,
                    x,
                    zero_vector(M, x),
                    inverse_retract(M, x, x, irm);
                    atol = eps(eltype(x)),
                )
            end
            zero_vector!(M, tv1, pts[1])
            @test isapprox(M, pts[1], tv1, zero_vector(M, pts[1]))
            log!(M, tv1, pts[1], pts[2])
            @test norm(M, pts[1], tv1) ≈ sqrt(inner(M, pts[1], tv1, tv1))

            @test isapprox(M, exp(M, pts[1], tv1, 1), pts[2])
            @test isapprox(M, exp(M, pts[1], tv1, 0), pts[1])

            @test distance(M, pts[1], pts[2]) ≈ norm(M, pts[1], tv1)

            @test mid_point(M, pts[1], pts[2]) == convert(T, [0.5, 0.5, 0.0])
            midp = allocate(pts[1])
            @test mid_point!(M, midp, pts[1], pts[2]) === midp
            @test midp == convert(T, [0.5, 0.5, 0.0])

            @testset "Geodesic interface test" begin
                @test isapprox(M, geodesic(M, pts[1], tv1)(0.0), pts[1])
                @test isapprox(M, geodesic(M, pts[1], tv1)(1.0), pts[2])
                @test isapprox(M, geodesic(M, pts[1], tv1, 1.0), pts[2])
                @test isapprox(M, geodesic(M, pts[1], tv1, 1.0 / 2), midp)
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2])(0.0), pts[1])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2])(1.0), pts[2])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2], 0.0), pts[1])
                @test isapprox(M, shortest_geodesic(M, pts[1], pts[2], 1.0), pts[2])
                @test all(
                    isapprox.(
                        Ref(M),
                        geodesic(M, pts[1], tv1, [0.0, 1.0 / 2, 1.0]),
                        [pts[1], midp, pts[2]],
                    ),
                )
                @test all(
                    isapprox.(
                        Ref(M),
                        shortest_geodesic(M, pts[1], pts[2], [0.0, 1.0 / 2, 1.0]),
                        [pts[1], midp, pts[2]],
                    ),
                )
            end

            @testset "basic linear algebra in tangent space" begin
                @test isapprox(
                    M,
                    pts[1],
                    0 * tv1,
                    zero_vector(M, pts[1]);
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
                @test a == ((T <: DefaultPoint) ? vec(X.value) : vec(X))
            end

            @testset "broadcasted linear algebra in tangent space" begin
                @test isapprox(M, pts[1], 3 * tv1, 2 .* tv1 .+ tv1)
                @test isapprox(M, pts[1], -tv1, tv1 .- 2 .* tv1)
                @test isapprox(M, pts[1], -tv1, .-tv1)
                v = similar(tv1)
                v .= 2 .* tv1 .+ tv1
                @test isapprox(M, pts[1], v, 3 * tv1)
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
                # test constructor
                @test default_vector_transport_method(M) == ParallelTransport()
                v1 = log(M, pts[1], pts[2])
                v2 = log(M, pts[1], pts[3])
                v1t1 = vector_transport_to(M, pts[1], v1, pts[3])
                v1t2 = zero(v1t1)
                vector_transport_to!(M, v1t2, pts[1], v1, v2, ProjectionTransport())
                v1t3 = vector_transport_direction(M, pts[1], v1, v2)
                @test ManifoldsBase.is_vector(M, pts[3], v1t1)
                @test ManifoldsBase.is_vector(M, pts[3], v1t3)
                @test isapprox(M, pts[3], v1t1, v1t3)
                # along a `Vector` of points
                c = [pts[1]]
                v1t4 = vector_transport_along(M, pts[1], v1, c)
                @test isapprox(M, pts[1], v1, v1t4)
                v1t5 = allocate(v1)
                vector_transport_along!(M, v1t5, pts[1], v1, c)
                @test isapprox(M, pts[1], v1, v1t5)
                # transport along more than one interims point
                @test vector_transport_along(M, pts[1], v1, pts[2:3]) == v1
                v1t6 = allocate(v1)
                vector_transport_along!(M, v1t6, pts[1], v1, pts[2:3])
                @test isapprox(M, pts[1], v1, v1t6)
                # along a custom type of points
                if T <: DefaultPoint
                    S = eltype(pts[1].value)
                    mat = reshape(pts[1].value, length(pts[1].value), 1)
                else
                    S = eltype(pts[1])
                    mat = reshape(pts[1], length(pts[1]), 1)
                end
                c2 = MatrixVectorTransport{S}(mat)
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
                @test vector_transport_to(
                    M,
                    pts[1],
                    v2,
                    pts[2],
                    ScaledVectorTransport(ParallelTransport()),
                ) == v2

                # along is also the identity
                c = [
                    mid_point(M, pts[1], pts[2]),
                    pts[2],
                    mid_point(M, pts[2], pts[3]),
                    pts[3],
                ]
                @test vector_transport_along(M, pts[1], v2, c, SchildsLadderTransport()) ==
                      v2
                @test vector_transport_along(M, pts[1], v2, c, PoleLadderTransport()) == v2
                @test vector_transport_along(M, pts[1], v2, c, ParallelTransport()) == v2
                # check mutating ones with defaults
                p = allocate(pts[1])
                ManifoldsBase.pole_ladder!(M, p, pts[1], pts[2], pts[3])
                # -log_p3 p == log_p1 p2
                @test isapprox(M, pts[3], -log(M, pts[3], p), log(M, pts[1], pts[2]))
                ManifoldsBase.schilds_ladder!(M, p, pts[1], pts[2], pts[3])
                @test isapprox(M, pts[3], log(M, pts[3], p), log(M, pts[1], pts[2]))

                @test repr(ParallelTransport()) == "ParallelTransport()"
                @test repr(ScaledVectorTransport(ParallelTransport())) ==
                      "ScaledVectorTransport(ParallelTransport())"
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

    @testset "mid_point on 0-index arrays" begin
        M = ManifoldsBase.DefaultManifold(1)
        p1 = fill(0.0)
        p2 = fill(1.0)
        @test isapprox(M, fill(0.5), mid_point(M, p1, p2))
    end

    @testset "Retraction" begin
        a = NLSolveInverseRetraction(ExponentialRetraction())
        @test a.retraction isa ExponentialRetraction

    end

    @testset "copy of points and vectors" begin
        M = ManifoldsBase.DefaultManifold(2)
        for (p, X) in (
            ([2.0, 3.0], [4.0, 5.0]),
            (DefaultPoint([2.0, 3.0]), DefaultTVector([4.0, 5.0])),
        )
            q = similar(p)
            copyto!(M, q, p)
            @test p == q
            r = copy(M, p)
            @test r == p
            Y = similar(X)
            copyto!(M, Y, p, X)
            @test Y == X
            Z = copy(M, p, X)
            @test Z == X
        end

        p1 = DefaultPoint([2.0, 3.0])
        p2 = copy(p1)
        @test (p1 == p2) && (p1 !== p2)
    end
    @testset "further vector and point automatic forwards" begin
        M = ManifoldsBase.DefaultManifold(3)
        p = DefaultPoint([1.0, 0.0, 0.0])
        q = DefaultPoint([0.0, 0.0, 0.0])
        X = DefaultTVector([0.0, 1.0, 0.0])
        Y = DefaultTVector([1.0, 0.0, 0.0])
        @test angle(M, p, X, Y) ≈ π / 2
        @test inverse_retract(M, p, q, LogarithmicInverseRetraction()) == -Y
        @test retract(M, q, Y, CustomDefinedRetraction()) == p
        @test retract(M, q, Y, ExponentialRetraction()) == p
        # rest not implemented - so they also fall back even onto mutating
        Z = similar(Y)
        r = similar(p)
        # test passthrough using the dummy implementations
        for retr in [PolarRetraction, ProjectionRetraction, QRRetraction, SoftmaxRetraction]
            @test retract(M, q, Y, retr()) == DefaultPoint(q.value + Y.value)
            @test retract!(M, r, q, Y, retr()) == DefaultPoint(q.value + Y.value)
        end
        @test retract(
            M,
            q,
            Y,
            ODEExponentialRetraction(PolarRetraction(), DefaultBasis()),
        ) == DefaultPoint(q.value + Y.value)
        @test retract!(
            M,
            r,
            q,
            Y,
            ODEExponentialRetraction(PolarRetraction(), DefaultBasis()),
        ) == DefaultPoint(q.value + Y.value)
        @test retract(M, q, Y, PadeRetraction(2)) == DefaultPoint(q.value + Y.value)
        @test retract!(M, r, q, Y, PadeRetraction(2)) == DefaultPoint(q.value + Y.value)
        @test retract!(M, r, q, Y, EmbeddedRetraction(ExponentialRetraction())) ==
              DefaultPoint(q.value + Y.value)
        @test retract(M, q, Y, EmbeddedRetraction(ExponentialRetraction())) ==
              DefaultPoint(q.value + Y.value)
        p2 = allocate(p, eltype(p.value), size(p.value))
        @test size(p2.value) == size(p.value)
        X2 = allocate(X, eltype(X.value), size(X.value))
        @test size(X2.value) == size(X.value)
        X3 = ManifoldsBase.allocate_result(M, log, p, q)
        @test log!(M, X3, p, q) == log(M, p, q)
        @test X3 == log(M, p, q)
        @test log!(M, X3, p, q) == log(M, p, q)
        @test X3 == log(M, p, q)
        @test inverse_retract(M, p, q, CustomDefinedInverseRetraction()) == -Y
        X4 = ManifoldsBase.allocate_result(M, inverse_retract, p, q)
        @test inverse_retract!(M, X4, p, q) == inverse_retract(M, p, q)
        @test X4 == inverse_retract(M, p, q)
        # rest not implemented but check passthrough
        for r in [
            PolarInverseRetraction,
            ProjectionInverseRetraction,
            QRInverseRetraction,
            SoftmaxInverseRetraction,
        ]
            @test inverse_retract(M, q, p, r()) == DefaultTVector(p.value - q.value)
            @test inverse_retract!(M, Z, q, p, r()) == DefaultTVector(p.value - q.value)
        end
        @test inverse_retract(
            M,
            q,
            p,
            EmbeddedInverseRetraction(LogarithmicInverseRetraction()),
        ) == DefaultTVector(p.value - q.value)
        @test inverse_retract(M, q, p, NLSolveInverseRetraction(ExponentialRetraction())) ==
              DefaultTVector(p.value - q.value)
        @test inverse_retract!(
            M,
            Z,
            q,
            p,
            EmbeddedInverseRetraction(LogarithmicInverseRetraction()),
        ) == DefaultTVector(p.value - q.value)
        @test inverse_retract!(
            M,
            Z,
            q,
            p,
            NLSolveInverseRetraction(ExponentialRetraction()),
        ) == DefaultTVector(p.value - q.value)
        c = ManifoldsBase.allocate_coordinates(M, p, Float64, manifold_dimension(M))
        @test c isa Vector
        @test length(c) == 3
        @test 2.0 \ X == DefaultTVector(2.0 \ X.value)
        @test X + Y == DefaultTVector(X.value + Y.value)
        @test +X == X
        @test (Y .= X) === Y
        # vector transport pass through
        @test vector_transport_to(M, p, X, q, ProjectionTransport()) == X
        @test vector_transport_direction(M, p, X, X, ProjectionTransport()) == X
        @test vector_transport_to!(M, Y, p, X, q, ProjectionTransport()) == X
        @test vector_transport_direction!(M, Y, p, X, X, ProjectionTransport()) == X
        @test vector_transport_along(M, p, X, X, ProjectionTransport()) == X
        @test vector_transport_along!(M, Z, p, X, X, ProjectionTransport()) == X
        @test vector_transport_to(M, p, X, :q, ProjectionTransport()) == X
        @test parallel_transport_to(M, p, X, q) == X
        @test parallel_transport_direction(M, p, X, X) == X
        @test parallel_transport_along(M, p, X, :c) == X
        @test parallel_transport_to!(M, Y, p, X, q) == X
        @test parallel_transport_direction!(M, Y, p, X, X) == X
        @test parallel_transport_along!(M, Y, p, X, :c) == X
    end
    @testset "DefaultManifold  and ONB" begin
        M = ManifoldsBase.DefaultManifold(3)
        p = [1.0, 0.0, 0.0]
        CB = get_basis(M, p, DefaultOrthonormalBasis())
        @test CB.data == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    end
    @testset "Show methods" begin
        @test repr(CayleyRetraction()) == "CayleyRetraction()"
        @test repr(PadeRetraction(2)) == "PadeRetraction(2)"
    end
    @testset "Further TestArrayRepresentation" begin
        M = ManifoldsBase.DefaultManifold(3)
        p = [1.0, 0.0, 0.0]
        X = [1.0, 0.0, 0.0]
        @test is_point(M, p, true)
        @test is_vector(M, p, X, true)
        pF = [1.0, 0.0]
        XF = [0.0, 0.0]
        m = ExponentialRetraction()
        @test_throws DomainError is_point(M, pF, true)
        @test_throws DomainError is_vector(M, p, XF, true)
        @test_throws DomainError is_vector(M, pF, XF, true; check_point = true)
        @test injectivity_radius(M) == Inf
        @test injectivity_radius(M, p) == Inf
        @test injectivity_radius(M, p, m) == Inf
        @test injectivity_radius(M, m) == Inf
    end
end
