using Test
using ManifoldsBase
using ManifoldsBase: DefaultManifold, connection, decorated_manifold, is_default_connection, is_default_metric, metric

using Random

using LinearAlgebra: I, dot, Diagonal, normalize
import ManifoldsBase: default_retraction_method, connection

struct TestDefaultManifold{N} <: AbstractManifold{ℝ} end
struct TestDefaultManifoldMetric <: AbstractMetric end
struct TestScaledDefaultManifoldMetric <: AbstractMetric end
struct TestRetraction <: AbstractRetractionMethod end
struct TestConnection <: AbstractAffineConnection end

ManifoldsBase.connection(::TestDefaultManifold) = TestConnection()
ManifoldsBase.default_retraction_method(::TestDefaultManifold) = TestRetraction()
function ManifoldsBase.default_retraction_method(
        ::MetricManifold{ℝ, <:TestDefaultManifold, <:TestDefaultManifoldMetric},
    )
    return TestRetraction()
end

ManifoldsBase.manifold_dimension(::TestDefaultManifold{N}) where {N} = N
function ManifoldsBase.get_coordinates_orthogonal(
        M::MetricManifold{ℝ, <:TestDefaultManifold, <:TestDefaultManifoldMetric},
        ::Any,
        X,
        ::ManifoldsBase.AbstractNumbers,
    )
    return 1 ./ [1.0:manifold_dimension(M)...] .* X
end
function ManifoldsBase.get_coordinates_orthogonal!(
        M::MetricManifold{ℝ, <:TestDefaultManifold, <:TestDefaultManifoldMetric},
        c,
        ::Any,
        X,
        ::ManifoldsBase.AbstractNumbers,
    )
    c .= 1 ./ [1.0:manifold_dimension(M)...] .* X
    return c
end
function ManifoldsBase.get_vector_orthonormal!(
        M::MetricManifold{ℝ, <:TestDefaultManifold, <:TestDefaultManifoldMetric},
        X,
        ::Any,
        c,
        ::ManifoldsBase.AbstractNumbers,
    )
    X .= [1.0:manifold_dimension(M)...] .* c
    return X
end
function ManifoldsBase.get_coordinates_orthogonal!(
        M::MetricManifold{ℝ, <:TestDefaultManifold, <:TestScaledDefaultManifoldMetric},
        c,
        ::Any,
        X,
    )
    c .= 1 ./ (2 .* [1.0:manifold_dimension(M)...]) .* X
    return c
end
function ManifoldsBase.get_vector_orthogonal!(
        M::MetricManifold{ℝ, <:TestDefaultManifold, <:TestScaledDefaultManifoldMetric},
        ::Any,
        c,
        ::ManifoldsBase.AbstractNumbers,
    )
    return 2 .* [1.0:manifold_dimension(M)...] .* c
end
function ManifoldsBase.get_vector_orthogonal!(
        M::MetricManifold{ℝ, <:TestDefaultManifold, <:TestScaledDefaultManifoldMetric},
        X,
        ::Any,
        c,
        ::ManifoldsBase.AbstractNumbers,
    )
    X .= 2 .* [1.0:manifold_dimension(M)...] .* c
    return X
end

struct BaseManifold{N} <: AbstractManifold{ℝ} end
struct BaseManifoldMetric{M} <: AbstractMetric end
struct DefaultBaseManifoldMetric <: AbstractMetric end
struct NotImplementedMetric <: AbstractMetric end

ManifoldsBase.manifold_dimension(::BaseManifold{N}) where {N} = N
ManifoldsBase.inner(::BaseManifold, p, X, Y) = 2 * dot(X, Y)
ManifoldsBase.representation_size(::BaseManifold{N}) where {N} = (N,)
function ManifoldsBase.rand!(rng::AbstractRNG, ::BaseManifold, p; kwargs...)
    return randn!(rng, p)
end
ManifoldsBase.get_embedding(::BaseManifold{N}) where {N} = ManifoldsBase.DefaultManifold(N)
ManifoldsBase.exp!(::BaseManifold, q, p, X) = q .= p + 2 * X
ManifoldsBase.exp_fused!(::BaseManifold, q, p, X, t::Number) = q .= p + 2 * t * X
ManifoldsBase.log!(::BaseManifold, Y, p, q) = Y .= (q - p) / 2
ManifoldsBase.project!(::BaseManifold, Y, p, X) = Y .= 2 .* X
ManifoldsBase.project!(::BaseManifold, q, p) = (q .= p)
ManifoldsBase.injectivity_radius(::BaseManifold) = Inf
ManifoldsBase.injectivity_radius(::BaseManifold, ::Any) = Inf
ManifoldsBase.injectivity_radius(::BaseManifold, ::AbstractRetractionMethod) = Inf
ManifoldsBase._injectivity_radius(::BaseManifold, ::ExponentialRetraction) = Inf
ManifoldsBase.injectivity_radius(::BaseManifold, ::Any, ::AbstractRetractionMethod) = Inf
ManifoldsBase._injectivity_radius(::BaseManifold, ::Any, ::ExponentialRetraction) = Inf
function ManifoldsBase.exp!(
        M::MetricManifold{ℝ, <:BaseManifold{N}, BaseManifoldMetric{N}},
        q,
        p,
        X,
    ) where {N}
    return exp!(base_manifold(M), q, p, X)
end
function ManifoldsBase.exp_fused!(
        M::MetricManifold{ℝ, <:BaseManifold{N}, BaseManifoldMetric{N}},
        q,
        p,
        X,
        t::Number,
    ) where {N}
    return ManifoldsBase.exp_fused!(base_manifold(M), q, p, X, t)
end
function ManifoldsBase.parallel_transport_to!(::BaseManifold, Y, p, X, q)
    return (Y .= X)
end
function ManifoldsBase.get_basis(
        ::BaseManifold{N},
        p,
        B::DefaultOrthonormalBasis{<:Any, ManifoldsBase.TangentSpaceType},
    ) where {N}
    return CachedBasis(B, [(Matrix{eltype(p)}(I, N, N)[:, i]) for i in 1:N])
end
function ManifoldsBase.get_coordinates_orthonormal!(
        ::BaseManifold,
        Y,
        p,
        X,
        ::ManifoldsBase.AbstractNumbers,
    )
    return Y .= X
end
function ManifoldsBase.get_vector_orthonormal!(
        ::BaseManifold,
        Y,
        p,
        X,
        ::ManifoldsBase.AbstractNumbers,
    )
    return Y .= X
end

ManifoldsBase.metric(::BaseManifold) = DefaultBaseManifoldMetric()

@testset "Metrics" begin
    # some tests failed due to insufficient accuracy for a particularly bad RNG state
    Random.seed!(42)
    @testset "Metric Basics" begin
        @test ManifoldsBase.EuclideanMetric() isa ManifoldsBase.AbstractMetric
        @test repr(MetricManifold(DefaultManifold(3), ManifoldsBase.EuclideanMetric())) ===
            "MetricManifold(DefaultManifold(3; field = ℝ), ManifoldsBase.EuclideanMetric())"
    end
    @testset "Connection Trait" begin
        E = DefaultManifold(3)
        M = ConnectionManifold(E, LeviCivitaConnection())
        @test is_default_connection(M)
        @test decorated_manifold(M) == E
        @test connection(M) == LeviCivitaConnection()
        @test is_default_connection(E, LeviCivitaConnection())
        @test !is_default_connection(TestDefaultManifold{3}(), LeviCivitaConnection())
        @test connection(TestDefaultManifold{3}()) == TestConnection()

        @test manifold_dimension(M) == manifold_dimension(E)
        @test representation_size(M) == representation_size(E)
    end

    @testset "solve_exp_ode error message" begin
        E = TestDefaultManifold{3}()
        g = TestDefaultManifoldMetric()
        M = MetricManifold(E, g)
        default_retraction_method(::TestDefaultManifold) = TestRetraction()
        @test TestDefaultManifoldMetric(E) === M
        @test g(E) === M
        p = [1.0, 2.0, 3.0]
        X = [2.0, 3.0, 4.0]
        q = similar(X)
        @test_throws MethodError exp(M, p, X)
        @test_throws MethodError ManifoldsBase.exp_fused(M, p, X, 1.0)
        @test_throws MethodError exp!(M, q, p, X)
        @test_throws MethodError ManifoldsBase.exp_fused!(M, q, p, X, 1.0)

        N = ConnectionManifold(E, LeviCivitaConnection())
        @test_throws MethodError exp(N, p, X)
        @test_throws MethodError ManifoldsBase.exp_fused(N, p, X, 1.0)
        @test_throws MethodError exp!(N, q, p, X)
        @test_throws MethodError ManifoldsBase.exp_fused!(N, q, p, X, 1.0)
    end

    @testset "Local Metric Error message" begin
        M = MetricManifold(BaseManifold{2}(), NotImplementedMetric())
        p = [3, 4]
        X = [3, 4]
        Y = [3, 4]
        Z = [3, 4]
        @test_throws MethodError inner(M, p, X, Y)
        @test_throws MethodError Weingarten(M, p, X, Y)
        @test_throws MethodError Weingarten!(M, Z, p, X, Y)
    end
    @testset "scaled DefaultManifold metric" begin
        n = 3
        E = TestDefaultManifold{n}()
        g = TestDefaultManifoldMetric()
        M = MetricManifold(E, g)
        @test repr(M) == "MetricManifold(TestDefaultManifold{3}(), TestDefaultManifoldMetric())"

        @test TestDefaultManifoldMetric()(E) === M
        @test TestDefaultManifoldMetric(E) === M
        @test connection(M) === ManifoldsBase.LeviCivitaConnection()

        @test manifold_dimension(M) == n
        @test representation_size(M) === nothing
        @test base_manifold(M) === E
        @test metric(M) === g

        p, X = randn(n), randn(n)

        @test ManifoldsBase.check_point(M, p) == ManifoldsBase.check_point(E, p)
        @test ManifoldsBase.check_vector(M, p, X) == ManifoldsBase.check_vector(E, p, X)
    end

    @testset "default_* functions" begin
        E = DefaultManifold(3)
        EM = MetricManifold(E, ManifoldsBase.EuclideanMetric())
        T = Vector{Float64}
        @test default_retraction_method(EM) === default_retraction_method(E)
        @test default_inverse_retraction_method(EM) === default_inverse_retraction_method(E)
        @test default_vector_transport_method(EM) === default_vector_transport_method(E)
        @test default_retraction_method(EM, T) == default_retraction_method(E, T)
        @test default_inverse_retraction_method(EM, T) == default_inverse_retraction_method(E, T)
        @test default_vector_transport_method(EM, T) == default_vector_transport_method(E, T)

        EC = ConnectionManifold(E, TestConnection())
        @test default_retraction_method(EC) === default_retraction_method(E)
        @test default_inverse_retraction_method(EC) === default_inverse_retraction_method(E)
        @test default_vector_transport_method(EC) === default_vector_transport_method(E)
        @test default_retraction_method(EC, T) == default_retraction_method(E, T)
        @test default_inverse_retraction_method(EC, T) == default_inverse_retraction_method(E, T)
        @test default_vector_transport_method(EC, T) == default_vector_transport_method(E, T)
    end

    @testset "is_metric_function" begin
        for f in [exp, inner, log]
            @test ManifoldsBase.is_metric_function(f)
        end
        for f in [manifold_dimension, representation_size]
            @test !ManifoldsBase.is_metric_function(f)
        end
    end

    @testset "Metric decorator" begin
        M = BaseManifold{3}()
        g = BaseManifoldMetric{3}()
        MM = MetricManifold(M, g)
        TP = Vector{Float64}

        @test DefaultBaseManifoldMetric(BaseManifold{3}()) ===
            MetricManifold(BaseManifold{3}(), DefaultBaseManifoldMetric())
        MT = DefaultBaseManifoldMetric()
        @test MT(BaseManifold{3}()) ===
            MetricManifold(BaseManifold{3}(), DefaultBaseManifoldMetric())

        g2 = DefaultBaseManifoldMetric()
        MM2 = MetricManifold(M, g2)

        @test is_default_metric(MM) == is_default_metric(base_manifold(MM), metric(MM))
        @test is_default_metric(MM2) == is_default_metric(base_manifold(MM2), metric(MM2))
        @test is_default_metric(MM2)

        @test get_embedding(MM) === get_embedding(M)
        @test get_embedding(MM, TP) === get_embedding(M, TP)

        @test convert(typeof(MM2), M) == MM2
        @test_throws ErrorException convert(typeof(MM), M)
        p = [0.1, 0.2, 0.4]
        X = [0.5, 0.7, 0.11]
        Y = [0.13, 0.17, 0.19]
        q = allocate(p)

        @test is_point(MM, rand(MM))
        @test is_point(MM, rand(Xoshiro(), MM))
        @test is_vector(MM, p, rand(MM; vector_at = p))
        rand!(MM, q)
        @test is_point(MM, q)
        @test embed(MM, p) == p
        @test embed!(MM, q, p) == q

        p2 = allocate(p)
        copyto!(MM, p2, p)
        p3 = allocate(p)
        copyto!(M, p3, p)
        @test p2 == p3
        X = zero_vector(MM, p)
        Y = allocate(X)
        copyto!(MM, Y, p, X)
        Y2 = allocate(X)
        copyto!(M, Y2, p, X)
        @test Y == Y2

        X = [0.5, 0.7, 0.11]

        @test inner(M, p, X, Y) == 2 * dot(X, Y)
        @test exp(M, p, X) == p + 2 * X
        @test ManifoldsBase.exp_fused(M, p, X, 0.5) == p + X
        @test exp(MM2, p, X) == exp(M, p, X)
        @test ManifoldsBase.exp_fused(MM2, p, X, 0.5) == ManifoldsBase.exp_fused(M, p, X, 0.5)
        @test exp!(MM, q, p, X) === exp!(M, q, p, X)
        @test ManifoldsBase.exp_fused!(MM, q, p, X, 0.5) ===
            ManifoldsBase.exp_fused!(M, q, p, X, 0.5)
        @test retract!(MM, q, p, X) === retract!(M, q, p, X)
        @test ManifoldsBase.retract_fused!(MM, q, p, X, 1) ===
            ManifoldsBase.retract_fused!(M, q, p, X, 1)
        @test ManifoldsBase.retract_fused(MM, p, X, 1) ==
            ManifoldsBase.retract_fused(M, p, X, 1)
        @test project(MM, p) == project(M, p)
        @test project(MM, p, X) == project(M, p, X)
        @test project!(MM, Y, p, X) === project!(M, Y, p, X)
        @test project!(MM, q, p) === project!(M, q, p)
        # without a definition for the metric from the embedding, no projection possible
        @test_throws MethodError log!(MM, Y, p, q) === project!(M, Y, p, q)
        @test_throws MethodError vector_transport_to!(MM, Y, p, X, q) ===
            vector_transport_to!(M, Y, p, X, q)
        # without DiffEq, these error
        @test_throws MethodError exp(MM, p, X, 1:3)
        # these always fall back anyways.
        @test zero_vector!(MM, X, p) === zero_vector!(M, X, p)

        @test default_approximation_method(MM, retract) === default_approximation_method(M, retract)

        @test injectivity_radius(MM, p) === injectivity_radius(M, p)
        @test injectivity_radius(MM) === injectivity_radius(M)
        @test injectivity_radius(MM, ProjectionRetraction()) ===
            injectivity_radius(M, ProjectionRetraction())
        @test injectivity_radius(MM, ExponentialRetraction()) ===
            injectivity_radius(M, ExponentialRetraction())
        @test injectivity_radius(MM) === injectivity_radius(M)

        @test is_point(MM, p) === is_point(M, p)
        @test is_vector(MM, p, X) === is_vector(M, p, X)

        @test inner(MM2, p, X, Y) === inner(M, p, X, Y)
        @test norm(MM2, p, X) === norm(M, p, X)
        @test distance(MM2, p, q) === distance(M, p, q)
        @test exp!(MM2, q, p, X) === exp!(M, q, p, X)
        @test exp(MM2, p, X) == exp(M, p, X)
        @test log!(MM2, X, p, q) === log!(M, X, p, q)
        @test log(MM2, p, q) == log(M, p, q)
        @test retract!(MM2, q, p, X) === retract!(M, q, p, X)
        for rm in [ExponentialRetraction(), EmbeddedRetraction(ExponentialRetraction())]
            @test retract(MM2, p, X, rm) == retract(M, p, X, rm)
            @test retract!(MM2, q, p, X, rm) === retract!(M, q, p, X, rm)
        end
        @test inverse_retract!(MM2, Y2, p, q) === inverse_retract!(M, Y2, q, p)
        for irm in [LogarithmicInverseRetraction(), EmbeddedInverseRetraction(LogarithmicInverseRetraction())]
            @test inverse_retract(MM2, p, q, irm) == inverse_retract(M, q, p, irm)
            @test inverse_retract!(MM2, Y2, p, q, irm) === inverse_retract!(M, Y2, q, p, irm)
        end
        @test ManifoldsBase.retract_fused!(MM2, q, p, X, 1) ===
            ManifoldsBase.retract_fused!(M, q, p, X, 1)
        @test_throws MethodError is_flat(MM)

        @test project!(MM2, q, p) === project!(M, q, p)
        @test project!(MM2, Y, p, X) === project!(M, Y, p, X)
        @test parallel_transport_to(MM2, p, X, q) == parallel_transport_to(M, q, X, p)
        @test parallel_transport_to!(MM2, Y, p, X, q) ==
            parallel_transport_to!(M, Y, q, X, p)
        @test project!(MM2, Y, p, X) === project!(M, Y, p, X)
        @test vector_transport_to(MM2, p, X, q) == vector_transport_to(M, q, X, p)
        @test vector_transport_to!(MM2, Y, p, X, q) == vector_transport_to!(M, Y, p, X, q)

        @test vector_transport_direction(MM2, p, X, Y) == vector_transport_to(M, p, X, Y)
        @test vector_transport_direction!(MM2, Y2, p, X, Y) == vector_transport_to!(M, Y2, p, X, Y)

        c = 2 * ones(3)
        m = ParallelTransport()
        @test zero_vector!(MM2, X, p) === zero_vector!(M, X, p)
        @test injectivity_radius(MM2, p) === injectivity_radius(M, p)
        @test injectivity_radius(MM2) === injectivity_radius(M)
        @test injectivity_radius(MM2, p, ExponentialRetraction()) ===
            injectivity_radius(M, p, ExponentialRetraction())
        @test injectivity_radius(MM2, ExponentialRetraction()) ===
            injectivity_radius(M, ExponentialRetraction())
        @test injectivity_radius(MM2, p, ProjectionRetraction()) ===
            injectivity_radius(M, p, ProjectionRetraction())
        @test injectivity_radius(MM2, ProjectionRetraction()) ===
            injectivity_radius(M, ProjectionRetraction())
        @test is_point(MM2, p) === is_point(M, p)
        @test is_vector(MM2, p, X) === is_vector(M, p, X)

        @test_throws MethodError log(MM, p, q)

        @test get_basis(M, p, DefaultOrthonormalBasis()).data ==
            get_basis(MM2, p, DefaultOrthonormalBasis()).data
        @test_throws MethodError get_basis(MM, p, DefaultOrthonormalBasis())

        @test_throws MethodError get_coordinates(MM, p, X, DefaultOrthonormalBasis())
        c2 = similar(c)
        @test_throws MethodError get_coordinates!(MM, c2, p, X, DefaultOrthonormalBasis())
        @test_throws MethodError get_vector(MM, p, c, DefaultOrthonormalBasis())
        @test_throws MethodError get_vector!(MM, Y2, p, c, DefaultOrthonormalBasis())
    end
end
