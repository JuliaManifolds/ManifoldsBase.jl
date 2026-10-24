using ManifoldsBase, Test
using ManifoldsBase: DefaultManifold, DiagonalizingBasisData, EuclideanMetric, FVector, PowerBasisData,
    ProductBasisData, TFVector, ℂ

# two objects created by f with separately allocated arrays are equal and have equal hashes,
# while an object created by g is different
function test_value_equality(f, g)
    a, b = f(), f()
    @test a == b
    @test hash(a) == hash(b)
    @test isequal(a, b)
    return @test a != g()
end

@testset "Equality by value" begin
    M = DefaultManifold(3)
    p = [1.0, 0.0, 0.0]
    q = [0.0, 1.0, 0.0]
    # a manifold that stores an array
    N(p) = TangentSpace(M, copy(p))

    @testset "manifolds" begin
        test_value_equality(() -> N(p), () -> N(q))
        test_value_equality(() -> MetricManifold(N(p), EuclideanMetric()), () -> MetricManifold(N(q), EuclideanMetric()))
        test_value_equality(() -> ConnectionManifold(N(p), LeviCivitaConnection()), () -> ConnectionManifold(N(q), LeviCivitaConnection()))
        test_value_equality(() -> EmbeddedManifold(N(p), M), () -> EmbeddedManifold(N(q), M))
        test_value_equality(() -> ProductManifold(N(p), M), () -> ProductManifold(N(q), M))
        test_value_equality(() -> PowerManifold(N(p), NestedPowerRepresentation(), 2), () -> PowerManifold(N(q), NestedPowerRepresentation(), 2))
        # the representation of a power manifold is a type parameter
        @test PowerManifold(M, NestedPowerRepresentation(), 2) != PowerManifold(M, NestedReplacingPowerRepresentation(), 2)
        test_value_equality(() -> ValidationManifold(N(p)), () -> ValidationManifold(N(q)))
        @test ValidationManifold(M) != ValidationManifold(M; ignore_contexts = [:Point])
    end

    @testset "bases" begin
        test_value_equality(() -> DiagonalizingOrthonormalBasis(copy(p)), () -> DiagonalizingOrthonormalBasis(q))
        # the number system is a type parameter
        @test DiagonalizingOrthonormalBasis(p) != DiagonalizingOrthonormalBasis(p, ℂ)
        test_value_equality(() -> get_basis(M, p, DefaultOrthonormalBasis()), () -> get_basis(M, p, DefaultOrthogonalBasis()))
        test_value_equality(
            () -> get_basis(M, p, DiagonalizingOrthonormalBasis(copy(q))),
            () -> get_basis(M, p, DiagonalizingOrthonormalBasis(p)),
        )
        test_value_equality(() -> DiagonalizingBasisData(copy(p), [0.0], [copy(q)]), () -> DiagonalizingBasisData(p, [1.0], [q]))
        test_value_equality(() -> PowerBasisData([copy(p)]), () -> PowerBasisData([q]))
        test_value_equality(() -> ProductBasisData((copy(p), 1)), () -> ProductBasisData((q, 1)))
        P = PowerManifold(M, NestedPowerRepresentation(), 2)
        test_value_equality(() -> get_basis(P, [p, q], DefaultOrthonormalBasis()), () -> get_basis(P, [p, q], DefaultOrthogonalBasis()))
        test_value_equality(() -> TFVector(copy(p), DefaultOrthonormalBasis()), () -> TFVector(q, DefaultOrthonormalBasis()))
        @test TFVector(p, DefaultOrthonormalBasis()) != CoTFVector(p, DefaultOrthonormalBasis())
    end

    @testset "methods" begin
        nl(X0) = NLSolveInverseRetraction(ExponentialRetraction(), copy(X0))
        test_value_equality(() -> nl(p), () -> nl(q))
        test_value_equality(() -> EmbeddedInverseRetraction(nl(p)), () -> EmbeddedInverseRetraction(nl(q)))
        test_value_equality(() -> StabilizedInverseRetraction(nl(p)), () -> StabilizedInverseRetraction(nl(q)))
        test_value_equality(() -> InverseProductRetraction(nl(p), nl(q)), () -> InverseProductRetraction(nl(q), nl(p)))
        test_value_equality(() -> InverseRetractionWithKeywords(nl(p); a = copy(p)), () -> InverseRetractionWithKeywords(nl(p); a = q))
        test_value_equality(
            () -> ShootingInverseRetraction(ExponentialRetraction(), nl(p), ParallelTransport(), 10, 1.0e-10, 100),
            () -> ShootingInverseRetraction(ExponentialRetraction(), nl(q), ParallelTransport(), 10, 1.0e-10, 100),
        )
        test_value_equality(() -> ApproximateLogarithmicInverseRetraction((; a = copy(p))), () -> ApproximateLogarithmicInverseRetraction((; a = q)))
        rk(a) = RetractionWithKeywords(ExponentialRetraction(); a = copy(a))
        test_value_equality(() -> rk(p), () -> rk(q))
        test_value_equality(() -> EmbeddedRetraction(rk(p)), () -> EmbeddedRetraction(rk(q)))
        test_value_equality(() -> StabilizedRetraction(rk(p)), () -> StabilizedRetraction(rk(q)))
        test_value_equality(() -> ProductRetraction(rk(p), rk(q)), () -> ProductRetraction(rk(q), rk(p)))
        test_value_equality(() -> ApproximateExponentialRetraction((; a = copy(p))), () -> ApproximateExponentialRetraction((; a = q)))
        vk(a) = VectorTransportWithKeywords(ParallelTransport(); a = copy(a))
        test_value_equality(() -> vk(p), () -> vk(q))
        test_value_equality(() -> EmbeddedVectorTransport(vk(p)), () -> EmbeddedVectorTransport(vk(q)))
        test_value_equality(() -> ScaledVectorTransport(vk(p)), () -> ScaledVectorTransport(vk(q)))
        test_value_equality(() -> ProductVectorTransport(vk(p), vk(q)), () -> ProductVectorTransport(vk(q), vk(p)))
        test_value_equality(() -> DifferentiatedRetractionVectorTransport(rk(p)), () -> DifferentiatedRetractionVectorTransport(rk(q)))
        test_value_equality(() -> PoleLadderTransport(rk(p), nl(p)), () -> PoleLadderTransport(rk(p), nl(q)))
        test_value_equality(() -> SchildsLadderTransport(rk(p), nl(p)), () -> SchildsLadderTransport(rk(q), nl(p)))
        test_value_equality(() -> VectorTransportDirection(vk(p), rk(p)), () -> VectorTransportDirection(vk(p), rk(q)))
        test_value_equality(() -> VectorTransportTo(vk(p), nl(p)), () -> VectorTransportTo(vk(q), nl(p)))
        test_value_equality(() -> ExtrinsicEstimation(GeodesicInterpolationWithinRadius(1.0)), () -> ExtrinsicEstimation(GeodesicInterpolationWithinRadius(2.0)))
    end

    @testset "points and vectors with forwarded fields" begin
        test_value_equality(() -> ValidationMPoint(copy(p)), () -> ValidationMPoint(q))
        test_value_equality(() -> ValidationTangentVector(copy(p)), () -> ValidationTangentVector(q))
        test_value_equality(() -> ManifoldsBase.Test.DefaultPoint(copy(p)), () -> ManifoldsBase.Test.DefaultPoint(q))
    end
end
