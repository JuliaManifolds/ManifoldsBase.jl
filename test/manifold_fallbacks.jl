using Test
using ManifoldsBase

struct NonManifold <: AbstractManifold{ManifoldsBase.ℝ} end

@testset "NotImplemented Errors" begin
    M = NonManifold()
    p = [1.0]
    q = similar(p)
    X = [2.0]
    Y = similar(X)
    for B in [
        VeeOrthogonalBasis(),
        DefaultBasis(),
        DefaultOrthogonalBasis(),
        DefaultOrthonormalBasis(),
        DiagonalizingOrthonormalBasis(X),
        CachedBasis(DefaultBasis(), X),
    ]
        if !(B isa CachedBasis)
            @test_throws MethodError get_basis(M, p, B)
            @test_throws MethodError get_vector(M, p, X, B)
            @test_throws MethodError get_vector!(M, Y, p, X, B)
        else
            @test get_basis(M, p, B) == B
            @test get_vector(M, p, X, B) == 4 # since we have 1 vector
            @test get_vector!(M, Y, p, X, B) == [4] # since Y is a vector
        end
        @test_throws MethodError get_coordinates(M, p, X, B)
        @test_throws MethodError get_coordinates!(M, Y, p, X, B)
    end
    @test_throws MethodError inverse_retract(M, p, q)
    @test_throws MethodError inverse_retract!(M, Y, p, q)
    @test_throws ErrorException PadeInverseRetraction(0)
    for IR in [
        CayleyInverseRetraction(),
        EmbeddedInverseRetraction(ProjectionInverseRetraction()),
        LogarithmicInverseRetraction(),
        PadeInverseRetraction(2),
        PolarInverseRetraction(),
        ProjectionInverseRetraction(),
        QRInverseRetraction(),
        NLSolveInverseRetraction(ExponentialRetraction()),
        SoftmaxInverseRetraction(),
    ]
        @test_throws MethodError inverse_retract(M, p, q, IR)
        @test_throws MethodError inverse_retract!(M, Y, p, q, IR)
    end
    @test_throws MethodError retract(M, p, X)
    @test_throws MethodError retract!(M, q, p, X)
    @test_throws ErrorException PadeRetraction(0)
    for R in [
        CayleyRetraction(),
        EmbeddedRetraction(ProjectionRetraction()),
        ExponentialRetraction(),
        ODEExponentialRetraction(ProjectionRetraction(), DefaultBasis()),
        ODEExponentialRetraction(ProjectionRetraction()),
        PadeRetraction(2),
        PolarRetraction(),
        ProjectionRetraction(),
        QRRetraction(),
        SoftmaxRetraction(),
        CayleyRetraction(),
        PadeRetraction(2),
    ]
        @test_throws MethodError retract(M, p, X, R)
        @test_throws MethodError retract(M, p, X, 1.0, R)
        @test_throws MethodError retract!(M, q, p, X, R)
        @test_throws MethodError retract!(M, q, p, X, 1.0, R)
    end
    for VT in [
        DifferentiatedRetractionVectorTransport(ExponentialRetraction()),
        ProjectionTransport(),
        ParallelTransport(),
    ]
        @test_throws MethodError vector_transport_along(M, p, X, :curve, VT)
        @test_throws MethodError vector_transport_along!(M, Y, p, X, :curve, VT)
        @test_throws MethodError vector_transport_along!(M, Y, p, X, [p], VT)
        @test_throws MethodError vector_transport_direction(M, p, X, X, VT)
        @test_throws MethodError vector_transport_direction!(M, Y, p, X, X, VT)
        @test_throws MethodError vector_transport_to(M, p, X, q, VT)
        @test_throws MethodError vector_transport_to!(M, Y, p, X, q, VT)
    end
end

@testset "Default Fallbacks and Error Messages" begin
    M = ManifoldsBase.DefaultManifold(3)
    p = [1.0, 0.0, 0.0]
    @test number_of_coordinates(M, ManifoldsBase.ℝ) == 3
    B = get_basis(M, p, DefaultBasis())
    @test_throws DomainError ODEExponentialRetraction(ProjectionRetraction(), B)
    @test_throws DomainError ODEExponentialRetraction(
        ExponentialRetraction(),
        DefaultBasis(),
    )
    @test_throws DomainError ODEExponentialRetraction(ExponentialRetraction(), B)
    @test_throws ErrorException PadeRetraction(0)
end

@testset "Estimation errors" begin
    @test_throws DomainError GeodesicInterpolationWithinRadius(-1)
end
