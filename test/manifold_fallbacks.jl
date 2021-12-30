using Test
using ManifoldsBase

struct NonManifold <: AbstractManifold{ManifoldsBase.ℝ} end

@testset "NotImplemented Errors" begin
    M = NonManifold()
    p = [1]
    q = similar(p)
    X = [2]
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
    for IR in [
        LogarithmicInverseRetraction(),
        PolarInverseRetraction(),
        ProjectionInverseRetraction(),
        QRInverseRetraction(),
    ]
        @test_throws MethodError inverse_retract(M, p, q, IR)
        @test_throws MethodError inverse_retract!(M, Y, p, q, IR)
    end
    @test_throws MethodError retract(M, p, X)
    @test_throws MethodError retract!(M, q, p, X)
    for R in
        [ExponentialRetraction(), PolarRetraction(), ProjectionRetraction(), QRRetraction()]
        @test_throws MethodError retract(M, p, X, R)
        @test_throws MethodError retract!(M, q, p, X, R)
    end
end

@testset "Default Fallbacks and Error Messages" begin
    M = ManifoldsBase.DefaultManifold(3)
    @test number_of_coordinates(M, ManifoldsBase.ℝ) == 3
end
