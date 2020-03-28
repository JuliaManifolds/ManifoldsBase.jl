using ManifoldsBase
using Test

struct ErrorTestManifold <: Manifold{ℝ} end

function ManifoldsBase.check_manifold_point(::ErrorTestManifold, x)
    if any(u -> u < 0, x)
        return DomainError(x, "<0")
    end
    return nothing
end
function ManifoldsBase.check_tangent_vector(M::ErrorTestManifold, x, v)
    mpe = check_manifold_point(M, x)
    mpe === nothing || return mpe
    if any(u -> u < 0, v)
        return DomainError(v, "<0")
    end
    return nothing
end

@testset "Domain errors" begin
    M = ErrorTestManifold()
    @test isa(check_manifold_point(M, [-1, 1]), DomainError)
    @test check_manifold_point(M, [1, 1]) === nothing
    @test !is_manifold_point(M, [-1, 1])
    @test is_manifold_point(M, [1, 1])
    @test_throws DomainError is_manifold_point(M, [-1, 1], true)

    @test isa(check_tangent_vector(M, [1, 1], [-1, 1]), DomainError)
    @test check_tangent_vector(M, [1, 1], [1, 1]) === nothing
    @test !is_tangent_vector(M, [1, 1], [-1, 1])
    @test is_tangent_vector(M, [1, 1], [1, 1])
    @test_throws DomainError is_tangent_vector(M, [1, 1], [-1, 1], true)
end
