using ManifoldsBase
using Test

struct ErrorTestManifold <: Manifold end

function ManifoldsBase.manifold_point_error(::ErrorTestManifold, x)
    if any(u -> u < 0, x)
        return DomainError(x, "<0")
    end
    return nothing
end
function ManifoldsBase.tangent_vector_error(M::ErrorTestManifold, x, v)
    mpe = manifold_point_error(M, x)
    mpe === nothing || return mpe
    if any(u -> u < 0, v)
        return DomainError(v, "<0")
    end
    return nothing
end

@testset "Domain errors" begin
    M = ErrorTestManifold()
    @test isa(manifold_point_error(M, [-1, 1]), DomainError)
    @test manifold_point_error(M, [1, 1]) === nothing
    @test !is_manifold_point(M, [-1, 1])
    @test is_manifold_point(M, [1, 1])
    @test_throws DomainError check_manifold_point(M, [-1, 1])

    @test isa(tangent_vector_error(M, [1, 1], [-1, 1]), DomainError)
    @test tangent_vector_error(M, [1, 1], [1, 1]) === nothing
    @test !is_tangent_vector(M, [1, 1], [-1, 1])
    @test is_tangent_vector(M, [1, 1], [1, 1])
    @test_throws DomainError check_tangent_vector(M, [1, 1], [-1, 1])
end
