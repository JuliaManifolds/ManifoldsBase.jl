using ManifoldsBase
using Test

struct ErrorTestManifold <: AbstractManifold{ℝ} end

function ManifoldsBase.check_point(::ErrorTestManifold, x; kwargs...)
    if any(u -> u < 0, x)
        return DomainError(x, "<0")
    end
    return nothing
end
function ManifoldsBase.check_vector(M::ErrorTestManifold, x, v; kwargs...)
    mpe = check_point(M, x)
    mpe === nothing || return mpe
    if any(u -> u < 0, v)
        return DomainError(v, "<0")
    end
    return nothing
end

@testset "Domain errors" begin
    M = ErrorTestManifold()
    @test isa(check_point(M, [-1, 1]), DomainError)
    @test check_point(M, [1, 1]) === nothing
    @test !is_point(M, [-1.0, 1.0])
    @test is_point(M, [1.0, 1.0])
    @test_throws DomainError is_point(M, [-1.0, 1.0], true)

    @test isa(check_vector(M, [1, 1], [-1, 1]), DomainError)
    @test check_vector(M, [1, 1], [1, 1]) === nothing
    @test !is_vector(M, [1.0, 1.0], [-1.0, 1.0])
    @test is_vector(M, [1.0, 1.0], [1.0, 1.0])
    @test_throws DomainError is_vector(M, [1.0, 1.0], [-1.0, 1.0], true)
end
