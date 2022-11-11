using ManifoldsBase
using ManifoldsBase: ℝ
using Test

struct ErrorTestManifold <: AbstractManifold{ℝ} end

function ManifoldsBase.check_size(::ErrorTestManifold, p)
    size(p) != (2,) && return DomainError(size(p), "size $p not (2,)")
    return nothing
end
function ManifoldsBase.check_size(::ErrorTestManifold, p, X)
    size(X) != (2,) && return DomainError(size(X), "size $X not (2,)")
    return nothing
end
function ManifoldsBase.check_point(::ErrorTestManifold, x)
    if any(u -> u < 0, x)
        return DomainError(x, "<0")
    end
    return nothing
end
function ManifoldsBase.check_vector(M::ErrorTestManifold, x, v)
    mpe = ManifoldsBase.check_point(M, x)
    mpe === nothing || return mpe
    if any(u -> u < 0, v)
        return DomainError(v, "<0")
    end
    return nothing
end

@testset "Domain errors" begin
    M = ErrorTestManifold()
    @test isa(ManifoldsBase.check_point(M, [-1, 1]), DomainError)
    @test isa(ManifoldsBase.check_size(M, [-1, 1, 1]), DomainError)
    @test isa(ManifoldsBase.check_size(M, [-1, 1], [1, 1, 1]), DomainError)
    @test ManifoldsBase.check_point(M, [1, 1]) === nothing
    @test !is_point(M, [-1, 1])
    @test !is_point(M, [1, 1, 1]) # checksize fails
    @test_throws DomainError is_point(M, [-1, 1, 1], true) # checksize errors
    @test_throws DomainError is_point(M, [-1, 1, 1], :error) # checksize errors
    cs = "DomainError with (3,)\nsize [-1, 1, 1] not (2,)"
    @test_logs (:info, cs) is_point(M, [-1,1,1], :info)
    @test_logs (:warn, cs) is_point(M, [-1,1,1], :warn)
    @test is_point(M, [1, 1])
    @test_throws DomainError is_point(M, [-1, 1], true)
    @test_throws DomainError is_point(M, [-1, 1], :error)
    ps = "DomainError with [-1, 1]\n<0"
    @test_logs (:info, ps)  is_point(M, [-1, 1], :info)
    @test_logs (:warn, ps)  is_point(M, [-1, 1], :warn)

    @test isa(ManifoldsBase.check_vector(M, [1, 1], [-1, 1]), DomainError)
    @test ManifoldsBase.check_vector(M, [1, 1], [1, 1]) === nothing
    @test !is_vector(M, [1, 1], [-1, 1])
    @test !is_vector(M, [1, 1], [1, 1, 1])
    @test_throws DomainError is_vector(M, [1, 1], [-1, 1, 1], true)
    @test_throws DomainError is_vector(M, [1, 1], [-1, 1, 1], :error)
    vs = "DomainError with (3,)\nsize [-1, 1, 1] not (2,)"
    @test_logs (:info, vs)  is_vector(M, [1, 1], [-1, 1, 1], :info)
    @test_logs (:warn, vs)  is_vector(M, [1, 1], [-1, 1, 1], :warn)
    @test !is_vector(M, [1, 1, 1], [1, 1, 1], false, true)
    @test_throws DomainError is_vector(M, [1, 1, 1], [1, 1], true, true)
    @test_throws DomainError is_vector(M, [1, 1, 1], [1, 1], :error, true)
    ps2 = "DomainError with (3,)\nsize [1, 1, 1] not (2,)"
    @test_logs (:info, ps2)  is_vector(M, [1, 1, 1], [1, 1], :info, true)
    @test_logs (:warn, ps2)  is_vector(M, [1, 1, 1], [1, 1], :warn, true)
    @test is_vector(M, [1, 1], [1, 1])
    @test_throws DomainError is_vector(M, [1, 1], [-1, 1], true)
end
