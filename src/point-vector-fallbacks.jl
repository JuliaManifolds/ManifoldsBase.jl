#
# Introduce default fallbacks for all basic functions on manifolds for
# (default) Manifold types
#
function angle(
    M::AbstractManifold,
    p::P,
    X::T,
    Y::T,
) where {P<:AbstractManifoldPoint,T<:TVector}
    return angle(M, p.value, X.value, Y.value)
end

function check_point(
    M::AbstractManifold,
    p::P,
    X::T;
    kwargs...,
) where {P<:AbstractManifoldPoint,T<:TVector}
    return check_point(M, p.value, X.value; kwargs...)
end

function distance(M::AbstractManifold, p::P, q::P) where {P<:AbstractManifoldPoint}
    return distance(M, p.value, q.value)
end

function embed!(M::AbstractManifold, q, p::P) where {P<:AbstractManifoldPoint}
    return embed!(M, q, p.value)
end

function embed!(
    M::AbstractManifold,
    Y,
    p::P,
    X::T,
) where {P<:AbstractManifoldPoint,T<:TVector}
    return embed!(M, Y, p.value, X.value)
end

function exp!(
    M::AbstractManifold,
    q::P,
    p::P,
    X::T,
) where {P<:AbstractManifoldPoint,T<:TVector}
    exp!(M, q.value, p.value, X.value)
    return q
end

function inner(
    M::AbstractManifold,
    p::P,
    X::T,
    Y::T,
) where {P<:AbstractManifoldPoint,T<:TVector}
    return inner(M, p.value, X.value, Y.value)
end

function isapprox(
    M::AbstractManifold,
    p::P,
    q::P;
    kwargs...,
) where {P<:AbstractManifoldPoint}
    return isapprox(M, p.value, q.value; kwargs...)
end

function isapprox(
    M::AbstractManifold,
    p::P,
    X::T,
    Y::T;
    kwargs...,
) where {P<:AbstractManifoldPoint,T<:TVector}
    return isapprox(M, p.value, X.value, Y.value; kwargs...)
end

function loq!(
    M::AbstractManifold,
    X::T,
    p::P,
    q::P,
) where {P<:AbstractManifoldPoint,T<:TVector}
    log!(M, X.value, p.value, q.value)
    return X
end
