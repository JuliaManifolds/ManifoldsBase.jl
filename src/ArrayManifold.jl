"""
    ArrayManifold{M<:Manifold} <: Manifold

A manifold to encapsulate manifolds working on array representations of `MPoints` and
`TVectors` in a transparent way, such that for these manifolds it's not necessary to
introduce explicit types for the points and tangent vectors, but they are
encapsulated/stripped automatically when needed.

This manifold is a decorator for a manifold, i.e. it decorates a manifold `M` with types
points, vectors, and covectors.
"""
struct ArrayManifold{M<:Manifold} <: Manifold
    manifold::M
end

"""
    ArrayMPoint <: MPoint

Represent a point on an [`ArrayManifold`](@ref), i.e. on a manifold where data can be
represented by arrays. The array is stored internally and semantically. This distinguished
the value from [`ArrayTVector`](@ref)s and [`ArrayCoTVector`](@ref)s.
"""
struct ArrayMPoint{V<:AbstractArray{<:Number}} <: MPoint
    value::V
end

"""
    ArrayTVector <: TVector

Represent a tangent vector to a point on an [`ArrayManifold`](@ref), i.e. on a manifold
where data can be represented by arrays. The array is stored internally and semantically.
This distinguished the value from [`ArrayMPoint`](@ref)s and [`ArrayCoTVector`](@ref)s.
"""
struct ArrayTVector{V<:AbstractArray{<:Number}} <: TVector
    value::V
end

"""
    ArrayCoTVector <: CoTVector

Represent a cotangent vector to a point on an [`ArrayManifold`](@ref), i.e. on a manifold
where data can be represented by arrays. The array is stored internally and semantically.
This distinguished the value from [`ArrayMPoint`](@ref)s and [`ArrayTVector`](@ref)s.
"""
struct ArrayCoTVector{V<:AbstractArray{<:Number}} <: TVector
    value::V
end

(+)(v1::ArrayCoTVector, v2::ArrayCoTVector) = ArrayCoTVector(v1.value + v2.value)
(-)(v1::ArrayCoTVector, v2::ArrayCoTVector) = ArrayCoTVector(v1.value - v2.value)
(-)(v::ArrayCoTVector) = ArrayCoTVector(-v.value)
(*)(a::Number, v::ArrayCoTVector) = ArrayCoTVector(a * v.value)

(+)(v1::ArrayTVector, v2::ArrayTVector) = ArrayTVector(v1.value + v2.value)
(-)(v1::ArrayTVector, v2::ArrayTVector) = ArrayTVector(v1.value - v2.value)
(-)(v::ArrayTVector) = ArrayTVector(-v.value)
(*)(a::Number, v::ArrayTVector) = ArrayTVector(a * v.value)


allocate(x::ArrayMPoint) = ArrayMPoint(allocate(x.value))
allocate(x::ArrayMPoint, ::Type{T}) where {T} = ArrayMPoint(allocate(x.value, T))
allocate(x::ArrayTVector) = ArrayTVector(allocate(x.value))
allocate(x::ArrayTVector, ::Type{T}) where {T} = ArrayTVector(allocate(x.value, T))

"""
    array_value(x)

Return the internal array value of a [`ArrayMPoint`](@ref), [`ArrayTVector`](@ref), or
[`ArrayCoTVector`](@ref) if the value `x` is encapsulated as such. Return `x` if it is
already an array.
"""
array_value(x::AbstractArray) = x
array_value(x::ArrayMPoint) = x.value
array_value(v::ArrayTVector) = v.value
array_value(v::ArrayCoTVector) = v.value

function base_manifold(M::ArrayManifold, depth::Val{N} = Val(-1)) where {N}
    return (N != 0) ? base_manifold(M.manifold, (N > 0) ? Val(N-1) : depth) : M
end

function check_manifold_point(M::ArrayManifold, x::MPoint; kwargs...)
    return check_manifold_point(M.manifold, array_value(x); kwargs...)
end

function check_tangent_vector(M::ArrayManifold, x::MPoint, v::TVector; kwargs...)
    return check_tangent_vector(M.manifold, array_value(x), array_value(v); kwargs...)
end

convert(::Type{V}, v::ArrayCoTVector{V}) where {V<:AbstractArray{<:Number}} = v.value
function convert(::Type{ArrayCoTVector{V}}, v::V) where {V<:AbstractArray{<:Number}}
    return ArrayCoTVector{V}(v)
end
convert(::Type{M}, m::ArrayManifold{M}) where {M<:Manifold} = m.manifold
convert(::Type{ArrayManifold{M}}, m::M) where {M<:Manifold} = ArrayManifold(m)
convert(::Type{V}, x::ArrayMPoint{V}) where {V<:AbstractArray{<:Number}} = x.value
convert(::Type{ArrayMPoint{V}}, x::V) where {V<:AbstractArray{<:Number}} = ArrayMPoint{V}(x)
convert(::Type{V}, v::ArrayTVector{V}) where {V<:AbstractArray{<:Number}} = v.value
function convert(::Type{ArrayTVector{V}}, v::V) where {V<:AbstractArray{<:Number}}
    return ArrayTVector{V}(v)
end

function copyto!(x::ArrayMPoint, y::ArrayMPoint)
    copyto!(x.value, y.value)
    return x
end
function copyto!(x::ArrayCoTVector, y::ArrayCoTVector)
    copyto!(x.value, y.value)
    return x
end
function copyto!(x::ArrayTVector, y::ArrayTVector)
    copyto!(x.value, y.value)
    return x
end

function distance(M::ArrayManifold, x, y; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_manifold_point(M, y, true; kwargs...)
    return distance(M.manifold, array_value(x), array_value(y))
end

function exp(M::ArrayManifold, x, v; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_tangent_vector(M, x, v, true; kwargs...)
    y = ArrayMPoint(exp(M.manifold, array_value(x), array_value(v)))
    is_manifold_point(M, y, true; kwargs...)
    return y
end

function exp!(M::ArrayManifold, y, x, v; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_tangent_vector(M, x, v, true; kwargs...)
    exp!(M.manifold, array_value(y), array_value(x), array_value(v))
    is_manifold_point(M, y, true; kwargs...)
    return y
end

injectivity_radius(M::ArrayManifold) = injectivity_radius(M.manifold)
function injectivity_radius(M::ArrayManifold, method::AbstractRetractionMethod)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(M::ArrayManifold, x; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(x))
end
function injectivity_radius(
    M::ArrayManifold,
    x,
    method::AbstractRetractionMethod;
    kwargs...,
)
    is_manifold_point(M, x, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(x), method)
end
function injectivity_radius(M::ArrayManifold, method::ExponentialRetraction)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(M::ArrayManifold, x, method::ExponentialRetraction; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(x), method)
end

function inner(M::ArrayManifold, x, v, w; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_tangent_vector(M, x, v, true; kwargs...)
    is_tangent_vector(M, x, w, true; kwargs...)
    return inner(M.manifold, array_value(x), array_value(v), array_value(w))
end

function isapprox(M::ArrayManifold, x, y; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_manifold_point(M, y, true; kwargs...)
    return isapprox(M.manifold, array_value(x), array_value(y); kwargs...)
end
function isapprox(M::ArrayManifold, x, v, w; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_tangent_vector(M, x, v, true; kwargs...)
    is_tangent_vector(M, x, w, true; kwargs...)
    return isapprox(M.manifold, array_value(x), array_value(v), array_value(w); kwargs...)
end

function log(M::ArrayManifold, x, y; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_manifold_point(M, y, true; kwargs...)
    v = ArrayTVector(log(M.manifold, array_value(x), array_value(y)))
    is_tangent_vector(M, x, v, true; kwargs...)
    return v
end

function log!(M::ArrayManifold, v, x, y; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    is_manifold_point(M, y, true; kwargs...)
    log!(M.manifold, array_value(v), array_value(x), array_value(y))
    is_tangent_vector(M, x, v, true; kwargs...)
    return v
end

manifold_dimension(M::ArrayManifold) = manifold_dimension(M.manifold)

number_eltype(::Type{ArrayMPoint{V}}) where {V} = number_eltype(V)
number_eltype(x::ArrayMPoint) = number_eltype(x.value)
number_eltype(::Type{ArrayCoTVector{V}}) where {V} = number_eltype(V)
number_eltype(x::ArrayCoTVector) = number_eltype(x.value)
number_eltype(::Type{ArrayTVector{V}}) where {V} = number_eltype(V)
number_eltype(x::ArrayTVector) = number_eltype(x.value)

function project_tangent!(M::ArrayManifold, w, x, v; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    project_tangent!(M.manifold, w.value, array_value(x), array_value(v))
    is_tangent_vector(M, x, w, true; kwargs...)
    return w
end

representation_size(M::ArrayManifold) = representation_size(M.manifold)

similar(x::ArrayMPoint) = ArrayMPoint(similar(x.value))
similar(x::ArrayMPoint, ::Type{T}) where {T} = ArrayMPoint(similar(x.value, T))
similar(x::ArrayCoTVector) = ArrayCoTVector(similar(x.value))
similar(x::ArrayCoTVector, ::Type{T}) where {T} = ArrayCoTVector(similar(x.value, T))
similar(x::ArrayTVector) = ArrayTVector(similar(x.value))
similar(x::ArrayTVector, ::Type{T}) where {T} = ArrayTVector(similar(x.value, T))

function vector_transport_along!(
    M::ArrayManifold,
    vto,
    x,
    v,
    c,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_tangent_vector(M, x, v, true; kwargs...)
    vector_transport_along!(
        M.manifold,
        array_value(vto),
        array_value(x),
        array_value(v),
        c,
        m,
    )
    is_tangent_vector(M, c(1), vto, true; kwargs...)
    return vto
end

function vector_transport_to!(
    M::ArrayManifold,
    vto,
    x,
    v,
    y,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_manifold_point(M, y, true; kwargs...)
    is_tangent_vector(M, x, v, true; kwargs...)
    vector_transport_to!(
        M.manifold,
        array_value(vto),
        array_value(x),
        array_value(v),
        array_value(y),
        m,
    )
    is_tangent_vector(M, y, vto, true; kwargs...)
    return vto
end

function vector_transport_to!(
    M::ArrayManifold,
    vto,
    x,
    v,
    y,
    m::ProjectionTransport;
    kwargs...,
)
    is_manifold_point(M, y, true; kwargs...)
    is_tangent_vector(M, x, v, true; kwargs...)
    vector_transport_to!(
        M.manifold,
        array_value(vto),
        array_value(x),
        array_value(v),
        array_value(y),
        m,
    )
    is_tangent_vector(M, y, vto, true; kwargs...)
    return vto
end

function zero_tangent_vector(M::ArrayManifold, x; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    w = zero_tangent_vector(M.manifold, array_value(x))
    is_tangent_vector(M, x, w, true; kwargs...)
    return w
end

function zero_tangent_vector!(M::ArrayManifold, v, x; kwargs...)
    is_manifold_point(M, x, true; kwargs...)
    zero_tangent_vector!(M.manifold, array_value(v), array_value(x); kwargs...)
    is_tangent_vector(M, x, v, true; kwargs...)
    return v
end
