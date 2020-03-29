"""
    ArrayManifold{M<:Manifold} <: Manifold

A manifold to encapsulate manifolds working on array representations of `MPoints` and
`TVectors` in a transparent way, such that for these manifolds it's not necessary to
introduce explicit types for the points and tangent vectors, but they are
encapsulated/stripped automatically when needed.

This manifold is a decorator for a manifold, i.e. it decorates a manifold `M` with types
points, vectors, and covectors.
"""
struct ArrayManifold{𝔽,M<:Manifold{𝔽}} <: AbstractDecoratorManifold{𝔽}
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

(+)(X::ArrayCoTVector, Y::ArrayCoTVector) = ArrayCoTVector(X.value + Y.value)
(-)(X::ArrayCoTVector, Y::ArrayCoTVector) = ArrayCoTVector(X.value - Y.value)
(-)(X::ArrayCoTVector) = ArrayCoTVector(-X.value)
(*)(a::Number, X::ArrayCoTVector) = ArrayCoTVector(a * X.value)

(+)(X::ArrayTVector, Y::ArrayTVector) = ArrayTVector(X.value + Y.value)
(-)(X::ArrayTVector, Y::ArrayTVector) = ArrayTVector(X.value - Y.value)
(-)(X::ArrayTVector) = ArrayTVector(-X.value)
(*)(a::Number, X::ArrayTVector) = ArrayTVector(a * X.value)


allocate(p::ArrayMPoint) = ArrayMPoint(allocate(p.value))
allocate(p::ArrayMPoint, ::Type{T}) where {T} = ArrayMPoint(allocate(p.value, T))
allocate(X::ArrayTVector) = ArrayTVector(allocate(X.value))
allocate(X::ArrayTVector, ::Type{T}) where {T} = ArrayTVector(allocate(X.value, T))

"""
    array_value(p)

Return the internal array value of an [`ArrayMPoint`](@ref), [`ArrayTVector`](@ref), or
[`ArrayCoTVector`](@ref) if the value `p` is encapsulated as such. Return `p` if it is
already an array.
"""
array_value(p::AbstractArray) = p
array_value(p::ArrayMPoint) = p.value
array_value(X::ArrayTVector) = X.value
array_value(ξ::ArrayCoTVector) = ξ.value

function check_manifold_point(M::ArrayManifold, p; kwargs...)
    return check_manifold_point(M.manifold, array_value(p); kwargs...)
end
function check_manifold_point(M::ArrayManifold, p::MPoint; kwargs...)
    return check_manifold_point(M.manifold, array_value(p); kwargs...)
end

function check_tangent_vector(M::ArrayManifold, p, X; kwargs...)
    return check_tangent_vector(M.manifold, array_value(p), array_value(X); kwargs...)
end
function check_tangent_vector(M::ArrayManifold, p::MPoint, X::TVector; kwargs...)
    return check_tangent_vector(M.manifold, array_value(p), array_value(X); kwargs...)
end

convert(::Type{V}, X::ArrayCoTVector{V}) where {V<:AbstractArray{<:Number}} = X.value
function convert(::Type{ArrayCoTVector{V}}, X::V) where {V<:AbstractArray{<:Number}}
    return ArrayCoTVector{V}(X)
end
convert(::Type{M}, m::ArrayManifold{𝔽,M}) where {𝔽,M<:Manifold{𝔽}} = m.manifold
convert(::Type{ArrayManifold{𝔽,M}}, m::M) where {𝔽,M<:Manifold{𝔽}} = ArrayManifold(m)
convert(::Type{V}, p::ArrayMPoint{V}) where {V<:AbstractArray{<:Number}} = p.value
convert(::Type{ArrayMPoint{V}}, x::V) where {V<:AbstractArray{<:Number}} = ArrayMPoint{V}(x)
convert(::Type{V}, X::ArrayTVector{V}) where {V<:AbstractArray{<:Number}} = X.value
function convert(::Type{ArrayTVector{V}}, X::V) where {V<:AbstractArray{<:Number}}
    return ArrayTVector{V}(X)
end

function copyto!(p::ArrayMPoint, q::ArrayMPoint)
    copyto!(p.value, q.value)
    return p
end
function copyto!(p::ArrayCoTVector, q::ArrayCoTVector)
    copyto!(p.value, q.value)
    return p
end
function copyto!(Y::ArrayTVector, X::ArrayTVector)
    copyto!(Y.value, X.value)
    return Y
end

function distance(M::ArrayManifold, p, q; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_manifold_point(M, q, true; kwargs...)
    return distance(M.manifold, array_value(p), array_value(q))
end

function exp(M::ArrayManifold, p, X; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    y = exp(M.manifold, array_value(p), array_value(X))
    is_manifold_point(M, y, true; kwargs...)
    return ArrayMPoint(y)
end

function exp!(M::ArrayManifold, q, p, X; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    exp!(M.manifold, array_value(q), array_value(p), array_value(X))
    is_manifold_point(M, q, true; kwargs...)
    return q
end

function get_basis(M::ArrayManifold, p, B::AbstractBasis; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    Ξ = get_basis(M.manifold, array_value(p), B)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    if N != manifold_dimension(M.manifold)
        throw(ErrorException(
            "For a basis of the tangent space at $(p) of $(M.manifold), $(manifold_dimension(M)) vectors are required, but get_basis $(B) computed $(N)"
        ))
    end
    # check that the vectors are linearly independent\
    bv_rank = rank(reduce(hcat, bvectors))
    if N != bv_rank
        throw(ErrorException(
            "For a basis of the tangent space at $(p) of $(M.manifold), $(manifold_dimension(M)) linearly independent vectors are required, but get_basis $(B) computed $(bv_rank)"
        ))
    end
    map(X -> is_tangent_vector(M, p, X, true; kwargs...), bvectors)
    return Ξ
end
function get_basis(
    M::ArrayManifold{𝔽},
    p,
    B::Union{AbstractOrthogonalBasis,CachedBasis{𝔽,<:AbstractOrthogonalBasis{𝔽}}};
    kwargs...,
) where {𝔽}
    is_manifold_point(M, p, true; kwargs...)
    Ξ = invoke(get_basis, Tuple{ArrayManifold,Any,AbstractBasis}, M, p, B; kwargs...)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    for i = 1:N
        for j = i+1:N
            dot_val = real(inner(M, p, bvectors[i], bvectors[j]))
            if !isapprox(dot_val, 0; atol = eps(eltype(p)))
                throw(ArgumentError("vectors number $i and $j are not orthonormal (inner product = $dot_val)"))
            end
        end
    end
    return Ξ
end
function get_basis(
    M::ArrayManifold{𝔽},
    p,
    B::Union{AbstractOrthonormalBasis{𝔽},CachedBasis{𝔽,<:AbstractOrthonormalBasis{𝔽}}};
    kwargs...,
) where {𝔽}
    is_manifold_point(M, p, true; kwargs...)
    Ξ = invoke(get_basis, Tuple{ArrayManifold,Any,AbstractOrthogonalBasis}, M, p, B; kwargs...)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    for i = 1:N
        Xi_norm = norm(M, p, bvectors[i])
        if !isapprox(Xi_norm, 1)
            throw(ArgumentError("vector number $i is not normalized (norm = $Xi_norm)"))
        end
    end
    return Ξ
end
for BT in DISAMBIGUATION_BASIS_TYPES
    if BT <: Union{AbstractOrthonormalBasis,CachedBasis{𝔽,<:AbstractOrthonormalBasis} where 𝔽}
        CT = AbstractOrthonormalBasis
    elseif BT <: Union{AbstractOrthogonalBasis,CachedBasis{𝔽,<:AbstractOrthogonalBasis} where 𝔽}
        CT = AbstractOrthogonalBasis
    else
        CT = AbstractBasis
    end
    eval(quote
        @invoke_maker 3 $CT get_basis(M::ArrayManifold, p, B::$BT; kwargs...)
    end)
end

function get_coordinates(M::ArrayManifold, p, X, B::AbstractBasis; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    return get_coordinates(M.manifold, p, X, B)
end
for BT in DISAMBIGUATION_BASIS_TYPES
    eval(quote
        @invoke_maker 4 AbstractBasis get_coordinates(M::ArrayManifold, p, X, B::$BT; kwargs...)
    end)
end

function get_coordinates!(M::ArrayManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    get_coordinates!(M.manifold, Y, p, X, B)
    return Y
end
for BT in DISAMBIGUATION_BASIS_TYPES
    eval(quote
        @invoke_maker 5 AbstractBasis get_coordinates!(M::ArrayManifold, Y, p, X, B::$BT; kwargs...)
    end)
end

function get_vector(M::ArrayManifold, p, X, B::AbstractBasis; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    size(X) == (manifold_dimension(M),) || error("Incorrect size of coefficient vector X")
    Y = get_vector(M.manifold, p, X, B)
    size(Y) == representation_size(M) || error("Incorrect size of tangent vector Y")
    return Y
end
for BT in DISAMBIGUATION_BASIS_TYPES
    eval(quote
        @invoke_maker 4 AbstractBasis get_vector(M::ArrayManifold, p, X, B::$BT; kwargs...)
    end)
end

function get_vector!(M::ArrayManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    size(X) == (manifold_dimension(M),) || error("Incorrect size of coefficient vector X")
    get_vector!(M.manifold, Y, p, X, B)
    size(Y) == representation_size(M) || error("Incorrect size of tangent vector Y")
    return Y
end
for BT in DISAMBIGUATION_BASIS_TYPES
    eval(quote
        @invoke_maker 5 AbstractBasis get_vector!(M::ArrayManifold, Y, p, X, B::$BT; kwargs...)
    end)
end

injectivity_radius(M::ArrayManifold) = injectivity_radius(M.manifold)
function injectivity_radius(M::ArrayManifold, method::AbstractRetractionMethod)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(M::ArrayManifold, p; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(p))
end
function injectivity_radius(
    M::ArrayManifold,
    p,
    method::AbstractRetractionMethod;
    kwargs...,
)
    is_manifold_point(M, p, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(p), method)
end
function injectivity_radius(M::ArrayManifold, method::ExponentialRetraction)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(M::ArrayManifold, p, method::ExponentialRetraction; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(p), method)
end

function inner(M::ArrayManifold, p, X, Y; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    is_tangent_vector(M, p, Y, true; kwargs...)
    return inner(M.manifold, array_value(p), array_value(X), array_value(Y))
end

function isapprox(M::ArrayManifold, p, q; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_manifold_point(M, q, true; kwargs...)
    return isapprox(M.manifold, array_value(p), array_value(q); kwargs...)
end
function isapprox(M::ArrayManifold, p, X, Y; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    is_tangent_vector(M, p, Y, true; kwargs...)
    return isapprox(M.manifold, array_value(p), array_value(X), array_value(Y); kwargs...)
end

function log(M::ArrayManifold, p, q; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_manifold_point(M, q, true; kwargs...)
    X = log(M.manifold, array_value(p), array_value(q))
    is_tangent_vector(M, p, X, true; kwargs...)
    return ArrayTVector(X)
end

function log!(M::ArrayManifold, X, p, q; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    is_manifold_point(M, q, true; kwargs...)
    log!(M.manifold, array_value(X), array_value(p), array_value(q))
    is_tangent_vector(M, p, X, true; kwargs...)
    return X
end

number_eltype(::Type{ArrayMPoint{V}}) where {V} = number_eltype(V)
number_eltype(p::ArrayMPoint) = number_eltype(p.value)
number_eltype(::Type{ArrayCoTVector{V}}) where {V} = number_eltype(V)
number_eltype(p::ArrayCoTVector) = number_eltype(p.value)
number_eltype(::Type{ArrayTVector{V}}) where {V} = number_eltype(V)
number_eltype(X::ArrayTVector) = number_eltype(X.value)

function project!(M::ArrayManifold, Y, p, X; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    project!(M.manifold, array_value(Y), array_value(p), array_value(X))
    is_tangent_vector(M, p, Y, true; kwargs...)
    return Y
end

similar(p::ArrayMPoint) = ArrayMPoint(similar(p.value))
similar(p::ArrayMPoint, ::Type{T}) where {T} = ArrayMPoint(similar(p.value, T))
similar(p::ArrayCoTVector) = ArrayCoTVector(similar(p.value))
similar(p::ArrayCoTVector, ::Type{T}) where {T} = ArrayCoTVector(similar(p.value, T))
similar(X::ArrayTVector) = ArrayTVector(similar(X.value))
similar(X::ArrayTVector, ::Type{T}) where {T} = ArrayTVector(similar(X.value, T))

function vector_transport_along!(
    M::ArrayManifold,
    Y,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_tangent_vector(M, p, X, true; kwargs...)
    vector_transport_along!(
        M.manifold,
        array_value(Y),
        array_value(p),
        array_value(X),
        c,
        m,
    )
    is_tangent_vector(M, c(1), Y, true; kwargs...)
    return Y
end

function vector_transport_to!(
    M::ArrayManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_manifold_point(M, q, true; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    vector_transport_to!(
        M.manifold,
        array_value(Y),
        array_value(p),
        array_value(X),
        array_value(q),
        m,
    )
    is_tangent_vector(M, q, Y, true; kwargs...)
    return Y
end

function vector_transport_to!(
    M::ArrayManifold,
    Y,
    p,
    X,
    q,
    m::ProjectionTransport;
    kwargs...,
)
    is_manifold_point(M, q, true; kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    vector_transport_to!(
        M.manifold,
        array_value(Y),
        array_value(p),
        array_value(X),
        array_value(q),
        m,
    )
    is_tangent_vector(M, q, Y, true; kwargs...)
    return Y
end

function zero_tangent_vector(M::ArrayManifold, p; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    w = zero_tangent_vector(M.manifold, array_value(p))
    is_tangent_vector(M, p, w, true; kwargs...)
    return w
end

function zero_tangent_vector!(M::ArrayManifold, X, p; kwargs...)
    is_manifold_point(M, p, true; kwargs...)
    zero_tangent_vector!(M.manifold, array_value(X), array_value(p); kwargs...)
    is_tangent_vector(M, p, X, true; kwargs...)
    return X
end
