"""
    ValidationManifold{𝔽,M<:AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}

A manifold to encapsulate manifolds working on array representations of [`AbstractManifoldPoint`](@ref)s
and [`TVector`](@ref)s in a transparent way, such that for these manifolds it's not
necessary to introduce explicit types for the points and tangent vectors, but they are
encapsulated/stripped automatically when needed.

This manifold is a decorator for a manifold, i.e. it decorates a [`AbstractManifold`](@ref) `M`
with types points, vectors, and covectors.

# Constructor

    ValidationManifold(M::AbstractManifold; error::Symbol = :error)

Generate the Validation manifold, where `error` is used as the symbol passed to all checks.
This `:error`s by default but could also be set to `:warn` for example
"""
struct ValidationManifold{𝔽,M<:AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}
    manifold::M
    mode::Symbol
end
function ValidationManifold(M::AbstractManifold; error::Symbol = :error)
    return ValidationManifold(M, error)
end

"""
    ValidationMPoint <: AbstractManifoldPoint

Represent a point on an [`ValidationManifold`](@ref), i.e. on a manifold where data can be
represented by arrays. The array is stored internally and semantically. This distinguished
the value from [`ValidationTVector`](@ref)s and [`ValidationCoTVector`](@ref)s.
"""
struct ValidationMPoint{V} <: AbstractManifoldPoint
    value::V
end

"""
    ValidationFibreVector{TType<:VectorSpaceType} <: AbstractFibreVector{TType}

Represent a tangent vector to a point on an [`ValidationManifold`](@ref), i.e. on a manifold
where data can be represented by arrays. The array is stored internally and semantically.
This distinguished the value from [`ValidationMPoint`](@ref)s vectors of other types.
"""
struct ValidationFibreVector{TType<:VectorSpaceType,V} <: AbstractFibreVector{TType}
    value::V
end
function ValidationFibreVector{TType}(value::V) where {TType,V}
    return ValidationFibreVector{TType,V}(value)
end

"""
    ValidationTVector = ValidationFibreVector{TangentSpaceType}

Represent a tangent vector to a point on an [`ValidationManifold`](@ref), i.e. on a manifold
where data can be represented by arrays. The array is stored internally and semantically.
This distinguished the value from [`ValidationMPoint`](@ref)s vectors of other types.
"""
const ValidationTVector = ValidationFibreVector{TangentSpaceType}

"""
    ValidationCoTVector = ValidationFibreVector{CotangentSpaceType}

Represent a cotangent vector to a point on an [`ValidationManifold`](@ref), i.e. on a manifold
where data can be represented by arrays. The array is stored internally and semantically.
This distinguished the value from [`ValidationMPoint`](@ref)s vectors of other types.
"""
const ValidationCoTVector = ValidationFibreVector{CotangentSpaceType}

@eval @manifold_vector_forwards ValidationFibreVector{TType} TType value

@eval @manifold_element_forwards ValidationMPoint value

@inline function active_traits(f, ::ValidationManifold, ::Any...)
    return merge_traits(IsExplicitDecorator())
end

"""
    array_value(p)

Return the internal array value of an [`ValidationMPoint`](@ref), [`ValidationTVector`](@ref), or
[`ValidationCoTVector`](@ref) if the value `p` is encapsulated as such. Return `p` if it is
already an array.
"""
array_value(p::AbstractArray) = p
array_value(p::ValidationMPoint) = p.value
array_value(X::ValidationFibreVector) = X.value

decorated_manifold(M::ValidationManifold) = M.manifold

convert(::Type{M}, m::ValidationManifold{𝔽,M}) where {𝔽,M<:AbstractManifold{𝔽}} = m.manifold
function convert(::Type{ValidationManifold{𝔽,M}}, m::M) where {𝔽,M<:AbstractManifold{𝔽}}
    return ValidationManifold(m)
end
convert(::Type{V}, p::ValidationMPoint{V}) where {V<:AbstractArray} = p.value
function convert(::Type{ValidationMPoint{V}}, x::V) where {V<:AbstractArray}
    return ValidationMPoint{V}(x)
end
function convert(
    ::Type{V},
    X::ValidationFibreVector{TType,V},
) where {TType,V<:AbstractArray}
    return X.value
end
function convert(
    ::Type{ValidationFibreVector{TType,V}},
    X::V,
) where {TType,V<:AbstractArray}
    return ValidationFibreVector{TType,V}(X)
end

function copyto!(M::ValidationManifold, q::ValidationMPoint, p::ValidationMPoint; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    copyto!(M.manifold, q.value, p.value)
    is_point(M, q; error = M.mode, kwargs...)
    return q
end
function copyto!(
    M::ValidationManifold,
    Y::ValidationFibreVector{TType},
    p::ValidationMPoint,
    X::ValidationFibreVector{TType};
    kwargs...,
) where {TType}
    is_point(M, p; error = M.mode, kwargs...)
    copyto!(M.manifold, Y.value, p.value, X.value)
    return p
end

function distance(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_point(M, q; error = M.mode, kwargs...)
    return distance(M.manifold, array_value(p), array_value(q))
end

function exp(M::ValidationManifold, p, X; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    y = exp(M.manifold, array_value(p), array_value(X))
    is_point(M, y; error = M.mode, kwargs...)
    return ValidationMPoint(y)
end

function exp!(M::ValidationManifold, q, p, X; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    exp!(M.manifold, array_value(q), array_value(p), array_value(X))
    is_point(M, q; error = M.mode, kwargs...)
    return q
end

function get_basis(M::ValidationManifold, p, B::AbstractBasis; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    Ξ = get_basis(M.manifold, array_value(p), B)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    if N != manifold_dimension(M.manifold)
        throw(
            ErrorException(
                "For a basis of the tangent space at $(p) of $(M.manifold), $(manifold_dimension(M)) vectors are required, but get_basis $(B) computed $(N)",
            ),
        )
    end
    # check that the vectors are linearly independent\
    bv_rank = rank(reduce(hcat, bvectors))
    if N != bv_rank
        throw(
            ErrorException(
                "For a basis of the tangent space at $(p) of $(M.manifold), $(manifold_dimension(M)) linearly independent vectors are required, but get_basis $(B) computed $(bv_rank)",
            ),
        )
    end
    map(X -> is_vector(M, p, X; error = M.mode, kwargs...), bvectors)
    return Ξ
end
function get_basis(
    M::ValidationManifold,
    p,
    B::Union{AbstractOrthogonalBasis,CachedBasis{𝔽,<:AbstractOrthogonalBasis{𝔽}}};
    kwargs...,
) where {𝔽}
    is_point(M, p; error = M.mode, kwargs...)
    Ξ = invoke(get_basis, Tuple{ValidationManifold,Any,AbstractBasis}, M, p, B; kwargs...)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    for i in 1:N
        for j in (i + 1):N
            dot_val = real(inner(M, p, bvectors[i], bvectors[j]))
            if !isapprox(dot_val, 0; atol = eps(eltype(p)))
                throw(
                    ArgumentError(
                        "vectors number $i and $j are not orthonormal (inner product = $dot_val)",
                    ),
                )
            end
        end
    end
    return Ξ
end
function get_basis(
    M::ValidationManifold,
    p,
    B::Union{AbstractOrthonormalBasis,CachedBasis{𝔽,<:AbstractOrthonormalBasis{𝔽}}};
    kwargs...,
) where {𝔽}
    is_point(M, p; error = M.mode, kwargs...)
    get_basis_invoke_types = Tuple{
        ValidationManifold,
        Any,
        Union{
            AbstractOrthogonalBasis,
            CachedBasis{𝔽2,<:AbstractOrthogonalBasis{𝔽2}},
        } where {𝔽2},
    }
    Ξ = invoke(get_basis, get_basis_invoke_types, M, p, B; kwargs...)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    for i in 1:N
        Xi_norm = norm(M, p, bvectors[i])
        if !isapprox(Xi_norm, 1)
            throw(ArgumentError("vector number $i is not normalized (norm = $Xi_norm)"))
        end
    end
    return Ξ
end

function get_coordinates(M::ValidationManifold, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p; error = :error, kwargs...)
    is_vector(M, p, X; error = :error, kwargs...)
    return get_coordinates(M.manifold, p, X, B)
end

function get_coordinates!(M::ValidationManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    get_coordinates!(M.manifold, Y, p, X, B)
    return Y
end

function get_vector(M::ValidationManifold, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    size(X) == (manifold_dimension(M),) || error("Incorrect size of coefficient vector X")
    Y = get_vector(M.manifold, p, X, B)
    size(Y) == representation_size(M) || error("Incorrect size of tangent vector Y")
    return Y
end

function get_vector!(M::ValidationManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    size(X) == (manifold_dimension(M),) || error("Incorrect size of coefficient vector X")
    get_vector!(M.manifold, Y, p, X, B)
    size(Y) == representation_size(M) || error("Incorrect size of tangent vector Y")
    return Y
end

injectivity_radius(M::ValidationManifold) = injectivity_radius(M.manifold)
function injectivity_radius(M::ValidationManifold, method::AbstractRetractionMethod)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(M::ValidationManifold, p; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    return injectivity_radius(M.manifold, array_value(p))
end
function injectivity_radius(
    M::ValidationManifold,
    p,
    method::AbstractRetractionMethod;
    kwargs...,
)
    is_point(M, p; error = M.mode, kwargs...)
    return injectivity_radius(M.manifold, array_value(p), method)
end

function inner(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    is_vector(M, p, Y; error = M.mode, kwargs...)
    return inner(M.manifold, array_value(p), array_value(X), array_value(Y))
end

function is_point(M::ValidationManifold, p; kw...)
    return is_point(M.manifold, array_value(p); kw...)
end
function is_vector(M::ValidationManifold, p, X, cbp::Bool = true; kw...)
    return is_vector(M.manifold, array_value(p), array_value(X), cbp; kw...)
end

function isapprox(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_point(M, q; error = M.mode, kwargs...)
    return isapprox(M.manifold, array_value(p), array_value(q); kwargs...)
end
function isapprox(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    is_vector(M, p, Y; error = M.mode, kwargs...)
    return isapprox(M.manifold, array_value(p), array_value(X), array_value(Y); kwargs...)
end

function log(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_point(M, q; error = M.mode, kwargs...)
    X = log(M.manifold, array_value(p), array_value(q))
    is_vector(M, p, X; error = M.mode, kwargs...)
    return ValidationTVector(X)
end

function log!(M::ValidationManifold, X, p, q; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    is_point(M, q; error = M.mode, kwargs...)
    log!(M.manifold, array_value(X), array_value(p), array_value(q))
    is_vector(M, p, X; error = M.mode, kwargs...)
    return X
end

function mid_point(M::ValidationManifold, p1, p2; kwargs...)
    is_point(M, p1; error = M.mode, kwargs...)
    is_point(M, p2; error = M.mode, kwargs...)
    q = mid_point(M.manifold, array_value(p1), array_value(p2))
    is_point(M, q; error = M.mode, kwargs...)
    return q
end

function mid_point!(M::ValidationManifold, q, p1, p2; kwargs...)
    is_point(M, p1; error = M.mode, kwargs...)
    is_point(M, p2; error = M.mode, kwargs...)
    mid_point!(M.manifold, array_value(q), array_value(p1), array_value(p2))
    is_point(M, q; error = M.mode, kwargs...)
    return q
end

number_eltype(::Type{ValidationMPoint{V}}) where {V} = number_eltype(V)
number_eltype(::Type{ValidationFibreVector{TType,V}}) where {TType,V} = number_eltype(V)

function project!(M::ValidationManifold, Y, p, X; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    project!(M.manifold, array_value(Y), array_value(p), array_value(X))
    is_vector(M, p, Y; error = M.mode, kwargs...)
    return Y
end

function vector_transport_along(
    M::ValidationManifold,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_vector(M, p, X; error = M.mode, kwargs...)
    Y = vector_transport_along(M.manifold, array_value(p), array_value(X), c, m)
    is_vector(M, c[end], Y; error = M.mode, kwargs...)
    return Y
end

function vector_transport_along!(
    M::ValidationManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_vector(M, p, X; error = M.mode, kwargs...)
    vector_transport_along!(
        M.manifold,
        array_value(Y),
        array_value(p),
        array_value(X),
        c,
        m,
    )
    is_vector(M, c[end], Y; error = M.mode, kwargs...)
    return Y
end

function vector_transport_to(
    M::ValidationManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_point(M, q; error = M.mode, kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    Y = vector_transport_to(M.manifold, array_value(p), array_value(X), array_value(q), m)
    is_vector(M, q, Y; error = M.mode, kwargs...)
    return Y
end
function vector_transport_to!(
    M::ValidationManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod;
    kwargs...,
)
    is_point(M, q; error = M.mode, kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    vector_transport_to!(
        M.manifold,
        array_value(Y),
        array_value(p),
        array_value(X),
        array_value(q),
        m,
    )
    is_vector(M, q, Y; error = M.mode, kwargs...)
    return Y
end

function zero_vector(M::ValidationManifold, p; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    w = zero_vector(M.manifold, array_value(p))
    is_vector(M, p, w; error = M.mode, kwargs...)
    return w
end

function zero_vector!(M::ValidationManifold, X, p; kwargs...)
    is_point(M, p; error = M.mode, kwargs...)
    zero_vector!(M.manifold, array_value(X), array_value(p); kwargs...)
    is_vector(M, p, X; error = M.mode, kwargs...)
    return X
end
