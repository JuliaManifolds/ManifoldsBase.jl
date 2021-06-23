"""
    ValidationManifold{𝔽,M<:AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}

A manifold to encapsulate manifolds working on array representations of [`AbstractManifoldPoint`](@ref)s
and [`TVector`](@ref)s in a transparent way, such that for these manifolds it's not
necessary to introduce explicit types for the points and tangent vectors, but they are
encapsulated/stripped automatically when needed.

This manifold is a decorator for a manifold, i.e. it decorates a [`AbstractManifold`](@ref) `M`
with types points, vectors, and covectors.
"""
struct ValidationManifold{𝔽,M<:AbstractManifold{𝔽}} <:
       AbstractDecoratorManifold{𝔽,DefaultDecoratorType}
    manifold::M
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


function (+)(X::ValidationFibreVector{TType}, Y::ValidationFibreVector{TType}) where {TType}
    return ValidationFibreVector{TType}(X.value + Y.value)
end
function (-)(X::ValidationFibreVector{TType}, Y::ValidationFibreVector{TType}) where {TType}
    return ValidationFibreVector{TType}(X.value - Y.value)
end
(-)(X::ValidationFibreVector{TType}) where {TType} = ValidationFibreVector{TType}(-X.value)
function (*)(a::Number, X::ValidationFibreVector{TType}) where {TType}
    return ValidationFibreVector{TType}(a * X.value)
end

allocate(p::ValidationMPoint) = ValidationMPoint(allocate(p.value))
allocate(p::ValidationMPoint, ::Type{T}) where {T} = ValidationMPoint(allocate(p.value, T))
function allocate(X::ValidationFibreVector{TType}) where {TType}
    return ValidationFibreVector{TType}(allocate(X.value))
end
function allocate(X::ValidationFibreVector{TType}, ::Type{T}) where {TType,T}
    return ValidationFibreVector{TType}(allocate(X.value, T))
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

function check_point(M::ValidationManifold, p; kwargs...)
    return check_point(M.manifold, array_value(p); kwargs...)
end
function check_point(M::ValidationManifold, p::AbstractManifoldPoint; kwargs...)
    return check_point(M.manifold, array_value(p); kwargs...)
end

function check_vector(M::ValidationManifold, p, X; kwargs...)
    return check_vector(M.manifold, array_value(p), array_value(X); kwargs...)
end
function check_vector(
    M::ValidationManifold,
    p::AbstractManifoldPoint,
    X::TVector;
    kwargs...,
)
    return check_vector(M.manifold, array_value(p), array_value(X); kwargs...)
end

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
    is_point(M, p, true; kwargs...)
    copyto!(M.manifold, q.value, p.value)
    is_point(M, q, true; kwargs...)
    return q
end
function copyto!(
    M::ValidationManifold,
    Y::ValidationFibreVector{TType},
    p::ValidationMPoint,
    X::ValidationFibreVector{TType};
    kwargs...,
) where {TType}
    is_point(M, p, true; kwargs...)
    copyto!(M.manifold, Y.value, p.value, X.value)
    return p
end

function distance(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p, true; kwargs...)
    is_point(M, q, true; kwargs...)
    return distance(M.manifold, array_value(p), array_value(q))
end

function exp(M::ValidationManifold, p, X; kwargs...)
    is_point(M, p, true; kwargs...)
    is_vector(M, p, X, true; kwargs...)
    y = exp(M.manifold, array_value(p), array_value(X))
    is_point(M, y, true; kwargs...)
    return ValidationMPoint(y)
end

function exp!(M::ValidationManifold, q, p, X; kwargs...)
    is_point(M, p, true; kwargs...)
    is_vector(M, p, X, true; kwargs...)
    exp!(M.manifold, array_value(q), array_value(p), array_value(X))
    is_point(M, q, true; kwargs...)
    return q
end

function get_basis(M::ValidationManifold, p, B::AbstractBasis; kwargs...)
    is_point(M, p, true; kwargs...)
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
    map(X -> is_vector(M, p, X, true; kwargs...), bvectors)
    return Ξ
end
function get_basis(
    M::ValidationManifold,
    p,
    B::Union{AbstractOrthogonalBasis,CachedBasis{𝔽,<:AbstractOrthogonalBasis{𝔽}}};
    kwargs...,
) where {𝔽}
    is_point(M, p, true; kwargs...)
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
    is_point(M, p, true; kwargs...)
    get_basis_invoke_types = Tuple{
        ValidationManifold,
        Any,
        Union{
            AbstractOrthogonalBasis,
            CachedBasis{𝔽,<:AbstractOrthogonalBasis{𝔽}},
        } where {𝔽},
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
for BT in DISAMBIGUATION_BASIS_TYPES
    if BT <:
       Union{AbstractOrthonormalBasis,CachedBasis{𝔽,<:AbstractOrthonormalBasis} where 𝔽}
        CT = AbstractOrthonormalBasis
    elseif BT <:
           Union{AbstractOrthogonalBasis,CachedBasis{𝔽,<:AbstractOrthogonalBasis} where 𝔽}
        CT = AbstractOrthogonalBasis
    else
        CT = AbstractBasis
    end
    eval(quote
        @invoke_maker 3 $CT get_basis(M::ValidationManifold, p, B::$BT; kwargs...)
    end)
end

function get_coordinates(M::ValidationManifold, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p, true; kwargs...)
    is_vector(M, p, X, true; kwargs...)
    return get_coordinates(M.manifold, p, X, B)
end
for BT in DISAMBIGUATION_BASIS_TYPES
    eval(
        quote
            @invoke_maker 4 AbstractBasis get_coordinates(
                M::ValidationManifold,
                p,
                X,
                B::$BT;
                kwargs...,
            )
        end,
    )
end

function get_coordinates!(M::ValidationManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_vector(M, p, X, true; kwargs...)
    get_coordinates!(M.manifold, Y, p, X, B)
    return Y
end
for BT in [DISAMBIGUATION_BASIS_TYPES..., DISAMBIGUATION_COTANGENT_BASIS_TYPES...]
    eval(
        quote
            @invoke_maker 5 AbstractBasis get_coordinates!(
                M::ValidationManifold,
                Y,
                p,
                X,
                B::$BT;
                kwargs...,
            )
        end,
    )
end

function get_vector(M::ValidationManifold, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p, true; kwargs...)
    size(X) == (manifold_dimension(M),) || error("Incorrect size of coefficient vector X")
    Y = get_vector(M.manifold, p, X, B)
    size(Y) == representation_size(M) || error("Incorrect size of tangent vector Y")
    return Y
end
for BT in DISAMBIGUATION_BASIS_TYPES
    eval(
        quote
            @invoke_maker 4 AbstractBasis get_vector(
                M::ValidationManifold,
                p,
                X,
                B::$BT;
                kwargs...,
            )
        end,
    )
end

function get_vector!(M::ValidationManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p, true; kwargs...)
    size(X) == (manifold_dimension(M),) || error("Incorrect size of coefficient vector X")
    get_vector!(M.manifold, Y, p, X, B)
    size(Y) == representation_size(M) || error("Incorrect size of tangent vector Y")
    return Y
end
for BT in [DISAMBIGUATION_BASIS_TYPES..., DISAMBIGUATION_COTANGENT_BASIS_TYPES...]
    eval(
        quote
            @invoke_maker 5 AbstractBasis get_vector!(
                M::ValidationManifold,
                Y,
                p,
                X,
                B::$BT;
                kwargs...,
            )
        end,
    )
end

injectivity_radius(M::ValidationManifold) = injectivity_radius(M.manifold)
function injectivity_radius(M::ValidationManifold, method::AbstractRetractionMethod)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(M::ValidationManifold, p; kwargs...)
    is_point(M, p, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(p))
end
function injectivity_radius(
    M::ValidationManifold,
    p,
    method::AbstractRetractionMethod;
    kwargs...,
)
    is_point(M, p, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(p), method)
end
function injectivity_radius(M::ValidationManifold, method::ExponentialRetraction)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(
    M::ValidationManifold,
    p,
    method::ExponentialRetraction;
    kwargs...,
)
    is_point(M, p, true; kwargs...)
    return injectivity_radius(M.manifold, array_value(p), method)
end

function inner(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p, true; kwargs...)
    is_vector(M, p, X, true; kwargs...)
    is_vector(M, p, Y, true; kwargs...)
    return inner(M.manifold, array_value(p), array_value(X), array_value(Y))
end

function isapprox(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p, true; kwargs...)
    is_point(M, q, true; kwargs...)
    return isapprox(M.manifold, array_value(p), array_value(q); kwargs...)
end
function isapprox(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p, true; kwargs...)
    is_vector(M, p, X, true; kwargs...)
    is_vector(M, p, Y, true; kwargs...)
    return isapprox(M.manifold, array_value(p), array_value(X), array_value(Y); kwargs...)
end

function log(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p, true; kwargs...)
    is_point(M, q, true; kwargs...)
    X = log(M.manifold, array_value(p), array_value(q))
    is_vector(M, p, X, true; kwargs...)
    return ValidationTVector(X)
end

function log!(M::ValidationManifold, X, p, q; kwargs...)
    is_point(M, p, true; kwargs...)
    is_point(M, q, true; kwargs...)
    log!(M.manifold, array_value(X), array_value(p), array_value(q))
    is_vector(M, p, X, true; kwargs...)
    return X
end

function mid_point(M::ValidationManifold, p1, p2; kwargs...)
    is_point(M, p1, true; kwargs...)
    is_point(M, p2, true; kwargs...)
    q = mid_point(M.manifold, array_value(p1), array_value(p2))
    is_point(M, q, true; kwargs...)
    return q
end

function mid_point!(M::ValidationManifold, q, p1, p2; kwargs...)
    is_point(M, p1, true; kwargs...)
    is_point(M, p2, true; kwargs...)
    mid_point!(M.manifold, array_value(q), array_value(p1), array_value(p2))
    is_point(M, q, true; kwargs...)
    return q
end

number_eltype(::Type{ValidationMPoint{V}}) where {V} = number_eltype(V)
number_eltype(p::ValidationMPoint) = number_eltype(p.value)
number_eltype(::Type{ValidationFibreVector{TType,V}}) where {TType,V} = number_eltype(V)
number_eltype(X::ValidationFibreVector) = number_eltype(X.value)

function project!(M::ValidationManifold, Y, p, X; kwargs...)
    is_point(M, p, true; kwargs...)
    project!(M.manifold, array_value(Y), array_value(p), array_value(X))
    is_vector(M, p, Y, true; kwargs...)
    return Y
end

similar(p::ValidationMPoint) = ValidationMPoint(similar(p.value))
similar(p::ValidationMPoint, ::Type{T}) where {T} = ValidationMPoint(similar(p.value, T))
function similar(X::ValidationFibreVector{TType}) where {TType}
    return ValidationFibreVector{TType}(similar(X.value))
end
function similar(X::ValidationFibreVector{TType}, ::Type{T}) where {TType,T}
    return ValidationFibreVector{TType}(similar(X.value, T))
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
    is_vector(M, p, X, true; kwargs...)
    vector_transport_along!(
        M.manifold,
        array_value(Y),
        array_value(p),
        array_value(X),
        c,
        m,
    )
    is_vector(M, c[end], Y, true; kwargs...)
    return Y
end
for VT in VECTOR_TRANSPORT_DISAMBIGUATION
    eval(
        quote
            @invoke_maker 6 AbstractVectorTransportMethod vector_transport_along!(
                M::ValidationManifold,
                vto,
                x,
                v,
                c::AbstractVector,
                B::$VT,
            )
        end,
    )
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
    is_point(M, q, true; kwargs...)
    is_vector(M, p, X, true; kwargs...)
    vector_transport_to!(
        M.manifold,
        array_value(Y),
        array_value(p),
        array_value(X),
        array_value(q),
        m,
    )
    is_vector(M, q, Y, true; kwargs...)
    return Y
end

for T in [
    PoleLadderTransport,
    ProjectionTransport,
    ScaledVectorTransport,
    SchildsLadderTransport,
]
    @eval begin
        function vector_transport_to!(M::ValidationManifold, Y, p, X, q, m::$T; kwargs...)
            is_point(M, q, true; kwargs...)
            is_vector(M, p, X, true; kwargs...)
            vector_transport_to!(
                M.manifold,
                array_value(Y),
                array_value(p),
                array_value(X),
                array_value(q),
                m,
            )
            is_vector(M, q, Y, true; kwargs...)
            return Y
        end
    end
end

function zero_vector(M::ValidationManifold, p; kwargs...)
    is_point(M, p, true; kwargs...)
    w = zero_vector(M.manifold, array_value(p))
    is_vector(M, p, w, true; kwargs...)
    return w
end

function zero_vector!(M::ValidationManifold, X, p; kwargs...)
    is_point(M, p, true; kwargs...)
    zero_vector!(M.manifold, array_value(X), array_value(p); kwargs...)
    is_vector(M, p, X, true; kwargs...)
    return X
end
