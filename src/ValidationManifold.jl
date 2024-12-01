"""
    ValidationManifold{𝔽,M<:AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}

A manifold to add tests to input and output values of functions defined in the interface.

Additionally the points and tangent vectors can also be encapsulated, cf.
[`ValidationMPoint`](@ref), [`ValidationTVector`](@ref), and [`ValidationCoTVector`](@ref).
These types can be used to see where some data is assumed to be from, when working on
manifolds where both points and tangent vectors are represented as (plain) arrays.

Using the `ignore_contexts` keyword allows to specify a single `Symbol` or a vector of `Symbols`
Of which contexts to ignore.

Current contexts are
* `:All`: disable all checks
* `:Point`: checks for points
* `:Vector`: checks for vectors
* `:Output`: checks for output
* `:Input`: checks for input variables

Using the `ignore_functions` keyword (dictionary) allows to disable/ignore certain checks
within single functions for this manifold.
The `key` of the dictionary has to be the `Function` to exclude something in.
The `value` is either a single symbol or a vector of symbols with the same meaning as the
`ignore_contexts` keyword, but limited to this function

# Examples
* `exp => :All` disables _all_ checks in the [`exp`](@ref) function
* `exp => :Point` excludes point checks in the [`exp`](@ref) function
* `exp => [:Point, :Vector]` excludes point and vector checks in the [`exp`](@ref) function

This manifold is a decorator for a manifold, i.e. it decorates a [`AbstractManifold`](@ref) `M`
with types points, vectors, and covectors.

# Fields

* `manifold::M`: The manifold to be decorated
* `mode::Symbol`: The mode to be used for error handling, either `:error` or `:warn`
* `ignore::Dict{<:Union{Symbol, Function},Symbol}`: A dictionary of disabled checks

# Constructor

    ValidationManifold(M::AbstractManifold; kwargs...)

Generate the Validation manifold

## Keyword arguments

* `error::Symbol=:error`: specify how errors in the validation should be reported.
  this is passed to [`is_point`](@ref) and [`is_vector`](@ref) as the `error` keyword argument.
  Available values are `:error`, `:warn`, `:info`, and `:none`. Every other value is treated as `:none`.
* `store_base_point::Bool=false`: specify whether or not to store the point `p` a tangent or cotangent vector
  is associated with. This can be useful for debugging purposes.
* `ignores=Dict{Union{Function,Symbol},Union{Symbol,Bool}()` a dictionary of disabled checks
"""
struct ValidationManifold{
    𝔽,
    M<:AbstractManifold{𝔽},
    D<:Dict{Function,Union{Symbol,Vector{Symbol}}},
    V<:AbstractVector{Symbol},
} <: AbstractDecoratorManifold{𝔽}
    manifold::M
    mode::Symbol
    store_base_point::Bool
    ignore_functions::D
    ignore_contexts::V
end
function ValidationManifold(
    M::AbstractManifold;
    error::Symbol = :error,
    store_base_point::Bool = false,
    ignore_functions::D = Dict{Function,Union{Symbol,Vector{Symbol}}}(),
    ignore_contexts::V = Vector{Symbol}(),
) where {D<:Dict{Function,Union{Symbol,Vector{Symbol}}},V<:AbstractVector{Symbol}}
    return ValidationManifold(M, error, store_base_point, ignore_functions, ignore_contexts)
end

"""
_vMc(M::ValidationManifold, f::Function, context::Symbol)
_vMc(M::ValidationManifold, f::Function, context::NTuple{N,Symbol}) where {N}

Return whether a check should be performed within `f` and the `context`(`s`) provided.

This function returns false and hence indicates not to check, when
* (one of the) `context`(`s`) is in the ignore list for `f` within `ignore_functions`
* (one of the) `context`(`s`) is in the `ignore_contexts` list

Otherwise the test is active.

!!! Note
   This function is internal and used very often, co it has a very short name;
    `_vMc` stands for "`ValidationManifold` check".
"""
function _vMc end

function _vMc(M::ValidationManifold, ::Nothing, context::Symbol)
    # Similarly for the global contexts
    (:All ∈ M.ignore_contexts) && return false
    (context ∈ M.ignore_contexts) && return false
end
function _vMc(M::ValidationManifold, f::Function, context::Symbol)
    if haskey(M.ignore_functions, f)
        # If :All is present -> deactivate
        !_VMc(M.ignore_functions[f], :All) && return false
        # If any of the provided contexts is present -> deactivate
        !_VMc(M.ignore_functions[f], context) && return false
    end
    !_vMc(M, nothing, context) && return false
    return true
end
function _vMc(M::ValidationManifold, f, contexts::NTuple{N,Symbol}) where {N}
    for c in contexts
        !_vMc(M, f, c) && return false
    end
    return true
end
# Sub tests: is any of a in b? Then return false
# If a is in or equal to b
_vMc(a::Symbol, b::Symbol) = !(a === b)
_vMc(a::Symbol, b::NTuple{N,Symbol} where {N}) = !(a ∈ b)
# If a is multiple, then test all of them
_vMc(a::NTuple{N,Symbol} where {N}, b::Symbol) = !(b ∈ a)
function _vMc(a::NTuple{N,Symbol} where {N}, b::NTuple{N,Symbol} where {N})
    for ai in a
        (ai ∈ b) && return false
    end
    return true
end

"""
    ValidationMPoint{P} <: AbstractManifoldPoint

Represent a point on an [`ValidationManifold`](@ref). The point is stored internally.

# Fields
* ` value::P`: the internally stored point on a manifold

# Constructor

        ValidationMPoint(value)

Create a point on the manifold with the value `value`.
"""
struct ValidationMPoint{P} <: AbstractManifoldPoint
    value::P
end

"""
    ValidationFibreVector{TType<:VectorSpaceType,V,P} <: AbstractFibreVector{TType}

Represent a tangent vector to a point on an [`ValidationManifold`](@ref).
The original vector of the manifold is stored internally. The corresponding base point
of the fibre can be stored as well.

The `TType` indicates the type of fibre, for example [`TangentSpaceType`](@ref) or [`CotangentSpaceType`](@ref).

# Fields

* `value::V`: the internally stored vector on the fibre
* `point::P`: the point the vector is associated with

# Constructor

        ValidationFibreVector{TType}(value, point=nothing)

"""
struct ValidationFibreVector{TType<:VectorSpaceType,V,P} <: AbstractFibreVector{TType}
    value::V
    point::P
end
function ValidationFibreVector{TType}(value::V, point::P = nothing) where {TType,V,P}
    return ValidationFibreVector{TType,V,P}(value, point)
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
    _value(p)

Return the internal value of an [`ValidationMPoint`](@ref), [`ValidationTVector`](@ref), or
[`ValidationCoTVector`](@ref) if the value `p` is encapsulated as such.
Return `p` if it is already an a (plain) value on a manifold.
"""
_value(p::AbstractArray) = p
_value(p::ValidationMPoint) = p.value
_value(X::ValidationFibreVector) = X.value

"""
    _msg(str,mode)

issue a message `str` according to the mode `mode` (as `@error`, `@warn`, `@info`).
"""
function _msg(str, mode)
    (mode === :error) && (@error str)
    (mode === :warn) && (@warn str)
    (mode === :info) && (@info str)
    return nothing
end

convert(::Type{M}, m::ValidationManifold{𝔽,M}) where {𝔽,M<:AbstractManifold{𝔽}} = m.manifold
function convert(::Type{ValidationManifold{𝔽,M}}, m::M) where {𝔽,M<:AbstractManifold{𝔽}}
    return ValidationManifold(m)
end
function convert(
    ::Type{V},
    p::ValidationMPoint{V},
) where {V<:Union{AbstractArray,AbstractManifoldPoint}}
    return p.value
end
function convert(::Type{ValidationMPoint{V}}, x::V) where {V<:AbstractArray}
    return ValidationMPoint{V}(x)
end

function convert(
    ::Type{V},
    X::ValidationFibreVector{TType,V,Nothing},
) where {TType,V<:Union{AbstractArray,AbstractFibreVector}}
    return X.value
end
function convert(::Type{ValidationFibreVector{TType,V,Nothing}}, X::V) where {TType,V}
    return ValidationFibreVector{TType}(X)
end

function copyto!(M::ValidationManifold, q::ValidationMPoint, p::ValidationMPoint; kwargs...)
    is_point(M, p; error = M.mode, within = copyto, context = (:Input,), kwargs...)
    copyto!(M.manifold, q.value, p.value)
    is_point(M, q; error = M.mode, within = copyto, context = (:Input,), kwargs...)
    return q
end
function copyto!(
    M::ValidationManifold,
    Y::ValidationFibreVector{TType},
    p::ValidationMPoint,
    X::ValidationFibreVector{TType};
    kwargs...,
) where {TType}
    is_point(M, p; error = M.mode, within = copyto, context = (:Input,), kwargs...)
    copyto!(M.manifold, Y.value, p.value, X.value)
    return p
end

decorated_manifold(M::ValidationManifold) = M.manifold

function distance(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; error = M.mode, within = distance, context = (:Input,), kwargs...)
    is_point(M, q; error = M.mode, within = distance, context = (:Input,), kwargs...)
    d = distance(M.manifold, _value(p), _value(q))
    (d < 0) && _vMc(M, distance, :Output) && _msg("Distance is negative: $d", M.mode)
    return d
end

function exp(M::ValidationManifold, p, X; kwargs...)
    is_point(M, p; error = M.mode, within = exp, context = (:Input,), kwargs...)
    is_vector(M, p, X; error = M.mode, within = exp, context = (:Input,), kwargs...)
    y = exp(M.manifold, _value(p), _value(X))
    is_point(M, y; error = M.mode, within = exp, context = (:Output,), kwargs...)
    return ValidationMPoint(y)
end

function exp!(M::ValidationManifold, q, p, X; kwargs...)
    is_point(M, p; error = M.mode, within = exp, context = (:Input,), kwargs...)
    is_vector(M, p, X; error = M.mode, within = exp, context = (:Input,), kwargs...)
    exp!(M.manifold, _value(q), _value(p), _value(X))
    is_point(M, q; error = M.mode, within = exp, context = (:Output,), kwargs...)
    return q
end

function get_basis(M::ValidationManifold, p, B::AbstractBasis; kwargs...)
    is_point(M, p; error = M.mode, within = get_basis, context = (:Input,), kwargs...)
    Ξ = get_basis(M.manifold, _value(p), B)
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
    map(
        X -> is_vector(
            M,
            p,
            X;
            error = M.mode,
            within = get_basis,
            context = (:Output,),
            kwargs...,
        ),
        bvectors,
    )
    return Ξ
end
function get_basis(
    M::ValidationManifold,
    p,
    B::Union{AbstractOrthogonalBasis,CachedBasis{𝔽,<:AbstractOrthogonalBasis{𝔽}} where {𝔽}};
    kwargs...,
)
    is_point(M, p; error = M.mode, within = get_basis, context = (:Input,), kwargs...)
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
    B::Union{
        AbstractOrthonormalBasis,
        <:CachedBasis{𝔽,<:AbstractOrthonormalBasis{𝔽}} where {𝔽},
    };
    kwargs...,
)
    is_point(M, p; error = M.mode, within = get_basis, context = (:Input,), kwargs...)
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
    is_point(M, p; error = :error, within = get_coordinates, context = (:Input,), kwargs...)
    is_vector(
        M,
        p,
        X;
        error = :error,
        within = get_coordinates,
        context = (:Input,),
        kwargs...,
    )
    return get_coordinates(M.manifold, p, X, B)
end

function get_coordinates!(M::ValidationManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_vector(
        M,
        p,
        X;
        error = M.mode,
        within = get_coordinates,
        context = (:Input,),
        kwargs...,
    )
    get_coordinates!(M.manifold, Y, p, X, B)
    return Y
end

function get_vector(M::ValidationManifold, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p; error = M.mode, within = get_vector, context = (:Input,), kwargs...)
    size(X) == (manifold_dimension(M),) || error("Incorrect size of coefficient vector X")
    Y = get_vector(M.manifold, p, X, B)
    size(Y) == representation_size(M) || error("Incorrect size of tangent vector Y")
    return Y
end

function get_vector!(M::ValidationManifold, Y, p, X, B::AbstractBasis; kwargs...)
    is_point(M, p; error = M.mode, within = get_vector, context = (:Input,), kwargs...)
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
    is_point(
        M,
        p;
        error = M.mode,
        within = injectivity_radius,
        context = (:Input,),
        kwargs...,
    )
    return injectivity_radius(M.manifold, _value(p))
end
function injectivity_radius(
    M::ValidationManifold,
    p,
    method::AbstractRetractionMethod;
    kwargs...,
)
    is_point(
        M,
        p;
        error = M.mode,
        within = injectivity_radius,
        context = (:Input,),
        kwargs...,
    )
    return injectivity_radius(M.manifold, _value(p), method)
end

function inner(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p; error = M.mode, within = inner, context = (:Input,), kwargs...)
    is_vector(M, p, X; error = M.mode, within = inner, context = (:Input,), kwargs...)
    is_vector(M, p, Y; error = M.mode, within = inner, context = (:Input,), kwargs...)
    return inner(M.manifold, _value(p), _value(X), _value(Y))
end

"""
    is_point(M::ValidationManifold, p; kwargs...)

perform [`is_point`](@ref) on a [`ValidationManifold`](@ref),
where two additional keywords can be used

* `within=nothing` to specify a function from within which this call was issued
* `context::NTuple{N,Symbol}=NTuple{0,Symbol}()` to specify one or more contexts, this
  call was issued in. The context `:Point` is added before checking whether the test
  should be performed

all other keywords are passed on.
"""
function is_point(
    M::ValidationManifold,
    p;
    within::Union{Nothing,<:Function} = nothing,
    context::Union{NTuple{N,Symbol} where N} = NTuple{0,Symbol}(),
    kwargs...,
)
    _vMc(M, within, (:Point, context...)) && return true
    return is_point(M.manifold, _value(p); kwargs...)
end

"""
    is_vector(M::ValidationManifold, p, X, cbp=true; kwargs...)

perform [`is_vector`](@ref) on a [`ValidationManifold`](@ref),
where two additional keywords can be used

* `within=nothing` to specify a function from within which this call was issued
* `context::NTuple{N,Symbol}=NTuple{0,Symbol}()` to specify one or more contexts, this
  call was issued in. The context `:Point` is added before checking whether the test
  should be performed

all other keywords are passed on.
"""
function is_vector(
    M::ValidationManifold,
    p,
    X,
    cbp::Bool = true;
    within::Union{Nothing,<:Function} = nothing,
    context::Union{NTuple{N,Symbol} where N} = NTuple{0,Symbol}(),
    kwargs...,
)
    _vMc(M, within, (:Point, context...)) && return true
    return is_vector(M.manifold, _value(p), _value(X), cbp; kwargs...)
end

function isapprox(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; error = M.mode, within = isapprox, context = (:Input,), kwargs...)
    is_point(M, q; error = M.mode, within = isapprox, context = (:Input,), kwargs...)
    return isapprox(M.manifold, _value(p), _value(q); kwargs...)
end
function isapprox(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p; error = M.mode, within = isapprox, context = (:Input,), kwargs...)
    is_vector(M, p, X; error = M.mode, within = isapprox, context = (:Input,), kwargs...)
    is_vector(M, p, Y; error = M.mode, within = isapprox, context = (:Input,), kwargs...)
    return isapprox(M.manifold, _value(p), _value(X), _value(Y); kwargs...)
end

function log(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; error = M.mode, within = log, context = (:Input,), kwargs...)
    is_point(M, q; error = M.mode, within = log, context = (:Input,), kwargs...)
    X = log(M.manifold, _value(p), _value(q))
    is_vector(M, p, X; error = M.mode, within = log, context = (:Output,), kwargs...)
    return ValidationTVector(X)
end

function log!(M::ValidationManifold, X, p, q; kwargs...)
    is_point(M, p; error = M.mode, within = log, context = (:Input,), kwargs...)
    is_point(M, q; error = M.mode, within = log, context = (:Input,), kwargs...)
    log!(M.manifold, _value(X), _value(p), _value(q))
    is_vector(M, p, X; error = M.mode, within = log, context = (:Output,), kwargs...)
    return X
end

function mid_point(M::ValidationManifold, p1, p2; kwargs...)
    is_point(M, p1; error = M.mode, within = mid_point, context = (:Input,), kwargs...)
    is_point(M, p2; error = M.mode, within = mid_point, context = (:Input,), kwargs...)
    q = mid_point(M.manifold, _value(p1), _value(p2))
    is_point(M, q; error = M.mode, within = mid_point, context = (:Output,), kwargs...)
    return q
end

function mid_point!(M::ValidationManifold, q, p1, p2; kwargs...)
    is_point(M, p1; error = M.mode, within = mid_point, context = (:Input,), kwargs...)
    is_point(M, p2; error = M.mode, within = mid_point, context = (:Input,), kwargs...)
    mid_point!(M.manifold, _value(q), _value(p1), _value(p2))
    is_point(M, q; error = M.mode, within = mid_point, context = (:Output,), kwargs...)
    return q
end

number_eltype(::Type{ValidationMPoint{V}}) where {V} = number_eltype(V)
number_eltype(::Type{ValidationFibreVector{TType,V,P}}) where {TType,V,P} = number_eltype(V)

function project!(M::ValidationManifold, Y, p, X; kwargs...)
    is_point(M, p; error = M.mode, within = project, context = (:Input,), kwargs...)
    project!(M.manifold, _value(Y), _value(p), _value(X))
    is_vector(M, p, Y; error = M.mode, within = project, context = (:Output,), kwargs...)
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
    is_vector(
        M,
        p,
        X;
        error = M.mode,
        within = vector_transport_along,
        context = (:Input,),
        kwargs...,
    )
    Y = vector_transport_along(M.manifold, _value(p), _value(X), c, m)
    is_vector(
        M,
        c[end],
        Y;
        error = M.mode,
        within = vector_transport_along,
        context = (:Output,),
        kwargs...,
    )
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
    is_vector(
        M,
        p,
        X;
        error = M.mode,
        within = vector_transport_along,
        context = (:Input,),
        kwargs...,
    )
    vector_transport_along!(M.manifold, _value(Y), _value(p), _value(X), c, m)
    is_vector(
        M,
        c[end],
        Y;
        error = M.mode,
        within = vector_transport_along,
        context = (:Output,),
        kwargs...,
    )
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
    is_point(
        M,
        q;
        error = M.mode,
        within = vector_transport_to,
        context = (:Input,),
        kwargs...,
    )
    is_vector(
        M,
        p,
        X;
        error = M.mode,
        within = vector_transport_to,
        context = (:Input,),
        kwargs...,
    )
    Y = vector_transport_to(M.manifold, _value(p), _value(X), _value(q), m)
    is_vector(
        M,
        q,
        Y;
        error = M.mode,
        within = vector_transport_to,
        context = (:Output,),
        kwargs...,
    )
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
    is_point(
        M,
        q;
        error = M.mode,
        within = vector_transport_to,
        context = (:Input,),
        kwargs...,
    )
    is_vector(
        M,
        p,
        X;
        error = M.mode,
        within = vector_transport_to,
        context = (:Input,),
        kwargs...,
    )
    vector_transport_to!(M.manifold, _value(Y), _value(p), _value(X), _value(q), m)
    is_vector(
        M,
        q,
        Y;
        error = M.mode,
        within = vector_transport_to,
        context = (:Output,),
        kwargs...,
    )
    return Y
end

function zero_vector(M::ValidationManifold, p; kwargs...)
    is_point(M, p; error = M.mode, within = zero_vector, context = (:Input,), kwargs...)
    w = zero_vector(M.manifold, _value(p))
    is_vector(
        M,
        p,
        w;
        error = M.mode,
        within = zero_vector,
        context = (:Output,),
        kwargs...,
    )
    return w
end

function zero_vector!(M::ValidationManifold, X, p; kwargs...)
    is_point(M, p; error = M.mode, within = zero_vector, context = (:Input,), kwargs...)
    zero_vector!(M.manifold, _value(X), _value(p); kwargs...)
    is_vector(
        M,
        p,
        X;
        error = M.mode,
        within = zero_vector,
        context = (:Output,),
        kwargs...,
    )
    return X
end
