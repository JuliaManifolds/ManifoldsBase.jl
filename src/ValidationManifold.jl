"""
    ValidationManifold{𝔽,M<:AbstractManifold{𝔽}} <: AbstractDecoratorManifold{𝔽}

A manifold to add tests to input and output values of functions defined in the interface.

Additionally the points and tangent vectors can also be encapsulated, cf.
[`ValidationMPoint`](@ref), [`ValidationTangentVector`](@ref), and [`ValidationCotangentVector`](@ref).
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
* `ignore_contexts::AbstractVector{Symbol}`: store contexts to be ignored of validation.
* `ignore_functions::Dict{<:Function,<:Union{Symbol,<:AbstractVector{Symbol}}`:
  store contexts to be ignored with in a function or its mutating variant.

# Constructors

    ValidationManifold(M::AbstractManifold; kwargs...)

Generate the Validation manifold

    ValidationManifold(M::AbstractManifold, V::ValidationManifold; kwargs...)

Generate the Validation manifold for `M` with the default values of `V`.


## Keyword arguments

* `error::Symbol=:error`: specify how errors in the validation should be reported.
  this is passed to [`is_point`](@ref) and [`is_vector`](@ref) as the `error` keyword argument.
  Available values are `:error`, `:warn`, `:info`, and `:none`. Every other value is treated as `:none`.
* `store_base_point::Bool=false`: specify whether or not to store the point `p` a tangent or cotangent vector
  is associated with. This can be useful for debugging purposes.
* `ignore_contexts = Vector{Symbol}()` a vector to indicate which validation contexts should not be performed.
* `ignore_functions=Dict{Function,Union{Symbol,Vector{Symbol}}}()` a dictionary to disable certain contexts within functions.
  The key here is the non-mutating function variant (if it exists). The contexts are thre same as in `ignore_contexts`.
"""
struct ValidationManifold{
        𝔽,
        M <: AbstractManifold{𝔽},
        D <: Dict{<:Function, <:Union{Symbol, <:AbstractVector{Symbol}}},
        V <: AbstractVector{Symbol},
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
        ignore_functions::D = Dict{Function, Union{Symbol, <:Vector{Symbol}}}(),
        ignore_contexts::V = Vector{Symbol}(),
    ) where {
        D <: Dict{<:Function, <:Union{Symbol, <:AbstractVector{Symbol}}},
        V <: AbstractVector{Symbol},
    }
    return ValidationManifold(M, error, store_base_point, ignore_functions, ignore_contexts)
end
function ValidationManifold(M::AbstractManifold, V::ValidationManifold; kwargs...)
    return ValidationManifold(
        M;
        error = V.mode,
        ignore_contexts = V.ignore_contexts,
        ignore_functions = V.ignore_functions,
        store_base_point = V.store_base_point,
        kwargs...,
    )
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
    return true
end
function _vMc(M::ValidationManifold, f::Function, context::Symbol)
    if haskey(M.ignore_functions, f)
        # If :All is present -> deactivate
        !_vMc(M.ignore_functions[f], :All) && return false
        # If any of the provided contexts is present -> deactivate
        !_vMc(M.ignore_functions[f], context) && return false
    end
    !_vMc(M, nothing, context) && return false
    return true
end
function _vMc(M::ValidationManifold, f, contexts::NTuple{N, Symbol}) where {N}
    for c in contexts
        !_vMc(M, f, c) && return false
    end
    return true
end
# Sub tests: is any of a in b? Then return false – b is from before always a symbol already
# If a and b are symbols, equality is checked
_vMc(a::Symbol, b::Symbol) = !(a === b)
# If a is a vector multiple, then return false if b appears in a
_vMc(a::Union{<:NTuple{N, Symbol} where {N}, <:AbstractVector{Symbol}}, b::Symbol) = !(b ∈ a)

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
struct ValidationFibreVector{TType <: VectorSpaceType, V, P} <: AbstractFibreVector{TType}
    value::V
    point::P
end
function ValidationFibreVector{TType}(value::V, point::P = nothing) where {TType, V, P}
    return ValidationFibreVector{TType, V, P}(value, point)
end

"""
    ValidationTangentVector = ValidationFibreVector{TangentSpaceType}

Represent a tangent vector to a point on an [`ValidationManifold`](@ref), i.e. on a manifold
where data can be represented by arrays. The array is stored internally and semantically.
This distinguished the value from [`ValidationMPoint`](@ref)s vectors of other types.
"""
const ValidationTangentVector = ValidationFibreVector{TangentSpaceType}

"""
    ValidationCotangentVector = ValidationFibreVector{CotangentSpaceType}

Represent a cotangent vector to a point on an [`ValidationManifold`](@ref), i.e. on a manifold
where data can be represented by arrays. The array is stored internally and semantically.
This distinguished the value from [`ValidationMPoint`](@ref)s vectors of other types.
"""
const ValidationCotangentVector = ValidationFibreVector{CotangentSpaceType}

@eval @manifold_vector_forwards ValidationFibreVector{TType} TType value

@eval @manifold_element_forwards ValidationMPoint value

@inline function active_traits(f, ::ValidationManifold, ::Any...)
    return merge_traits(IsExplicitDecorator())
end

"""
    internal_value(p)

Return the internal value of an [`ValidationMPoint`](@ref), [`ValidationTangentVector`](@ref), or
[`ValidationCotangentVector`](@ref) if the value `p` is encapsulated as such.
Return `p` if it is already an a (plain) value on a manifold.
"""
internal_value(p) = p
internal_value(p::ValidationMPoint) = p.value
internal_value(X::ValidationFibreVector) = X.value

_update_basepoint!(::ValidationManifold, X, p) = X
function _update_basepoint!(
        ::ValidationManifold,
        X::ValidationTangentVector{P, Nothing},
        p,
    ) where {P}
    return X
end
function _update_basepoint!(
        M::ValidationManifold,
        X::ValidationTangentVector{P, V},
        p,
    ) where {P, V}
    copyto!(M.manifold, X.point, p)
    return X
end


"""
    _msg(str; error=:None, within::Union{Nothing,<:Function} = nothing,
    context::Union{NTuple{N,Symbol} where N} = NTuple{0,Symbol}())

issue a message `str` according to the mode `mode` (as `@error`, `@warn`, `@info`).
"""
function _msg(
        M::ValidationManifold,
        str;
        error = M.mode,
        within::Union{Nothing, <:Function} = nothing,
        context::Union{NTuple{N, Symbol} where {N}} = NTuple{0, Symbol}(),
    )
    !_vMc(M, within, context) && return nothing
    (error === :error) && (throw(ErrorException(str)))
    (error === :warn) && (@warn str)
    (error === :info) && (@info str)
    return nothing
end
function _msg(
        M::ValidationManifold,
        err::Union{DomainError, ArgumentError, ErrorException};
        error = M.mode,
        within::Union{Nothing, <:Function} = nothing,
        context::Union{NTuple{N, Symbol} where {N}} = NTuple{0, Symbol}(),
    )
    !_vMc(M, within, context) && return nothing
    (error === :error) && (throw(err))
    (error === :warn) && (@warn "$err")
    (error === :info) && (@info "$err")
    return nothing
end

convert(::Type{M}, m::ValidationManifold{𝔽, M}) where {𝔽, M <: AbstractManifold{𝔽}} = m.manifold
function convert(::Type{<:ValidationManifold{𝔽, M}}, m::M) where {𝔽, M <: AbstractManifold{𝔽}}
    return ValidationManifold(m)
end
function convert(
        ::Type{V},
        p::ValidationMPoint{V},
    ) where {V <: Union{AbstractArray, AbstractManifoldPoint}}
    return p.value
end
function convert(::Type{ValidationMPoint{V}}, x::V) where {V <: AbstractArray}
    return ValidationMPoint{V}(x)
end

function convert(
        ::Type{V},
        X::ValidationFibreVector{TType, V, Nothing},
    ) where {TType, V <: Union{AbstractArray, AbstractFibreVector}}
    return X.value
end
function convert(::Type{ValidationFibreVector{TType, V, Nothing}}, X::V) where {TType, V}
    return ValidationFibreVector{TType}(X)
end

function copyto!(M::ValidationManifold, q::ValidationMPoint, p::ValidationMPoint; kwargs...)
    is_point(M, p; within = copyto!, context = (:Input,), kwargs...)
    copyto!(M.manifold, q.value, p.value)
    is_point(M, q; within = copyto!, context = (:Input,), kwargs...)
    return q
end
function copyto!(
        M::ValidationManifold,
        Y::ValidationFibreVector{TType},
        p::ValidationMPoint,
        X::ValidationFibreVector{TType};
        kwargs...,
    ) where {TType}
    is_point(M, p; within = copyto!, context = (:Input,), kwargs...)
    copyto!(M.manifold, Y.value, p.value, X.value)
    return p
end

decorated_manifold(M::ValidationManifold) = M.manifold

function distance(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; within = distance, context = (:Input,), kwargs...)
    is_point(M, q; within = distance, context = (:Input,), kwargs...)
    d = distance(M.manifold, internal_value(p), internal_value(q))
    (d < 0) && _msg(
        M,
        DomainError(d, "Distance is negative.");
        within = distance,
        context = (:Output,),
    )
    return d
end

function embed(M::ValidationManifold, p; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    q = embed(M.manifold, internal_value(p))
    MEV = ValidationManifold(get_embedding(M.manifold), M)
    is_point(MEV, q; error = MEV.mode, within = embed, context = (:Output,), kwargs...)
    return ValidationMPoint(q)
end
function embed(M::ValidationManifold, p, X; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = embed, context = (:Input,), kwargs...)
    Y = embed(M.manifold, internal_value(p), internal_value(X))
    MEV = ValidationManifold(get_embedding(M.manifold), M)
    q = embed(M.manifold, internal_value(p))
    is_vector(MEV, q, Y; within = embed, context = (:Output,), kwargs...)
    return ValidationTangentVector(Y, M.store_base_point ? q : nothing)
end

function embed!(M::ValidationManifold, q, p; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    embed!(M.manifold, internal_value(q), internal_value(p))
    MEV = ValidationManifold(get_embedding(M.manifold), M)
    is_point(MEV, q; error = MEV.mode, within = embed, context = (:Output,), kwargs...)
    return q
end
function embed!(M::ValidationManifold, Y, p, X; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = embed, context = (:Input,), kwargs...)
    embed!(M.manifold, internal_value(Y), internal_value(p), internal_value(X))
    MEV = ValidationManifold(get_embedding(M.manifold), M)
    q = embed(M.manifold, internal_value(p))
    _update_basepoint!(M, Y, q)
    is_vector(MEV, q, Y; within = embed, context = (:Output,), kwargs...)
    return Y
end

function embed_project(M::ValidationManifold, p; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    q = embed_project(M.manifold, internal_value(p))
    is_point(M, q; within = embed, context = (:Output,), kwargs...)
    return ValidationMPoint(q)
end
function embed_project(M::ValidationManifold, p, X; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = embed, context = (:Input,), kwargs...)
    Y = embed_project(M.manifold, internal_value(p), internal_value(X))
    is_vector(M, p, Y; within = embed, context = (:Output,), kwargs...)
    return ValidationTangentVector(Y, M.store_base_point ? copy(M, p) : nothing)
end

function embed_project!(M::ValidationManifold, q, p; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    embed_project!(M.manifold, internal_value(q), internal_value(p))
    is_point(M, q; within = embed, context = (:Output,), kwargs...)
    return q
end
function embed_project!(M::ValidationManifold, Y, p, X; kwargs...)
    is_point(M, p; within = embed, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = embed, context = (:Input,), kwargs...)
    embed_project!(M.manifold, internal_value(Y), internal_value(p), internal_value(X))
    _update_basepoint!(M, Y, p)
    is_vector(M, p, Y; within = embed, context = (:Output,), kwargs...)
    return Y
end

function exp(M::ValidationManifold, p, X; kwargs...)
    is_point(M, p; within = exp, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = exp, context = (:Input,), kwargs...)
    y = exp(M.manifold, internal_value(p), internal_value(X))
    is_point(M, y; within = exp, context = (:Output,), kwargs...)
    return ValidationMPoint(y)
end

function exp!(M::ValidationManifold, q, p, X; kwargs...)
    is_point(M, p; within = exp, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = exp, context = (:Input,), kwargs...)
    exp!(M.manifold, internal_value(q), internal_value(p), internal_value(X))
    is_point(M, q; within = exp, context = (:Output,), kwargs...)
    return q
end

function get_basis(M::ValidationManifold, p, B::AbstractBasis; kwargs...)
    is_point(M, p; within = get_basis, context = (:Input,), kwargs...)
    Ξ = get_basis(M.manifold, internal_value(p), B)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    if N != manifold_dimension(M.manifold)
        _msg(
            M,
            "For a basis of the tangent space at $(p) of $(M.manifold), $(manifold_dimension(M)) vectors are required, but get_basis $(B) computed $(N)";
            within = get_basis,
            context = (:Output,),
        )
    end
    # check that the vectors are linearly independent\
    bv_rank = rank(reduce(hcat, bvectors))
    if N != bv_rank
        _msg(
            M,
            "For a basis of the tangent space at $(p) of $(M.manifold), $(manifold_dimension(M)) linearly independent vectors are required, but get_basis $(B) computed $(bv_rank)";
            within = get_basis,
            context = (:Output,),
        )
    end
    map(
        X -> is_vector(M, p, X; within = get_basis, context = (:Output,), kwargs...),
        bvectors,
    )
    return Ξ
end
function get_basis(
        M::ValidationManifold,
        p,
        B::Union{AbstractOrthogonalBasis, CachedBasis{𝔽, <:AbstractOrthogonalBasis{𝔽}} where {𝔽}};
        kwargs...,
    )
    is_point(M, p; within = get_basis, context = (:Input,), kwargs...)
    Ξ = invoke(get_basis, Tuple{ValidationManifold, Any, AbstractBasis}, M, p, B; kwargs...)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    for i in 1:N
        for j in (i + 1):N
            dot_val = real(inner(M, p, bvectors[i], bvectors[j]))
            if !isapprox(dot_val, 0; atol = eps(eltype(p)))
                _msg(
                    M,
                    ArgumentError(
                        "vectors number $i and $j are not orthonormal (inner product = $dot_val)",
                    );
                    within = get_basis,
                    context = (:Output,),
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
            <:CachedBasis{𝔽, <:AbstractOrthonormalBasis{𝔽}} where {𝔽},
        };
        kwargs...,
    )
    is_point(M, p; within = get_basis, context = (:Input,), kwargs...)
    get_basis_invoke_types = Tuple{
        ValidationManifold,
        Any,
        Union{
            AbstractOrthogonalBasis,
            CachedBasis{𝔽2, <:AbstractOrthogonalBasis{𝔽2}},
        } where {𝔽2},
    }
    Ξ = invoke(get_basis, get_basis_invoke_types, M, p, B; kwargs...)
    bvectors = get_vectors(M, p, Ξ)
    N = length(bvectors)
    for i in 1:N
        Xi_norm = norm(M, p, bvectors[i])
        if !isapprox(Xi_norm, 1)
            _msg(
                M,
                ArgumentError("vector number $i is not normalized (norm = $Xi_norm)");
                within = get_basis,
                context = (:Output,),
            )
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

function get_coordinates!(M::ValidationManifold, c, p, X, B::AbstractBasis; kwargs...)
    is_vector(M, p, X; within = get_coordinates, context = (:Input,), kwargs...)
    get_coordinates!(M.manifold, c, internal_value(p), internal_value(X), B)
    return c
end

function get_vector(M::ValidationManifold, p, c, B::AbstractBasis; kwargs...)
    is_point(M, p; within = get_vector, context = (:Input,), kwargs...)
    if size(c) !== (manifold_dimension(M),)
        _msg(
            M,
            ArgumentError(
                "Incorrect size of coefficient vector X ($(size(c))), expected $(manifold_dimension(M)).",
            );
            within = get_vector,
            context = (:Input,),
        )
    end
    Y = get_vector(M.manifold, internal_value(p), internal_value(c), B)
    is_vector(M, p, Y; within = get_vector, context = (:Output,), kwargs...)
    return Y
end

function get_vector!(M::ValidationManifold, Y, p, c, B::AbstractBasis; kwargs...)
    is_point(M, p; within = get_vector, context = (:Input,), kwargs...)
    if size(c) !== (manifold_dimension(M),)
        _msg(
            M,
            ArgumentError(
                "Incorrect size of coefficient vector X ($(size(c))), expected $(manifold_dimension(M)).",
            );
            within = get_vector,
            context = (:Input,),
        )
    end
    get_vector!(M.manifold, internal_value(Y), internal_value(p), internal_value(c), B)
    _update_basepoint!(M, Y, p)
    is_vector(M, p, Y; within = get_vector, context = (:Output,), kwargs...)
    return Y
end

injectivity_radius(M::ValidationManifold) = injectivity_radius(M.manifold)
function injectivity_radius(M::ValidationManifold, method::AbstractRetractionMethod)
    return injectivity_radius(M.manifold, method)
end
function injectivity_radius(M::ValidationManifold, p; kwargs...)
    is_point(M, p; within = injectivity_radius, context = (:Input,), kwargs...)
    return injectivity_radius(M.manifold, internal_value(p))
end
function injectivity_radius(
        M::ValidationManifold,
        p,
        method::AbstractRetractionMethod;
        kwargs...,
    )
    is_point(M, p; within = injectivity_radius, context = (:Input,), kwargs...)
    return injectivity_radius(M.manifold, internal_value(p), method)
end

function inner(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p; within = inner, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = inner, context = (:Input,), kwargs...)
    is_vector(M, p, Y; within = inner, context = (:Input,), kwargs...)
    return inner(M.manifold, internal_value(p), internal_value(X), internal_value(Y))
end

"""
    is_point(M::ValidationManifold, p; kwargs...)

Perform [`is_point`](@ref) on a [`ValidationManifold`](@ref),
where two additional keywords can be used

* `within=nothing` to specify a function from within which this call was issued
* `context::NTuple{N,Symbol} where N=()` to specify one or more contexts, this
  call was issued in. The context `:Point` is added before checking whether the test
  should be performed

all other keywords are passed on.
"""
function is_point(
        M::ValidationManifold,
        p;
        error::Symbol = M.mode,
        within::Union{Nothing, Function} = nothing,
        context::NTuple{N, Symbol} where {N} = (),
        kwargs...,
    )
    !_vMc(M, within, (:Point, context...)) && return true
    return is_point(M.manifold, internal_value(p); error = error, kwargs...)
end

"""
    is_vector(M::ValidationManifold, p, X, cbp=true; kwargs...)

perform [`is_vector`](@ref) on a [`ValidationManifold`](@ref),
where two additional keywords can be used

* `within=nothing` to specify a function from within which this call was issued
* `context::NTuple{N,Symbol} where N=()` to specify one or more contexts, this
  call was issued in. The context `:Point` is added before checking whether the test
  should be performed

all other keywords are passed on.
"""
function is_vector(
        M::ValidationManifold,
        p,
        X,
        cbp::Bool = true;
        error::Symbol = M.mode,
        within::Union{Nothing, Function} = nothing,
        context::NTuple{N, Symbol} where {N} = (),
        kwargs...,
    )
    !_vMc(M, within, (:Vector, context...)) && return true
    return is_vector(
        M.manifold,
        internal_value(p),
        internal_value(X),
        cbp;
        error = error,
        kwargs...,
    )
end

function isapprox(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; within = isapprox, context = (:Input,), kwargs...)
    is_point(M, q; within = isapprox, context = (:Input,), kwargs...)
    return isapprox(M.manifold, internal_value(p), internal_value(q); kwargs...)
end
function isapprox(M::ValidationManifold, p, X, Y; kwargs...)
    is_point(M, p; within = isapprox, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = isapprox, context = (:Input,), kwargs...)
    is_vector(M, p, Y; within = isapprox, context = (:Input,), kwargs...)
    return isapprox(
        M.manifold,
        internal_value(p),
        internal_value(X),
        internal_value(Y);
        kwargs...,
    )
end

function log(M::ValidationManifold, p, q; kwargs...)
    is_point(M, p; within = log, context = (:Input,), kwargs...)
    is_point(M, q; within = log, context = (:Input,), kwargs...)
    X = log(M.manifold, internal_value(p), internal_value(q))
    is_vector(M, p, X; within = log, context = (:Output,), kwargs...)
    return ValidationTangentVector(X, M.store_base_point ? copy(M, p) : nothing)
end

function log!(M::ValidationManifold, X, p, q; kwargs...)
    is_point(M, p; within = log, context = (:Input,), kwargs...)
    is_point(M, q; within = log, context = (:Input,), kwargs...)
    log!(M.manifold, internal_value(X), internal_value(p), internal_value(q))
    _update_basepoint!(M, X, p)
    is_vector(M, p, X; within = log, context = (:Output,), kwargs...)
    return X
end

function mid_point(M::ValidationManifold, p1, p2; kwargs...)
    is_point(M, p1; within = mid_point, context = (:Input,), kwargs...)
    is_point(M, p2; within = mid_point, context = (:Input,), kwargs...)
    q = mid_point(M.manifold, internal_value(p1), internal_value(p2))
    is_point(M, q; within = mid_point, context = (:Output,), kwargs...)
    return q
end

function mid_point!(M::ValidationManifold, q, p1, p2; kwargs...)
    is_point(M, p1; within = mid_point, context = (:Input,), kwargs...)
    is_point(M, p2; within = mid_point, context = (:Input,), kwargs...)
    mid_point!(M.manifold, internal_value(q), internal_value(p1), internal_value(p2))
    is_point(M, q; within = mid_point, context = (:Output,), kwargs...)
    return q
end

function norm(M::ValidationManifold, p, X; kwargs...)
    is_point(M, p; within = norm, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = norm, context = (:Input,), kwargs...)
    n = norm(M.manifold, internal_value(p), internal_value(X))
    (n < 0) &&
        _msg(M, DomainError(n, "Norm is negative."); within = norm, context = (:Output,))
    return n
end

number_eltype(::Type{ValidationMPoint{V}}) where {V} = number_eltype(V)
number_eltype(::Type{ValidationFibreVector{TType, V, P}}) where {TType, V, P} = number_eltype(V)

function project!(M::ValidationManifold, Y, p, X; kwargs...)
    is_point(M, p; within = project, context = (:Input,), kwargs...)
    project!(M.manifold, internal_value(Y), internal_value(p), internal_value(X))
    _update_basepoint!(M, Y, p)
    is_vector(M, p, Y; within = project, context = (:Output,), kwargs...)
    return Y
end

function rand(M::ValidationManifold; vector_at = nothing, kwargs...)
    if vector_at !== nothing
        is_point(M, vector_at; within = rand, context = (:Input,), kwargs...)
    end
    pX = rand(M.manifold; vector_at = vector_at, kwargs...)
    if vector_at !== nothing
        is_vector(M, vector_at, pX; within = rand, context = (:Output,), kwargs...)
    else
        is_point(M, pX; within = rand, context = (:Output,), kwargs...)
    end
    return pX
end

function riemann_tensor(M::ValidationManifold, p, X, Y, Z; kwargs...)
    is_point(M, p; within = riemann_tensor, context = (:Input,), kwargs...)
    for W in (X, Y, Z)
        is_vector(M, p, W; within = riemann_tensor, context = (:Input,), kwargs...)
    end
    W = riemann_tensor(
        M.manifold,
        internal_value(p),
        internal_value(X),
        internal_value(Y),
        internal_value(Z),
    )
    is_vector(M, p, W; within = riemann_tensor, context = (:Output,), kwargs...)
    return ValidationTangentVector(W, M.store_base_point ? copy(M, p) : nothing)
end

function riemann_tensor!(M::ValidationManifold, W, p, X, Y, Z; kwargs...)
    is_point(M, p; within = riemann_tensor, context = (:Input,), kwargs...)
    for W in (X, Y, Z)
        is_vector(M, p, W; within = riemann_tensor, context = (:Input,), kwargs...)
    end
    riemann_tensor!(
        M.manifold,
        internal_value(W),
        internal_value(p),
        internal_value(X),
        internal_value(Y),
        internal_value(Z),
    )
    _update_basepoint!(M, W, p)
    is_vector(M, p, W; within = riemann_tensor, context = (:Output,), kwargs...)
    return W
end

function show(io::IO, M::ValidationManifold)
    s = """
    ValidationManifold of $(M.manifold)
        * mode = :$(M.mode)
        * store_base_point = $(M.store_base_point)
    """
    if length(M.ignore_contexts) > 0
        s *= "    * ignore_context = $(M.ignore_contexts)\n"
    end
    if length(M.ignore_functions) > 0
        s *= "    * ignore_functions = $(M.ignore_functions)"
    end
    return print(io, s)
end

function vector_transport_to(
        M::ValidationManifold,
        p,
        X,
        q,
        m::AbstractVectorTransportMethod;
        kwargs...,
    )
    is_point(M, q; within = vector_transport_to, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = vector_transport_to, context = (:Input,), kwargs...)
    Y = vector_transport_to(
        M.manifold,
        internal_value(p),
        internal_value(X),
        internal_value(q),
        m,
    )
    is_vector(M, q, Y; within = vector_transport_to, context = (:Output,), kwargs...)
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
    is_point(M, q; within = vector_transport_to, context = (:Input,), kwargs...)
    is_vector(M, p, X; within = vector_transport_to, context = (:Input,), kwargs...)
    vector_transport_to!(
        M.manifold,
        internal_value(Y),
        internal_value(p),
        internal_value(X),
        internal_value(q),
        m,
    )
    _update_basepoint!(M, Y, q)
    is_vector(M, q, Y; within = vector_transport_to, context = (:Output,), kwargs...)
    return Y
end

function zero_vector(M::ValidationManifold, p; kwargs...)
    is_point(M, p; within = zero_vector, context = (:Input,), kwargs...)
    w = zero_vector(M.manifold, internal_value(p))
    is_vector(M, p, w; within = zero_vector, context = (:Output,), kwargs...)
    return w
end

function zero_vector!(M::ValidationManifold, X, p; kwargs...)
    is_point(M, p; within = zero_vector, context = (:Input,), kwargs...)
    zero_vector!(M.manifold, internal_value(X), internal_value(p); kwargs...)
    _update_basepoint!(M, X, p)
    is_vector(M, p, X; within = zero_vector, context = (:Output,), kwargs...)
    return X
end
