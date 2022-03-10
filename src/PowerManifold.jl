"""
    AbstractPowerRepresentation

An abstract representation type of points and tangent vectors on a power manifold.
"""
abstract type AbstractPowerRepresentation end

"""
    NestedPowerRepresentation

Representation of points and tangent vectors on a power manifold using arrays
of size equal to `TSize` of a [`PowerManifold`](@ref).
Each element of such array stores a single point or tangent vector.

For modifying operations, each element of the outer array is modified in-place, differently
than in [`NestedReplacingPowerRepresentation`](@ref).
"""
struct NestedPowerRepresentation <: AbstractPowerRepresentation end

"""
    NestedReplacingPowerRepresentation

Representation of points and tangent vectors on a power manifold using arrays
of size equal to `TSize` of a [`PowerManifold`](@ref).
Each element of such array stores a single point or tangent vector.

For modifying operations, each element of the outer array is replaced using non-modifying
operations, differently than for [`NestedReplacingPowerRepresentation`](@ref).
"""
struct NestedReplacingPowerRepresentation <: AbstractPowerRepresentation end

@doc raw"""
    AbstractPowerManifold{𝔽,M,TPR} <: AbstractManifold{𝔽}

An abstract [`AbstractManifold`](@ref) to represent manifolds that are build as powers
of another [`AbstractManifold`](@ref) `M` with representation type `TPR`, a subtype of
[`AbstractPowerRepresentation`](@ref).
"""
abstract type AbstractPowerManifold{
    𝔽,
    M<:AbstractManifold{𝔽},
    TPR<:AbstractPowerRepresentation,
} <: AbstractManifold{𝔽} end

@doc raw"""
    PowerManifold{𝔽,TM<:AbstractManifold,TSize<:Tuple,TPR<:AbstractPowerRepresentation} <: AbstractPowerManifold{𝔽,TM}

The power manifold $\mathcal M^{n_1× n_2 × … × n_d}$ with power geometry
 `TSize` statically defines the number of elements along each axis.

For example, a manifold-valued time series would be represented by a power manifold with
$d$ equal to 1 and $n_1$ equal to the number of samples. A manifold-valued image
(for example in diffusion tensor imaging) would be represented by a two-axis power
manifold ($d=2$) with $n_1$ and $n_2$ equal to width and height of the image.

While the size of the manifold is static, points on the power manifold
would not be represented by statically-sized arrays.

# Constructor

    PowerManifold(M::PowerManifold, N_1, N_2, ..., N_d)
    PowerManifold(M::AbstractManifold, NestedPowerRepresentation(), N_1, N_2, ..., N_d)
    M^(N_1, N_2, ..., N_d)

Generate the power manifold $M^{N_1 × N_2 × … × N_d}$.
By default, a [`PowerManifold`](@ref} is expanded further, i.e. for `M=PowerManifold(N,3)`
`PowerManifold(M,2)` is equivalend to `PowerManifold(N,3,2)`. Points are then 3×2 matrices
of points on `N`.
Providing a [`NestedPowerRepresentation`](@ref) as the second argument to the constructor
can be used to nest manifold, i.e. `PowerManifold(M,NestedPowerRepresentation(),2)`
represents vectors of length 2 whose elements are vectors of length 3 of points on N
in a nested array representation.

Since there is no default [`AbstractPowerRepresentation`](@ref) within this interface, the
`^` operator is only available for `PowerManifold`s and concatenates dimensions.
"""
struct PowerManifold{𝔽,TM<:AbstractManifold{𝔽},TSize,TPR<:AbstractPowerRepresentation} <:
       AbstractPowerManifold{𝔽,TM,TPR}
    manifold::TM
end

function PowerManifold(
    M::AbstractManifold{𝔽},
    ::TPR,
    size::Integer...,
) where {𝔽,TPR<:AbstractPowerRepresentation}
    return PowerManifold{𝔽,typeof(M),Tuple{size...},TPR}(M)
end
function PowerManifold(
    M::PowerManifold{𝔽,TM,TSize,TPR},
    size::Integer...,
) where {𝔽,TM<:AbstractManifold{𝔽},TSize,TPR<:AbstractPowerRepresentation}
    return PowerManifold{𝔽,TM,Tuple{TSize.parameters...,size...},TPR}(M.manifold)
end
function PowerManifold(
    M::PowerManifold{𝔽,TM,TSize},
    ::TPR,
    size::Integer...,
) where {𝔽,TM<:AbstractManifold{𝔽},TSize,TPR<:AbstractPowerRepresentation}
    return PowerManifold{𝔽,TM,Tuple{TSize.parameters...,size...},TPR}(M.manifold)
end
function PowerManifold(
    M::PowerManifold{𝔽,TM,TSize},
    ::TPR,
    size::Integer...,
) where {
    𝔽,
    TM<:AbstractManifold{𝔽},
    TSize,
    TPR<:Union{NestedPowerRepresentation,NestedReplacingPowerRepresentation},
}
    return PowerManifold{𝔽,PowerManifold{𝔽,TM,TSize},Tuple{size...},TPR}(M)
end

"""
    PowerBasisData{TB<:AbstractArray}

Data storage for an array of basis data.
"""
struct PowerBasisData{TB<:AbstractArray}
    bases::TB
end

const PowerManifoldNested =
    AbstractPowerManifold{𝔽,<:AbstractManifold{𝔽},NestedPowerRepresentation} where {𝔽}

const PowerManifoldNestedReplacing = AbstractPowerManifold{
    𝔽,
    <:AbstractManifold{𝔽},
    NestedReplacingPowerRepresentation,
} where {𝔽}

_access_nested(x, i::Int) = x[i]
_access_nested(x, i::Tuple) = x[i...]

function Base.:^(
    M::PowerManifold{
        𝔽,
        TM,
        TSize,
        <:Union{NestedPowerRepresentation,NestedReplacingPowerRepresentation},
    },
    size::Integer...,
) where {𝔽,TM<:AbstractManifold{𝔽},TSize}
    return PowerManifold(M, size...)
end

function allocate_result(M::PowerManifoldNested, f, x...)
    if representation_size(M.manifold) === ()
        return allocate(x[1])
    else
        return [
            allocate_result(M.manifold, f, map(y -> _access_nested(y, i), x)...) for
            i in get_iterator(M)
        ]
    end
end
# avoid ambituities - though usually not used
function allocate_result(
    M::PowerManifoldNested,
    f::typeof(get_coordinates),
    p,
    X,
    B::AbstractBasis,
)
    return invoke(
        allocate_result,
        Tuple{AbstractManifold,typeof(get_coordinates),Any,Any,AbstractBasis},
        M,
        f,
        p,
        X,
        B,
    )
end
function allocate_result(::PowerManifoldNestedReplacing, f, x...)
    return copy(x[1])
end
# the following is not used but necessary to avoid ambiguities
function allocate_result(
    M::PowerManifoldNestedReplacing,
    f::typeof(get_coordinates),
    p,
    X,
    B::AbstractBasis,
)
    return invoke(
        allocate_result,
        Tuple{AbstractManifold,typeof(get_coordinates),Any,Any,AbstractBasis},
        M,
        f,
        p,
        X,
        B,
    )
end
function allocate_result(M::PowerManifoldNested, f::typeof(get_vector), p, X)
    return [allocate_result(M.manifold, f, _access_nested(p, i)) for i in get_iterator(M)]
end
function allocate_result(::PowerManifoldNestedReplacing, ::typeof(get_vector), p, X)
    return copy(p)
end
function allocation_promotion_function(M::AbstractPowerManifold, f, args::Tuple)
    return allocation_promotion_function(M.manifold, f, args)
end


"""
    check_point(M::AbstractPowerManifold, p; kwargs...)

Check whether `p` is a valid point on an [`AbstractPowerManifold`](@ref) `M`,
i.e. each element of `p` has to be a valid point on the base manifold.
If `p` is not a point on `M` a [`CompositeManifoldError`](@ref) consisting of all error messages of the
components, for which the tests fail is returned.

The tolerance for the last test can be set using the `kwargs...`.
"""
function check_point(M::AbstractPowerManifold, p; kwargs...)
    rep_size = representation_size(M.manifold)
    e = [
        (i, check_point(M.manifold, _read(M, rep_size, p, i); kwargs...)) for
        i in get_iterator(M)
    ]
    errors = filter((x) -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

"""
    check_power_size(M, p)
    check_power_size(M, p, X)

Check whether p hase the right size to represent points on M generically, i.e. just cheking the overall sizes, not the individual ones per manifold
"""
function check_power_size(M::AbstractPowerManifold, p)
    d = prod(representation_size(M.manifold)) * prod(power_dimensions(M))
    (d != length(p)) && return DomainError(
        length(p),
        "The point $p can not be a point on $M, since its number of elements does not match the required overall representation size ($d)",
    )
    return nothing
end
function check_power_size(M::Union{PowerManifoldNested,PowerManifoldNestedReplacing}, p)
    d = prod(power_dimensions(M))
    (d != length(p)) && return DomainError(
        length(p),
        "The point $p can not be a point on $M, since its number of elements does not match the power dimensions ($d)",
    )
    return nothing
end

function check_power_size(M::AbstractPowerManifold, p, X)
    d = prod(representation_size(M.manifold)) * prod(power_dimensions(M))
    (d != length(X)) && return DomainError(
        length(X),
        "The tangent vector $X can not belong to a trangent space at on $M, since its number of elements does not match the required overall representation size ($d)",
    )
    return nothing
end
function check_power_size(M::Union{PowerManifoldNested,PowerManifoldNestedReplacing}, p, X)
    d = prod(power_dimensions(M))
    (d != length(X)) && return DomainError(
        length(X),
        "The point $p can not be a point on $M, since its number of elements does not match the power dimensions ($d)",
    )
    return nothing
end

function check_size(M::AbstractPowerManifold, p)
    cps = check_power_size(M, p)
    (cps === nothing) || return cps
    rep_size = representation_size(M.manifold)
    e = [(i, check_size(M.manifold, _read(M, rep_size, p, i))) for i in get_iterator(M)]
    errors = filter((x) -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

function check_size(M::AbstractPowerManifold, p, X)
    cps = check_power_size(M, p, X)
    (cps === nothing) || return cps
    rep_size = representation_size(M.manifold)
    e = [
        (i, check_size(M.manifold, _read(M, rep_size, p, i), _read(M, rep_size, X, i);)) for i in get_iterator(M)
    ]
    errors = filter((x) -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

"""
    check_vector(M::AbstractPowerManifold, p, X; kwargs... )

Check whether `X` is a tangent vector to `p` an the [`AbstractPowerManifold`](@ref)
`M`, i.e. atfer [`check_point`](@ref)`(M, p)`, and all projections to
base manifolds must be respective tangent vectors.
If `X` is not a tangent vector to `p` on `M` a [`CompositeManifoldError`](@ref) consisting of all error
messages of the components, for which the tests fail is returned.

The tolerance for the last test can be set using the `kwargs...`.
"""
function check_vector(M::AbstractPowerManifold, p, X; kwargs...)
    rep_size = representation_size(M.manifold)
    e = [
        (
            i,
            check_vector(
                M.manifold,
                _read(M, rep_size, p, i),
                _read(M, rep_size, X, i);
                kwargs...,
            ),
        ) for i in get_iterator(M)
    ]
    errors = filter((x) -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

@doc raw"""
    copyto!(M::PowerManifoldNested, q, p)

Copy the values elementwise, i.e. call `copyto!(M.manifold, b, a)` for all elements `a` and
`b` of `p` and `q`, respectively.
"""
function copyto!(M::PowerManifoldNested, q, p)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        copyto!(M.manifold, _write(M, rep_size, q, i), _read(M, rep_size, p, i))
    end
    return q
end

@doc raw"""
    copyto!(M::PowerManifoldNested, Y, p, X)

Copy the values elementwise, i.e. call `copyto!(M.manifold, B, a, A)` for all elements
`A`, `a` and `B` of `X`, `p`, and `Y`, respectively.
"""
function copyto!(M::PowerManifoldNested, Y, p, X)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        copyto!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
        )
    end
    return Y
end

@doc raw"""
    distance(M::AbstractPowerManifold, p, q)

Compute the distance between `q` and `p` on an [`AbstractPowerManifold`](@ref),
i.e. from the element wise distances the Forbenius norm is computed.
"""
function distance(M::AbstractPowerManifold, p, q)
    sum_squares = zero(number_eltype(p))
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        sum_squares +=
            distance(M.manifold, _read(M, rep_size, p, i), _read(M, rep_size, q, i))^2
    end
    return sqrt(sum_squares)
end

@doc raw"""
    exp(M::AbstractPowerManifold, p, X)

Compute the exponential map from `p` in direction `X` on the [`AbstractPowerManifold`](@ref) `M`,
which can be computed using the base manifolds exponential map elementwise.
"""
exp(::AbstractPowerManifold, ::Any...)

function exp!(M::AbstractPowerManifold, q, p, X)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        exp!(
            M.manifold,
            _write(M, rep_size, q, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
        )
    end
    return q
end
function exp!(M::PowerManifoldNestedReplacing, q, p, X)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        q[i...] = exp(M.manifold, _read(M, rep_size, p, i), _read(M, rep_size, X, i))
    end
    return q
end

function get_basis(M::AbstractPowerManifold, p, B::AbstractBasis)
    rep_size = representation_size(M.manifold)
    vs = [get_basis(M.manifold, _read(M, rep_size, p, i), B) for i in get_iterator(M)]
    return CachedBasis(B, PowerBasisData(vs))
end
function get_basis(M::AbstractPowerManifold, p, B::DiagonalizingOrthonormalBasis)
    rep_size = representation_size(M.manifold)
    vs = [
        get_basis(
            M.manifold,
            _read(M, rep_size, p, i),
            DiagonalizingOrthonormalBasis(_read(M, rep_size, B.frame_direction, i)),
        ) for i in get_iterator(M)
    ]
    return CachedBasis(B, PowerBasisData(vs))
end

"""
    get_component(M::AbstractPowerManifold, p, idx...)

Get the component of a point `p` on an [`AbstractPowerManifold`](@ref) `M` at index `idx`.
"""
function get_component(M::AbstractPowerManifold, p, idx...)
    rep_size = representation_size(M.manifold)
    return _read(M, rep_size, p, idx)
end

function get_coordinates(M::AbstractPowerManifold, p, X, B::AbstractBasis)
    rep_size = representation_size(M.manifold)
    vs = [
        get_coordinates(M.manifold, _read(M, rep_size, p, i), _read(M, rep_size, X, i), B) for i in get_iterator(M)
    ]
    return reduce(vcat, reshape(vs, length(vs)))
end
function get_coordinates(
    M::AbstractPowerManifold,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis,<:PowerBasisData},
) where {𝔽}
    rep_size = representation_size(M.manifold)
    vs = [
        get_coordinates(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _access_nested(B.data.bases, i),
        ) for i in get_iterator(M)
    ]
    return reduce(vcat, reshape(vs, length(vs)))
end

function get_coordinates!(M::AbstractPowerManifold, c, p, X, B::AbstractBasis)
    rep_size = representation_size(M.manifold)
    dim = manifold_dimension(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        # TODO: this view is really suboptimal when `dim` can be statically determined
        get_coordinates!(
            M.manifold,
            view(c, v_iter:(v_iter + dim - 1)),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            B,
        )
        v_iter += dim
    end
    return c
end
function get_coordinates!(
    M::AbstractPowerManifold,
    c,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis,<:PowerBasisData},
) where {𝔽}
    rep_size = representation_size(M.manifold)
    dim = manifold_dimension(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        # TODO: this view is really suboptimal when `dim` can be statically determined
        get_coordinates!(
            M.manifold,
            view(c, v_iter:(v_iter + dim - 1)),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _access_nested(B.data.bases, i),
        )
        v_iter += dim
    end
    return c
end

get_iterator(::PowerManifold{𝔽,<:AbstractManifold{𝔽},Tuple{N}}) where {𝔽,N} = Base.OneTo(N)
@generated function get_iterator(
    ::PowerManifold{𝔽,<:AbstractManifold{𝔽},SizeTuple},
) where {𝔽,SizeTuple}
    size_tuple = size_to_tuple(SizeTuple)
    return Base.product(map(Base.OneTo, size_tuple)...)
end

function get_vector(
    M::AbstractPowerManifold,
    p,
    c,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:PowerBasisData},
) where {𝔽}
    Y = allocate_result(M, get_vector, p, c)
    return get_vector!(M, Y, p, c, B)
end
function get_vector!(
    M::AbstractPowerManifold,
    Y,
    p,
    c,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:PowerBasisData},
) where {𝔽}
    dim = manifold_dimension(M.manifold)
    rep_size = representation_size(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        get_vector!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            c[v_iter:(v_iter + dim - 1)],
            _access_nested(B.data.bases, i),
        )
        v_iter += dim
    end
    return Y
end
function get_vector!(
    M::PowerManifoldNestedReplacing,
    Y,
    p,
    c,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:PowerBasisData},
) where {𝔽}
    dim = manifold_dimension(M.manifold)
    rep_size = representation_size(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        Y[i...] = get_vector(
            M.manifold,
            _read(M, rep_size, p, i),
            c[v_iter:(v_iter + dim - 1)],
            _access_nested(B.data.bases, i),
        )
        v_iter += dim
    end
    return Y
end
function get_vector(M::AbstractPowerManifold, p, c, B::AbstractBasis)
    Y = allocate_result(M, get_vector, p, c)
    return get_vector!(M, Y, p, c, B)
end
function get_vector!(M::AbstractPowerManifold, Y, p, c, B::AbstractBasis)
    dim = manifold_dimension(M.manifold)
    rep_size = representation_size(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        get_vector!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            c[v_iter:(v_iter + dim - 1)],
            B,
        )
        v_iter += dim
    end
    return Y
end
function get_vector!(M::PowerManifoldNestedReplacing, Y, p, c, B::AbstractBasis)
    dim = manifold_dimension(M.manifold)
    rep_size = representation_size(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        Y[i...] = get_vector(
            M.manifold,
            _read(M, rep_size, p, i),
            c[v_iter:(v_iter + dim - 1)],
            B,
        )
        v_iter += dim
    end
    return Y
end

"""
    getindex(p, M::AbstractPowerManifold, i::Union{Integer,Colon,AbstractVector}...)
    p[M::AbstractPowerManifold, i...]

Access the element(s) at index `[i...]` of a point `p` on an [`AbstractPowerManifold`](@ref)
`M` by linear or multidimensional indexing.
See also [Array Indexing](https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing-1) in Julia.
"""
Base.@propagate_inbounds function Base.getindex(
    p::AbstractArray,
    M::AbstractPowerManifold,
    I::Union{Integer,Colon,AbstractVector}...,
)
    return get_component(M, p, I...)
end

@doc raw"""
    injectivity_radius(M::AbstractPowerManifold[, p])

the injectivity radius on an [`AbstractPowerManifold`](@ref) is for the global case
equal to the one of its base manifold. For a given point `p` it's equal to the
minimum of all radii in the array entries.
"""
function injectivity_radius(M::AbstractPowerManifold, p)
    radius = 0.0
    initialized = false
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        cur_rad = injectivity_radius(M.manifold, _read(M, rep_size, p, i))
        if initialized
            radius = min(cur_rad, radius)
        else
            radius = cur_rad
            initialized = true
        end
    end
    return radius
end
injectivity_radius(M::AbstractPowerManifold) = injectivity_radius(M.manifold)
function injectivity_radius(M::AbstractPowerManifold, ::AbstractRetractionMethod)
    return injectivity_radius(M)
end

@doc raw"""
    inner(M::AbstractPowerManifold, p, X, Y)

Compute the inner product of `X` and `Y` from the tangent space at `p` on an
[`AbstractPowerManifold`](@ref) `M`, i.e. for each arrays entry the tangent
vector entries from `X` and `Y` are in the tangent space of the corresponding
element from `p`.
The inner product is then the sum of the elementwise inner products.
"""
function inner(M::AbstractPowerManifold, p, X, Y)
    result = zero(number_eltype(X))
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        result += inner(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, Y, i),
        )
    end
    return result
end

function Base.isapprox(M::AbstractPowerManifold, p, q; kwargs...)
    result = true
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        result &= isapprox(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i);
            kwargs...,
        )
    end
    return result
end
function Base.isapprox(M::AbstractPowerManifold, p, X, Y; kwargs...)
    result = true
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        result &= isapprox(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, Y, i);
            kwargs...,
        )
    end
    return result
end

@doc raw"""
    inverse_retract(M::AbstractPowerManifold, p, q, m::AbstractInverseRetractionMethod)

Compute the inverse retraction from `p` with respect to `q` on an [`AbstractPowerManifold`](@ref) `M`
using an [`AbstractInverseRetractionMethod`](@ref).
Then this method is performed elementwise, so the inverse
retraction method has to be one that is available on the base [`AbstractManifold`](@ref).
"""
inverse_retract(::AbstractPowerManifold, ::Any...)

function inverse_retract!(
    M::AbstractPowerManifold,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = LogarithmicInverseRetraction(),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        inverse_retract!(
            M.manifold,
            _write(M, rep_size, X, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i),
            m,
        )
    end
    return X
end
function inverse_retract!(
    M::PowerManifoldNestedReplacing,
    X,
    p,
    q,
    m::AbstractInverseRetractionMethod = LogarithmicInverseRetraction(),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        X[i...] = inverse_retract(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i),
            m,
        )
    end
    return X
end

@doc raw"""
    log(M::AbstractPowerManifold, p, q)

Compute the logarithmic map from `p` to `q` on the [`AbstractPowerManifold`](@ref) `M`,
which can be computed using the base manifolds logarithmic map elementwise.
"""
log(::AbstractPowerManifold, ::Any...)

function log!(M::AbstractPowerManifold, X, p, q)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        log!(
            M.manifold,
            _write(M, rep_size, X, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i),
        )
    end
    return X
end
function log!(M::PowerManifoldNestedReplacing, X, p, q)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        X[i...] = log(M.manifold, _read(M, rep_size, p, i), _read(M, rep_size, q, i))
    end
    return X
end

@doc raw"""
    manifold_dimension(M::PowerManifold)

Returns the manifold-dimension of an [`PowerManifold`](@ref) `M`
$=\mathcal N = (\mathcal M)^{n_1,…,n_d}$, i.e. with $n=(n_1,…,n_d)$ the array
size of the power manifold and $d_{\mathcal M}$ the dimension of the base manifold
$\mathcal M$, the manifold is of dimension

````math
\dim(\mathcal N) = \dim(\mathcal M)\prod_{i=1}^d n_i = n_1n_2\cdot…\cdot n_d \dim(\mathcal M).
````
"""
function manifold_dimension(M::PowerManifold{𝔽,<:AbstractManifold,TSize}) where {𝔽,TSize}
    return manifold_dimension(M.manifold) * prod(size_to_tuple(TSize))
end

function mid_point!(M::AbstractPowerManifold, q, p1, p2)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        mid_point!(
            M.manifold,
            _write(M, rep_size, q, i),
            _read(M, rep_size, p1, i),
            _read(M, rep_size, p2, i),
        )
    end
    return q
end
function mid_point!(M::PowerManifoldNestedReplacing, q, p1, p2)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        q[i...] =
            mid_point(M.manifold, _read(M, rep_size, p1, i), _read(M, rep_size, p2, i))
    end
    return q
end

@doc raw"""
    norm(M::AbstractPowerManifold, p, X)

Compute the norm of `X` from the tangent space of `p` on an
[`AbstractPowerManifold`](@ref) `M`, i.e. from the element wise norms the
Frobenius norm is computed.
"""
function LinearAlgebra.norm(M::AbstractPowerManifold, p, X)
    sum_squares = zero(number_eltype(X))
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        sum_squares +=
            norm(M.manifold, _read(M, rep_size, p, i), _read(M, rep_size, X, i))^2
    end
    return sqrt(sum_squares)
end


function parallel_transport_direction!(M::AbstractPowerManifold, Y, p, X, d)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        parallel_transport_direction!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, d, i),
        )
    end
    return Y
end
function parallel_transport_direction(M::AbstractPowerManifold, p, X, d)
    Y = allocate_result(M, vector_transport_direction, p, X, d)
    return parallel_transport_direction!(M, Y, p, X, d)
end

function parallel_transport_direction!(M::PowerManifoldNestedReplacing, Y, p, X, d)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        Y[i...] = parallel_transport_direction(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, d, i),
        )
    end
    return Y
end
function parallel_transport_direction(M::PowerManifoldNestedReplacing, p, X, d)
    Y = allocate_result(M, parallel_transport_direction, p, X, d)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        Y[i...] = parallel_transport_direction(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, d, i),
        )
    end
    return Y
end

function parallel_transport_to(M::AbstractPowerManifold, p, X, q)
    Y = allocate_result(M, vector_transport_to, p, X)
    return parallel_transport_to!(M, Y, p, X, q)
end
function parallel_transport_to!(M::AbstractPowerManifold, Y, p, X, q)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        vector_transport_to!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, q, i),
        )
    end
    return Y
end
function parallel_transport_to!(M::PowerManifoldNestedReplacing, Y, p, X, q)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        Y[i...] = parallel_transport_to(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, q, i),
        )
    end
    return Y
end

@doc raw"""
    power_dimensions(M::PowerManifold)

return the power of `M`,
"""
function power_dimensions(::PowerManifold{𝔽,<:AbstractManifold,TSize}) where {𝔽,TSize}
    return size_to_tuple(TSize)
end

@doc raw"""
    project(M::AbstractPowerManifold, p)

Project the point `p` from the embedding onto the [`AbstractPowerManifold`](@ref) `M`
by projecting all components.
"""
project(::AbstractPowerManifold, ::Any)

function project!(M::AbstractPowerManifold, q, p)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        project!(M.manifold, _write(M, rep_size, q, i), _read(M, rep_size, p, i))
    end
    return q
end
function project!(M::PowerManifoldNestedReplacing, q, p)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        q[i...] = project(M.manifold, _read(M, rep_size, p, i))
    end
    return q
end

@doc raw"""
    project(M::AbstractPowerManifold, p, X)

Project the point `X` onto the tangent space at `p` on the
[`AbstractPowerManifold`](@ref) `M` by projecting all components.
"""
project(::AbstractPowerManifold, ::Any, ::Any)

function project!(M::AbstractPowerManifold, Z, q, Y)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        project!(
            M.manifold,
            _write(M, rep_size, Z, i),
            _read(M, rep_size, q, i),
            _read(M, rep_size, Y, i),
        )
    end
    return Z
end
function project!(M::PowerManifoldNestedReplacing, Z, q, Y)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        q[i...] = project(M.manifold, _read(M, rep_size, q, i), _read(M, rep_size, Y, i))
    end
    return Z
end

Base.@propagate_inbounds @inline function _read(
    M::AbstractPowerManifold,
    rep_size::Tuple,
    x::AbstractArray,
    i::Int,
)
    return _read(M, rep_size, x, (i,))
end

Base.@propagate_inbounds @inline function _read(
    ::Union{PowerManifoldNested,PowerManifoldNestedReplacing},
    rep_size::Tuple,
    x::AbstractArray,
    i::Tuple,
)
    return x[i...]
end

@generated function rep_size_to_colons(rep_size::Tuple)
    N = length(rep_size.parameters)
    return ntuple(i -> Colon(), N)
end


@doc raw"""
    retract(M::AbstractPowerManifold, p, X, method::AbstractRetractionMethod)

Compute the retraction from `p` with tangent vector `X` on an [`AbstractPowerManifold`](@ref) `M`
using a [`AbstractRetractionMethod`](@ref).
Then this method is performed elementwise, so the retraction
method has to be one that is available on the base [`AbstractManifold`](@ref).
"""
retract(::AbstractPowerManifold, ::Any...)

function retract!(
    M::AbstractPowerManifold,
    q,
    p,
    X,
    m::AbstractRetractionMethod = ExponentialRetraction(),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        retract!(
            M.manifold,
            _write(M, rep_size, q, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            m,
        )
    end
    return q
end
function retract!(
    M::PowerManifoldNestedReplacing,
    q,
    p,
    X,
    m::AbstractRetractionMethod = ExponentialRetraction(),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        q[i...] = retract(M.manifold, _read(M, rep_size, p, i), _read(M, rep_size, X, i), m)
    end
    return q
end

"""
    set_component!(M::AbstractPowerManifold, q, p, idx...)

Set the component of a point `q` on an [`AbstractPowerManifold`](@ref) `M` at index `idx`
to `p`, which itself is a point on the [`AbstractManifold`](@ref) the power manifold is build on.
"""
function set_component!(M::AbstractPowerManifold, q, p, idx...)
    rep_size = representation_size(M.manifold)
    return copyto!(_write(M, rep_size, q, idx), p)
end
function set_component!(::PowerManifoldNestedReplacing, q, p, idx...)
    return q[idx...] = p
end
"""
    setindex!(q, p, M::AbstractPowerManifold, i::Union{Integer,Colon,AbstractVector}...)
    q[M::AbstractPowerManifold, i...] = p

Set the element(s) at index `[i...]` of a point `q` on an [`AbstractPowerManifold`](@ref)
`M` by linear or multidimensional indexing to `q`.
See also [Array Indexing](https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing-1) in Julia.
"""
Base.@propagate_inbounds function Base.setindex!(
    q::AbstractArray,
    p,
    M::AbstractPowerManifold,
    I::Union{Integer,Colon,AbstractVector}...,
)
    return set_component!(M, q, p, I...)
end

function Base.show(io::IO, M::PowerManifold{𝔽,TM,TSize,TPR}) where {𝔽,TM,TSize,TPR}
    return print(
        io,
        "PowerManifold($(M.manifold), $(TPR()), $(join(TSize.parameters, ", ")))",
    )
end

function Base.show(
    io::IO,
    mime::MIME"text/plain",
    B::CachedBasis{𝔽,T,D},
) where {T<:AbstractBasis,D<:PowerBasisData,𝔽}
    println(io, "$(T()) for a power manifold")
    for i in Base.product(map(Base.OneTo, size(B.data.bases))...)
        println(io, "Basis for component $i:")
        show(io, mime, _access_nested(B.data.bases, i))
        println(io)
    end
    return nothing
end

function vector_transport_direction!(
    M::AbstractPowerManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        vector_transport_direction!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, d, i),
            m,
        )
    end
    return Y
end
function vector_transport_direction(
    M::AbstractPowerManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    Y = allocate_result(M, vector_transport_direction, p, X, d)
    return vector_transport_direction!(M, Y, p, X, d, m)
end

function vector_transport_direction!(
    M::PowerManifoldNestedReplacing,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        Y[i...] = vector_transport_direction(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, d, i),
            m,
        )
    end
    return Y
end
function vector_transport_direction(
    M::PowerManifoldNestedReplacing,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    Y = allocate_result(M, vector_transport_direction, p, X, d)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        Y[i...] = vector_transport_direction(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, d, i),
            m,
        )
    end
    return Y
end

@doc raw"""
    vector_transport_to(M::AbstractPowerManifold, p, X, q, method::AbstractVectorTransportMethod)

Compute the vector transport the tangent vector `X`at `p` to `q` on the
[`PowerManifold`](@ref) `M` using an [`AbstractVectorTransportMethod`](@ref) `m`.
This method is performed elementwise, i.e. the method `m` has to be implemented on the
base manifold.
"""
vector_transport_to(
    ::AbstractPowerManifold,
    ::Any,
    ::Any,
    ::Any,
    ::AbstractVectorTransportMethod,
)
function vector_transport_to(
    M::AbstractPowerManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    Y = allocate_result(M, vector_transport_to, p, X)
    return vector_transport_to!(M, Y, p, X, q, m)
end
function vector_transport_to!(
    M::AbstractPowerManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        vector_transport_to!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, q, i),
            m,
        )
    end
    return Y
end
function vector_transport_to!(
    M::PowerManifoldNestedReplacing,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        Y[i...] = vector_transport_to(
            M.manifold,
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, q, i),
            m,
        )
    end
    return Y
end

"""
    view(p, M::AbstractPowerManifold, i::Union{Integer,Colon,AbstractVector}...)

Get the view of the element(s) at index `[i...]` of a point `p` on an
[`AbstractPowerManifold`](@ref) `M` by linear or multidimensional indexing.
"""
function Base.view(
    p::AbstractArray,
    M::AbstractPowerManifold,
    I::Union{Integer,Colon,AbstractVector}...,
)
    rep_size = representation_size(M.manifold)
    return _write(M, rep_size, p, I...)
end

@inline function _write(M::AbstractPowerManifold, rep_size::Tuple, x::AbstractArray, i::Int)
    return _write(M, rep_size, x, (i,))
end

@inline function _write(
    ::AbstractPowerManifold,
    rep_size::Tuple,
    x::AbstractArray,
    i::Tuple,
)
    return view(x[i...], rep_size_to_colons(rep_size)...)
end
@inline function _write(::PowerManifoldNested, ::Tuple{}, x::AbstractArray, i::Tuple)
    return view(x, i...)
end

function zero_vector!(M::AbstractPowerManifold, X, p)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        zero_vector!(M.manifold, _write(M, rep_size, X, i), _read(M, rep_size, p, i))
    end
    return X
end
function zero_vector!(M::PowerManifoldNestedReplacing, X, p)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        X[i...] = zero_vector(M.manifold, _read(M, rep_size, p, i))
    end
    return X
end
