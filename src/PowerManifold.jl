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

[`GraphManifold`](@ref) uses this representation.
"""
struct NestedPowerRepresentation <: AbstractPowerRepresentation end

@doc raw"""
    AbstractPowerManifold{𝔽,M,TPR} <: Manifold{𝔽}

An abstract [`Manifold`](@ref) to represent manifolds that are build as powers
of another [`Manifold`](@ref) `M` with representation type `TPR`, a subtype of
[`AbstractPowerRepresentation`](@ref).
"""
abstract type AbstractPowerManifold{𝔽,M<:Manifold{𝔽},TPR<:AbstractPowerRepresentation} <:
              Manifold{𝔽} end

@doc raw"""
    PowerManifold{𝔽,TM<:Manifold,TSize<:Tuple,TPR<:AbstractPowerRepresentation} <: AbstractPowerManifold{𝔽,TM}

The power manifold $\mathcal M^{n_1× n_2 × … × n_d}$ with power geometry
 `TSize` statically defines the number of elements along each axis.

For example, a manifold-valued time series would be represented by a power manifold with
$d$ equal to 1 and $n_1$ equal to the number of samples. A manifold-valued image
(for example in diffusion tensor imaging) would be represented by a two-axis power
manifold ($d=2$) with $n_1$ and $n_2$ equal to width and height of the image.

While the size of the manifold is static, points on the power manifold
would not be represented by statically-sized arrays. Operations on small
power manifolds might be faster if they are represented as [`ProductManifold`](@ref).

# Constructor

    PowerManifold(M, N_1, N_2, ..., N_d)
    PowerManifold(M, NestedPowerRepresentation(), N_1, N_2, ..., N_d)
    M^(N_1, N_2, ..., N_d)

Generate the power manifold $M^{N_1 × N_2 × … × N_d}$.
By default, the [`ArrayPowerRepresentation`](@ref) of points
and tangent vectors is used, although a different one, for example
[`NestedPowerRepresentation`](@ref), can be given as the second argument to the
constructor.
When `M` is a `PowerManifold` (not any [`AbstractPowerManifold`](@ref)) itself, given
dimensions will be appended to the dimensions already present, for example
`PowerManifold(PowerManifold(Sphere(2), 2), 3)` is equivalent to
`PowerManifold(Sphere(2), 2, 3)`. This feature preserves the representation of the inner
power manifold (unless it's explicitly overridden).
If you specify `NestedPowerRepresentation()`, the sizes are not concatenated but you end up
with a nested power manifold within a power manifold.
"""
struct PowerManifold{𝔽,TM<:Manifold{𝔽},TSize,TPR<:AbstractPowerRepresentation} <:
       AbstractPowerManifold{𝔽,TM,TPR}
    manifold::TM
end

function PowerManifold(
    M::Manifold{𝔽},
    ::TPR,
    size::Integer...,
) where {𝔽,TPR<:AbstractPowerRepresentation}
    return PowerManifold{𝔽,typeof(M),Tuple{size...},TPR}(M)
end
function PowerManifold(
    M::PowerManifold{𝔽,TM,TSize,TPR},
    size::Integer...,
) where {𝔽,TM<:Manifold{𝔽},TSize,TPR<:AbstractPowerRepresentation}
    return PowerManifold{𝔽,TM,Tuple{TSize.parameters...,size...},TPR}(M.manifold)
end
function PowerManifold(
    M::PowerManifold{𝔽,TM,TSize},
    ::TPR,
    size::Integer...,
) where {𝔽,TM<:Manifold{𝔽},TSize,TPR<:AbstractPowerRepresentation}
    return PowerManifold{𝔽,TM,Tuple{TSize.parameters...,size...},TPR}(M.manifold)
end
function PowerManifold(
    M::PowerManifold{𝔽,TM,TSize},
    ::NestedPowerRepresentation,
    size::Integer...,
) where {𝔽,TM<:Manifold{𝔽},TSize}
    return PowerManifold{
        𝔽,
        PowerManifold{𝔽,TM,TSize},
        Tuple{size...},
        NestedPowerRepresentation,
    }(
        M,
    )
end

"""
    PowerRetraction{TR<:AbstractRetractionMethod} <: AbstractRetractionMethod

The `PowerRetraction` avoids ambiguities between dispatching on the [`AbstractPowerManifold`](@ref)
and dispatching on the [`AbstractRetractionMethod`](@ref) and encapsulates this.
This container should only be used in rare cases outside of this package. Usually a
subtype of the [`AbstractPowerManifold`](@ref) should define a way how to treat
its [`AbstractRetractionMethod`](@ref)s.

# Constructor

    PowerRetraction(retraction::AbstractRetractionMethod)
"""
struct PowerRetraction{TR<:AbstractRetractionMethod} <: AbstractRetractionMethod
    retraction::TR
end

"""
    InversePowerRetraction{TR<:AbstractInverseRetractionMethod} <: AbstractInverseRetractionMethod

The `InversePowerRetraction` avoids ambiguities between dispatching on the [`AbstractPowerManifold`](@ref)
and dispatching on the [`AbstractInverseRetractionMethod`](@ref) and encapsulates this.
This container should only be used in rare cases outside of this package. Usually a
subtype of the [`AbstractPowerManifold`](@ref) should define a way how to treat
its [`AbstractRetractionMethod`](@ref)s.

# Constructor

    InversePowerRetraction(inverse_retractions::AbstractInverseRetractionMethod...)
"""
struct InversePowerRetraction{TR<:AbstractInverseRetractionMethod} <:
       AbstractInverseRetractionMethod
    inverse_retraction::TR
end

"""
    PowerVectorTransport{TR<:AbstractVectorTransportMethod} <:
       AbstractVectorTransportMethod

The `PowerVectorTransport` avoids ambiguities between dispatching on the [`AbstractPowerManifold`](@ref)
and dispatching on the [`AbstractVectorTransportMethod`](@ref) and encapsulates this.
This container should only be used in rare cases outside of this package. Usually a
subtype of the [`AbstractPowerManifold`](@ref) should define a way how to treat
its [`AbstractVectorTransportMethod`](@ref)s.

# Constructor

    PowerVectorTransport(method::AbstractVectorTransportMethod)
"""
struct PowerVectorTransport{TR<:AbstractVectorTransportMethod} <:
       AbstractVectorTransportMethod
    method::TR
end

"""
    PowerBasisData{TB<:AbstractArray}

Data storage for an array of basis data.
"""
struct PowerBasisData{TB<:AbstractArray}
    bases::TB
end

const PowerManifoldNested =
    AbstractPowerManifold{𝔽,<:Manifold{𝔽},NestedPowerRepresentation} where {𝔽}

_access_nested(x, i::Int) = x[i]
_access_nested(x, i::Tuple) = x[i...]

function allocate_result(M::PowerManifoldNested, f, x...)
    return [
        allocate_result(M.manifold, f, map(y -> _access_nested(y, i), x)...)
        for i in get_iterator(M)
    ]
end

function allocate_result(M::PowerManifoldNested, f::typeof(get_vector), p, X)
    return [allocate_result(M.manifold, f, _access_nested(p, i)) for i in get_iterator(M)]
end
function allocate_result(
    M::PowerManifoldNested,
    f::typeof(get_coordinates),
    p,
    X,
    B::AbstractBasis,
)
    return invoke(
        allocate_result,
        Tuple{Manifold,typeof(get_coordinates),Any,Any,typeof(B)},
        M,
        f,
        p,
        X,
        B,
    )
end

function allocation_promotion_function(M::AbstractPowerManifold, f, args::Tuple)
    return allocation_promotion_function(M.manifold, f, args)
end


"""
    check_manifold_point(M::AbstractProductManifold, p; kwargs...)

Check whether `p` is a valid point on an [`AbstractPowerManifold`](@ref) `M`,
i.e. each element of `p` has to be a valid point on the base manifold.
If `p` is not a point on `M` a `CompositeManifoldError` consisting of all error messages of the
components, for which the tests fail is returned.

The tolerance for the last test can be set using the `kwargs...`.
"""
function check_manifold_point(M::AbstractPowerManifold, p; kwargs...)
    rep_size = representation_size(M.manifold)
    e = [
        (i, check_manifold_point(M.manifold, _read(M, rep_size, p, i); kwargs...))
        for i in get_iterator(M)
    ]
    errors = filter((x) -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

"""
    check_tangent_vector(M::AbstractPowerManifold, p, X; check_base_point = true, kwargs... )

Check whether `X` is a tangent vector to `p` an the [`AbstractPowerManifold`](@ref)
`M`, i.e. atfer [`check_manifold_point`](@ref)`(M, p)`, and all projections to
base manifolds must be respective tangent vectors.
The optional parameter `check_base_point` indicates, whether to call [`check_manifold_point`](@ref)  for `p`.
If `X` is not a tangent vector to `p` on `M` a `CompositeManifoldError` consisting of all error
messages of the components, for which the tests fail is returned.


The tolerance for the last test can be set using the `kwargs...`.
"""
function check_tangent_vector(
    M::AbstractPowerManifold,
    p,
    X;
    check_base_point = true,
    kwargs...,
)
    if check_base_point
        mpe = check_manifold_point(M, p)
        mpe === nothing || return mpe
    end
    rep_size = representation_size(M.manifold)
    e = [
        (
            i,
            check_tangent_vector(
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
for BT in ManifoldsBase.DISAMBIGUATION_BASIS_TYPES
    if BT == DiagonalizingOrthonormalBasis
        continue
    end
    eval(quote
        @invoke_maker 3 AbstractBasis get_basis(M::AbstractPowerManifold, p, B::$BT)
    end)
end

"""
    get_component(M::AbstractPowerManifold, p, idx...)

Get the component of a point `p` on an [`AbstractPowerManifold`](@ref) `M` at index `idx`.
"""
function get_component(M::AbstractPowerManifold, p, idx...)
    rep_size = representation_size(M.manifold)
    return _read(M, rep_size, p, idx)
end

function get_coordinates(M::AbstractPowerManifold, p, X, B::DefaultOrthonormalBasis)
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

function get_coordinates!(M::AbstractPowerManifold, Y, p, X, B::DefaultOrthonormalBasis)
    rep_size = representation_size(M.manifold)
    dim = manifold_dimension(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        # TODO: this view is really suboptimal when `dim` can be statically determined
        get_coordinates!(
            M.manifold,
            view(Y, v_iter:(v_iter + dim - 1)),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            B,
        )
        v_iter += dim
    end
    return Y
end
function get_coordinates!(
    M::AbstractPowerManifold,
    Y,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis,<:PowerBasisData},
) where {𝔽}
    rep_size = representation_size(M.manifold)
    dim = manifold_dimension(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        get_coordinates!(
            M.manifold,
            view(Y, v_iter:(v_iter + dim - 1)),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _access_nested(B.data.bases, i),
        )
        v_iter += dim
    end
    return Y
end

get_iterator(M::PowerManifold{𝔽,<:Manifold{𝔽},Tuple{N}}) where {𝔽,N} = Base.OneTo(N)
@generated function get_iterator(
    M::PowerManifold{𝔽,<:Manifold{𝔽},SizeTuple},
) where {𝔽,SizeTuple}
    size_tuple = size_to_tuple(SizeTuple)
    return Base.product(map(Base.OneTo, size_tuple)...)
end

function get_vector!(
    M::PowerManifold,
    Y,
    p,
    X,
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
            X[v_iter:(v_iter + dim - 1)],
            _access_nested(B.data.bases, i),
        )
        v_iter += dim
    end
    return Y
end
function get_vector!(M::AbstractPowerManifold, Y, p, X, B::DefaultOrthonormalBasis)
    dim = manifold_dimension(M.manifold)
    rep_size = representation_size(M.manifold)
    v_iter = 1
    for i in get_iterator(M)
        get_vector!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            X[v_iter:(v_iter + dim - 1)],
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
Base.@propagate_inbounds function Base.getindex(
    p::AbstractArray,
    M::AbstractPowerManifold,
    I::Integer...,
)
    return collect(get_component(M, p, I...))
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
eval(
    quote
        @invoke_maker 1 Manifold injectivity_radius(
            M::AbstractPowerManifold,
            rm::AbstractRetractionMethod,
        )
    end,
)
eval(
    quote
        @invoke_maker 1 Manifold injectivity_radius(
            M::AbstractPowerManifold,
            rm::ExponentialRetraction,
        )
    end,
)

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
    inverse_retract(M::AbstractPowerManifold, p, q, m::InversePowerRetraction)

Compute the inverse retraction from `p` with respect to `q` on an [`AbstractPowerManifold`](@ref) `M`
using an [`InversePowerRetraction`](@ref), which by default encapsulates a inverse retraction
of the base manifold. Then this method is performed elementwise, so the encapsulated inverse
retraction method has to be one that is available on the base [`Manifold`](@ref).
"""
inverse_retract(::AbstractPowerManifold, ::Any...)

function inverse_retract!(M::AbstractPowerManifold, X, p, q, method::InversePowerRetraction)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        inverse_retract!(
            M.manifold,
            _write(M, rep_size, X, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i),
            method.inverse_retraction,
        )
    end
    return X
end
# log and power have to be explicitly stated to avoid an ambiguity in the third case with AbstractPower
@invoke_maker 5 AbstractInverseRetractionMethod inverse_retract!(
    M::AbstractPowerManifold,
    X,
    q,
    p,
    m::LogarithmicInverseRetraction,
)
function inverse_retract!(
    M::AbstractPowerManifold,
    X,
    q,
    p,
    m::AbstractInverseRetractionMethod,
)
    return inverse_retract!(M, X, q, p, InversePowerRetraction(m))
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
function manifold_dimension(M::PowerManifold{𝔽,<:Manifold,TSize}) where {𝔽,TSize}
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

@doc raw"""
    power_dimensions(M::PowerManifold)

return the power of `M`,
"""
function power_dimensions(::PowerManifold{𝔽,<:Manifold,TSize}) where {𝔽,TSize}
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

Base.@propagate_inbounds @inline function _read(
    M::AbstractPowerManifold,
    rep_size::Tuple,
    x::AbstractArray,
    i::Int,
)
    return _read(M, rep_size, x, (i,))
end

Base.@propagate_inbounds @inline function _read(
    ::PowerManifoldNested,
    rep_size::Tuple,
    x::AbstractArray,
    i::Tuple,
)
    return view(x[i...], rep_size_to_colons(rep_size)...)
end

@generated function rep_size_to_colons(rep_size::Tuple)
    N = length(rep_size.parameters)
    return ntuple(i -> Colon(), N)
end


@doc raw"""
    retract(M::AbstractPowerManifold, p, X, method::PowerRetraction)

Compute the retraction from `p` with tangent vector `X` on an [`AbstractPowerManifold`](@ref) `M`
using a [`PowerRetraction`](@ref), which by default encapsulates a retraction of the
base manifold. Then this method is performed elementwise, so the encapsulated retraction
method has to be one that is available on the base [`Manifold`](@ref).
"""
retract(::AbstractPowerManifold, ::Any...)

function retract!(M::AbstractPowerManifold, q, p, X, method::PowerRetraction)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        retract!(
            M.manifold,
            _write(M, rep_size, q, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            method.retraction,
        )
    end
    return q
end
# exp and power have to be explicitly stated, since the third case otherwise introduces and ambiguity.
@invoke_maker 5 AbstractRetractionMethod retract!(
    M::AbstractPowerManifold,
    q,
    p,
    X,
    m::ExponentialRetraction,
)
function retract!(M::AbstractPowerManifold, q, p, X, m::AbstractRetractionMethod)
    return retract!(M, q, p, X, PowerRetraction(m))
end


"""
    set_component!(M::AbstractPowerManifold, q, p, idx...)

Set the component of a point `q` on an [`AbstractPowerManifold`](@ref) `M` at index `idx`
to `p`, which itself is a point on the [`Manifold`](@ref) the power manifold is build on.
"""
function set_component!(M::AbstractPowerManifold, q, p, idx...)
    rep_size = representation_size(M.manifold)
    return copyto!(_write(M, rep_size, q, idx), p)
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

function vector_transport_direction(M::AbstractPowerManifold, p, X, d)
    return vector_transport_direction(M, p, X, d, PowerVectorTransport(ParallelTransport()))
end

function vector_transport_direction!(M::AbstractPowerManifold, Y, p, X, d)
    return vector_transport_direction!(
        M,
        Y,
        p,
        X,
        d,
        PowerVectorTransport(ParallelTransport()),
    )
end
function vector_transport_direction!(
    M::AbstractPowerManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod,
)
    return vector_transport_direction!(M, Y, p, X, d, PowerVectorTransport(m))
end
function vector_transport_direction!(
    M::AbstractPowerManifold,
    Y,
    p,
    X,
    d,
    m::PowerVectorTransport,
)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        vector_transport_direction!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, d, i),
            m.method,
        )
    end
    return Y
end

@doc raw"""
    vector_transport_to(M::AbstractPowerManifold, p, X, q, method::PowerVectorTransport)

Compute the vector transport the tangent vector `X`at `p` to `q` on the
[`PowerManifold`](@ref) `M` using an [`PowerVectorTransport`](@ref) `m`.
This method is performed elementwise, i.e. the method `m` has to be implemented on the
base manifold.
"""
vector_transport_to(::AbstractPowerManifold, ::Any, ::Any, ::Any, ::PowerVectorTransport)
function vector_transport_to(M::AbstractPowerManifold, p, X, q)
    return vector_transport_to(M, p, X, q, PowerVectorTransport(ParallelTransport()))
end

function vector_transport_to!(M::AbstractPowerManifold, Y, p, X, q)
    return vector_transport_to!(M, Y, p, X, q, PowerVectorTransport(ParallelTransport()))
end
function vector_transport_to!(M::AbstractPowerManifold, Y, p, X, q, m::PowerVectorTransport)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        vector_transport_to!(
            M.manifold,
            _write(M, rep_size, Y, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, X, i),
            _read(M, rep_size, q, i),
            m.method,
        )
    end
    return Y
end

function vector_transport_to!(
    M::AbstractPowerManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod,
)
    return vector_transport_to!(M, Y, p, X, q, PowerVectorTransport(m))
end
for VT in ManifoldsBase.VECTOR_TRANSPORT_DISAMBIGUATION
    eval(
        quote
            @invoke_maker 6 AbstractVectorTransportMethod vector_transport_to!(
                M::AbstractPowerManifold,
                Y,
                p,
                X,
                q,
                B::$VT,
            )
        end,
    )
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

@inline function _write(M::PowerManifoldNested, rep_size::Tuple, x::AbstractArray, i::Tuple)
    return view(x[i...], rep_size_to_colons(rep_size)...)
end
