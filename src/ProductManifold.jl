@doc raw"""
    ProductManifold{𝔽,TM<:Tuple} <: AbstractManifold{𝔽}

Product manifold $M_1 × M_2 × …  × M_n$ with product geometry.

# Constructor

    ProductManifold(M_1, M_2, ..., M_n)

generates the product manifold $M_1 × M_2 × … × M_n$.
Alternatively, the same manifold can be contructed using the `×` operator:
`M_1 × M_2 × M_3`.
"""
struct ProductManifold{𝔽,TM<:Tuple} <: AbstractDecoratorManifold{𝔽}
    manifolds::TM
end

function ProductManifold(manifolds::AbstractManifold...)
    𝔽 = ManifoldsBase._unify_number_systems((number_system.(manifolds))...)
    return ProductManifold{𝔽,typeof(manifolds)}(manifolds)
end

"""
    getindex(M::ProductManifold, i)
    M[i]

access the `i`th manifold component from the [`ProductManifold`](@ref) `M`.
"""
@inline Base.getindex(M::ProductManifold, i::Integer) = M.manifolds[i]

ProductManifold() = throw(MethodError("No method matching ProductManifold()."))

const PRODUCT_BASIS_LIST = [
    VeeOrthogonalBasis,
    DefaultBasis,
    DefaultBasis{<:Any,TangentSpaceType},
    DefaultOrthogonalBasis,
    DefaultOrthogonalBasis{<:Any,TangentSpaceType},
    DefaultOrthonormalBasis,
    DefaultOrthonormalBasis{<:Any,TangentSpaceType},
    ProjectedOrthonormalBasis{:gram_schmidt,ℝ},
    ProjectedOrthonormalBasis{:svd,ℝ},
]

"""
    ProductBasisData

A typed tuple to store tuples of data of stored/precomputed bases for a [`ProductManifold`](@ref).
"""
struct ProductBasisData{T<:Tuple}
    parts::T
end

const PRODUCT_BASIS_LIST_CACHED = [CachedBasis]

"""
    ProductMetric <: AbstractMetric

A type to represent the product of metrics for a [`ProductManifold`](@ref).
"""
struct ProductMetric <: AbstractMetric end

"""
    ProductRetraction(retractions::AbstractRetractionMethod...)

Product retraction of `retractions`. Works on [`ProductManifold`](@ref).
"""
struct ProductRetraction{TR<:Tuple} <: AbstractRetractionMethod
    retractions::TR
end

function ProductRetraction(retractions::AbstractRetractionMethod...)
    return ProductRetraction{typeof(retractions)}(retractions)
end

"""
    InverseProductRetraction(retractions::AbstractInverseRetractionMethod...)

Product inverse retraction of `inverse retractions`. Works on [`ProductManifold`](@ref).
"""
struct InverseProductRetraction{TR<:Tuple} <: AbstractInverseRetractionMethod
    inverse_retractions::TR
end

function InverseProductRetraction(inverse_retractions::AbstractInverseRetractionMethod...)
    return InverseProductRetraction{typeof(inverse_retractions)}(inverse_retractions)
end

function allocation_promotion_function(M::ProductManifold, f, args::Tuple)
    apfs = map(MM -> allocation_promotion_function(MM, f, args), M.manifolds)
    return reduce(combine_allocation_promotion_functions, apfs)
end

"""
    ProductVectorTransport(methods::AbstractVectorTransportMethod...)

Product vector transport type of `methods`. Works on [`ProductManifold`](@ref).
"""
struct ProductVectorTransport{TR<:Tuple} <: AbstractVectorTransportMethod
    methods::TR
end

function ProductVectorTransport(methods::AbstractVectorTransportMethod...)
    return ProductVectorTransport{typeof(methods)}(methods)
end

"""
    change_representer(M::ProductManifold, ::AbstractMetric, p, X)

Since the metric on a product manifold decouples, the change of a representer can be done elementwise
"""
change_representer(::ProductManifold, ::AbstractMetric, ::Any, ::Any)

function change_representer!(M::ProductManifold, Y, G::AbstractMetric, p, X)
    map(
        (m, y, P, x) -> change_representer!(m, y, G, P, x),
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
    )
    return Y
end

"""
    change_metric(M::ProductManifold, ::AbstractMetric, p, X)

Since the metric on a product manifold decouples, the change of metric can be done elementwise.
"""
change_metric(::ProductManifold, ::AbstractMetric, ::Any, ::Any)

function change_metric!(M::ProductManifold, Y, G::AbstractMetric, p, X)
    map(
        (m, y, P, x) -> change_metric!(m, y, G, P, x),
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
    )
    return Y
end

"""
    check_point(M::ProductManifold, p; kwargs...)

Check whether `p` is a valid point on the [`ProductManifold`](@ref) `M`.
If `p` is not a point on `M` a [`CompositeManifoldError`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/functions.html#ManifoldsBase.CompositeManifoldError).consisting of all error messages of the
components, for which the tests fail is returned.

The tolerance for the last test can be set using the `kwargs...`.
"""
function check_point(M::ProductManifold, p; kwargs...)
    try
        submanifold_components(M, p)
    catch e
        return DomainError("Point $p does not support submanifold_components")
    end
    ts = ziptuples(Tuple(1:length(M.manifolds)), M.manifolds, submanifold_components(M, p))
    e = [(t[1], check_point(t[2:end]...; kwargs...)) for t in ts]
    errors = filter((x) -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

"""
    check_size(M::ProductManifold, p; kwargs...)

Check whether `p` is of valid size on the [`ProductManifold`](@ref) `M`.
If `p` has components of wrong size a [`CompositeManifoldError`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/functions.html#ManifoldsBase.CompositeManifoldError).consisting of all error messages of the
components, for which the tests fail is returned.

The tolerance for the last test can be set using the `kwargs...`.
"""
function check_size(M::ProductManifold, p)
    try
        submanifold_components(M, p)
    catch e
        return DomainError("Point $p does not support submanifold_components")
    end
    ts = ziptuples(Tuple(1:length(M.manifolds)), M.manifolds, submanifold_components(M, p))
    e = [(t[1], check_size(t[2:end]...)) for t in ts]
    errors = filter((x) -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

function check_size(M::ProductManifold, p, X)
    try
        submanifold_components(M, X)
    catch e
        return DomainError("Vector $X does not support submanifold_components")
    end
    ts = ziptuples(
        Tuple(1:length(M.manifolds)),
        M.manifolds,
        submanifold_components(M, p),
        submanifold_components(M, X),
    )
    e = [(t[1], check_size(t[2:end]...)) for t in ts]
    errors = filter(x -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end
"""
    check_vector(M::ProductManifold, p, X; kwargs... )

Check whether `X` is a tangent vector to `p` on the [`ProductManifold`](@ref)
`M`, i.e. all projections to base manifolds must be respective tangent vectors.
If `X` is not a tangent vector to `p` on `M` a [`CompositeManifoldError`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/functions.html#ManifoldsBase.CompositeManifoldError).consisting
of all error messages of the components, for which the tests fail is returned.

The tolerance for the last test can be set using the `kwargs...`.
"""
function check_vector(M::ProductManifold, p, X; kwargs...)
    try
        submanifold_components(M, X)
    catch e
        return DomainError("Vector $X does not support submanifold_components")
    end
    ts = ziptuples(
        Tuple(1:length(M.manifolds)),
        M.manifolds,
        submanifold_components(M, p),
        submanifold_components(M, X),
    )
    e = [(t[1], check_vector(t[2:end]...; kwargs...)) for t in ts]
    errors = filter(x -> !(x[2] === nothing), e)
    cerr = [ComponentManifoldError(er...) for er in errors]
    (length(errors) > 1) && return CompositeManifoldError(cerr)
    (length(errors) == 1) && return cerr[1]
    return nothing
end

function default_retraction_method(M::ProductManifold)
    return ProductRetraction(map(default_retraction_method, M.manifolds)...)
end

function default_inverse_retraction_method(M::ProductManifold)
    return InverseProductRetraction(map(default_inverse_retraction_method, M.manifolds)...)
end

function default_vector_transport_method(M::ProductManifold)
    return ProductVectorTransport(map(default_vector_transport_method, M.manifolds)...)
end

@doc raw"""
    distance(M::ProductManifold, p, q)

Compute the distance between two points `p` and `q` on the [`ProductManifold`](@ref) `M`, which is
the 2-norm of the elementwise distances on the internal manifolds that build `M`.
"""
function distance(M::ProductManifold, p, q)
    return sqrt(
        sum(
            map(
                distance,
                M.manifolds,
                submanifold_components(M, p),
                submanifold_components(M, q),
            ) .^ 2,
        ),
    )
end

@doc raw"""
    exp(M::ProductManifold, p, X)

compute the exponential map from `p` in the direction of `X` on the [`ProductManifold`](@ref) `M`,
which is the elementwise exponential map on the internal manifolds that build `M`.
"""
exp(::ProductManifold, ::Any...)

function exp!(M::ProductManifold, q, p, X)
    map(
        exp!,
        M.manifolds,
        submanifold_components(M, q),
        submanifold_components(M, p),
        submanifold_components(M, X),
    )
    return q
end
function exp!(M::ProductManifold, q, p, X, t::Number)
    map(
        (N, qc, pc, Xc) -> exp!(N, qc, pc, Xc, t),
        M.manifolds,
        submanifold_components(M, q),
        submanifold_components(M, p),
        submanifold_components(M, X),
    )
    return q
end

function get_basis(M::ProductManifold, p, B::AbstractBasis)
    parts = map(t -> get_basis(t..., B), ziptuples(M.manifolds, submanifold_components(p)))
    return CachedBasis(B, ProductBasisData(parts))
end
function get_basis(M::ProductManifold, p, B::CachedBasis)
    return invoke(get_basis, Tuple{AbstractManifold,Any,CachedBasis}, M, p, B)
end
function get_basis(M::ProductManifold, p, B::DiagonalizingOrthonormalBasis)
    vs = map(
        ziptuples(
            M.manifolds,
            submanifold_components(p),
            submanifold_components(B.frame_direction),
        ),
    ) do t
        return get_basis(t[1], t[2], DiagonalizingOrthonormalBasis(t[3]))
    end
    return CachedBasis(B, ProductBasisData(vs))
end

"""
    get_component(M::ProductManifold, p, i)

Get the `i`th component of a point `p` on a [`ProductManifold`](@ref) `M`.
"""
function get_component(M::ProductManifold, p, i)
    return submanifold_component(M, p, i)
end

function get_coordinates(M::ProductManifold, p, X, B::AbstractBasis)
    reps = map(
        t -> get_coordinates(t..., B),
        ziptuples(M.manifolds, submanifold_components(M, p), submanifold_components(M, X)),
    )
    return vcat(reps...)
end
function get_coordinates(
    M::ProductManifold,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:ProductBasisData},
) where {𝔽}
    reps = map(
        get_coordinates,
        M.manifolds,
        submanifold_components(M, p),
        submanifold_components(M, X),
        B.data.parts,
    )
    return vcat(reps...)
end

function get_coordinates!(M::ProductManifold, Xⁱ, p, X, B::AbstractBasis)
    dim = manifold_dimension(M)
    @assert length(Xⁱ) == dim
    i = one(dim)
    ts = ziptuples(M.manifolds, submanifold_components(M, p), submanifold_components(M, X))
    for t in ts
        SM = first(t)
        dim = manifold_dimension(SM)
        tXⁱ = @inbounds view(Xⁱ, i:(i + dim - 1))
        get_coordinates!(SM, tXⁱ, Base.tail(t)..., B)
        i += dim
    end
    return Xⁱ
end
function get_coordinates!(
    M::ProductManifold,
    Xⁱ,
    p,
    X,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:ProductBasisData},
) where {𝔽}
    dim = manifold_dimension(M)
    @assert length(Xⁱ) == dim
    i = one(dim)
    ts = ziptuples(
        M.manifolds,
        submanifold_components(M, p),
        submanifold_components(M, X),
        B.data.parts,
    )
    for t in ts
        SM = first(t)
        dim = manifold_dimension(SM)
        tXⁱ = @inbounds view(Xⁱ, i:(i + dim - 1))
        get_coordinates!(SM, tXⁱ, Base.tail(t)...)
        i += dim
    end
    return Xⁱ
end

function _get_dim_ranges(dims::NTuple{N,Any}) where {N}
    dims_acc = accumulate(+, (1, dims...))
    return ntuple(i -> (dims_acc[i]:(dims_acc[i] + dims[i] - 1)), Val(N))
end

function get_vector!(M::ProductManifold, X, p, Xⁱ, B::AbstractBasis)
    dims = map(manifold_dimension, M.manifolds)
    @assert length(Xⁱ) == sum(dims)
    dim_ranges = _get_dim_ranges(dims)
    tXⁱ = map(dr -> (@inbounds view(Xⁱ, dr)), dim_ranges)
    ts = ziptuples(
        M.manifolds,
        submanifold_components(M, X),
        submanifold_components(M, p),
        tXⁱ,
    )
    map(ts) do t
        return get_vector!(t..., B)
    end
    return X
end
function get_vector!(
    M::ProductManifold,
    X,
    p,
    Xⁱ,
    B::CachedBasis{𝔽,<:AbstractBasis{𝔽},<:ProductBasisData},
) where {𝔽}
    dims = map(manifold_dimension, M.manifolds)
    @assert length(Xⁱ) == sum(dims)
    dim_ranges = _get_dim_ranges(dims)
    tXⁱ = map(dr -> (@inbounds view(Xⁱ, dr)), dim_ranges)
    ts = ziptuples(
        M.manifolds,
        submanifold_components(M, X),
        submanifold_components(M, p),
        tXⁱ,
        B.data.parts,
    )
    map(ts) do t
        return get_vector!(t...)
    end
    return X
end

@doc raw"""
    injectivity_radius(M::ProductManifold)
    injectivity_radius(M::ProductManifold, x)

Compute the injectivity radius on the [`ProductManifold`](@ref), which is the
minimum of the factor manifolds.
"""
injectivity_radius(::ProductManifold, ::Any...)
function injectivity_radius(M::ProductManifold, p)
    return min(map(injectivity_radius, M.manifolds, submanifold_components(M, p))...)
end
function injectivity_radius(M::ProductManifold, p, m::AbstractRetractionMethod)
    return min(
        map(
            (lM, lp) -> injectivity_radius(lM, lp, m),
            M.manifolds,
            submanifold_components(M, p),
        )...,
    )
end
function injectivity_radius(M::ProductManifold, p, m::ProductRetraction)
    return min(
        map(
            (lM, lp, lm) -> injectivity_radius(lM, lp, lm),
            M.manifolds,
            submanifold_components(M, p),
            m.retractions,
        )...,
    )
end
injectivity_radius(M::ProductManifold) = min(map(injectivity_radius, M.manifolds)...)
function injectivity_radius(M::ProductManifold, m::AbstractRetractionMethod)
    return min(map(manif -> injectivity_radius(manif, m), M.manifolds)...)
end
function injectivity_radius(M::ProductManifold, m::ProductRetraction)
    return min(map((lM, lm) -> injectivity_radius(lM, lm), M.manifolds, m.retractions)...)
end

@doc raw"""
    inner(M::ProductManifold, p, X, Y)

compute the inner product of two tangent vectors `X`, `Y` from the tangent space
at `p` on the [`ProductManifold`](@ref) `M`, which is just the sum of the
internal manifolds that build `M`.
"""
function inner(M::ProductManifold, p, X, Y)
    subproducts = map(
        inner,
        M.manifolds,
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, Y),
    )
    return sum(subproducts)
end

@doc raw"""
    inverse_retract(M::ProductManifold, p, q, m::InverseProductRetraction)

Compute the inverse retraction from `p` with respect to `q` on the [`ProductManifold`](@ref)
`M` using an [`InverseProductRetraction`](@ref), which by default encapsulates a inverse
retraction for each manifold of the product. Then this method is performed elementwise,
so the encapsulated inverse retraction methods have to be available per factor.
"""
inverse_retract(::ProductManifold, ::Any, ::Any, ::Any, ::InverseProductRetraction)

function inverse_retract!(M::ProductManifold, Y, p, q, method::InverseProductRetraction)
    map(
        inverse_retract!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, q),
        method.inverse_retractions,
    )
    return Y
end

function _isapprox(M::ProductManifold, p, q; kwargs...)
    return all(
        t -> isapprox(t...; kwargs...),
        ziptuples(M.manifolds, submanifold_components(M, p), submanifold_components(M, q)),
    )
end
function _isapprox(M::ProductManifold, p, X, Y; kwargs...)
    return all(
        t -> isapprox(t...; kwargs...),
        ziptuples(
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
            submanifold_components(M, Y),
        ),
    )
end

"""
    is_flat(::ProductManifold)

Return true if and only if all component manifolds of [`ProductManifold`](@ref) `M` are flat.
"""
function is_flat(M::ProductManifold)
    return all(is_flat, M.manifolds)
end

@doc raw"""
    log(M::ProductManifold, p, q)

Compute the logarithmic map from `p` to `q` on the [`ProductManifold`](@ref) `M`,
which can be computed using the logarithmic maps of the manifolds elementwise.
"""
log(::ProductManifold, ::Any...)

function log!(M::ProductManifold, X, p, q)
    map(
        log!,
        M.manifolds,
        submanifold_components(M, X),
        submanifold_components(M, p),
        submanifold_components(M, q),
    )
    return X
end

@doc raw"""
    manifold_dimension(M::ProductManifold)

Return the manifold dimension of the [`ProductManifold`](@ref), which is the sum of the
manifold dimensions the product is made of.
"""
manifold_dimension(M::ProductManifold) = mapreduce(manifold_dimension, +, M.manifolds)

function mid_point!(M::ProductManifold, q, p1, p2)
    map(
        mid_point!,
        M.manifolds,
        submanifold_components(M, q),
        submanifold_components(M, p1),
        submanifold_components(M, p2),
    )
    return q
end

@doc raw"""
    norm(M::ProductManifold, p, X)

Compute the norm of `X` from the tangent space of `p` on the [`ProductManifold`](@ref),
i.e. from the element wise norms the 2-norm is computed.
"""
function LinearAlgebra.norm(M::ProductManifold, p, X)
    norms_squared = (
        map(
            norm,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, X),
        ) .^ 2
    )
    return sqrt(sum(norms_squared))
end

"""
    number_of_components(M::ProductManifold{<:NTuple{N,Any}}) where {N}

Calculate the number of manifolds multiplied in the given [`ProductManifold`](@ref) `M`.
"""
number_of_components(::ProductManifold{𝔽,<:NTuple{N,Any}}) where {𝔽,N} = N

function parallel_transport_direction!(M::ProductManifold, Y, p, X, d)
    map(
        parallel_transport_direction!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, d),
    )
    return Y
end

function parallel_transport_to!(M::ProductManifold, Y, p, X, q)
    map(
        parallel_transport_to!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, q),
    )
    return Y
end

function project!(M::ProductManifold, q, p)
    map(project!, M.manifolds, submanifold_components(M, q), submanifold_components(M, p))
    return q
end

function project!(M::ProductManifold, Y, p, X)
    map(
        project!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
    )
    return Y
end

function Random.rand!(
    M::ProductManifold,
    pX;
    vector_at = nothing,
    parts_kwargs = map(_ -> (;), M.manifolds),
)
    return rand!(
        Random.default_rng(),
        M,
        pX;
        vector_at = vector_at,
        parts_kwargs = parts_kwargs,
    )
end
function Random.rand!(
    rng::AbstractRNG,
    M::ProductManifold,
    pX;
    vector_at = nothing,
    parts_kwargs = map(_ -> (;), M.manifolds),
)
    if vector_at === nothing
        map(
            (N, q, kwargs) -> rand!(rng, N, q; kwargs...),
            M.manifolds,
            submanifold_components(M, pX),
            parts_kwargs,
        )
    else
        map(
            (N, X, p, kwargs) -> rand!(rng, N, X; vector_at = p, kwargs...),
            M.manifolds,
            submanifold_components(M, pX),
            submanifold_components(M, vector_at),
            parts_kwargs,
        )
    end
    return pX
end

@doc raw"""
    retract(M::ProductManifold, p, X, m::ProductRetraction)

Compute the retraction from `p` with tangent vector `X` on the [`ProductManifold`](@ref) `M`
using an [`ProductRetraction`](@ref), which by default encapsulates retractions of the
base manifolds. Then this method is performed elementwise, so the encapsulated retractions
method has to be one that is available on the manifolds.
"""
retract(::ProductManifold, ::Any...)

function _retract!(M::ProductManifold, q, p, X, t::Number, method::ProductRetraction)
    map(
        (N, qc, pc, Xc, rm) -> retract!(N, qc, pc, Xc, t, rm),
        M.manifolds,
        submanifold_components(M, q),
        submanifold_components(M, p),
        submanifold_components(M, X),
        method.retractions,
    )
    return q
end

function representation_size(M::ProductManifold)
    return (mapreduce(m -> prod(representation_size(m)), +, M.manifolds),)
end

@doc raw"""
    riemann_tensor(M::ProductManifold, p, X, Y, Z)

Compute the Riemann tensor at point from `p` with tangent vectors `X`, `Y` and `Z` on
the [`ProductManifold`](@ref) `M`.
"""
riemann_tensor(M::ProductManifold, p, X, Y, X)

function riemann_tensor!(M::ProductManifold, Xresult, p, X, Y, Z)
    map(
        riemann_tensor!,
        M.manifolds,
        submanifold_components(M, Xresult),
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, Y),
        submanifold_components(M, Z),
    )
    return Xresult
end

"""
    select_from_tuple(t::NTuple{N, Any}, positions::Val{P})

Selects elements of tuple `t` at positions specified by the second argument.
For example `select_from_tuple(("a", "b", "c"), Val((3, 1, 1)))` returns
`("c", "a", "a")`.
"""
@generated function select_from_tuple(t::NTuple{N,Any}, positions::Val{P}) where {N,P}
    for k in P
        (k < 0 || k > N) && error("positions must be between 1 and $N")
    end
    return Expr(:tuple, [Expr(:ref, :t, k) for k in P]...)
end

"""
    set_component!(M::ProductManifold, q, p, i)

Set the `i`th component of a point `q` on a [`ProductManifold`](@ref) `M` to `p`, where `p` is a point on the [`AbstractManifold`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/types.html#ManifoldsBase.AbstractManifold)  this factor of the product manifold consists of.
"""
function set_component!(M::ProductManifold, q, p, i)
    return copyto!(submanifold_component(M, q, i), p)
end

function _show_submanifold(io::IO, M::AbstractManifold; pre = "")
    sx = sprint(show, "text/plain", M, context = io, sizehint = 0)
    if occursin('\n', sx)
        sx = sprint(show, M, context = io, sizehint = 0)
    end
    sx = replace(sx, '\n' => "\n$(pre)")
    print(io, pre, sx)
    return nothing
end

function _show_submanifold_range(io::IO, Ms, range; pre = "")
    for i in range
        M = Ms[i]
        print(io, '\n')
        _show_submanifold(io, M; pre = pre)
    end
    return nothing
end

function _show_product_manifold_no_header(io::IO, M)
    n = length(M.manifolds)
    sz = displaysize(io)
    screen_height, screen_width = sz[1] - 4, sz[2]
    half_height = div(screen_height, 2)
    inds = 1:n
    pre = " "
    if n > screen_height
        inds = [1:half_height; (n - div(screen_height - 1, 2) + 1):n]
    end
    if n ≤ screen_height
        _show_submanifold_range(io, M.manifolds, 1:n; pre = pre)
    else
        _show_submanifold_range(io, M.manifolds, 1:half_height; pre = pre)
        print(io, "\n$(pre)⋮")
        _show_submanifold_range(
            io,
            M.manifolds,
            (n - div(screen_height - 1, 2) + 1):n;
            pre = pre,
        )
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", M::ProductManifold)
    n = length(M.manifolds)
    print(io, "ProductManifold with $(n) submanifold$(n == 1 ? "" : "s"):")
    return _show_product_manifold_no_header(io, M)
end

function Base.show(io::IO, M::ProductManifold)
    return print(io, "ProductManifold(", join(M.manifolds, ", "), ")")
end

function Base.show(
    io::IO,
    mime::MIME"text/plain",
    B::CachedBasis{𝔽,T,D},
) where {𝔽,T<:AbstractBasis{𝔽},D<:ProductBasisData}
    println(io, "$(T) for a product manifold")
    for (i, cb) in enumerate(B.data.parts)
        println(io, "Basis for component $i:")
        show(io, mime, cb)
        println(io)
    end
    return nothing
end

"""
    submanifold(M::ProductManifold, i::Integer)

Extract the `i`th factor of the product manifold `M`.
"""
submanifold(M::ProductManifold, i::Integer) = M.manifolds[i]

"""
    submanifold(M::ProductManifold, i::Val)
    submanifold(M::ProductManifold, i::AbstractVector)

Extract the factor of the product manifold `M` indicated by indices in `i`.
For example, for `i` equal to `Val((1, 3))` the product manifold constructed
from the first and the third factor is returned.

The version with `AbstractVector` is not type-stable, for better preformance use `Val`.
"""
function submanifold(M::ProductManifold, i::Val)
    return ProductManifold(select_from_tuple(M.manifolds, i)...)
end
submanifold(M::ProductManifold, i::AbstractVector) = submanifold(M, Val(tuple(i...)))

function vector_transport_direction!(
    M::ProductManifold,
    Y,
    p,
    X,
    d,
    m::ProductVectorTransport,
)
    map(
        vector_transport_direction!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, d),
        m.methods,
    )
    return Y
end

@doc raw"""
    vector_transport_to(M::ProductManifold, p, X, q, m::ProductVectorTransport)

Compute the vector transport the tangent vector `X`at `p` to `q` on the
[`ProductManifold`](@ref) `M` using an [`ProductVectorTransport`](@ref) `m`.
This method is performed elementwise, i.e. the method `m` has to be implemented on the
base manifold.
"""
vector_transport_to(::ProductManifold, ::Any, ::Any, ::Any, ::ProductVectorTransport)

function vector_transport_to!(M::ProductManifold, Y, p, X, q, m::ProductVectorTransport)
    map(
        vector_transport_to!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, q),
        m.methods,
    )
    return Y
end
function vector_transport_to!(M::ProductManifold, Y, p, X, q, m::ParallelTransport)
    map(
        (iM, iY, ip, iX, id) -> vector_transport_to!(iM, iY, ip, iX, id, m),
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, q),
    ),
    return Y
end

@doc raw"""
    Y = Weingarten(M::ProductManifold, p, X, V)
    Weingarten!(M::ProductManifold, Y, p, X, V)

Since the metric decouples, also the computation of the Weingarten map
``\mathcal W_p`` can be computed elementwise on the single elements of the [`ProductManifold`](@ref) `M`.
"""
Weingarten(::ProductManifold, p, X, V)

function Weingarten!(M::ProductManifold, Y, p, X, V)
    map(
        Weingarten!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, X),
        submanifold_components(M, V),
    )
    return Y
end

function zero_vector!(M::ProductManifold, X, p)
    map(
        zero_vector!,
        M.manifolds,
        submanifold_components(M, X),
        submanifold_components(M, p),
    )
    return X
end




@doc raw"""
    submanifold_component(M::AbstractManifold, p, i::Integer)
    submanifold_component(M::AbstractManifold, p, ::Val{i}) where {i}
    submanifold_component(p, i::Integer)
    submanifold_component(p, ::Val{i}) where {i}

Project the product array `p` on `M` to its `i`th component. A new array is returned.
"""
submanifold_component(::Any...)
@inline function submanifold_component(M::AbstractManifold, p, i::Integer)
    return submanifold_component(M, p, Val(i))
end
@inline submanifold_component(M::AbstractManifold, p, i::Val) = submanifold_component(p, i)
@inline submanifold_component(p, i::Integer) = submanifold_component(p, Val(i))

@doc raw"""
    submanifold_components(M::AbstractManifold, p)
    submanifold_components(p)

Get the projected components of `p` on the submanifolds of `M`. The components are returned in a Tuple.
"""
submanifold_components(::Any...)
@inline submanifold_components(::AbstractManifold, p) = submanifold_components(p)

"""
    ziptuples(a, b[, c[, d[, e]]])

Zips tuples `a`, `b`, and remaining in a fast, type-stable way. If they have different
lengths, the result is trimmed to the length of the shorter tuple.
"""
@generated function ziptuples(a::NTuple{N,Any}, b::NTuple{M,Any}) where {N,M}
    ex = Expr(:tuple)
    for i in 1:min(N, M)
        push!(ex.args, :((a[$i], b[$i])))
    end
    return ex
end
@generated function ziptuples(
    a::NTuple{N,Any},
    b::NTuple{M,Any},
    c::NTuple{L,Any},
) where {N,M,L}
    ex = Expr(:tuple)
    for i in 1:min(N, M, L)
        push!(ex.args, :((a[$i], b[$i], c[$i])))
    end
    return ex
end
@generated function ziptuples(
    a::NTuple{N,Any},
    b::NTuple{M,Any},
    c::NTuple{L,Any},
    d::NTuple{K,Any},
) where {N,M,L,K}
    ex = Expr(:tuple)
    for i in 1:min(N, M, L, K)
        push!(ex.args, :((a[$i], b[$i], c[$i], d[$i])))
    end
    return ex
end
@generated function ziptuples(
    a::NTuple{N,Any},
    b::NTuple{M,Any},
    c::NTuple{L,Any},
    d::NTuple{K,Any},
    e::NTuple{J,Any},
) where {N,M,L,K,J}
    ex = Expr(:tuple)
    for i in 1:min(N, M, L, K, J)
        push!(ex.args, :((a[$i], b[$i], c[$i], d[$i], e[$i])))
    end
    return ex
end
