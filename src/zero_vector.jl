"""
    ZeroVector <: AbstractTangentVector

A type to represent the zero vector. This can be used generically
when it is beneficial to not allocate memory in certain places of computations.
"""
struct ZeroVector <: AbstractTangentVector end

# A new constructor
@doc raw"""
    zero_vector(M::AbstractManifold, p, allocate::Bool = true)

Return the tangent vector from the tangent space ``T_p\mathcal M`` at `p` on the
[`AbstractManifold`](@ref) `M`, that represents the zero vector.
Setting `allocate = false` avoids allocation of memory for this vector.
"""
function zero_vector(M::AbstractManifold, p, allocate::Bool)
    return allocate ? zero_vector(M, p) : ZeroVector()
end
zero_vector!(M, Y::ZeroVector, p) = Y

#
#
# From the main ManifoldsBase.jl file

# overwrite functions where tangent vectors appear
# TODO: Is this all we need for allocation?
allocate(::ZeroVector) = ZeroVector()

# TODO: angle? Would currently error for an actual zero vector anyways

are_linearly_independent(M::AbstractManifold, p, ::ZeroVector, Y; kwargs...) = false
are_linearly_independent(M::AbstractManifold, p, X, ::ZeroVector; kwargs...) = false
are_linearly_independent(M::AbstractManifold, p, ::ZeroVector, ::ZeroVector; kwargs...) = false

function check_approx(M::AbstractManifold, p, ::ZeroVector, Y; kwargs...)
    v = norm(M, p, Y)
    return isapprox(v, 0; kwargs...) ? nothing : ApproximatelyError(v, "The tangent vector $Y in the tangent space at $p on $M is not (approximately) the zero vector.")
end
function check_approx(M::AbstractManifold, p, X, ::ZeroVector; kwargs...)
    v = norm(M, p, X)
    return isapprox(v, 0; kwargs...) ? nothing : ApproximatelyError(v, "The tangent vector $X in the tangent space at $p on $M is not (approximately) the zero vector.")
end
check_approx(::AbstractManifold, p, ::ZeroVector, ::ZeroVector; kwargs...) = nothing

# The zero vector is always in the tangent space, so both checks of vector and size
# always report no error (return `nothing`)
check_vector(::AbstractManifold, p, ::ZeroVector) = nothing
check_size(::AbstractManifold, p, ::ZeroVector) = nothing
copy(M::AbstractManifold, p, X::ZeroVector) = ZeroVector()

copyto!(M::AbstractManifold, Y, p, ::ZeroVector) = zero_vector!(M, Y, p)
copyto!(::AbstractManifold, Y::ZeroVector, p, ::ZeroVector) = Y
embed(::AbstractManifold, p, ::ZeroVector) = ZeroVector()
embed!(M::AbstractManifold, Y, p, ::ZeroVector) = zero_vector!(get_embedding(M), Y, embed(M, p))
embed!(::AbstractManifold, Y::ZeroVector, p, ::ZeroVector) = Y

inner(::AbstractManifold, p, ::ZeroVector, Y) = zero(eltype(Y))
inner(::AbstractManifold, p, X, ::ZeroVector) = zero(eltype(X))
inner(::AbstractManifold, p, ::ZeroVector, ::ZeroVector) = zero(eltype(p))

norm(::AbstractManifold, p, ::ZeroVector) = real(zero(eltype(p)))

# TODO: project! ?
# TODO: Riemannian_tensor? That would be 8 cases
# TODO: sectional_curvature?
# TODO: Weingarten!? (also 8 cases)
# With a zero vector both do not make that much sense.
# TODO: For the following in general: Allocating variants necessary? – or are the inplace enough?

# Bases – here it really seems easier to maybe leave allocation to allocate result first?
function get_coordinates!(
        ::AbstractManifold, c, p, ::ZeroVector, ::AbstractBasis = default_basis(M, typeof(p))
    )
    return fill!(c, 0)
end

# exp - TODO: allocating ones necessary?
exp(M::AbstractManifold, p, ::ZeroVector) = copy(M, p)
exp!(M::AbstractManifold, q, p, ::ZeroVector) = copy!(M, q, p)
exp_fused(M::AbstractManifold, p, ::ZeroVector, ::Number = 1) = copy(M, p)
exp_fused!(M::AbstractManifold, q, p, ::ZeroVector, ::Number = 1) = copyto!(M, q, p)

# Metric
change_metric(::AbstractManifold, G::AbstractMetric, p, X::ZeroVector) = ZeroVector()
change_metric!(M::AbstractManifold, Y, G::AbstractMetric, p, X::ZeroVector) = zero_vector!(M, Y, p)
change_representer(::AbstractManifold, G::AbstractMetric, p, X::ZeroVector) = ZeroVector()
change_representer!(M::AbstractManifold, Y, G::AbstractMetric, p, X::ZeroVector) = zero_vector!(M, Y, p)
# TODO: to avoid ambiguities, probably also necessary on MetricManifold ?
# * exp/exp!/exp_fused/exp/fused!/get_coordinates/get_coordinates!inner on MetricManifold
# * parallel transports / project?
# * vector transports
zero_vector(M::MetricManifold, p, allocate::Bool) = zero_vector(M.manifold, p, allocate)

# projections
project(::AbstractManifold, p, ::ZeroVector) = ZeroVector()
project!(M::AbstractManifold, Y, p, ::ZeroVector) = zero_vector!(M, Y, p)

# parallel transport
parallel_transport_direction(M::AbstractManifold, p, ::ZeroVector, d; kwargs...) = ZeroVector()
parallel_transport_direction!(M::AbstractManifold, Y, p, ::ZeroVector, d; kwargs...) = zero_vector!(M, Y, exp(M, p, d))
parallel_transport_to(M::AbstractManifold, p, ::ZeroVector, q; kwargs...) = ZeroVector()
parallel_transport_to!(M::AbstractManifold, Y, p, ::ZeroVector, q; kwargs...) = zero_vector!(M, Y, q)

# retract
function retract(
        M::AbstractManifold, p, ::ZeroVector,
        ::AbstractRetractionMethod = default_retraction_method(M, typeof(p));
        kwargs...
    )
    return copy(M, p)
end
function retract!(
        M::AbstractManifold, q, p, ::ZeroVector,
        ::AbstractRetractionMethod = default_retraction_method(M, typeof(p));
        kwargs...
    )
    return copyto!(M, q, p)
end
function retract_fused(
        M::AbstractManifold, p, ::ZeroVector, ::Number,
        ::AbstractRetractionMethod = default_retraction_method(M, typeof(p)),
    )
    return copy(M, p)
end
function retract_fused!(
        ::AbstractManifold, q, p, ::ZeroVector, ::Number,
        ::AbstractRetractionMethod = default_retraction_method(M, typeof(p));
        kwargs...
    )
    return copyto!(M, q, p)
end

# vector transports
function vector_transport_direction(
        M::AbstractManifold, p, ::ZeroVector, d,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)); kwargs...
    )
    return ZeroVector()
end
function vector_transport_direction!(
        M::AbstractManifold, Y, p, ::ZeroVector, d,
        ::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)); kwargs...
    )
    r = default_retraction_method(M, typeof(p))
    return zero_vector!(M, Y, retract(M, p, d, r))
end
function vector_transport_to(
        M::AbstractManifold, p, ::ZeroVector, q,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    return ZeroVector()
end
function vector_transport_to!(
        M::AbstractManifold, Y, p, ::ZeroVector, q,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    return zero_vector!(M, Y, q)
end

# Vector spaces

Base.:+(::ZeroVector, Y) = Y
Base.:+(X, ::ZeroVector) = X
Base.:+(::ZeroVector, ::ZeroVector) = ZeroVector()

Base.:-(::ZeroVector, Y) = -Y
Base.:-(X, ::ZeroVector) = X
Base.:-(::ZeroVector, ::ZeroVector) = ZeroVector()

Base.:-(::ZeroVector) = ZeroVector()

Base.:*(α::Number, ::ZeroVector) = ZeroVector()
