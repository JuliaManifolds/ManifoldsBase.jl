"""
    ShootingInverseRetraction <: ApproximateInverseRetraction

Approximating the inverse of a retraction using the shooting method.

This implementation of the shooting method works by using another inverse retraction to form
the first guess of the vector. This guess is updated by shooting the vector, guessing the
vector pointing from the shooting result to the target point, and transporting this vector
update back to the initial point on a discretized grid. This process is repeated until the
norm of the vector update falls below a specified tolerance or the maximum number of
iterations is reached.

# Fields
- `retraction::AbstractRetractionMethod`: The retraction whose inverse is approximated.
- `initial_inverse_retraction::AbstractInverseRetractionMethod`: The inverse retraction used
    to form the initial guess of the vector.
- `vector_transport::AbstractVectorTransportMethod`: The vector transport used to transport
    the initial guess of the vector.
- `num_transport_points::Int`: The number of discretization points used for vector
    transport in the shooting method. 2 is the minimum number of points, including just the
    endpoints.
- `tolerance::Real`: The tolerance for the shooting method.
- `max_iterations::Int`: The maximum number of iterations for the shooting method.
"""
struct ShootingInverseRetraction{
    R<:AbstractRetractionMethod,
    IR<:AbstractInverseRetractionMethod,
    VT<:AbstractVectorTransportMethod,
    T<:Real,
} <: ApproximateInverseRetraction
    retraction::R
    initial_inverse_retraction::IR
    vector_transport::VT
    num_transport_points::Int
    tolerance::T
    max_iterations::Int
end

function _inverse_retract!(M::AbstractManifold, X, p, q, m::ShootingInverseRetraction)
    return inverse_retract_shooting!(M, X, p, q, m)
end

"""
    inverse_retract_shooting!(M::AbstractManifold, X, p, q, m::ShootingInverseRetraction)

Approximate the inverse of a retraction using the shooting method.
"""
function inverse_retract_shooting!(
    M::AbstractManifold,
    X,
    p,
    q,
    m::ShootingInverseRetraction,
)
    inverse_retract!(M, X, p, q, m.initial_inverse_retraction)
    gap = norm(M, p, X)
    gap < m.tolerance && return X
    T = real(Base.promote_type(number_eltype(X), number_eltype(p), number_eltype(q)))
    transport_grid = range(one(T), zero(T); length = m.num_transport_points)[2:(end - 1)]
    ΔX = allocate(X)
    ΔXnew = tX = allocate(ΔX)
    retr_tX = allocate_result(M, retract, p, X)
    if m.num_transport_points > 2
        retr_tX_new = allocate_result(M, retract, p, X)
    end
    iteration = 1
    while (gap > m.tolerance) && (iteration < m.max_iterations)
        retract!(M, retr_tX, p, X, m.retraction)
        inverse_retract!(M, ΔX, retr_tX, q, m.initial_inverse_retraction)
        gap = norm(M, retr_tX, ΔX)
        for t in transport_grid
            tX .= t .* X
            retract!(M, retr_tX_new, p, tX, m.retraction)
            vector_transport_to!(M, ΔXnew, retr_tX, ΔX, retr_tX_new, m.vector_transport)
            # realias storage
            retr_tX, retr_tX_new, ΔX, ΔXnew, tX = retr_tX_new, retr_tX, ΔXnew, ΔX, ΔX
        end
        vector_transport_to!(M, ΔXnew, retr_tX, ΔX, p, m.vector_transport)
        X .+= ΔXnew
        iteration += 1
    end
    return X
end
