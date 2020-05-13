
"""
    AbstractVectorTransportMethod

Abstract type for methods for transporting vectors.
"""
abstract type AbstractVectorTransportMethod end

"""
    ParallelTransport <: AbstractVectorTransportMethod

Specify to use parallel transport as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref).
"""
struct ParallelTransport <: AbstractVectorTransportMethod end

"""
    ProjectionTransport <: AbstractVectorTransportMethod

Specify to use projection onto tangent space as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref). See [`project`](@ref) for details.
"""
struct ProjectionTransport <: AbstractVectorTransportMethod end


@doc raw"""
    PoleLadderTransport <: AbstractVectorTransportMethod

Specify to use [`pole_ladder`](@ref) as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref), i.e.

Let $X\in T_p\mathcal M$ be a tangent vector at $p\in\mathcal M$ and $q\in\mathcal M$ the
point to transport to. Then $x = \exp_pX$ is used to call
`y = [`pole_ladder`](@ref)`(M,p,x,q)` and the resulting vector is obtained by computing
$Y = -\log_qy$.

The [`PoleLadderTransport`](@ref) posesses two advantages compared to
[`SchildsLadderTransport`](@ref):
* it is cheaper to evaluate, if you want to transport several vectors, since the
  mid point $c$ then stays unchanged.
* while both methods are exact if the curvature is zero, pole ladder is even exact in
  symmetric Riemannian manifolds[^Pennec2018]

The pole ladder was was proposed in [^LorenziPennec2014]. Its name stems from the fact that
it resembles a pole ladder when applied to a sequence of points usccessively.

# Constructor
````julia
PoleLadderTransport(
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
````
Construct the classical pole ladder that employs exp and log, i.e. as proposed
in[^LorenziPennec2014]. For an even cheaper transport the inner operations can be
changed to an [`AbstractRetractionMethod`](@ref) `retraction` and an
[`AbstractInverseRetractionMethod`](@ref) `inverse_retraction`, respectively.

[^LorenziPennec2014]:
    > Lorenzi, M. and Pennec, X: Efficient parallel transport of deformations in time
    > series of images: From Schild’s to pole ladder.
    > Journal of Mathematical Imaging and Vision (2014), 50(1), pp. 5–17
    > doi [10.1007/s10851-013-0470-3](https://doi.org/10.1007/s10851-013-0470-3),
    > hal: [hal-00870489](https://hal.inria.fr/hal-00870489)
[^Pennec2018]:
    > Pennec, X: Parallel Transport with Pole Ladder: a Third Order Scheme in Affine
    > Connection Spaces which is Exact in Affine Symmetric Spaces.
    > arXiv: [1805.11436](https://arxiv.org/abs/1805.11436)

"""
struct PoleLadderTransport{
    RT<:AbstractRetractionMethod,
    IRT<:AbstractInverseRetractionMethod,
} <: AbstractVectorTransportMethod
    retraction::RT
    inverse_retraction::IRT
    function PoleLadderTransport(
        retraction = ExponentialRetraction(),
        inverse_retraction = LogarithmicInverseRetraction(),
    )
        return new{typeof(retraction),typeof(inverse_retraction)}(
            retraction,
            inverse_retraction,
        )
    end
end

@doc raw"""
    SchildsLadderTransport <: AbstractVectorTransportMethod

Specify to use [`schilds_ladder`](@ref) as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref), i.e.

Let $X\in T_p\mathcal M$ be a tangent vector at $p\in\mathcal M$ and $q\in\mathcal M$ the
point to transport to. Then

````math
P^{\mathrm{S}}_{q\gets p}(X) = \log_q\bigl( \retr_p( 2\retr_p^{-1}c)\bigr),
````
where $c$ is the mid point between $q$ and $d=\exp_pX$.

This method employs the internal function [`schilds_ladder`](@ref)`(M,p,d,q)` that avoids
leaving the manifold.

The name stems from the image of this paralleltogram in a repeated application yielding the
image of a ladder. The approximation was proposed in [^EhlersPiraniSchild1972].

# Constructor
````julia
SchildsLadderTransport(
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
````
Construct the classical Schilds ladder that employs exp and log, i.e. as proposed
in[^EhlersPiraniSchild1972]. For an even cheaper transport these inner operations can be
changed to an [`AbstractRetractionMethod`](@ref) `retraction` and an
[`AbstractInverseRetractionMethod`](@ref) `inverse_retraction`, respectively.

[^EhlersPiraniSchild1972]:
    > Ehlers, J., Pirani, F.A.E., Schild, A.: The geometry of free fall and light
    > propagation. In: O’Raifeartaigh, L. (ed.) General Relativity: Papers in Honour of
    > J. L. Synge, pp. 63–84. Clarendon Press, Oxford (1972).
    > reprint doi: [10.1007/s10714-012-1353-4](https://doi.org/10.1007/s10714-012-1353-4)

"""
struct SchildsLadderTransport{
    RT<:AbstractRetractionMethod,
    IRT<:AbstractInverseRetractionMethod,
} <: AbstractVectorTransportMethod
    retraction::RT
    inverse_retraction::IRT
    function SchildsLadderTransport(
        retraction = ExponentialRetraction(),
        inverse_retraction = LogarithmicInverseRetraction(),
    )
        return new{typeof(retraction),typeof(inverse_retraction)}(
            retraction,
            inverse_retraction,
        )
    end
end

@doc raw"""
    pole_ladder(
        M,
        p,
        d,
        q,
        c = mid_point(M, p, q);
        retraction=ExponentialRetraction(),
        inverse_retraction=LogarithmicInverseRetraction()
    )

Compute an inner step of the pole ladder, that can be used as a [`vector_transport_to`](@ref).
Let $c = \gamma_{p,q}(\frac{1}{2})$ mid point between `p` and `q`, then the pole ladder is
given by

````math
    \operatorname{Pl}(p,d,q) = \operatorname{retr}_d (2\operatorname{retr}_d^{-1}c)
````

Where the classical pole ladder employs $\operatorname{retr}_d=\exp_d$
and $\operatorname{retr}_d^{-1}=\log_d$ but for an even cheaper transport these can be set
to different [`AbstractRetractionMethod`](@ref) and [`AbstractInverseRetractionMethod`](@ref).

When you have $X=log_pd$ and $Y = -\log_q \operatorname{Pl}(p,d,q)$,
you will obtain the [`PoleLadderTransport`](@ref). When performing multiple steps, this
method avoidsd the switching to the tangent space. Keep in mind that after $n$ successive
steps the tangent vector reads $Y_n = (-1)^n\log_q \operatorname{Pl}(p_{n-1},d_{n-1},p_n)$.

It is cheaper to evaluate than [`schilds_ladder`](@ref), sinc if you want to form multiple
ladder steps between `p` and `q`, but with different `d`, there is just one evaluation of a geodesic
each., since the center `c` can be reused.
"""
function pole_ladder(
    M,
    p,
    d,
    q,
    c = mid_point(M, p, q);
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
    return retract(M, d, 2 * inverse_retract(M, d, c, inverse_retraction), retraction)
end
@doc raw"""
    pole_ladder(
        M,
        pl,
        p,
        d,
        q,
        c = mid_point(M, p, q),
        X = allocate_result_type(M, log, d, c);
        retraction = ExponentialRetraction(),
        inverse_retraction = LogarithmicInverseRetraction()
    )

Compute the [`pole_ladder`](@ref), i.e. the result is saved in `pl`.
`X` is used for storing intermediate inverse retraction.
"""
function pole_ladder!(
    M,
    pl,
    p,
    d,
    q,
    c = mid_point(M, p, q),
    X = allocate_result(M, log, d, c);
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
    inverse_retract!(M, X, d, c, inverse_retraction)
    X *= 2
    return retract!(M, pl, d, X, retraction)
end

@doc raw"""
    schilds_ladder(
        M,
        p,
        d,
        q,
        c = mid_point(M, q, d);
        retraction = ExponentialRetraction(),
        inverse_retraction = LogarithmicInverseRetraction()
    )

Perform an inner step of schilds ladder, which can be used as a
[`vector_transport_to`](@ref), see [`SchildsLadderTransport`](@ref).
Let $c = \gamma_{q,d}(\frac{1}{2})$ denote the mid point
on the shortest geodesic connecting $q$ and the point $d$. Then Schild's ladder reads as

````math
\operatorname{Sl}(p,d,q) = \operatorname{retr}_x( 2\operatorname{retr}_x^{-1} c
````

Where the classical Schilds ladder employs $\operatorname{retr}_d=\exp_d$
and $\operatorname{retr}_d^{-1}=\log_d$ but for an even cheaper transport these can be set
to different [`AbstractRetractionMethod`](@ref) and [`AbstractInverseRetractionMethod`](@ref).

In consistency with [`pole_ladder`](@ref) you can change the way the mid point is computed
using the optional parameter `c`, but note that here it's the mid point between `q` and d`.

When you have $X=log_pd$ and $Y = \log_q \operatorname{Sl}(p,d,q)$,
you will obtain the [`PoleLadderTransport`](@ref).
Then the approximation to the transported vector is given by $\log_q\operatorname{Sl}(p,d,q)$.

When performing multiple steps, this method avoidsd the switching to the tangent space.
Hence after $n$ successive steps the tangent vector reads
$Y_n = \log_q \operatorname{Pl}(p_{n-1},d_{n-1},p_n)$.
"""
function schilds_ladder(
    M,
    p,
    d,
    q,
    c = mid_point(M, q, d);
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
    return retract(M, p, 2 * inverse_retract(M, p, c, inverse_retraction), retraction)
end
@doc raw"""
    schilds_ladder!(
        M,
        sl
        p,
        d,
        q,
        c = mid_point(M, q, d),
        X = allocate_result_type(M, log, d, c);
        retraction = ExponentialRetraction(),
        inverse_retraction = LogarithmicInverseRetraction()
    )

Compute [`schilds_ladder`](@ref) and return the value in the parameter `sl`.
If the required mid point `c` was computed before, it can be passed using `c`,
and the allocation of new memory can be avoided providing a tangent vector `X`
for the interims result.
"""
function schilds_ladder!(
    M,
    sl,
    p,
    d,
    q,
    c = mid_point(M, q, d),
    X = allocate_result(M, log, d, c);
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
    inverse_retract!(M, X, p, c, inverse_retraction)
    X *= 2
    return retract!(M, sl, p, X, retraction)
end

"""
    vector_transport_along(M::Manifold, p, X, c)
    vector_transport_along(M::Manifold, p, X, c, method::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`Manifold`](@ref) `M`
along the curve `c` such that `c(0)` is equal to `p` to the point `c(1)` using the `method`,
which defaults to [`ParallelTransport`](@ref).
"""
function vector_transport_along(M::Manifold, p, X, c)
    return vector_transport_along(M, p, X, c, ParallelTransport())
end
function vector_transport_along(M::Manifold, p, X, c, m::AbstractVectorTransportMethod)
    Y = allocate_result(M, vector_transport_along, X, p)
    vector_transport_along!(M, Y, p, X, c, m)
    return Y
end

"""
    vector_transport_along!(M::Manifold, Y, p, X, c)
    vector_transport_along!(M::Manifold, Y, p, X, c, method::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`Manifold`](@ref) `M`
along the curve `c` such that `c(0)` is equal to `p` to the point `c(1)` using the `method`,
which defaults to [`ParallelTransport`](@ref). The result is saved to `Y`.
"""
function vector_transport_along!(M::Manifold, Y, p, X, c)
    return vector_transport_along!(M, Y, p, X, c, ParallelTransport())
end
function vector_transport_along!(
    M::Manifold,
    Y,
    p,
    X,
    c,
    method::AbstractVectorTransportMethod,
)
    return error(manifold_function_not_implemented_message(
        M,
        vector_transport_along!,
        M,
        Y,
        p,
        X,
        c,
        method,
    ))
end
@doc raw"""
    vector_transport_along!(
        M::Manifold,
        Y,
        p,
        X,
        c::AbstractVector,
        method::AbstractVectorTransportMethod
    ) where {T}

Compute the vector transport along a discretized curve `c` using an
[`AbstractVectorTransportMethod`](@ref) `method` succesively along the sampled curve.
"""
function vector_transport_along!(
    M::Manifold,
    Y,
    p,
    X,
    c::AbstractVector,
    method::AbstractVectorTransportMethod,
)
    n = length(c)
    if n == 0
        copyto!(Y, X)
    else
        # we shouldn't assume that vector_transport_to! works when both input and output
        # vectors are the same object
        Y2 = allocate(X)
        vector_transport_to!(M, Y2, p, X, c[1], method)
        for i in 1:(length(c) - 1)
            vector_transport_to!(M, Y, c[i], Y2, c[i + 1], method)
            copyto!(Y2, Y)
        end
    end
    return Y
end
@doc raw"""
    function vector_transport_along!(
        M::Manifold,
        Y,
        p,
        X,
        c::AbstractVector,
        method::PoleLadderTransport
    )

Compute the vector transport along a discretized curve using
[`PoleLadderTransport`](@ref) succesively along the sampled curve.
This method is avoiding additional allocations as well as inner exp/log by performing all
ladder steps on the manifold and only computing one tangent vector in the end.
"""
function vector_transport_along!(
    M::Manifold,
    Y,
    p,
    X,
    c::AbstractVector,
    method::PoleLadderTransport,
)
    clen = length(c)
    if clen == 0
        copyto!(Y, X)
    else
        d = retract(M, p, X, method.retraction)
        m = mid_point(M, p, c[1])
        pole_ladder!(
            M,
            d,
            p,
            d,
            c[1],
            m,
            Y;
            retraction = method.retraction,
            inverse_retraction = method.inverse_retraction,
        )
        for i in 1:(clen - 1)
            # precompute mid point inplace
            ci = c[i]
            cip1 = c[i + 1]
            mid_point!(M, m, ci, cip1)
            # compute new ladder point
            pole_ladder!(
                M,
                d,
                ci,
                d,
                cip1,
                m,
                Y;
                retraction = method.retraction,
                inverse_retraction = method.inverse_retraction,
            )
        end
        inverse_retract!(M, Y, c[clen], d, method.inverse_retraction)
        Y *= (-1)^clen
    end
    return Y
end
@doc raw"""
    vector_transport_along!(
        M::Manifold,
        Y,
        p,
        X,
        c::AbstractVector,
        method::SchildsLadderTransport
    )

Compute the vector transport along a discretized curve using
[`SchildsLadderTransport`](@ref) succesively along the sampled curve.
This method is avoiding additional allocations as well as inner exp/log by performing all
ladder steps on the manifold and only computing one tangent vector in the end.
"""
function vector_transport_along!(
    M::Manifold,
    Y,
    p,
    X,
    c::AbstractVector,
    method::SchildsLadderTransport,
)
    clen = length(c)
    if clen == 0
        copyto!(Y, X)
    else
        d = retract(M, p, X, method.retraction)
        m = mid_point(M, c[1], d)
        schilds_ladder!(
            M,
            d,
            p,
            d,
            c[1],
            m,
            Y;
            retraction = method.retraction,
            inverse_retraction = method.inverse_retraction,
        )
        for i in 1:(clen - 1)
            ci = c[i]
            cip1 = c[i + 1]
            # precompute mid point inplace
            mid_point!(M, m, cip1, d)
            # compute new ladder point
            schilds_ladder!(
                M,
                d,
                ci,
                d,
                cip1,
                m,
                Y;
                retraction = method.retraction,
                inverse_retraction = method.inverse_retraction,
            )
        end
        inverse_retract!(M, Y, c[clen], d, method.inverse_retraction)
    end
    return Y
end


"""
    vector_transport_direction(M::Manifold, p, X, d)
    vector_transport_direction(M::Manifold, p, X, d, method::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`Manifold`](@ref) `M`
in the direction indicated by the tangent vector `d` at `p`. By default, [`exp`](@ref) and
[`vector_transport_to!`](@ref) are used with the `method`, which defaults
to [`ParallelTransport`](@ref).
"""
function vector_transport_direction(M::Manifold, p, X, d)
    return vector_transport_direction(M, p, X, d, ParallelTransport())
end
function vector_transport_direction(
    M::Manifold,
    p,
    X,
    d,
    method::AbstractVectorTransportMethod,
)
    Y = allocate_result(M, vector_transport_direction, X, p, d)
    vector_transport_direction!(M, Y, p, X, d, method)
    return Y
end

"""
    vector_transport_direction!(M::Manifold, Y, p, X, d)
    vector_transport_direction!(M::Manifold, Y, p, X, d, method::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`Manifold`](@ref) `M`
in the direction indicated by the tangent vector `d` at `p`. By default, [`exp`](@ref) and
[`vector_transport_to!`](@ref) are used with the `method`, which defaults
to [`ParallelTransport`](@ref). The result is saved to `Y`.
"""
function vector_transport_direction!(M::Manifold, Y, p, X, d)
    return vector_transport_direction!(M, Y, p, X, d, ParallelTransport())
end
function vector_transport_direction!(
    M::Manifold,
    Y,
    p,
    X,
    d,
    method::AbstractVectorTransportMethod,
)
    y = exp(M, p, d)
    return vector_transport_to!(M, Y, p, X, y, method)
end

"""
    vector_transport_to(M::Manifold, p, X, q)
    vector_transport_to(M::Manifold, p, X, q, method::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`Manifold`](@ref) `M`
along the [`shortest_geodesic`](@ref) to the tangent space at another point `q`.
By default, the [`AbstractVectorTransportMethod`](@ref) `method` is
[`ParallelTransport`](@ref).
"""
function vector_transport_to(M::Manifold, p, X, q)
    return vector_transport_to(M, p, X, q, ParallelTransport())
end
function vector_transport_to(M::Manifold, p, X, q, method::AbstractVectorTransportMethod)
    Y = allocate_result(M, vector_transport_to, X, p, q)
    vector_transport_to!(M, Y, p, X, q, method)
    return Y
end

"""
    vector_transport_to!(M::Manifold, Y, p, X, q)
    vector_transport_to!(M::Manifold, Y, p, X, q, method::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`Manifold`](@ref) `M`
along the [`shortest_geodesic`](@ref) to the tangent space at another point `q`.
By default, the [`AbstractVectorTransportMethod`](@ref) `method` is
[`ParallelTransport`](@ref). The result is saved to `Y`.
"""
function vector_transport_to!(M::Manifold, Y, p, q, X)
    return vector_transport_to!(M, Y, p, q, X, ParallelTransport())
end
"""
    vector_transport_to!(M::Manifold, Y, p, X, q, method::ProjectionTransport)

Transport a vector `X` from the tangent space at `p` on the [`Manifold`](@ref) `M` by
interpreting it as an element of the embedding and then projecting it onto the tangent space
at `q`. This function needs to be separately implemented for each manifold because
projection [`project`](@ref project(M::Manifold, p, X)) may also change vector
representation (if it's different than in the embedding) and it is assumed that the vector
`X` already has the correct representation for `M`.
"""
vector_transport_to!(M::Manifold, Y, p, X, q, m::ProjectionTransport)

@doc raw"""
    vector_transport_to!(M::Manifold, Y, p, X, q, method::PoleLadderTransport)

Perform a vector transport by using [`PoleLadderTransport`](@ref).
"""
function vector_transport_to!(M::Manifold, Y, p, X, q, m::PoleLadderTransport)
    inverse_retract!(
        M,
        Y,
        q,
        pole_ladder(
            M,
            p,
            retract(M, p, X, m.retraction),
            q;
            retraction = m.retraction,
            inverse_retraction = m.inverse_retraction,
        ),
        m.inverse_retraction,
    )
    copyto!(Y, -Y)
    return Y
end

@doc raw"""
    vector_transport_to!(M::Manifold, Y, p, X, q, method::SchildsLadderTransport)

Perform a vector transport by using [`SchildsLadderTransport`](@ref).
"""
function vector_transport_to!(M::Manifold, Y, p, X, q, m::SchildsLadderTransport)
    return inverse_retract!(
        M,
        Y,
        q,
        schilds_ladder(
            M,
            p,
            retract(M, p, X, m.retraction),
            q;
            retraction = m.retraction,
            inverse_retraction = m.inverse_retraction,
        ),
        m.inverse_retraction,
    )
end

function vector_transport_to!(
    M::Manifold,
    Y,
    p,
    X,
    q,
    method::AbstractVectorTransportMethod,
)
    return error(manifold_function_not_implemented_message(
        M,
        vector_transport_to!,
        Y,
        p,
        X,
        q,
        method,
    ))
end

const VECTOR_TRANSPORT_DISAMBIGUATION = [PoleLadderTransport, SchildsLadderTransport]
