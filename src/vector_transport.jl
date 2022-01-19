
"""
    AbstractVectorTransportMethod

Abstract type for methods for transporting vectors. Such vector transports are not
necessarily linear.

# See also

[`AbstractLinearVectorTransportMethod`](@ref)
"""
abstract type AbstractVectorTransportMethod end

"""
    AbstractLinearVectorTransportMethod <: AbstractVectorTransportMethod

Abstract type for linear methods for transporting vectors, that is transport of a linear
combination of vectors is a linear combination of transported vectors.
"""
abstract type AbstractLinearVectorTransportMethod <: AbstractVectorTransportMethod end

@doc raw"""
    DifferentiatedRetractionVectorTransport{R<:AbstractRetractionMethod} <:
        AbstractVectorTransportMethod

A type to specify a vector transport that is given by differentiating a retraction.
This can be introduced in two ways. Let ``\mathcal M`` be a Riemannian manifold,
``p∈\mathcal M`` a point, and ``X,Y∈ T_p\mathcal M`` denote two tangent vectors at ``p``.

Given a retraction (cf. [`AbstractRetractionMethod`](@ref)) ``\operatorname{retr}``,
the vector transport of `X` in direction `Y` (cf. [`vector_transport_direction`](@ref))
by differentiation this retraction, is given by

```math
\mathcal T^{\operatorname{retr}}_{p,Y}X
= D_Y\operatorname{retr}_p(Y)[X]
= \frac{\mathrm{d}}{\mathrm{d}t}\operatorname{retr}_p(Y+tX)\Bigr|_{t=0}.
```
see [^AbsilMahonySepulchre2008], Section 8.1.2 for more details.

This can be phrased similarly as a [`vector_transport_to`](@ref) by introducing
``q=\operatorname{retr}_pX`` and defining

```math
\mathcal T^{\operatorname{retr}}_{q \gets p}X = \mathcal T^{\operatorname{retr}}_{p,Y}X
```

which in practice usually requires the [`inverse_retract`](@ref) to exists in order to
compute ``Y = \operatorname{retr}_p^{-1}q``.

# Constructor

    DifferentiatedRetractionVectorTransport(m::AbstractRetractionMethod)
"""
struct DifferentiatedRetractionVectorTransport{R<:AbstractRetractionMethod} <:
       AbstractLinearVectorTransportMethod
    retraction::R
end

@doc raw"""
    ParallelTransport <: AbstractVectorTransportMethod

Compute the vector transport by parallel transport, see
[`parallel_transport_to`](@ref)
"""
struct ParallelTransport <: AbstractLinearVectorTransportMethod end

"""
    ProjectionTransport <: AbstractVectorTransportMethod

Specify to use projection onto tangent space as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref). See [`project`](@ref) for details.
"""
struct ProjectionTransport <: AbstractLinearVectorTransportMethod end


@doc raw"""
    PoleLadderTransport <: AbstractVectorTransportMethod

Specify to use [`pole_ladder`](@ref) as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref), i.e.

Let $X∈ T_p\mathcal M$ be a tangent vector at $p∈\mathcal M$ and $q∈\mathcal M$ the
point to transport to. Then $x = \exp_pX$ is used to call
`y = `[`pole_ladder`](@ref)`(M, p, x, q)` and the resulting vector is obtained by computing
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
} <: AbstractLinearVectorTransportMethod
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
    ScaledVectorTransport{T} <: AbstractVectorTransportMethod

Introduce a scaled variant of any [`AbstractVectorTransportMethod`](@ref) `T`,
as introduced in [^SatoIwai2013] for some ``X∈ T_p\mathcal M`` as

```math
    \mathcal T^{\mathrm{S}}(X) = \frac{\lVert X\rVert_p}{\lVert \mathcal T(X)\rVert_q}\mathcal T(X).
```

Note that the resulting point `q` has to be known, i.e. for [`vector_transport_direction`](@ref)
the curve or more precisely its end point has to be known (via an exponential map or a
retraction). Therefore a default implementation is only provided for the [`vector_transport_to`](@ref)

# Constructor

    ScaledVectorTransport(m::AbstractVectorTransportMethod)

[^SatoIwai2013]:
    > Sato, H., Iwai, T.: _A new, globally convergent Riemannian conjugate gradient method_,
    > Optimization, 2013, Volume 64(4), pp. 1011–1031.
    > doi: [10.1080/02331934.2013.836650](https://doi.org/10.1080/02331934.2013.836650),
    > arXiv: [1302.0125](https://arxiv.org/abs/1302.0125).
"""
struct ScaledVectorTransport{T<:AbstractVectorTransportMethod} <:
       AbstractVectorTransportMethod
    method::T
end

@doc raw"""
    SchildsLadderTransport <: AbstractVectorTransportMethod

Specify to use [`schilds_ladder`](@ref) as vector transport method within
[`vector_transport_to`](@ref), [`vector_transport_direction`](@ref), or
[`vector_transport_along`](@ref), i.e.

Let $X∈ T_p\mathcal M$ be a tangent vector at $p∈\mathcal M$ and $q∈\mathcal M$ the
point to transport to. Then

````math
P^{\mathrm{S}}_{q\gets p}(X) =
    \log_q\bigl( \operatorname{retr}_p ( 2\operatorname{retr}_p^{-1}c ) \bigr),
````

where $c$ is the mid point between $q$ and $d=\exp_pX$.

This method employs the internal function [`schilds_ladder`](@ref)`(M, p, d, q)` that avoids
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
} <: AbstractLinearVectorTransportMethod
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

"""
    default_vector_transport_method(M::AbstractManifold)

The [`AbstractVectorTransportMethod`](@ref) that is used when calling
[`vector_transport_along`](@ref), [`vector_transport_to`](@ref), or
[`vector_transport_direction`](@ref) without specifying the vector transport method.
By default, this is [`ParallelTransport`](@ref).
"""
function default_vector_transport_method(::AbstractManifold)
    return ParallelTransport()
end

@doc raw"""
    pole_ladder(
        M,
        p,
        d,
        q,
        c = mid_point(M, p, q);
        retraction=default_retraction_method(M),
        inverse_retraction=default_inverse_retraction_method(M)
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
method avoids the switching to the tangent space. Keep in mind that after $n$ successive
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
    retraction = default_retraction_method(M),
    inverse_retraction = default_inverse_retraction_method(M),
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
        retraction = default_retraction_method(M),
        inverse_retraction = default_inverse_retraction_method(M),
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
    retraction = default_retraction_method(M),
    inverse_retraction = default_inverse_retraction_method(M),
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
        retraction = default_retraction_method(M),
        inverse_retraction = default_inverse_retraction_method(M),
    )

Perform an inner step of schilds ladder, which can be used as a
[`vector_transport_to`](@ref), see [`SchildsLadderTransport`](@ref).
Let $c = \gamma_{q,d}(\frac{1}{2})$ denote the mid point
on the shortest geodesic connecting $q$ and the point $d$. Then Schild's ladder reads as

````math
\operatorname{Sl}(p,d,q) = \operatorname{retr}_x( 2\operatorname{retr}_p^{-1} c)
````

Where the classical Schilds ladder employs $\operatorname{retr}_d=\exp_d$
and $\operatorname{retr}_d^{-1}=\log_d$ but for an even cheaper transport these can be set
to different [`AbstractRetractionMethod`](@ref) and [`AbstractInverseRetractionMethod`](@ref).

In consistency with [`pole_ladder`](@ref) you can change the way the mid point is computed
using the optional parameter `c`, but note that here it's the mid point between `q` and `d`.

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
    retraction = default_retraction_method(M),
    inverse_retraction = default_inverse_retraction_method(M),
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
        retraction = default_retraction_method(M),
        inverse_retraction = default_inverse_retraction_method(M),
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
    retraction = default_retraction_method(M),
    inverse_retraction = default_inverse_retraction_method(M),
)
    inverse_retract!(M, X, p, c, inverse_retraction)
    X *= 2
    return retract!(M, sl, p, X, retraction)
end

function show(io::IO, ::ParallelTransport)
    return print(io, "ParallelTransport()")
end
function show(io::IO, m::ScaledVectorTransport)
    return print(io, "ScaledVectorTransport($(m.method))")
end

"""
    vector_transport_along(M::AbstractManifold, p, X, c)
    vector_transport_along(M::AbstractManifold, p, X, c, m::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`AbstractManifold`](@ref) `M`
along the curve represented by `c` using the `method`, which defaults to
[`default_vector_transport_method`](@ref)`(M)`.
"""
function vector_transport_along(
    M::AbstractManifold,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return _vector_transport_along(M, p, X, c, m)
end
function _vector_transport_along(M::AbstractManifold, p, X, c, ::ParallelTransport)
    return parallel_transport_along(M, p, X, c)
end
function _vector_transport_along(
    M::AbstractManifold,
    p,
    X,
    c,
    m::DifferentiatedRetractionVectorTransport,
)
    return vector_transport_along_diff(M, p, X, c, m.retraction)
end
function vector_transport_along_diff(
    M::AbstractManifold,
    p,
    X,
    c,
    m::AbstractRetractionMethod,
)
    Y = allocate_result(M, vector_transport_along, X, p)
    return vector_transport_along_diff!(M, Y, p, X, c, m)
end
function _vector_transport_along(M::AbstractManifold, p, X, c, ::ProjectionTransport)
    return vector_transport_along_project(M, p, X, c)
end
function vector_transport_along_project(M::AbstractManifold, p, X, c)
    Y = allocate_result(M, vector_transport_along, X, p)
    return vector_transport_along_project!(M, Y, p, X, c)
end
function _vector_transport_along(M::AbstractManifold, p, X, c, m::PoleLadderTransport)
    Y = allocate_result(M, vector_transport_along, X, p)
    return _vector_transport_along!(M, Y, p, X, c, m)
end
function _vector_transport_along(M::AbstractManifold, p, X, c, m::SchildsLadderTransport)
    Y = allocate_result(M, vector_transport_along, X, p)
    return _vector_transport_along!(M, Y, p, X, c, m)
end

"""
    vector_transport_along!(M::AbstractManifold, Y, p, X, c)
    vector_transport_along!(M::AbstractManifold, Y, p, X, c, m::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`AbstractManifold`](@ref) `M`
along the curve represented by `c` using the `method`, which defaults to
[`default_vector_transport_method`](@ref)`(M)`. The result is saved to `Y`.
"""
function vector_transport_along!(
    M::AbstractManifold,
    Y,
    p,
    X,
    c,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
)
    return _vector_transport_along!(M, Y, p, X, c, m)
end
function _vector_transport_along!(M::AbstractManifold, Y, p, X, c, ::ParallelTransport)
    return parallel_transport_along!(M, Y, p, X, c)
end
function _vector_transport_along!(
    M::AbstractManifold,
    Y,
    p,
    X,
    c,
    m::DifferentiatedRetractionVectorTransport,
)
    return vector_transport_along_diff!(M, Y, p, X, c, m.retraction)
end
function vector_transport_along_diff! end
function _vector_transport_along!(M::AbstractManifold, Y, p, X, c, ::ProjectionTransport)
    return vector_transport_along_project!(M, Y, p, X, c)
end
function vector_transport_along_project! end

@doc raw"""
    vector_transport_along!(
        M::AbstractManifold,
        Y,
        p,
        X,
        c::AbstractVector,
        m::AbstractVectorTransportMethod
    ) where {T}

Compute the vector transport along a discretized curve `c` using an
[`AbstractVectorTransportMethod`](@ref) `method` succesively along the sampled curve.
"""
function vector_transport_along!(
    M::AbstractManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::AbstractVectorTransportMethod,
)
    n = length(c)
    if n == 0
        copyto!(Y, X)
    else
        vector_transport_to!(M, Y, p, X, c[1], m)
        for i in 1:(length(c) - 1)
            vector_transport_to!(M, Y, c[i], Y, c[i + 1], m)
        end
    end
    return Y
end

@doc raw"""
    function vector_transport_along!(
        M::AbstractManifold,
        Y,
        p,
        X,
        c::AbstractVector,
        m::PoleLadderTransport
    )

Compute the vector transport along a discretized curve using
[`PoleLadderTransport`](@ref) succesively along the sampled curve.
This method is avoiding additional allocations as well as inner exp/log by performing all
ladder steps on the manifold and only computing one tangent vector in the end.
"""
function _vector_transport_along!(
    M::AbstractManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::PoleLadderTransport,
)
    clen = length(c)
    if clen == 0
        copyto!(Y, X)
    else
        d = retract(M, p, X, m.retraction)
        mp = mid_point(M, p, c[1])
        pole_ladder!(
            M,
            d,
            p,
            d,
            c[1],
            mp,
            Y;
            retraction = m.retraction,
            inverse_retraction = m.inverse_retraction,
        )
        for i in 1:(clen - 1)
            # precompute mid point inplace
            ci = c[i]
            cip1 = c[i + 1]
            mid_point!(M, mp, ci, cip1)
            # compute new ladder point
            pole_ladder!(
                M,
                d,
                ci,
                d,
                cip1,
                mp,
                Y;
                retraction = m.retraction,
                inverse_retraction = m.inverse_retraction,
            )
        end
        inverse_retract!(M, Y, c[clen], d, m.inverse_retraction)
        Y *= (-1)^clen
    end
    return Y
end
@doc raw"""
    vector_transport_along!(
        M::AbstractManifold,
        Y,
        p,
        X,
        c::AbstractVector,
        m::SchildsLadderTransport
    )

Compute the vector transport along a discretized curve using
[`SchildsLadderTransport`](@ref) succesively along the sampled curve.
This method is avoiding additional allocations as well as inner exp/log by performing all
ladder steps on the manifold and only computing one tangent vector in the end.
"""
function _vector_transport_along!(
    M::AbstractManifold,
    Y,
    p,
    X,
    c::AbstractVector,
    m::SchildsLadderTransport,
)
    clen = length(c)
    if clen == 0
        copyto!(Y, X)
    else
        d = retract(M, p, X, m.retraction)
        mp = mid_point(M, c[1], d)
        schilds_ladder!(
            M,
            d,
            p,
            d,
            c[1],
            mp,
            Y;
            retraction = m.retraction,
            inverse_retraction = m.inverse_retraction,
        )
        for i in 1:(clen - 1)
            ci = c[i]
            cip1 = c[i + 1]
            # precompute mid point inplace
            mid_point!(M, mp, cip1, d)
            # compute new ladder point
            schilds_ladder!(
                M,
                d,
                ci,
                d,
                cip1,
                mp,
                Y;
                retraction = m.retraction,
                inverse_retraction = m.inverse_retraction,
            )
        end
        inverse_retract!(M, Y, c[clen], d, m.inverse_retraction)
    end
    return Y
end

@doc raw"""
    vector_transport_direction(M::AbstractManifold, p, X, d)
    vector_transport_direction(M::AbstractManifold, p, X, d, m::AbstractVectorTransportMethod)
    vector_transport_direction(M::AbstractManifold, p, X, d, m::AbstractVectorTransportMethod, r::AbstractRetractionMethod)

Given an [`AbstractManifold`](@ref) ``\mathcal M`` the vector transport is a generalization of the
[`parallel_transport_direction`](@ref) that identifies vectors from different tangent spaces.

More precisely using [^AbsilMahonySepulchre2008], Def. 8.1.1, a vector transport
``T_{p,d}: T_p\mathcal M \to T_q\mathcal M``, ``p∈ \mathcal M``, ``Y∈ T_p\mathcal M`` is a smooth mapping
associated to a retraction ``\operatorname{retr}_p(Y) = q`` such that

1. (associated retraction) ``\mathcal T_{p,d}X ∈ T_q\mathcal M`` if and only if ``q = \operatorname{retr}_p(d)``.
2. (consistency) ``\mathcal T_{p,0_p}X = X`` for all ``X∈T_p\mathcal M``
3. (linearity) ``\mathcal T_{p,d}(αX+βY) = α\mathcal T_{p,d}X + β\mathcal T_{p,d}Y``

For the [`AbstractVectorTransportMethod`](@ref) we might even omit the third point,
but the [`AbstractLinearVectorTransportMethod`](@ref)s are linear.

# Input Parameters
* `M` a manifold
* `p` indicating the tangent space of
* `X` the tangent vector to be transported
* `d` indicating a transport direction (and distance through its length)
* `m` an [`AbstractVectorTransportMethod`](@ref), by default [`default_vector_transport_method`](@ref), so usually [`ParallelTransport`](@ref)
* `r` a [`AbstractRetractionMethod`](@ref) by default [`default_retraction_method`](@ref), so usually [`ExponentialRetraction`](@ref)

Instead of spcifying a start direction `d` one can equivalently also specify a target tanget space
``T_q\mathcal M``, see [`vector_transport_to`](@ref).
This is what this method defaults to, if no implementation if given.

[^AbsilMahonySepulchre2008]:
    > Absil, P.-A., Mahony, R. and Sepulchre R.,
    > _Optimization Algorithms on Matrix Manifolds_
    > Princeton University Press, 2008,
    > doi: [10.1515/9781400830244](https://doi.org/10.1515/9781400830244)
    > [open access](http://press.princeton.edu/chapters/absil/)
"""
function vector_transport_direction(
    M::AbstractManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
    r::AbstractRetractionMethod = default_retraction_method(M),
)
    return _vector_transport_direction(M, p, X, d, m, r)
end
function _vector_transport_direction(
    M::AbstractManifold,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
    r::AbstractRetractionMethod = default_retraction_method(M),
)
    return vector_transport_to(M, p, X, retract(M, p, d, r), m, r)
end
function _vector_transport_direction(
    M::AbstractManifold,
    p,
    X,
    d,
    ::ParallelTransport,
    ::ExponentialRetraction,
)
    return parallel_transport_direction(M, p, X, d)
end


"""
    vector_transport_direction!(M::AbstractManifold, Y, p, X, d)
    vector_transport_direction!(M::AbstractManifold, Y, p, X, d, m::AbstractVectorTransportMethod)
    vector_transport_direction!(M::AbstractManifold, Y, p, X, d, m::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`AbstractManifold`](@ref) `M`
in the direction indicated by the tangent vector `d` at `p`. By default, [`retract`](@ref) and
[`vector_transport_to!`](@ref) are used with the `m` and `r`, which default
to [`default_vector_transport_method`](@ref)`(M)` and [`default_retraction_method`](@ref)`(M)`, respectively.
The result is saved to `Y`.

See [`vector_transport_direction`](@ref) for more details.
"""
function vector_transport_direction!(
    M::AbstractManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
    r::AbstractRetractionMethod = default_retraction_method(M),
)
    return _vector_transport_direction!(M, Y, p, X, d, m, r)
end


function _vector_transport_direction!(
    M::AbstractManifold,
    Y,
    p,
    X,
    d,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
    r::AbstractRetractionMethod = default_retraction_method(M),
)
    return vector_transport_to!(M, Y, p, X, retract(M, p, d, r), m)
end
function _vector_transport_direction!(
    M::AbstractManifold,
    Y,
    p,
    X,
    d,
    ::ParallelTransport,
    ::ExponentialRetraction,
)
    return parallel_transport_direction!(M, Y, p, X, d)
end

@doc raw"""
    vector_transport_to(M::AbstractManifold, p, X, q)
    vector_transport_to(M::AbstractManifold, p, X, q, m::AbstractVectorTransportMethod)
    vector_transport_to(M::AbstractManifold, p, X, q, m::AbstractVectorTransportMethod, r::AbstractRetractionMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`AbstractManifold`](@ref) `M`
along a curve specified by the [`AbstractRetractionMethod`](@ref) `r` to the tangent space at another point `q`.
By `m` is the [`default_vector_transport_method`](@ref)`(M)` and `r` the [`default_retraction_method`](@ref).

Note that some methods do not require the [`AbstractRetractionMethod`](@ref), for example [`ProjectionTransport`](@ref),
while for others, it has to match the type like for the [`DifferentiatedRetractionVectorTransport`](@ref)
or the [`parallel_transport_to`](@ref) which is always meant along the [`shortest_geodesic`](@ref).

This method is equivalent to using ``d = \operatorname{retr}^{-1}_p(q)`` in [`vector_transport_direction`](@ref)`(M, p, X, q, m, r)`,
where you can find the formal definition.
However, this method form is – to the best of our knowledge – moreoften implemented.
"""
function vector_transport_to(
    M::AbstractManifold,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
    r::AbstractRetractionMethod = default_retraction_method(M),
)
    return _vector_transport_to(M, p, X, q, m, r)
end
function _vector_transport_to(
    M::AbstractManifold,
    p,
    X,
    q,
    ::ParallelTransport,
    ::ExponentialRetraction,
)
    return parallel_transport_to(M, p, X, q)
end
function _vector_transport_to(
    M::AbstractManifold,
    p,
    X,
    q,
    ::DifferentiatedRetractionVectorTransport{R},
    r::R,
) where {R<:AbstractRetractionMethod}
    return vector_transport_to_diff(M, p, X, q, r)
end
function vector_transport_to_diff(M::AbstractManifold, p, X, q, r)
    Y = allocate_result(M, vector_transport_to, X, p)
    return vector_transport_to_diff!(M, Y, p, X, q, r)
end
function _vector_transport_to(
    M::AbstractManifold,
    p,
    X,
    q,
    ::ProjectionTransport,
    ::AbstractRetractionMethod,
)
    return vector_transport_to_project(M, p, X, q)
end
function vector_transport_to_project(M::AbstractManifold, p, X, q)
    Y = allocate_result(M, vector_transport_to, X, p)
    return vector_transport_to_project!(M, Y, p, X, q)
end

"""
    vector_transport_to!(M::AbstractManifold, Y, p, X, q)
    vector_transport_to!(M::AbstractManifold, Y, p, X, q, m::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`AbstractManifold`](@ref) `M`
to `q` using the [`AbstractVectorTransportMethod`](@ref) `m` and the [`AbstractRetractionMethod`](@ref) `r`.

The result is computed in `Y`.
See [`vector_transport_to`](@ref) for more details.
"""
function vector_transport_to!(
    M::AbstractManifold,
    Y,
    p,
    X,
    q,
    m::AbstractVectorTransportMethod = default_vector_transport_method(M),
    r::AbstractRetractionMethod = default_retraction_method(M),
)
    return _vector_transport_to!(M, Y, p, X, q, m, r)
end
function _vector_transport_to!(
    M::AbstractManifold,
    Y,
    p,
    X,
    q,
    ::ParallelTransport,
    ::ExponentialRetraction,
)
    return parallel_transport_to!(M, Y, p, X, q)
end
function _vector_transport_to!(
    M::AbstractManifold,
    Y,
    p,
    X,
    q,
    ::DifferentiatedRetractionVectorTransport{R},
    r::R,
) where {R}
    return vector_transport_to_diff!(M, Y, p, X, q, r)
end
function vector_transport_to_diff! end
function _vector_transport_to!(
    M::AbstractManifold,
    Y,
    p,
    X,
    q,
    ::ProjectionTransport,
    ::AbstractRetractionMethod,
)
    return vector_transport_to_project!(M, Y, p, X, q)
end
function vector_transport_to_project! end

function _vector_transport_to(
    M::AbstractManifold,
    p,
    X,
    c,
    m::PoleLadderTransport{R},
    r::R,
) where {R}
    Y = allocate_result(M, vector_transport_to, X, p)
    return _vector_transport_to!(M, Y, p, X, c, m, r)
end

function _vector_transport_to!(
    M::AbstractManifold,
    Y,
    p,
    X,
    q,
    m::PoleLadderTransport{R},
    ::R,
) where {R}
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

function _vector_transport_to(
    M::AbstractManifold,
    p,
    X,
    c,
    m::ScaledVectorTransport,
    r::AbstractRetractionMethod,
)
    Y = allocate_result(M, vector_transport_to, X, p)
    return _vector_transport_to!(M, Y, p, X, c, m, r)
end
function _vector_transport_to!(
    M::AbstractManifold,
    Y,
    p,
    X,
    q,
    m::ScaledVectorTransport,
    r::AbstractRetractionMethod,
)
    vector_transport_to!(M, Y, p, X, q, m.method, r)
    Y .*= norm(M, p, X) / norm(M, q, Y)
    return Y
end

function _vector_transport_to(
    M::AbstractManifold,
    p,
    X,
    c,
    m::SchildsLadderTransport{R},
    r::R,
) where {R}
    Y = allocate_result(M, vector_transport_to, X, p)
    return _vector_transport_to!(M, Y, p, X, c, m, r)
end

function _vector_transport_to!(
    M::AbstractManifold,
    Y,
    p,
    X,
    q,
    m::SchildsLadderTransport{R},
    ::R,
) where {R}
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
