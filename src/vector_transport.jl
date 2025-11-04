"""
    AbstractVectorTransportMethod <: AbstractApproximationMethod

Abstract type for methods for transporting vectors. Such vector transports are not
necessarily linear.

# See also

[`AbstractLinearVectorTransportMethod`](@ref)
"""
abstract type AbstractVectorTransportMethod <: AbstractApproximationMethod end

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
see [AbsilMahonySepulchre:2008](@cite), Section 8.1.2 for more details.

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
struct DifferentiatedRetractionVectorTransport{R <: AbstractRetractionMethod} <:
    AbstractLinearVectorTransportMethod
    retraction::R
end

@doc raw"""
    EmbeddedVectorTransport{T<:AbstractVectorTransportMethod} <: AbstractVectorTransportMethod

Compute a vector transport by using the vector transport of type `T` in the embedding and projecting the result.

# Constructor

    EmbeddedVectorTransport(vt::AbstractVectorTransportMethod)

Generate the vector transport with vector transport `vt` to use in the embedding.
"""
struct EmbeddedVectorTransport{T <: AbstractVectorTransportMethod} <:
    AbstractVectorTransportMethod
    vector_transport::T
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
[`vector_transport_to`](@ref) or [`vector_transport_direction`](@ref).
See [`project`](@ref) for details.
"""
struct ProjectionTransport <: AbstractLinearVectorTransportMethod end


@doc raw"""
    PoleLadderTransport <: AbstractVectorTransportMethod

Specify to use [`pole_ladder`](@ref) as vector transport method within
[`vector_transport_to`](@ref) or [`vector_transport_direction`](@ref), i.e.

Let $X∈ T_p\mathcal M$ be a tangent vector at $p∈\mathcal M$ and $q∈\mathcal M$ the
point to transport to. Then $x = \exp_pX$ is used to call
`y = `[`pole_ladder`](@ref)`(M, p, x, q)` and the resulting vector is obtained by computing
$Y = -\log_qy$.

The [`PoleLadderTransport`](@ref) posesses two advantages compared to
[`SchildsLadderTransport`](@ref):
* it is cheaper to evaluate, if you want to transport several vectors, since the
  mid point $c$ then stays unchanged.
* while both methods are exact if the curvature is zero, pole ladder is even exact in
  symmetric Riemannian manifolds [Pennec:2018](@cite)

The pole ladder was was proposed in [LorenziPennec:2013](@cite). Its name stems from the fact that
it resembles a pole ladder when applied to a sequence of points usccessively.

# Constructor
````julia
PoleLadderTransport(
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
````
Construct the classical pole ladder that employs exp and log, i.e. as proposed
in[LorenziPennec:2013](@cite). For an even cheaper transport the inner operations can be
changed to an [`AbstractRetractionMethod`](@ref) `retraction` and an
[`AbstractInverseRetractionMethod`](@ref) `inverse_retraction`, respectively.
"""
struct PoleLadderTransport{
        RT <: AbstractRetractionMethod, IRT <: AbstractInverseRetractionMethod,
    } <: AbstractLinearVectorTransportMethod
    retraction::RT
    inverse_retraction::IRT
    function PoleLadderTransport(
            retraction = ExponentialRetraction(),
            inverse_retraction = LogarithmicInverseRetraction(),
        )
        return new{typeof(retraction), typeof(inverse_retraction)}(
            retraction,
            inverse_retraction,
        )
    end
end

@doc raw"""
    ScaledVectorTransport{T} <: AbstractVectorTransportMethod

Introduce a scaled variant of any [`AbstractVectorTransportMethod`](@ref) `T`,
as introduced in [SatoIwai:2013](@cite) for some ``X∈ T_p\mathcal M`` as

```math
    \mathcal T^{\mathrm{S}}(X) = \frac{\lVert X\rVert_p}{\lVert \mathcal T(X)\rVert_q}\mathcal T(X).
```

Note that the resulting point `q` has to be known, i.e. for [`vector_transport_direction`](@ref)
the curve or more precisely its end point has to be known (via an exponential map or a
retraction). Therefore a default implementation is only provided for the [`vector_transport_to`](@ref)

# Constructor

    ScaledVectorTransport(m::AbstractVectorTransportMethod)
"""
struct ScaledVectorTransport{T <: AbstractVectorTransportMethod} <:
    AbstractVectorTransportMethod
    method::T
end

@doc raw"""
    SchildsLadderTransport <: AbstractVectorTransportMethod

Specify to use [`schilds_ladder`](@ref) as vector transport method within
[`vector_transport_to`](@ref) or [`vector_transport_direction`](@ref), i.e.

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
image of a ladder. The approximation was proposed in [EhlersPiraniSchild:1972](@cite).

# Constructor
````julia
SchildsLadderTransport(
    retraction = ExponentialRetraction(),
    inverse_retraction = LogarithmicInverseRetraction(),
)
````
Construct the classical Schilds ladder that employs exp and log, i.e. as proposed
in [EhlersPiraniSchild:1972](@cite). For an even cheaper transport these inner operations can be
changed to an [`AbstractRetractionMethod`](@ref) `retraction` and an
[`AbstractInverseRetractionMethod`](@ref) `inverse_retraction`, respectively.
"""
struct SchildsLadderTransport{
        RT <: AbstractRetractionMethod, IRT <: AbstractInverseRetractionMethod,
    } <: AbstractLinearVectorTransportMethod
    retraction::RT
    inverse_retraction::IRT
    function SchildsLadderTransport(
            retraction = ExponentialRetraction(),
            inverse_retraction = LogarithmicInverseRetraction(),
        )
        return new{typeof(retraction), typeof(inverse_retraction)}(
            retraction,
            inverse_retraction,
        )
    end
end

@doc raw"""
    VectorTransportDirection{VM<:AbstractVectorTransportMethod,RM<:AbstractRetractionMethod}
        <: AbstractVectorTransportMethod

Specify a [`vector_transport_direction`](@ref) using a [`AbstractVectorTransportMethod`](@ref)
with explicitly using the [`AbstractRetractionMethod`](@ref) to determine the point in
the specified direction where to transsport to.
Note that you only need this for the non-default (non-implicit) second retraction method
associated to a vector transport, i.e. when a first implementation assumed
an implicit associated retraction.
"""
struct VectorTransportDirection{
        VM <: AbstractVectorTransportMethod, RM <: AbstractRetractionMethod,
    } <: AbstractVectorTransportMethod
    retraction::RM
    vector_transport::VM
    function VectorTransportDirection(
            vector_transport = ParallelTransport(),
            retraction = ExponentialRetraction(),
        )
        return new{typeof(vector_transport), typeof(retraction)}(
            retraction,
            vector_transport,
        )
    end
end

@doc raw"""
    VectorTransportTo{VM<:AbstractVectorTransportMethod,RM<:AbstractRetractionMethod}
        <: AbstractVectorTransportMethod

Specify a [`vector_transport_to`](@ref) using a [`AbstractVectorTransportMethod`](@ref)
with explicitly using the [`AbstractInverseRetractionMethod`](@ref) to determine the direction
that transports from  in `p`to `q`.
Note that you only need this for the non-default (non-implicit) second retraction method
associated to a vector transport, i.e. when a first implementation assumed
an implicit associated retraction.
"""
struct VectorTransportTo{
        VM <: AbstractVectorTransportMethod, IM <: AbstractInverseRetractionMethod,
    } <: AbstractVectorTransportMethod
    inverse_retraction::IM
    vector_transport::VM
    function VectorTransportTo(
            vector_transport = ParallelTransport(),
            inverse_retraction = LogarithmicInverseRetraction(),
        )
        return new{typeof(vector_transport), typeof(inverse_retraction)}(
            inverse_retraction,
            vector_transport,
        )
    end
end

"""
    VectorTransportWithKeywords{V<:AbstractVectorTransportMethod, K} <: AbstractVectorTransportMethod

Since vector transports might have keywords, this type is a way to set them as an own type to be
used as a specific vector transport.
Another reason for this type is that we dispatch on the vector transport first and only the
last layer would be implemented with keywords, so this way they can be passed down.

## Fields

* `vector_transport` the vector transport that is decorated with keywords
* `kwargs` the keyword arguments

Note that you can nest this type. Then the most outer specification of a keyword is used.

## Constructor

    VectorTransportWithKeywords(m::T; kwargs...) where {T <: AbstractVectorTransportMethod}

Specify the subtype `T <: `[`AbstractVectorTransportMethod`](@ref) to have keywords `kwargs...`.
"""
struct VectorTransportWithKeywords{T <: AbstractVectorTransportMethod, K} <:
    AbstractVectorTransportMethod
    vector_transport::T
    kwargs::K
end
function VectorTransportWithKeywords(
        m::T;
        kwargs...,
    ) where {T <: AbstractVectorTransportMethod}
    return VectorTransportWithKeywords{T, typeof(kwargs)}(m, kwargs)
end

"""
    default_vector_transport_method(M::AbstractManifold)
    default_vector_transport_method(M::AbstractManifold, ::Type{T}) where {T}

The [`AbstractVectorTransportMethod`](@ref) that is used when calling [`vector_transport_to`](@ref)
or [`vector_transport_direction`](@ref) without specifying the vector transport method.
By default, this is [`ParallelTransport`](@ref).

This method can also be specified more precisely with a point type `T`, for the case
that on a `M` there are two different representations of points, which provide
different vector transport methods.
"""
function default_vector_transport_method(::AbstractManifold)
    return ParallelTransport()
end
function default_vector_transport_method(M::AbstractManifold, ::Type{T}) where {T}
    return default_vector_transport_method(M)
end


@doc raw"""
    pole_ladder(
        M, p, d, q, c = mid_point(M, p, q);
        retraction=default_retraction_method(M, typeof(p)),
        inverse_retraction=default_inverse_retraction_method(M, typeof(p))
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
        M, p, d, q, c = mid_point(M, p, q);
        retraction = default_retraction_method(M, typeof(p)),
        inverse_retraction = default_inverse_retraction_method(M, typeof(p)),
    )
    return retract(M, d, 2 * inverse_retract(M, d, c, inverse_retraction), retraction)
end
@doc raw"""
    pole_ladder(
        M, pl, p, d, q, c = mid_point(M, p, q), X = allocate_result_type(M, log, d, c);
        retraction = default_retraction_method(M, typeof(p)),
        inverse_retraction = default_inverse_retraction_method(M, typeof(p)),
    )

Compute the [`pole_ladder`](@ref), i.e. the result is saved in `pl`.
`X` is used for storing intermediate inverse retraction.
"""
function pole_ladder!(
        M, pl, p, d, q, c = mid_point(M, p, q), X = allocate_result(M, log, d, c);
        retraction = default_retraction_method(M, typeof(p)),
        inverse_retraction = default_inverse_retraction_method(M, typeof(p)),
    )
    inverse_retract!(M, X, d, c, inverse_retraction)
    X *= 2
    return retract!(M, pl, d, X, retraction)
end

@doc raw"""
    schilds_ladder(
        M, p, d, q, c = mid_point(M, q, d);
        retraction = default_retraction_method(M, typeof(p)),
        inverse_retraction = default_inverse_retraction_method(M, typeof(p)),
    )

Perform an inner step of schilds ladder, which can be used as a
[`vector_transport_to`](@ref), see [`SchildsLadderTransport`](@ref).
Let $c = \gamma_{q,d}(\frac{1}{2})$ denote the mid point
on the shortest geodesic connecting $q$ and the point $d$. Then Schild's ladder reads as

````math
\operatorname{Sl}(p,d,q) = \operatorname{retr}_p( 2\operatorname{retr}_p^{-1} c)
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
        M, p, d, q, c = mid_point(M, q, d);
        retraction = default_retraction_method(M, typeof(p)),
        inverse_retraction = default_inverse_retraction_method(M, typeof(p)),
    )
    return retract(M, p, 2 * inverse_retract(M, p, c, inverse_retraction), retraction)
end
@doc raw"""
    schilds_ladder!(M, sl, p, d, q, c = mid_point(M, q, d),
        X = allocate_result_type(M, log, d, c);
        retraction = default_retraction_method(M, typeof(p)),
        inverse_retraction = default_inverse_retraction_method(M, typeof(p)),
    )

Compute [`schilds_ladder`](@ref) and return the value in the parameter `sl`.
If the required mid point `c` was computed before, it can be passed using `c`,
and the allocation of new memory can be avoided providing a tangent vector `X`
for the interims result.
"""
function schilds_ladder!(
        M, sl, p, d, q, c = mid_point(M, q, d), X = allocate_result(M, log, d, c);
        retraction = default_retraction_method(M, typeof(p)),
        inverse_retraction = default_inverse_retraction_method(M, typeof(p)),
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

@doc raw"""
    vector_transport_direction(M::AbstractManifold, p, X, d)
    vector_transport_direction(M::AbstractManifold, p, X, d, m::AbstractVectorTransportMethod)

Given an [`AbstractManifold`](@ref) ``\mathcal M`` the vector transport is a generalization of the
[`parallel_transport_direction`](@ref) that identifies vectors from different tangent spaces.

More precisely using [AbsilMahonySepulchre:2008](@cite), Def. 8.1.1, a vector transport
``T_{p,d}: T_p\mathcal M \to T_q\mathcal M``, ``p∈ \mathcal M``, ``Y∈ T_p\mathcal M`` is a smooth mapping
associated to a retraction ``\operatorname{retr}_p(Y) = q`` such that

1. (associated retraction) ``\mathcal T_{p,d}X ∈ T_q\mathcal M`` if and only if ``q = \operatorname{retr}_p(d)``.
2. (consistency) ``\mathcal T_{p,0_p}X = X`` for all ``X∈T_p\mathcal M``
3. (linearity) ``\mathcal T_{p,d}(αX+βY) = α\mathcal T_{p,d}X + β\mathcal T_{p,d}Y``

For the [`AbstractVectorTransportMethod`](@ref) we might even omit the third point.
The [`AbstractLinearVectorTransportMethod`](@ref)s are linear.

# Input Parameters
* `M` a manifold
* `p` indicating the tangent space of
* `X` the tangent vector to be transported
* `d` indicating a transport direction (and distance through its length)
* `m` an [`AbstractVectorTransportMethod`](@ref), by default [`default_vector_transport_method`](@ref), so usually [`ParallelTransport`](@ref)

Usually this method requires a [`AbstractRetractionMethod`](@ref) as well.
By default this is assumed to be the [`default_retraction_method`](@ref) or
implicitly given (and documented) for a vector transport.
To explicitly distinguish different retractions for a vector transport,
see [`VectorTransportDirection`](@ref).

Instead of spcifying a start direction `d` one can equivalently also specify a target tanget space
``T_q\mathcal M``, see [`vector_transport_to`](@ref).
By default [`vector_transport_direction`](@ref) falls back to using [`vector_transport_to`](@ref),
using the [`default_retraction_method`](@ref) on `M`.
"""
function vector_transport_direction(
        M::AbstractManifold, p, X, d,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    return _vector_transport_direction(M, p, X, d, m)
end
function _vector_transport_direction(
        M::AbstractManifold, p, X, d,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p));
        kwargs...,
    )
    # allocate first
    Y = allocate_result(M, vector_transport_direction, X, p, d)
    return vector_transport_direction!(M, Y, p, X, d, m; kwargs...)
end
function _vector_transport_direction(
        M::AbstractManifold,
        p,
        X,
        d,
        m::VectorTransportDirection;
        kwargs...,
    )
    mv =
        length(kwargs) > 0 ? VectorTransportWithKeywords(m.vector_transport; kwargs...) :
        m.vector_transport
    mr = m.retraction
    return vector_transport_to(M, p, X, retract(M, p, d, mr), mv)
end
function _vector_transport_direction(
        M::AbstractManifold,
        p,
        X,
        d,
        m::VectorTransportWithKeywords;
        kwargs...,
    )
    return _vector_transport_direction(
        M,
        p,
        X,
        d,
        m.vector_transport;
        kwargs...,
        m.kwargs...,
    )
end

"""
    vector_transport_direction!(M::AbstractManifold, Y, p, X, d)
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
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p));
        kwargs...,
    )
    return _vector_transport_direction!(M, Y, p, X, d, m; kwargs...)
end
function _vector_transport_direction!(
        M::AbstractManifold,
        Y,
        p,
        X,
        d,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p));
        kwargs...,
    )
    r = default_retraction_method(M, typeof(p))
    v = length(kwargs) > 0 ? VectorTransportWithKeywords(m; kwargs...) : m
    return vector_transport_to!(M, Y, p, X, retract(M, p, d, r), v)
end
function _vector_transport_direction!(
        M::AbstractManifold,
        Y,
        p,
        X,
        d,
        m::VectorTransportDirection;
        kwargs...,
    )
    v =
        length(kwargs) > 0 ? VectorTransportWithKeywords(m.vector_transport; kwargs...) :
        m.vector_transport
    return vector_transport_to!(M, Y, p, X, retract(M, p, d, m.retraction), v)
end
function _vector_transport_direction!(
        M::AbstractManifold,
        Y,
        p,
        X,
        d,
        m::DifferentiatedRetractionVectorTransport;
        kwargs...,
    )
    return vector_transport_direction_diff!(M, Y, p, X, d, m.retraction; kwargs...)
end
function _vector_transport_direction!(
        M::AbstractManifold,
        Y,
        p,
        X,
        d,
        m::EmbeddedVectorTransport;
        kwargs...,
    )
    return vector_transport_direction_embedded!(
        M,
        Y,
        p,
        X,
        d,
        m.vector_transport;
        kwargs...,
    )
end
@doc raw"""
    vector_transport_direction_diff!(M::AbstractManifold, Y, p, X, d, m::AbstractRetractionMethod)

Compute the vector transport of `X` from ``T_p\mathcal M`` into the direction `d`
using the differential of the [`AbstractRetractionMethod`](@ref) `m` in place of `Y`.
"""
vector_transport_direction_diff!(M, Y, p, X, d, m)

function vector_transport_direction_diff! end

function _vector_transport_direction!(
        M::AbstractManifold,
        Y,
        p,
        X,
        d,
        ::ParallelTransport;
        kwargs...,
    )
    return parallel_transport_direction!(M, Y, p, X, d; kwargs...)
end

function _vector_transport_direction!(
        M::AbstractManifold,
        Y,
        p,
        X,
        d,
        m::VectorTransportWithKeywords;
        kwargs...,
    )
    return _vector_transport_direction!(
        M,
        Y,
        p,
        X,
        d,
        m.vector_transport;
        kwargs...,
        m.kwargs...,
    )
end

@doc raw"""
    vector_transport_direction_embedded!(M::AbstractManifold, Y, p, X, d, m::AbstractVectorTransportMethod)

Compute the vector transport of `X` from ``T_p\mathcal M`` into the direction `d`
using the [`AbstractRetractionMethod`](@ref) `m` in the embedding.

The default implementataion requires one allocation for the points and tangent vectors in the
embedding and the resulting point, but the final projection is performed in place of `Y`
"""
function vector_transport_direction_embedded!(
        M::AbstractManifold,
        Y,
        p::P,
        X,
        d,
        m::AbstractVectorTransportMethod,
    ) where {P}
    p_e = embed(M, p)
    d_e = embed(M, d)
    X_e = embed(M, p, X)
    Y_e = vector_transport_direction(get_embedding(M, P), p_e, X_e, d_e, m)
    q = exp(M, p, d)
    return project!(M, Y, q, Y_e)
end

@doc raw"""
    vector_transport_to(M::AbstractManifold, p, X, q)
    vector_transport_to(M::AbstractManifold, p, X, q, m::AbstractVectorTransportMethod)
    vector_transport_to(M::AbstractManifold, p, X, q, m::AbstractVectorTransportMethod)

Transport a vector `X` from the tangent space at a point `p` on the [`AbstractManifold`](@ref) `M`
along a curve implicitly given by an [`AbstractRetractionMethod`](@ref) associated to `m`.
By default `m` is the [`default_vector_transport_method`](@ref)`(M)`.
To explicitly specify a (different) retraction to the implicitly assumeed retraction, see [`VectorTransportTo`](@ref).
Note that some vector transport methods might also carry their own retraction they are associated to,
like the  [`DifferentiatedRetractionVectorTransport`](@ref) and some are even independent of the retraction, for example the [`ProjectionTransport`](@ref).

This method is equivalent to using ``d = \operatorname{retr}^{-1}_p(q)`` in [`vector_transport_direction`](@ref)`(M, p, X, q, m, r)`,
where you can find the formal definition. This is the fallback for [`VectorTransportTo`](@ref).
"""
function vector_transport_to(
        M::AbstractManifold,
        p,
        X,
        q,
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    return _vector_transport_to(M, p, X, q, m)
end
function _vector_transport_to(M::AbstractManifold, p, X, q, m::VectorTransportTo; kwargs...)
    d = inverse_retract(M, p, q, m.inverse_retraction)
    v =
        length(kwargs) > 0 ? VectorTransportWithKeywords(m.vector_transport; kwargs...) :
        m.vector_transport
    return vector_transport_direction(M, p, X, d, v)
end
function _vector_transport_to(M::AbstractManifold, p, X, q, ::ParallelTransport; kwargs...)
    return parallel_transport_to(M, p, X, q; kwargs...)
end
function _vector_transport_to(
        M::AbstractManifold,
        p,
        X,
        q,
        m::AbstractVectorTransportMethod;
        kwargs...,
    )
    Y = allocate_result(M, vector_transport_to, X, p)
    return vector_transport_to!(M, Y, p, X, q, m; kwargs...)
end
function _vector_transport_to(
        M::AbstractManifold,
        p,
        X,
        q,
        m::VectorTransportWithKeywords;
        kwargs...,
    )
    return _vector_transport_to(M, p, X, q, m.vector_transport; kwargs..., m.kwargs...)
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
        m::AbstractVectorTransportMethod = default_vector_transport_method(M, typeof(p)),
    )
    return _vector_transport_to!(M, Y, p, X, q, m)
end
function _vector_transport_to!(
        M::AbstractManifold,
        Y,
        p,
        X,
        q,
        m::VectorTransportTo;
        kwargs...,
    )
    d = inverse_retract(M, p, q, m.inverse_retraction)
    return _vector_transport_direction!(M, Y, p, X, d, m.vector_transport; kwargs...)
end
function _vector_transport_to!(
        M::AbstractManifold,
        Y,
        p,
        X,
        q,
        ::ParallelTransport;
        kwargs...,
    )
    return parallel_transport_to!(M, Y, p, X, q; kwargs...)
end
function _vector_transport_to!(
        M::AbstractManifold,
        Y,
        p,
        X,
        q,
        m::DifferentiatedRetractionVectorTransport;
        kwargs...,
    )
    return vector_transport_to_diff!(M, Y, p, X, q, m.retraction; kwargs...)
end

function _vector_transport_to!(
        M::AbstractManifold,
        Y,
        p,
        X,
        q,
        m::EmbeddedVectorTransport;
        kwargs...,
    )
    return vector_transport_to_embedded!(M, Y, p, X, q, m.vector_transport; kwargs...)
end
function _vector_transport_to!(
        M::AbstractManifold,
        Y,
        p,
        X,
        q,
        ::ProjectionTransport;
        kwargs...,
    )
    return vector_transport_to_project!(M, Y, p, X, q; kwargs...)
end
function _vector_transport_to!(
        M::AbstractManifold,
        Y,
        p,
        X,
        q,
        m::PoleLadderTransport;
        kwargs...,
    )
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
            kwargs...,
        ),
        m.inverse_retraction,
    )
    copyto!(Y, -Y)
    return Y
end
function _vector_transport_to!(
        M::AbstractManifold,
        Y,
        p,
        X,
        q,
        m::ScaledVectorTransport;
        kwargs...,
    )
    v = length(kwargs) > 0 ? VectorTransportWithKeywords(m.method; kwargs...) : m.method
    vector_transport_to!(M, Y, p, X, q, v)
    Y .*= norm(M, p, X) / norm(M, q, Y)
    return Y
end
function _vector_transport_to!(
        M::AbstractManifold, Y, p, X, q, m::VectorTransportWithKeywords; kwargs...,
    )
    return _vector_transport_to!(M, Y, p, X, q, m.vector_transport; kwargs..., m.kwargs...)
end
function _vector_transport_to!(
        M::AbstractManifold, Y, p, X, q, m::SchildsLadderTransport;
        kwargs...,
    )
    return inverse_retract!(
        M, Y, q,
        schilds_ladder(
            M, p, retract(M, p, X, m.retraction), q;
            retraction = m.retraction,
            inverse_retraction = m.inverse_retraction,
            kwargs...,
        ),
        m.inverse_retraction,
    )
end

@doc raw"""
    vector_transport_to_diff(M::AbstractManifold, p, X, q, r)

Compute a vector transport by using a [`DifferentiatedRetractionVectorTransport`](@ref) `r` in place of `Y`.
"""
vector_transport_to_diff!(M::AbstractManifold, Y, p, X, q, r)

function vector_transport_to_diff! end

@doc raw"""
    vector_transport_to_embedded!(M::AbstractManifold, Y, p, X, q, m::AbstractRetractionMethod)

Compute the vector transport of `X` from ``T_p\mathcal M`` to the point `q`
using the  of the [`AbstractRetractionMethod`](@ref) `m` in th embedding.

The default implementataion requires one allocation for the points and tangent vectors in the
embedding and the resulting point, but the final projection is performed in place of `Y`
"""
function vector_transport_to_embedded!(M::AbstractManifold, Y, p::P, X, q, m) where {P}
    p_e = embed(M, p)
    X_e = embed(M, p, X)
    q_e = embed(M, q)
    Y_e = vector_transport_to(get_embedding(M, P), p_e, X_e, q_e, m)
    return project!(M, Y, q, Y_e)
end

@doc raw"""
    vector_transport_to_project!(M::AbstractManifold, Y, p, X, q)

Compute a vector transport by projecting ``X\in T_p\mathcal M`` onto the tangent
space ``T_q\mathcal M`` at ``q`` in place of `Y`.
"""
function vector_transport_to_project!(M::AbstractManifold, Y, p, X, q; kwargs...)
    # Note that we have to use embed (not embed!) since we do not have memory to store this embedded value in
    return project!(M, Y, q, embed(M, p, X); kwargs...)
end

# default estimation fallbacks with and without the T
function default_approximation_method(
        M::AbstractManifold,
        ::typeof(vector_transport_direction),
    )
    return default_vector_transport_method(M)
end
function default_approximation_method(M::AbstractManifold, ::typeof(vector_transport_to))
    return default_vector_transport_method(M)
end
function default_approximation_method(
        M::AbstractManifold,
        ::typeof(vector_transport_direction),
        T,
    )
    return default_vector_transport_method(M, T)
end
function default_approximation_method(M::AbstractManifold, ::typeof(vector_transport_to), T)
    return default_vector_transport_method(M, T)
end
