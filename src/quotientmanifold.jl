_doc_canonical_project = raw"""
    canonical_project(M::AbstractManifold, p)
    canonical_project!(M::AbstractManifold, q, p)

Compute the canonical projection ``π`` on a quotient manifold ``\mathcal M``.
The canonical (or natural) projection ``π`` from the total space ``\mathcal N``
onto ``\mathcal M`` given by

```math
    π = π_{\mathcal N, \mathcal M} : \mathcal N → \mathcal M, p ↦ π_{\mathcal N, \mathcal M}(p) = [p].
```

in other words, this function implicitly assumes, that the total space ``\mathcal N`` is given.
"""
@doc "$(_doc_canonical_project)"
function canonical_project(M::AbstractManifold, p)
    q = allocate_result(M, canonical_project, p)
    return canonical_project!(M, q, p)
end

function canonical_project! end

@doc "$(_doc_canonical_project)"
canonical_project!(M::AbstractManifold, q, p)

_doc_diff_canonical_project = raw"""
    diff_canonical_project(M::AbstractManifold, p, X)
    diff_canonical_project!(M::AbstractManifold, Y, p, X)

Compute the differential of the canonical projection ``π`` on a quotient manifold
``\mathcal M``.
The canonical (or natural) projection ``π`` from the total space ``\mathcal N``
onto ``\mathcal M``, such that its differential

```math
 Dπ(p) : T_p\mathcal N → T_{π(p)}\mathcal M
```

where again the total space might be implicitly assumed.
"""

@doc "$(_doc_diff_canonical_project)"
function diff_canonical_project(M::AbstractManifold, p, X)
    q = allocate_result(M, diff_canonical_project, p, X)
    return diff_canonical_project!(M, q, p, X)
end

function diff_canonical_project! end
@doc "$(_doc_diff_canonical_project)"
diff_canonical_project!(M::AbstractManifold, q, p)

_doc_horizontal_lift = raw"""
    horizontal_lift(N::AbstractManifold, q, X)
    horizontal_lift!(N::AbstractManifold, Y, q, X)

Given a point `q` in total space of the quotient manifold `N` such that ``p=π(q)`` is a point on
a quotient manifold `M` (implicitly given for the first case) and a tangent vector `X` this
method computes a tangent vector `Y` on the horizontal space of ``T_q\mathcal N``,
i.e. the subspace that is orthogonal to the kernel of ``Dπ(q)``.
"""

@doc "$(_doc_horizontal_lift)"
function horizontal_lift(N::AbstractManifold, q, X)
    Y = zero_vector(N, q)
    return horizontal_lift!(N, Y, q, X)
end

function horizontal_lift! end
@doc "$(_doc_horizontal_lift)"
horizontal_lift!(N::AbstractManifold, Y, q, X)

_doc_horizontal_lift = raw"""
    horizontal_component(M::AbstractManifold, p, X)
    horizontal_component!(M::AbstractManifold, Y, p, X)

Compute the horizontal component of tangent vector `X` at point `p`
in the total space of quotient manifold `N`.

This is often written as the space ``\mathrm{Hor}_p^π\mathcal N``.
"""

@doc "$(_doc_horizontal_lift)"
function horizontal_component(N::AbstractManifold, p, X)
    Y = allocate_result(N, horizontal_component, X, p)
    return horizontal_component!(N, Y, p, X)
end

function horizontal_component! end
@doc "$(_doc_horizontal_lift)"
horizontal_component!(N::AbstractManifold, Y, p, X)

function get_total_space end
@doc raw"""
    get_total_space(M::AbstractManifold)

Return the total space of a quotient manifold.
"""
get_total_space(::AbstractManifold)

_doc_vertical_component = raw"""
    vertical_component(N::AbstractManifold, p, X)
    vertical_component!(N::AbstractManifold, Y, p, X)

Compute the vertical component of tangent vector `X` at point `p`
in the total space of quotient manifold `N`.

This is often written as the space ``\mathrm{ver}_p^π\mathcal N``.
"""

@doc "$(_doc_vertical_component)"
function vertical_component(N::AbstractManifold, p, X)
    return X - horizontal_component(N, p, X)
end

@doc "$(_doc_vertical_component)"
function vertical_component!(N::AbstractManifold, Y, p, X)
    horizontal_component!(N, Y, p, X)
    Y .*= -1
    Y .+= X
    return Y
end
