
"""
    AbstractHomogeneousSpace{𝔽, GO, A <: ActionDirection, T <: AbstractDecoratorType} <: AbstractDecoratorManifold{𝔽,T}

An abstract homogeneous space, i.e. an [`AbstractManifold`](@ref) `M`
together with an [`AbstractGroupManifold`](@ref) `G` (to be precise its group action)
or a [`AbstractGroupOperation`](@ref) `O` that acts transitively on the Manifold.

Note that `M` equal to `G` yields the case of just a [`GroupManifold`](@ref) `G`.
"""
abstract type AbstractHomogeneousSpace{𝔽, O, A<: ActionDirection, T<:AbstractDecoratorType} <:
              AbstractDecoratorManifold{𝔽,T} end

"""
    base_group(H::AbstractHomogeneousSpace)

The group that acts on the [`AbstractHomogeneousSpace`](@ref) `H`.
"""
base_group(::AbstractHomogeneousSpace)

"""
    g_manifold(A::AbstractGroupAction)

The manifold the action `A` acts upon.
"""
base_manifold(::AbstractHomogeneousSpace)

allocate_result(H::AbstractHomogeneousSpace, f, p...) = allocate_result(base_manifold(H), f, p...)

@doc raw"""
    HomogeneousSpace{𝔽,O,T<:AbstractDecoratorType} <:
              Abstract{𝔽,T} end}
"""
struct HomogeneousSpace{𝔽, GO<:Union{AbstractGroupManifold{𝔽,O},AbstractGroupAction}, M<: AbstractManifold, T<:AbstractDecoratorType} <: AbstractHomogeneousSpace{𝔽,O,T}
    groupOrAction::G
    manifold::M
end

base_manifold(H::HomogeneousSpace) = H.manifold

base_group(H::HomogeneousSpace{𝔽,GO<:AbstractGroupManifold{𝔽}}) where 𝔽 = H.groupOrAction

@doc raw"""
    apply(H::AbstractHomogeneousSpace, a, p)

Apply action `a` to the point `p` using map $τ_a$, specified by `A`.
Unless otherwise specified, the right action is defined in terms of the left action:

````math
\mathrm{R}_a = \mathrm{L}_{a^{-1}}
````
"""
function apply(H::AbstractHomogeneousSpace, a, p)
    q = allocate_result(H, apply, p, a)
    apply!(A, q, a, p)
    return q
end

"""
    apply!(A::AbstractGroupAction, q, a, p)

Apply action `a` to the point `p` with the rule specified by `A`.
The result is saved in `q`.
"""
apply!(A::AbstractHomogeneousSpace{LeftAction}, q, a, p)
    return error(
        "apply! not implemented for action $(typeof(A)) and points $(typeof(q)), $(typeof(p)) and $(typeof(a)).",
    )
end
# TODO Missing: Left/Right action distinction – maybe the group shoiuld
#
#
# --- stopped refacvtoring until the difference between action and operation is clear!
# (then maybe actiondirection shoudl be called operationdirection?)
function apply!(A::AbstractHomogeneousSpace, q, a, p)
    ainv = inv(base_group(A), a)
    apply!(switch_direction(A), q, ainv, p)
    return q
end

"""
    inverse_apply(A::AbstractGroupAction, a, p)

Apply inverse of action `a` to the point `p`. The action is specified by `A`.
"""
function inverse_apply(A::AbstractGroupAction, a, p)
    q = allocate_result(A, inverse_apply, p, a)
    inverse_apply!(A, q, a, p)
    return q
end

"""
    inverse_apply!(A::AbstractGroupAction, q, a, p)

Apply inverse of action `a` to the point `p` with the rule specified by `A`.
The result is saved in `q`.
"""
function inverse_apply!(A::AbstractGroupAction, q, a, p)
    inva = inv(base_group(A), a)
    apply!(A, q, inva, p)
    return q
end

@doc raw"""
    apply_diff(A::AbstractGroupAction, a, p, X)

For group point $p ∈ \mathcal M$ and tangent vector $X ∈ T_p \mathcal M$, compute the action
on $X$ of the differential of the action of $a ∈ \mathcal{G}$, specified by rule `A`.
Written as $(\mathrm{d}τ_a)_p$, with the specified left or right convention, the
differential transports vectors

````math
(\mathrm{d}τ_a)_p : T_p \mathcal M → T_{τ_a p} \mathcal M
````
"""
function apply_diff(A::AbstractGroupAction, a, p, X)
    return error(
        "apply_diff not implemented for action $(typeof(A)), points $(typeof(a)) and $(typeof(p)), and vector $(typeof(X))",
    )
end

function apply_diff!(A::AbstractGroupAction, Y, a, p, X)
    return error(
        "apply_diff! not implemented for action $(typeof(A)), points $(typeof(a)) and $(typeof(p)), vectors $(typeof(Y)) and $(typeof(X))",
    )
end

@doc raw"""
    inverse_apply_diff(A::AbstractGroupAction, a, p, X)

For group point $p ∈ \mathcal M$ and tangent vector $X ∈ T_p \mathcal M$, compute the action
on $X$ of the differential of the inverse action of $a ∈ \mathcal{G}$, specified by rule
`A`. Written as $(\mathrm{d}τ_a^{-1})_p$, with the specified left or right convention,
the differential transports vectors

````math
(\mathrm{d}τ_a^{-1})_p : T_p \mathcal M → T_{τ_a^{-1} p} \mathcal M
````
"""
function inverse_apply_diff(A::AbstractGroupAction, a, p, X)
    return apply_diff(A, inv(base_group(A), a), p, X)
end

function inverse_apply_diff!(A::AbstractGroupAction, Y, a, p, X)
    return apply_diff!(A, Y, inv(base_group(A), a), p, X)
end

compose(A::AbstractGroupAction{LeftAction}, a, b) = compose(base_group(A), a, b)
compose(A::AbstractGroupAction{RightAction}, a, b) = compose(base_group(A), b, a)

compose!(A::AbstractGroupAction{LeftAction}, q, a, b) = compose!(base_group(A), q, a, b)
compose!(A::AbstractGroupAction{RightAction}, q, a, b) = compose!(base_group(A), q, b, a)

@doc raw"""
    optimal_alignment(A::AbstractGroupAction, p, q)

Calculate an action element $a$ of action `A` that acts upon `p` to produce
the element closest to `q` in the metric of the G-manifold:
```math
\arg\min_{a ∈ \mathcal{G}} d_{\mathcal M}(τ_a p, q)
```
where $\mathcal{G}$ is the group that acts on the G-manifold $\mathcal M$.
"""
function optimal_alignment(A::AbstractGroupAction, p, q)
    return error(
        "optimal_alignment not implemented for $(typeof(A)) and points $(typeof(p)) and $(typeof(q)).",
    )
end

"""
    optimal_alignment!(A::AbstractGroupAction, x, p, q)

Calculate an action element of action `A` that acts upon `p` to produce the element closest
to `q`.
The result is written to `x`.
"""
function optimal_alignment!(A::AbstractGroupAction, x, p, q)
    return copyto!(x, optimal_alignment(A, p, q))
end


@doc raw"""
    GroupOperationAction(group::AbstractGroupManifold, AD::ActionDirection = LeftAction())

Action of a group upon itself via left or right translation.
"""
struct GroupOperationAction{G,AD} <: AbstractGroupAction{AD}
    group::G
end

function GroupOperationAction(
    G::AbstractGroupManifold,
    ::TAD=LeftAction(),
) where {TAD<:ActionDirection}
    return GroupOperationAction{typeof(G),TAD}(G)
end

function Base.show(io::IO, A::GroupOperationAction)
    return print(io, "GroupOperationAction($(A.group), $(direction(A)))")
end

base_group(A::GroupOperationAction) = A.group

g_manifold(A::GroupOperationAction) = A.group

function switch_direction(A::GroupOperationAction)
    return GroupOperationAction(A.group, switch_direction(direction(A)))
end

apply(A::GroupOperationAction, a, p) = translate(A.group, a, p, direction(A))

apply!(A::GroupOperationAction, q, a, p) = translate!(A.group, q, a, p, direction(A))

function inverse_apply(A::GroupOperationAction, a, p)
    return inverse_translate(A.group, a, p, direction(A))
end

function inverse_apply!(A::GroupOperationAction, q, a, p)
    return inverse_translate!(A.group, q, a, p, direction(A))
end

function apply_diff(A::GroupOperationAction, a, p, X)
    return translate_diff(A.group, a, p, X, direction(A))
end

function apply_diff!(A::GroupOperationAction, Y, a, p, X)
    return translate_diff!(A.group, Y, a, p, X, direction(A))
end

function inverse_apply_diff(A::GroupOperationAction, a, p, X)
    return inverse_translate_diff(A.group, a, p, X, direction(A))
end

function inverse_apply_diff!(A::GroupOperationAction, Y, a, p, X)
    return inverse_translate_diff!(A.group, Y, a, p, X, direction(A))
end

function optimal_alignment(A::GroupOperationAction, p, q)
    return inverse_apply(switch_direction(A), p, q)
end

function optimal_alignment!(A::GroupOperationAction, x, p, q)
    return inverse_apply!(switch_direction(A), x, p, q)
end