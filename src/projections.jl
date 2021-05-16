"""
    project(M::Manifold, p)

Project point `p` from the ambient space of the [`Manifold`](@ref) `M` to `M`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Additionally, the projection includes changing data representation, if applicable,
i.e. if the points on `M` are not represented in the same array data, the data is changed
accordingly.

See also: [`EmbeddedManifold`](@ref), [`embed`](@ref embed(M::Manifold, p))
"""
function project(M::Manifold, p)
    q = allocate_result(M, project, p)
    project!(M, q, p)
    return q
end

"""
    project!(M::Manifold, q, p)

Project point `p` from the ambient space onto the [`Manifold`](@ref) `M`. The result is
storedin `q`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Additionally, the projection includes changing data representation, if applicable,
i.e. if the points on `M` are not represented in the same array data, the data is changed
accordingly.

See also: [`EmbeddedManifold`](@ref), [`embed!`](@ref embed!(M::Manifold, q, p))
"""
function project!(M::Manifold, q, p)
    return error(manifold_function_not_implemented_message(M, project!, q, p))
end

"""
    project(M::Manifold, p, X)

Project ambient space representation of a vector `X` to a tangent vector at point `p` on
the [`Manifold`](@ref) `M`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `project` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. after projection the change to the
Lie algebra is perfomed, too.

See also: [`EmbeddedManifold`](@ref), [`embed`](@ref embed(M::Manifold, p, X))
"""
function project(M::Manifold, p, X)
    # Note that the order is switched,
    # since the allocation by default takes the type of the first.
    Y = allocate_result(M, project, X, p)
    project!(M, Y, p, X)
    return Y
end

"""
    project!(M::Manifold, Y, p, X)

Project ambient space representation of a vector `X` to a tangent vector at point `p` on
the [`Manifold`](@ref) `M`. The result is saved in vector `Y`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `project!` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. after projection the change to the
Lie algebra is perfomed, too.

See also: [`EmbeddedManifold`](@ref), [`embed!`](@ref embed!(M::Manifold, Y, p, X))
"""
function project!(M::Manifold, Y, p, X)
    return error(manifold_function_not_implemented_message(M, project!, Y, p, X))
end
