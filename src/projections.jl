"""
    project(M::AbstractManifold, p)

Project point `p` from the ambient space of the [`AbstractManifold`](@ref) `M` to `M`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Additionally, the projection includes changing data representation, if applicable,
i.e. if the points on `M` are not represented in the same array data, the data is changed
accordingly.

See also: [`EmbeddedManifold`](@ref), [`embed`](@ref embed(M::AbstractManifold, p))
"""
function project(M::AbstractManifold, p)
    q = allocate_result_point(M, project, p)
    project!(M, q, p)
    return q
end

"""
    project!(M::AbstractManifold, q, p)

Project point `p` from the ambient space onto the [`AbstractManifold`](@ref) `M`. The result is
storedin `q`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Additionally, the projection includes changing data representation, if applicable,
i.e. if the points on `M` are not represented in the same array data, the data is changed
accordingly.

See also: [`EmbeddedManifold`](@ref), [`embed!`](@ref embed!(M::AbstractManifold, q, p))
"""
function project!(M::AbstractManifold, q, p)
    return error(manifold_function_not_implemented_message(M, project!, q, p))
end

"""
    project(M::AbstractManifold, p, X)

Project ambient space representation of a vector `X` to a tangent vector at point `p` on
the [`AbstractManifold`](@ref) `M`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `project` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. after projection the change to the
Lie algebra is perfomed, too.

See also: [`EmbeddedManifold`](@ref), [`embed`](@ref embed(M::AbstractManifold, p, X))
"""
function project(M::AbstractManifold, p, X)
    Y = allocate_result_vector(M, project, p, X)
    project!(M, Y, p, X)
    return Y
end

"""
    project!(M::AbstractManifold, Y, p, X)

Project ambient space representation of a vector `X` to a tangent vector at point `p` on
the [`AbstractManifold`](@ref) `M`. The result is saved in vector `Y`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given.
Additionally, `project!` includes changing data representation, if applicable, i.e.
if the tangents on `M` are not represented in the same way as points on the embedding,
the representation is changed accordingly. This is the case for example for Lie groups,
when tangent vectors are represented in the Lie algebra. after projection the change to the
Lie algebra is perfomed, too.

See also: [`EmbeddedManifold`](@ref), [`embed!`](@ref embed!(M::AbstractManifold, Y, p, X))
"""
function project!(M::AbstractManifold, Y, p, X)
    return error(manifold_function_not_implemented_message(M, project!, Y, p, X))
end
