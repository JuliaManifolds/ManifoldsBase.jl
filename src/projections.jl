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
    local q
    try
        q = allocate_result_embedding(M, project, p)
    catch e
        # because we want `project` to default to identity
        q = allocate_result(M, project, p)
    end
    project!(M, q, p)
    return q
end

"""
    project!(M::AbstractManifold, q, p)

Project point `p` from the ambient space onto the [`AbstractManifold`](@ref) `M`. The result is
stored in `q`.
This method is only available for manifolds where implicitly an embedding or ambient space
is given. Additionally, the projection includes changing data representation, if applicable,
i.e. if the points on `M` are not represented in the same array data, the data is changed
accordingly.

See also: [`EmbeddedManifold`](@ref), [`embed!`](@ref embed!(M::AbstractManifold, q, p))
"""
project!(::AbstractManifold, q, p)

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
Lie algebra is performed, too.

See also: [`EmbeddedManifold`](@ref), [`embed`](@ref embed(M::AbstractManifold, p, X))
"""
function project(M::AbstractManifold, p, X)
    # Note that the order is switched,
    # since the allocation by default takes the type of the first.
    local Y
    try
        Y = allocate_result_embedding(M, project, X, p)
    catch e
        # because we want `project` to default to identity
        Y = allocate_result(M, project, X, p)
    end
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
project!(::AbstractManifold, Y, p, X)
