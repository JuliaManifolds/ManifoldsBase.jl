"""
    EmbeddedManifold{𝔽, MT <: AbstractManifold, NT <: AbstractManifold} <: AbstractDecoratorManifold{𝔽}

A type to represent an explicit embedding of a [`AbstractManifold`](@ref) `M` of type `MT` embedded
into a manifold `N` of type `NT`.
By default, an embedded manifold is set to be embedded, but neither isometrically embedded
nor a submanifold.

!!! note
    This type is not required if a manifold `M` is to be embedded in one specific manifold `N`.
     One can then just implement [`embed!`](@ref) and [`project!`](@ref).
    You can further pass functions to the embedding, for example, when it is an isometric embedding,
    by using an [`AbstractDecoratorManifold`](@ref).
    Only for a second –maybe considered non-default–
    embedding, this type should be considered in order to dispatch on different embed
    and project methods for different embeddings `N`.

# Fields

* `manifold` the manifold that is an embedded manifold
* `embedding` a second manifold, the first one is embedded into

# Constructor

    EmbeddedManifold(M, N)

Generate the `EmbeddedManifold` of the [`AbstractManifold`](@ref) `M` into the
[`AbstractManifold`](@ref) `N`.
"""
struct EmbeddedManifold{𝔽, MT <: AbstractManifold{𝔽}, NT <: AbstractManifold} <:
    AbstractDecoratorManifold{𝔽}
    manifold::MT
    embedding::NT
end

"""
    get_embedding_type(::EmbeddedManifold)

Specify the type of embedding. This by default returns [`EmbeddedManifoldType`](@ref)`()`.
It can be further specified by dispatching in the parameters of the [`EmbeddedManifold`](@ref).
"""
get_embedding_type(::EmbeddedManifold) = EmbeddedManifoldType()

function allocate_result(M::EmbeddedManifold, f::typeof(project), x...)
    T = allocate_result_type(M, f, x)
    return allocate(M, x[1], T, representation_size(base_manifold(M), T))
end

"""
    decorated_manifold(M::EmbeddedManifold, d::Val{N} = Val(-1))

Return the manifold of `M` that is decorated with its embedding. For this specific
type the internally stored enhanced manifold `M.manifold` is returned.

See also [`base_manifold`](@ref), where this is used to (potentially) completely undecorate the manifold.
"""
decorated_manifold(M::EmbeddedManifold) = M.manifold


"""
    get_embedding(M::EmbeddedManifold)

Return the embedding [`EmbeddedManifold`](@ref) `N` of `M`, if it exists.
"""
function get_embedding(M::EmbeddedManifold)
    return M.embedding
end

function project(M::EmbeddedManifold, p)
    q = allocate_result(M, project, p)
    project!(M, q, p)
    return q
end

function show(
        io::IO, M::EmbeddedManifold{𝔽, MT, NT},
    ) where {𝔽, MT <: AbstractManifold{𝔽}, NT <: AbstractManifold}
    return print(io, "EmbeddedManifold($(M.manifold), $(M.embedding))")
end
