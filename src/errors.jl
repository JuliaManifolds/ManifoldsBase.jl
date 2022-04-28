@doc """
    CompnentError{I,E} <: Exception

Store an error that occured in a component, where the additional `index` is stored.

# Fields

* `index` index where the error occured`
* `error` error that occured.
"""
struct ComponentManifoldError{I,E} <: Exception where {I,E<:Exception}
    index::I
    error::E
end
function ComponentManifoldError(i::I, e::E) where {I,E<:Exception}
    return ComponentManifoldError{I,E}(i, e)
end

@doc """
    CompositeManifoldError{T} <: Exception

A composite type to collect a set of errors that occured. Mainly used in conjunction
with [`ComponentManifoldError`](@ref) to store a set of errors that occured.

# Fields
* `errors` a `Vector` of `<:Exceptions`.
"""
struct CompositeManifoldError{T} <: Exception where {T<:Exception}
    errors::Vector{T}
end
CompositeManifoldError() = CompositeManifoldError{Exception}(Exception[])
function CompositeManifoldError(errors::Vector{T}) where {T<:Exception}
    return CompositeManifoldError{T}(errors)
end

isempty(c::CompositeManifoldError) = isempty(c.errors)
length(c::CompositeManifoldError) = length(c.errors)

function Base.show(io::IO, ex::ComponentManifoldError)
    return print(io, "ComponentManifoldError($(ex.index), $(ex.error))")
end
function Base.show(io::IO, ex::CompositeManifoldError)
    print(io, "CompositeManifoldError(")
    if !isempty(ex)
        print(io, "[")
        start = true
        for e in ex.errors
            show(io, e)
            print(io, ", ")
        end
        print(io, "]")
    end
    return print(io, ")")
end

function Base.showerror(io::IO, ex::ComponentManifoldError)
    print(io, "At #$(ex.index): ")
    return showerror(io, ex.error)
end

function Base.showerror(io::IO, ex::CompositeManifoldError)
    return if !isempty(ex)
        print(io, "CompositeManifoldError: ")
        showerror(io, ex.errors[1])
        remaining = length(ex) - 1
        if remaining > 0
            print(io, string("\n\n...and ", remaining, " more error(s).\n"))
        end
    else
        print(io, "CompositeManifoldError()\n")
    end
end

"""
    OutOfInjectivityRadiusError

An error thrown when a function (for example [`log`](@ref)arithmic map or
[`inverse_retract`](@ref)) is given arguments outside of its [`injectivity_radius`](@ref).
"""
struct OutOfInjectivityRadiusError <: Exception end

"""
    ManifoldDomainError{<:Exception} <: Exception

An error to represent a nested (Domain) error on a manifold, for example
if a point or tangent vector is invalid because its representation in some
embedding is already invalid.
"""
struct ManifoldDomainError{E} <: Exception where {E<:Exception}
    outer_text::String
    error::E
end

function Base.showerror(io::IO, ex::ManifoldDomainError)
    print(io, ex.outer_text)
    return print(io, ex.error)
end
