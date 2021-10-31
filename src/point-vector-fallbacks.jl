
"""
    manifold_thing_forwards(T::Symbol, field::Symbol)

Introduce basic fallbacks for type `T` that represents points or vectors for a manifold.
Fallbacks will work by forwarding to field passed in the second argument.

List of forwarded functions:
* [`allocate`](@ref),
* [`copy`](@ref),
* [`copyto!`](@ref),
* [`number_eltype`](@ref) (only for values, not the type itself),
* `similar`,
* `==`.
"""
macro manifold_thing_forwards(T, Twhere, field::Symbol)
    return esc(
        quote
            allocate(p::$T) where {$Twhere} = $T(allocate(p.$field))
            allocate(p::$T, ::Type{P}) where {P,$Twhere} = $T(allocate(p.$field, P))
            function allocate(p::$T, ::Type{P}, dims::Tuple) where {P,$Twhere}
                return $T(allocate(p.$field, P, dims))
            end

            @inline Base.copy(p::$T) where {$Twhere} = $T(copy(p.$field))

            function Base.copyto!(q::$T, p::$T) where {$Twhere}
                copyto!(q.$field, p.$field)
                return q
            end

            number_eltype(p::$T) where {$Twhere} = typeof(one(eltype(p.$field)))

            Base.similar(p::$T) where {$Twhere} = $T(similar(p.$field))
            Base.similar(p::$T, ::Type{P}) where {P,$Twhere} = $T(similar(p.$field, P))

            Base.:(==)(p::$T, q::$T) where {$Twhere} = (p.$field == q.$field)
        end,
    )
end


"""
    default_manifold_fallbacks(TM::Symbol, TP::Symbol, TV::Symbol, pfield::Symbol, vfield::Symbol)

Introduce default fallbacks for all basic functions on manifolds, for manifold of type `TM`,
points of type `TP`, tangent vectors of type `TV`, with forwarding to fields `pfield` and
`vfield` for, respectively, point and tangent vector functions
"""
macro default_manifold_fallbacks(
    TM::Symbol,
    TP::Symbol,
    TV::Symbol,
    pfield::Symbol,
    vfield::Symbol,
)
    return esc(quote
        function angle(M::$TM, p::$TP, X::$TV, Y::$TV)
            return angle(M, p.$pfield, X.$vfield, Y.$vfield)
        end

        function check_point(M::$TM, p::$TP, X::$TV; kwargs...)
            return check_point(M, p.$pfield, X.$vfield; kwargs...)
        end

        function distance(M::$TM, p::$TP, q::$TP)
            return distance(M, p.$pfield, q.$pfield)
        end

        function embed!(M::$TM, q, p::$TP)
            return embed!(M, q, p.$pfield)
        end

        function embed!(M::$TM, Y, p::$TP, X::$TV)
            return embed!(M, Y, p.$pfield, X.$vfield)
        end

        function exp!(M::$TM, q::$TP, p::$TP, X::$TV)
            exp!(M, q.$pfield, p.$pfield, X.$vfield)
            return q
        end

        function inner(M::$TM, p::$TP, X::$TV, Y::$TV)
            return inner(M, p.$pfield, X.$vfield, Y.$vfield)
        end

        function isapprox(M::$TM, p::$TP, q::$TP; kwargs...)
            return isapprox(M, p.$pfield, q.$pfield; kwargs...)
        end

        function isapprox(M::$TM, p::$TP, X::$TV, Y::$TV; kwargs...)
            return isapprox(M, p.$pfield, X.$vfield, Y.$vfield; kwargs...)
        end

        function log!(M::$TM, X::$TV, p::$TP, q::$TP)
            log!(M, X.$vfield, p.$pfield, q.$pfield)
            return X
        end

    end)
end



@doc raw"""
    manifold_vector_forwards(T, Twhere, field::Symbol)

Introduce basic fallbacks for type `T` that represents vectors from a vector bundle for a
manifold. `Twhere` is put into `where` clause of each method. Fallbacks work by forwarding
to field passed as `field`.

List of forwarded functions:
* basic arithmetic (`*`, `/`, `\`, `+`, `-`),
* all things from [`@manifold_thing_forwards`](@ref),
* broadcasting support.

# example

    @eval @manifold_vector_forwards ValidationFibreVector{TType} TType value
"""
macro manifold_vector_forwards(T, Twhere, field::Symbol)
    return esc(
        quote
            Base.:*(X::$T, s::Number) where {$Twhere} = $T(X.$field * s)
            Base.:*(s::Number, X::$T) where {$Twhere} = $T(s * X.$field)
            Base.:/(X::$T, s::Number) where {$Twhere} = $T(X.$field / s)
            Base.:\(s::Number, X::$T) where {$Twhere} = $T(s \ X.$field)
            Base.:+(X::$T, Y::$T) where {$Twhere} = $T(X.$field + Y.$field)
            Base.:-(X::$T, Y::$T) where {$Twhere} = $T(X.$field - Y.$field)
            Base.:-(X::$T) where {$Twhere} = $T(-X.$field)
            Base.:+(X::$T) where {$Twhere} = $T(X.$field)

            @eval @manifold_thing_forwards $T $Twhere $field

            Base.axes(p::$T) where {$Twhere} = axes(p.$field)

            function Broadcast.BroadcastStyle(::Type{<:$T}) where {$Twhere}
                return Broadcast.Style{$T}()
            end
            function Broadcast.BroadcastStyle(
                ::Broadcast.AbstractArrayStyle{0},
                b::Broadcast.Style{$T},
            ) where {$Twhere}
                return b
            end

            function Broadcast.instantiate(
                bc::Broadcast.Broadcasted{Broadcast.Style{$T},Nothing},
            ) where {$Twhere}
                return bc
            end
            function Broadcast.instantiate(
                bc::Broadcast.Broadcasted{Broadcast.Style{$T}},
            ) where {$Twhere}
                Broadcast.check_broadcast_axes(bc.axes, bc.args...)
                return bc
            end

            Broadcast.broadcastable(X::$T) where {$Twhere} = X

            @inline function Base.copy(
                bc::Broadcast.Broadcasted{Broadcast.Style{$T}},
            ) where {$Twhere}
                return $T(Broadcast._broadcast_getindex(bc, 1))
            end

            Base.@propagate_inbounds function Broadcast._broadcast_getindex(
                X::$T,
                I,
            ) where {$Twhere}
                return X.$field
            end

            @inline function Base.copyto!(
                dest::$T,
                bc::Broadcast.Broadcasted{Broadcast.Style{$T}},
            ) where {$Twhere}
                axes(dest) == axes(bc) || Broadcast.throwdm(axes(dest), axes(bc))
                # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
                if bc.f === identity && bc.args isa Tuple{T} # only a single input argument to broadcast!
                    A = bc.args[1]
                    if axes(dest) == axes(A)
                        return copyto!(dest, A)
                    end
                end
                bc′ = Broadcast.preprocess(dest, bc)
                # Performance may vary depending on whether `@inbounds` is placed outside the
                # for loop or not. (cf. https://github.com/JuliaLang/julia/issues/38086)
                copyto!(dest.$field, bc′[1])
                return dest
            end
        end,
    )
end
