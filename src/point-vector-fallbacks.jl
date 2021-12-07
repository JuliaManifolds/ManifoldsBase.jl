
"""
    manifold_element_forwards(T, field::Symbol)
    manifold_element_forwards(T, Twhere, field::Symbol)

Introduce basic fallbacks for type `T` (which can be a subtype of `Twhere`) that represents
points or vectors for a manifold.
Fallbacks will work by forwarding to the field passed in `field``

List of forwarded functions:
* [`allocate`](@ref),
* [`copy`](@ref),
* [`copyto!`](@ref),
* [`number_eltype`](@ref) (only for values, not the type itself),
* `similar`,
* `==`.
"""
macro manifold_element_forwards(T, field::Symbol)
    return esc(quote
        @manifold_element_forwards ($T) _ ($field)
    end)
end
macro manifold_element_forwards(T, Twhere, field::Symbol)
    return esc(
        quote
            ManifoldsBase.allocate(p::$T) where {$Twhere} = $T(allocate(p.$field))
            function ManifoldsBase.allocate(p::$T, ::Type{P}) where {P,$Twhere}
                return $T(allocate(p.$field, P))
            end
            function ManifoldsBase.allocate(
                p::$T,
                ::Type{P},
                dims::Tuple,
            ) where {P,$Twhere}
                return $T(allocate(p.$field, P, dims))
            end

            @inline Base.copy(p::$T) where {$Twhere} = $T(copy(p.$field))

            function Base.copyto!(q::$T, p::$T) where {$Twhere}
                copyto!(q.$field, p.$field)
                return q
            end

            function ManifoldsBase.number_eltype(p::$T) where {$Twhere}
                return typeof(one(eltype(p.$field)))
            end

            Base.similar(p::$T) where {$Twhere} = $T(similar(p.$field))
            Base.similar(p::$T, ::Type{P}) where {P,$Twhere} = $T(similar(p.$field, P))

            Base.:(==)(p::$T, q::$T) where {$Twhere} = (p.$field == q.$field)
        end,
    )
end


"""
    default_manifold_fallbacks(TM, TP, TV, pfield::Symbol, vfield::Symbol)

Introduce default fallbacks for all basic functions on manifolds, for manifold of type `TM`,
points of type `TP`, tangent vectors of type `TV`, with forwarding to fields `pfield` and
`vfield` for point and tangent vector functions, respectively.
"""
macro default_manifold_fallbacks(TM, TP, TV, pfield::Symbol, vfield::Symbol)
    block = quote
        function ManifoldsBase.allocate_coordinates(M::$TM, p::$TP, T, n::Int)
            return ManifoldsBase.allocate_coordinates(M, p.$pfield, T, n)
        end

        function ManifoldsBase.angle(M::$TM, p::$TP, X::$TV, Y::$TV)
            return angle(M, p.$pfield, X.$vfield, Y.$vfield)
        end

        function ManifoldsBase.check_point(M::$TM, p::$TP; kwargs...)
            return check_point(M, p.$pfield; kwargs...)
        end

        function ManifoldsBase.check_vector(M::$TM, p::$TP, X::$TV; kwargs...)
            return check_vector(M, p.$pfield, X.$vfield; kwargs...)
        end

        function ManifoldsBase.distance(M::$TM, p::$TP, q::$TP)
            return distance(M, p.$pfield, q.$pfield)
        end

        function ManifoldsBase.embed!(M::$TM, q::$TP, p::$TP)
            return embed!(M, q.$pfield, p.$pfield)
        end

        function ManifoldsBase.embed!(M::$TM, Y::$TV, p::$TP, X::$TV)
            return embed!(M, Y.$vfield, p.$pfield, X.$vfield)
        end

        function ManifoldsBase.exp!(M::$TM, q::$TP, p::$TP, X::$TV)
            exp!(M, q.$pfield, p.$pfield, X.$vfield)
            return q
        end

        function ManifoldsBase.inner(M::$TM, p::$TP, X::$TV, Y::$TV)
            return inner(M, p.$pfield, X.$vfield, Y.$vfield)
        end

        function ManifoldsBase.inverse_retract!(
            M::$TM,
            X::$TV,
            p::$TP,
            q::$TP,
            m::AbstractInverseRetractionMethod,
        )
            inverse_retract!(M, X.$vfield, p.$pfield, q.$pfield, m)
            return X
        end
        function ManifoldsBase.inverse_retract!(
            M::$TM,
            X::$TV,
            p::$TP,
            q::$TP,
            m::LogarithmicInverseRetraction,
        )
            inverse_retract!(M, X.$vfield, p.$pfield, q.$pfield, m)
            return X
        end

        function ManifoldsBase.isapprox(M::$TM, p::$TP, q::$TP; kwargs...)
            return isapprox(M, p.$pfield, q.$pfield; kwargs...)
        end

        function ManifoldsBase.isapprox(M::$TM, p::$TP, X::$TV, Y::$TV; kwargs...)
            return isapprox(M, p.$pfield, X.$vfield, Y.$vfield; kwargs...)
        end

        function ManifoldsBase.allocate_result(::$TM, ::typeof(log), p::$TP, ::$TP)
            a = allocate(p.$vfield)
            return $TV(a)
        end
        function ManifoldsBase.allocate_result(
            ::$TM,
            ::typeof(inverse_retract),
            p::$TP,
            ::$TP,
        )
            a = allocate(p.$vfield)
            return $TV(a)
        end

        function ManifoldsBase.log!(M::$TM, X::$TV, p::$TP, q::$TP)
            log!(M, X.$vfield, p.$pfield, q.$pfield)
            return X
        end

        function ManifoldsBase.norm(M::$TM, p::$TP, X::$TV)
            return norm(M, p.$pfield, X.$vfield)
        end

        function ManifoldsBase.retract!(
            M::$TM,
            q::$TP,
            p::$TP,
            X::$TV,
            m::AbstractRetractionMethod,
        )
            retract!(M, q.$pfield, p.$pfield, X.$vfield, m)
            return X
        end
        function ManifoldsBase.retract!(
            M::$TM,
            q::$TP,
            p::$TP,
            X::$TV,
            m::ExponentialRetraction,
        )
            retract!(M, q.$pfield, p.$pfield, X.$vfield, m)
            return X
        end

        function ManifoldsBase.vector_transport_along!(M::$TM, Y::$TV, p::$TP, X::$TV, c)
            vector_transport_along!(M, Y.$vfield, p.$pfield, X.$vfield, c)
            return Y
        end
        function ManifoldsBase.vector_transport_direction!(
            M::$TM,
            Y::$TV,
            p::$TP,
            X::$TV,
            d::$TV,
            m,
        )
            vector_transport_direction!(M, Y.$vfield, p.$pfield, X.$vfield, d.$vfield, m)
            return Y
        end
        function ManifoldsBase.vector_transport_to!(
            M::$TM,
            Y::$TV,
            p::$TP,
            X::$TV,
            q::$TP,
            m,
        )
            vector_transport_to!(M, Y.$vfield, p.$pfield, X.$vfield, q.$pfield, m)
            return Y
        end

        function ManifoldsBase.zero_vector(M::$TM, p::$TP)
            return $TV(zero_vector(M, p.$pfield))
        end

        function ManifoldsBase.zero_vector!(M::$TM, X::$TV, p::$TP)
            zero_vector!(M, X.$vfield, p.$pfield)
            return X
        end
    end

    for BT in [
        ManifoldsBase.DISAMBIGUATION_BASIS_TYPES...,
        ManifoldsBase.DISAMBIGUATION_COTANGENT_BASIS_TYPES...,
    ]
        push!(
            block.args,
            quote
                function ManifoldsBase.get_coordinates!(M::$TM, Y, p::$TP, X::$TV, B::$BT)
                    return get_coordinates!(M, Y, p.$pfield, X.$vfield, B)
                end

                function ManifoldsBase.get_vector(M::$TM, p::$TP, X, B::$BT)
                    return $TV(get_vector(M, p.$pfield, X, B))
                end

                function ManifoldsBase.get_vector!(M::$TM, Y::$TV, p::$TP, X, B::$BT)
                    return get_vector!(M, Y.$vfield, p.$pfield, X, B)
                end
            end,
        )
    end

    for VTM in [ParallelTransport, VECTOR_TRANSPORT_DISAMBIGUATION...]
        push!(
            block.args,
            quote
                function ManifoldsBase.vector_transport_direction!(
                    M::$TM,
                    Y::$TV,
                    p::$TP,
                    X::$TV,
                    d::$TV,
                    m::$VTM,
                )
                    vector_transport_direction!(
                        M,
                        Y.$vfield,
                        p.$pfield,
                        X.$vfield,
                        d.$vfield,
                        m,
                    )
                    return Y
                end
                function ManifoldsBase.vector_transport_to!(
                    M::$TM,
                    Y::$TV,
                    p::$TP,
                    X::$TV,
                    q::$TP,
                    m::$VTM,
                )
                    vector_transport_to!(M, Y.$vfield, p.$pfield, X.$vfield, q.$pfield, m)
                    return Y
                end
            end,
        )
    end

    return esc(block)
end



@doc raw"""
    manifold_vector_forwards(T, field::Symbol)
    manifold_vector_forwards(T, Twhere, field::Symbol)

Introduce basic fallbacks for type `T` that represents vectors from a vector bundle for a
manifold. `Twhere` is put into `where` clause of each method. Fallbacks work by forwarding
to field passed as `field`.

List of forwarded functions:
* basic arithmetic (`*`, `/`, `\`, `+`, `-`),
* all things from [`@manifold_element_forwards`](@ref),
* broadcasting support.

# example

    @eval @manifold_vector_forwards ValidationFibreVector{TType} TType value
"""
macro manifold_vector_forwards(T, field::Symbol)
    return esc(quote
        @manifold_vector_forwards ($T) _ ($field)
    end)
end
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
            Base.zero(X::$T) where {$Twhere} = $T(zero(X.$field))

            @eval @manifold_element_forwards $T $Twhere $field

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
