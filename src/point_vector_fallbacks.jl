
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
* `size`,
* `==`.
"""
macro manifold_element_forwards(T, field::Symbol)
    return esc(quote
        ManifoldsBase.@manifold_element_forwards ($T) _unused ($field)
    end)
end
macro manifold_element_forwards(T, Twhere, field::Symbol)
    TWT = if Twhere === :_unused
        T
    else
        :($T where {$Twhere})
    end
    code = quote
        function ManifoldsBase.number_eltype(p::$TWT)
            return typeof(one(eltype(p.$field)))
        end

        Base.size(p::$TWT) = size(p.$field)
    end
    if Twhere === :_unused
        push!(
            code.args,
            quote
                ManifoldsBase.allocate(p::$T) = $T(allocate(p.$field))
                function ManifoldsBase.allocate(p::$T, ::Type{P}) where {P}
                    return $T(allocate(p.$field, P))
                end
                function ManifoldsBase.allocate(p::$T, ::Type{P}, dims::Tuple) where {P}
                    return $T(allocate(p.$field, P, dims))
                end
                @inline Base.copy(p::$T) = $T(copy(p.$field))
                function Base.copyto!(q::$T, p::$T)
                    copyto!(q.$field, p.$field)
                    return q
                end
                Base.similar(p::$T) = $T(similar(p.$field))
                Base.similar(p::$T, ::Type{P}) where {P} = $T(similar(p.$field, P))
                Base.:(==)(p::$T, q::$T) = (p.$field == q.$field)
            end,
        )
    else
        push!(
            code.args,
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
                Base.similar(p::$T) where {$Twhere} = $T(similar(p.$field))
                Base.similar(p::$T, ::Type{P}) where {P,$Twhere} = $T(similar(p.$field, P))
                Base.:(==)(p::$T, q::$T) where {$Twhere} = (p.$field == q.$field)
            end,
        )
    end

    return esc(code)
end


"""
    default_manifold_fallbacks(TM, TP, TV, pfield::Symbol, vfield::Symbol)

Introduce default fallbacks for all basic functions on manifolds, for manifold of type `TM`,
points of type `TP`, tangent vectors of type `TV`, with forwarding to fields `pfield` and
`vfield` for point and tangent vector functions, respectively.
"""
macro default_manifold_fallbacks(TM, TP, TV, pfield::Symbol, vfield::Symbol)
    block = quote
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
        function ManifoldsBase.allocate_coordinates(M::$TM, p::$TP, T, n::Int)
            return ManifoldsBase.allocate_coordinates(M, p.$pfield, T, n)
        end

        function ManifoldsBase.angle(M::$TM, p::$TP, X::$TV, Y::$TV)
            return ManifoldsBase.angle(M, p.$pfield, X.$vfield, Y.$vfield)
        end

        function ManifoldsBase.check_point(M::$TM, p::$TP; kwargs...)
            return ManifoldsBase.check_point(M, p.$pfield; kwargs...)
        end

        function ManifoldsBase.check_vector(M::$TM, p::$TP, X::$TV; kwargs...)
            return ManifoldsBase.check_vector(M, p.$pfield, X.$vfield; kwargs...)
        end

        function ManifoldsBase.check_approx(M::$TM, p::$TP, q::$TP; kwargs...)
            return ManifoldsBase.check_approx(M, p.$pfield, q.$pfield; kwargs...)
        end
        function ManifoldsBase.check_approx(M::$TM, p::$TP, X::$TV, Y::$TV; kwargs...)
            return ManifoldsBase.check_approx(M, p.$pfield, X.$vfield, Y.$vfield; kwargs...)
        end

        function ManifoldsBase.distance(M::$TM, p::$TP, q::$TP)
            return ManifoldsBase.distance(M, p.$pfield, q.$pfield)
        end

        function ManifoldsBase.embed!(M::$TM, q::$TP, p::$TP)
            ManifoldsBase.embed!(M, q.$pfield, p.$pfield)
            return q
        end

        function ManifoldsBase.embed!(M::$TM, Y::$TV, p::$TP, X::$TV)
            ManifoldsBase.embed!(M, Y.$vfield, p.$pfield, X.$vfield)
            return Y
        end

        function ManifoldsBase.exp!(M::$TM, q::$TP, p::$TP, X::$TV)
            ManifoldsBase.exp!(M, q.$pfield, p.$pfield, X.$vfield)
            return q
        end

        function ManifoldsBase.exp_fused!(M::$TM, q::$TP, p::$TP, X::$TV, t::Number)
            ManifoldsBase.exp_fused!(M, q.$pfield, p.$pfield, X.$vfield, t)
            return q
        end

        function ManifoldsBase.inner(M::$TM, p::$TP, X::$TV, Y::$TV)
            return ManifoldsBase.inner(M, p.$pfield, X.$vfield, Y.$vfield)
        end

        function ManifoldsBase.inverse_retract!(
            M::$TM,
            X::$TV,
            p::$TP,
            q::$TP,
            m::LogarithmicInverseRetraction,
        )
            ManifoldsBase.inverse_retract!(M, X.$vfield, p.$pfield, q.$pfield, m)
            return X
        end

        function ManifoldsBase.isapprox(M::$TM, p::$TP, q::$TP; kwargs...)
            return ManifoldsBase.isapprox(M, p.$pfield, q.$pfield; kwargs...)
        end

        function ManifoldsBase.isapprox(M::$TM, p::$TP, X::$TV, Y::$TV; kwargs...)
            return ManifoldsBase.isapprox(M, p.$pfield, X.$vfield, Y.$vfield; kwargs...)
        end

        function ManifoldsBase.log!(M::$TM, X::$TV, p::$TP, q::$TP)
            ManifoldsBase.log!(M, X.$vfield, p.$pfield, q.$pfield)
            return X
        end

        function ManifoldsBase.norm(M::$TM, p::$TP, X::$TV)
            return ManifoldsBase.norm(M, p.$pfield, X.$vfield)
        end

        function ManifoldsBase.retract!(
            M::$TM,
            q::$TP,
            p::$TP,
            X::$TV,
            m::ExponentialRetraction,
        )
            ManifoldsBase.retract!(M, q.$pfield, p.$pfield, X.$vfield, m)
            return q
        end

        function ManifoldsBase.retract_fused!(
            M::$TM,
            q::$TP,
            p::$TP,
            X::$TV,
            t::Number,
            m::ExponentialRetraction,
        )
            ManifoldsBase.retract_fused!(M, q.$pfield, p.$pfield, X.$vfield, t, m)
            return q
        end

        function ManifoldsBase.zero_vector(M::$TM, p::$TP)
            return $TV(zero_vector(M, p.$pfield))
        end

        function ManifoldsBase.zero_vector!(M::$TM, X::$TV, p::$TP)
            ManifoldsBase.zero_vector!(M, X.$vfield, p.$pfield)
            return X
        end
    end
    for f_postfix in [:default, :orthogonal, :orthonormal, :vee, :cached, :diagonalizing]
        ca = Symbol("get_coordinates_$(f_postfix)")
        cm = Symbol("get_coordinates_$(f_postfix)!")
        va = Symbol("get_vector_$(f_postfix)")
        vm = Symbol("get_vector_$(f_postfix)!")
        B_types = if f_postfix in [:default, :orthogonal, :orthonormal, :vee]
            [:AbstractNumbers, :RealNumbers, :ComplexNumbers]
        elseif f_postfix === :cached
            [:CachedBasis]
        elseif f_postfix === :diagonalizing
            [:DiagonalizingOrthonormalBasis]
        else
            [:Any]
        end
        for B_type in B_types
            push!(
                block.args,
                quote
                    function ManifoldsBase.$ca(M::$TM, p::$TP, X::$TV, B::$B_type)
                        return ManifoldsBase.$ca(M, p.$pfield, X.$vfield, B)
                    end
                    function ManifoldsBase.$cm(M::$TM, Y, p::$TP, X::$TV, B::$B_type)
                        ManifoldsBase.$cm(M, Y, p.$pfield, X.$vfield, B)
                        return Y
                    end
                    function ManifoldsBase.$va(M::$TM, p::$TP, X, B::$B_type)
                        return $TV(ManifoldsBase.$va(M, p.$pfield, X, B))
                    end
                    function ManifoldsBase.$vm(M::$TM, Y::$TV, p::$TP, X, B::$B_type)
                        ManifoldsBase.$vm(M, Y.$vfield, p.$pfield, X, B)
                        return Y
                    end
                end,
            )
        end
    end
    for f_postfix in [:polar, :project, :qr, :softmax, :sasaki]
        rm = Symbol("retract_$(f_postfix)!")
        push!(block.args, quote
            function ManifoldsBase.$rm(M::$TM, q, p::$TP, X::$TV)
                ManifoldsBase.$rm(M, q.$pfield, p.$pfield, X.$vfield)
                return q
            end
        end)
        rmf = Symbol("retract_$(f_postfix)_fused!")
        push!(
            block.args,
            quote
                function ManifoldsBase.$rmf(M::$TM, q, p::$TP, X::$TV, t::Number)
                    ManifoldsBase.$rmf(M, q.$pfield, p.$pfield, X.$vfield, t)
                    return q
                end
            end,
        )
    end
    push!(
        block.args,
        quote
            function ManifoldsBase.retract_exp_ode!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                m::AbstractRetractionMethod,
                B::ManifoldsBase.AbstractBasis,
            )
                ManifoldsBase.retract_exp_ode!(M, q.$pfield, p.$pfield, X.$vfield, m, B)
                return q
            end
            function ManifoldsBase.retract_exp_ode_fused!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                t::Number,
                m::AbstractRetractionMethod,
                B::ManifoldsBase.AbstractBasis,
            )
                ManifoldsBase.retract_exp_ode_fused!(
                    M,
                    q.$pfield,
                    p.$pfield,
                    X.$vfield,
                    t,
                    m,
                    B,
                )
                return q
            end
            function ManifoldsBase.retract_pade!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                m::PadeRetraction,
            )
                ManifoldsBase.retract_pade!(M, q.$pfield, p.$pfield, X.$vfield, m)
                return q
            end
            function ManifoldsBase.retract_pade_fused!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                t::Number,
                m::PadeRetraction,
            )
                ManifoldsBase.retract_pade_fused!(M, q.$pfield, p.$pfield, X.$vfield, t, m)
                return q
            end
            function ManifoldsBase.retract_embedded!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                m::AbstractRetractionMethod,
            )
                ManifoldsBase.retract_embedded!(M, q.$pfield, p.$pfield, X.$vfield, m)
                return q
            end
            function ManifoldsBase.retract_embedded_fused!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                t::Number,
                m::AbstractRetractionMethod,
            )
                ManifoldsBase.retract_embedded_fused!(
                    M,
                    q.$pfield,
                    p.$pfield,
                    X.$vfield,
                    t,
                    m,
                )
                return q
            end
            function ManifoldsBase.retract_sasaki!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                m::SasakiRetraction,
            )
                ManifoldsBase.retract_sasaki!(M, q.$pfield, p.$pfield, X.$vfield, m)
                return q
            end
            function ManifoldsBase.retract_sasaki_fused!(
                M::$TM,
                q::$TP,
                p::$TP,
                X::$TV,
                t::Number,
                m::SasakiRetraction,
            )
                ManifoldsBase.retract_sasaki_fused!(
                    M,
                    q.$pfield,
                    p.$pfield,
                    X.$vfield,
                    t,
                    m,
                )
                return q
            end
        end,
    )
    for f_postfix in [:polar, :project, :qr, :softmax]
        rm = Symbol("inverse_retract_$(f_postfix)!")
        push!(block.args, quote
            function ManifoldsBase.$rm(M::$TM, Y::$TV, p::$TP, q::$TP)
                ManifoldsBase.$rm(M, Y.$vfield, p.$pfield, q.$pfield)
                return Y
            end
        end)
    end
    push!(
        block.args,
        quote
            function ManifoldsBase.inverse_retract_embedded!(
                M::$TM,
                X::$TV,
                p::$TP,
                q::$TP,
                m::AbstractInverseRetractionMethod,
            )
                ManifoldsBase.inverse_retract_embedded!(
                    M,
                    X.$vfield,
                    p.$pfield,
                    q.$pfield,
                    m,
                )
                return X
            end

            function ManifoldsBase.inverse_retract_nlsolve!(
                M::$TM,
                X::$TV,
                p::$TP,
                q::$TP,
                m::NLSolveInverseRetraction,
            )
                ManifoldsBase.inverse_retract_nlsolve!(
                    M,
                    X.$vfield,
                    p.$pfield,
                    q.$pfield,
                    m,
                )
                return X
            end
        end,
    )
    # forward vector transports

    for sub in [:project, :diff, :embedded]
        # project & diff
        vttm = Symbol("vector_transport_to_$(sub)!")
        push!(
            block.args,
            quote
                function ManifoldsBase.$vttm(M::$TM, Y::$TV, p::$TP, X::$TV, q::$TP)
                    ManifoldsBase.$vttm(M, Y.$vfield, p.$pfield, X.$vfield, q.$pfield)
                    return Y
                end
            end,
        )
    end
    # parallel transports
    push!(
        block.args,
        quote
            function ManifoldsBase.parallel_transport_direction(
                M::$TM,
                p::$TP,
                X::$TV,
                d::$TV,
            )
                return $TV(
                    ManifoldsBase.parallel_transport_direction(
                        M,
                        p.$pfield,
                        X.$vfield,
                        d.$vfield,
                    ),
                )
            end
            function ManifoldsBase.parallel_transport_direction!(
                M::$TM,
                Y::$TV,
                p::$TP,
                X::$TV,
                d::$TV,
            )
                ManifoldsBase.parallel_transport_direction!(
                    M,
                    Y.$vfield,
                    p.$pfield,
                    X.$vfield,
                    d.$vfield,
                )
                return Y
            end
            function ManifoldsBase.parallel_transport_to(M::$TM, p::$TP, X::$TV, q::$TP)
                return $TV(
                    ManifoldsBase.parallel_transport_to(M, p.$pfield, X.$vfield, q.$pfield),
                )
            end
            function ManifoldsBase.parallel_transport_to!(
                M::$TM,
                Y::$TV,
                p::$TP,
                X::$TV,
                q::$TP,
            )
                ManifoldsBase.parallel_transport_to!(
                    M,
                    Y.$vfield,
                    p.$pfield,
                    X.$vfield,
                    q.$pfield,
                )
                return Y
            end
        end,
    )
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
        ManifoldsBase.@manifold_vector_forwards ($T) _unused ($field)
    end)
end
macro manifold_vector_forwards(T, Twhere, field::Symbol)
    TWT = if Twhere === :_unused
        T
    else
        :($T where {$Twhere})
    end
    code = quote
        @eval ManifoldsBase.@manifold_element_forwards $T $Twhere $field

        Base.axes(p::$TWT) = axes(p.$field)

        Broadcast.broadcastable(X::$TWT) = X

        Base.@propagate_inbounds function Broadcast._broadcast_getindex(X::$TWT, I)
            return X.$field
        end
    end

    broadcast_copyto_code = quote
        axes(dest) == axes(bc) || Broadcast.throwdm(axes(dest), axes(bc))
        # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
        if bc.f === identity && bc.args isa Tuple{$T} # only a single input argument to broadcast!
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

    if Twhere === :_unused
        push!(
            code.args,
            quote
                Base.:*(X::$T, s::Number) = $T(X.$field * s)
                Base.:*(s::Number, X::$T) = $T(s * X.$field)
                Base.:/(X::$T, s::Number) = $T(X.$field / s)
                Base.:\(s::Number, X::$T) = $T(s \ X.$field)
                Base.:-(X::$T) = $T(-X.$field)
                Base.:+(X::$T) = $T(X.$field)
                Base.zero(X::$T) = $T(zero(X.$field))
                Base.:+(X::$T, Y::$T) = $T(X.$field + Y.$field)
                Base.:-(X::$T, Y::$T) = $T(X.$field - Y.$field)

                function Broadcast.BroadcastStyle(::Type{<:$T})
                    return Broadcast.Style{$T}()
                end

                function Broadcast.BroadcastStyle(
                    ::Broadcast.AbstractArrayStyle{0},
                    b::Broadcast.Style{$T},
                )
                    return b
                end

                function Broadcast.instantiate(
                    bc::Broadcast.Broadcasted{Broadcast.Style{$T},Nothing},
                )
                    return bc
                end
                function Broadcast.instantiate(
                    bc::Broadcast.Broadcasted{Broadcast.Style{$T}},
                )
                    Broadcast.check_broadcast_axes(bc.axes, bc.args...)
                    return bc
                end

                @inline function Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{$T}})
                    return $T(Broadcast._broadcast_getindex(bc, 1))
                end

                @inline function Base.copyto!(
                    dest::$T,
                    bc::Broadcast.Broadcasted{Broadcast.Style{$T}},
                )
                    return $broadcast_copyto_code
                end
            end,
        )
    else
        push!(
            code.args,
            quote
                Base.:*(X::$T, s::Number) where {$Twhere} = $T(X.$field * s)
                Base.:*(s::Number, X::$T) where {$Twhere} = $T(s * X.$field)
                Base.:/(X::$T, s::Number) where {$Twhere} = $T(X.$field / s)
                Base.:\(s::Number, X::$T) where {$Twhere} = $T(s \ X.$field)
                Base.:-(X::$T) where {$Twhere} = $T(-X.$field)
                Base.:+(X::$T) where {$Twhere} = $T(X.$field)
                Base.zero(X::$T) where {$Twhere} = $T(zero(X.$field))
                Base.:+(X::$T, Y::$T) where {$Twhere} = $T(X.$field + Y.$field)
                Base.:-(X::$T, Y::$T) where {$Twhere} = $T(X.$field - Y.$field)

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

                @inline function Base.copy(
                    bc::Broadcast.Broadcasted{Broadcast.Style{$T}},
                ) where {$Twhere}
                    return $T(Broadcast._broadcast_getindex(bc, 1))
                end

                @inline function Base.copyto!(
                    dest::$T,
                    bc::Broadcast.Broadcasted{Broadcast.Style{$T}},
                ) where {$Twhere}
                    return $broadcast_copyto_code
                end
            end,
        )
    end
    return esc(code)
end
