module ManifoldsBaseQuaternionsExt

if isdefined(Base, :get_extension)
    using ManifoldsBase
    using ManifoldsBase: ℍ, QuaternionNumbers
    using Quaternions
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..ManifoldsBase
    using ..ManifoldsBase: ℍ, QuaternionNumbers
    using ..Quaternions
end

@inline function ManifoldsBase.allocate_result_type(
    ::AbstractManifold{ℍ},
    f::TF,
    args::Tuple{},
) where {TF}
    return QuaternionF64
end

@inline function ManifoldsBase.coordinate_eltype(
    ::AbstractManifold,
    p,
    𝔽::QuaternionNumbers,
)
    return quat(number_eltype(p))
end

end
