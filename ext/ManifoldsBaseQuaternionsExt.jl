module ManifoldsBaseQuaternionsExt

if isdefined(Base, :get_extension)
    using ManifoldsBase
    using ManifoldsBase: ℍ
    using Quaternions
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..ManifoldsBase
    using ..ManifoldsBase: ℍ
    using ..Quaternions
end

@inline function ManifoldsBase.allocate_result_type(
    ::AbstractManifold{ℍ},
    f::TF,
    args::Tuple{},
) where {TF}
    return QuaternionF64
end

end
