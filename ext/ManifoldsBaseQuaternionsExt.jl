module ManifoldsBaseQuaternionsExt

using ManifoldsBase
using ManifoldsBase: ℍ, QuaternionNumbers
using Quaternions

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
