"""
    ManifoldsBaseTestSuite.jl

A small module to provide common dummy types and defaults for testing.

! Note
    Whenever you rewrite tests maybe think about moving common definitions here.
"""
module ManifoldsBaseTestSuite
using ManifoldsBase, Test


# Validation Specials


struct CustomValidationManifoldRetraction <: ManifoldsBase.AbstractRetractionMethod end

function ManifoldsBase.injectivity_radius(
    ::ManifoldsBase.DefaultManifold,
    ::CustomValidationManifoldRetraction,
)
    return 10.0
end
function ManifoldsBase.injectivity_radius(
    ::ManifoldsBase.DefaultManifold,
    p,
    ::CustomValidationManifoldRetraction,
)
    return 11.0
end

struct ValidationDummyManifold <: ManifoldsBase.AbstractManifold{ℝ} end
ManifoldsBase.check_point(::ValidationDummyManifold, p) = nothing
ManifoldsBase.check_vector(::ValidationDummyManifold, p, X) = nothing
ManifoldsBase.project!(::ValidationDummyManifold, Y, p, X) = (Y .= 2 .* X)
ManifoldsBase.project!(::ValidationDummyManifold, q, p) = (q .= 2 .* p)
ManifoldsBase.distance(::ValidationDummyManifold, p, q) = -1.0
ManifoldsBase.norm(::ValidationDummyManifold, p, v) = -1.0
ManifoldsBase.get_embedding(::ValidationDummyManifold) = ManifoldsBase.DefaultManifold(3)
end
