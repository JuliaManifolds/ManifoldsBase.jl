"""
    @fieldwise_equality T

Define `==` and `hash` for objects of type `T`, which compare all fields by value with `==`
and hash all fields, respectively.

Types that store mutable data, like a matrix of a metric or an array of basis vectors,
or that wrap other objects which do, would otherwise be compared by identity with `===`.
Type parameters that are not determined by the fields have to be compared separately.
"""
macro fieldwise_equality(T)
    return esc(
        quote
            Base.:(==)(a::$T, b::$T) = ManifoldsBase._fieldwise_equal(a, b)
            Base.hash(a::$T, h::UInt) = ManifoldsBase._fieldwise_hash(a, h)
        end
    )
end

function _fieldwise_equal(a, b)
    return all(name -> getfield(a, name) == getfield(b, name), fieldnames(typeof(a)))
end

function _fieldwise_hash(a, h::UInt)
    fields = map(name -> getfield(a, name), fieldnames(typeof(a)))
    return foldr(hash, (nameof(typeof(a)), fields...); init = h)
end

# manifolds that wrap other manifolds or store data
@fieldwise_equality ConnectionManifold
@fieldwise_equality EmbeddedManifold
@fieldwise_equality MetricManifold
@fieldwise_equality ProductManifold
@fieldwise_equality ValidationManifold

function Base.:(==)(
        M::PowerManifold{𝔽, <:AbstractManifold, <:Any, TPR},
        N::PowerManifold{𝔽, <:AbstractManifold, <:Any, TPR},
    ) where {𝔽, TPR}
    return _fieldwise_equal(M, N)
end
function Base.hash(M::PowerManifold{𝔽, <:AbstractManifold, <:Any, TPR}, h::UInt) where {𝔽, TPR}
    return _fieldwise_hash(M, hash(TPR, h))
end

# bases and their data
function Base.:(==)(
        B::DiagonalizingOrthonormalBasis{𝔽}, C::DiagonalizingOrthonormalBasis{𝔽}
    ) where {𝔽}
    return B.frame_direction == C.frame_direction
end
function Base.hash(B::DiagonalizingOrthonormalBasis{𝔽}, h::UInt) where {𝔽}
    return _fieldwise_hash(B, hash(𝔽, h))
end
function Base.:(==)(B::CachedBasis{𝔽, TB}, C::CachedBasis{𝔽, TB}) where {𝔽, TB}
    return B.data == C.data
end
function Base.hash(B::CachedBasis{𝔽, TB}, h::UInt) where {𝔽, TB}
    return _fieldwise_hash(B, hash(TB, h))
end
@fieldwise_equality DiagonalizingBasisData
@fieldwise_equality PowerBasisData
@fieldwise_equality ProductBasisData

# vectors with respect to a basis
@fieldwise_equality FVector

# retractions, inverse retractions, vector transports and approximation methods
# that wrap other methods or store data
@fieldwise_equality ApproximateExponentialRetraction
@fieldwise_equality ApproximateLogarithmicInverseRetraction
@fieldwise_equality DifferentiatedRetractionVectorTransport
@fieldwise_equality EmbeddedInverseRetraction
@fieldwise_equality EmbeddedRetraction
@fieldwise_equality EmbeddedVectorTransport
@fieldwise_equality ExtrinsicEstimation
@fieldwise_equality InverseProductRetraction
@fieldwise_equality InverseRetractionWithKeywords
@fieldwise_equality NLSolveInverseRetraction
@fieldwise_equality PoleLadderTransport
@fieldwise_equality ProductRetraction
@fieldwise_equality ProductVectorTransport
@fieldwise_equality RetractionWithKeywords
@fieldwise_equality ScaledVectorTransport
@fieldwise_equality SchildsLadderTransport
@fieldwise_equality ShootingInverseRetraction
@fieldwise_equality StabilizedInverseRetraction
@fieldwise_equality StabilizedRetraction
@fieldwise_equality VectorTransportDirection
@fieldwise_equality VectorTransportTo
@fieldwise_equality VectorTransportWithKeywords
