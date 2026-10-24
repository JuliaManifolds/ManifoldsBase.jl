"""
    @fieldwise_equality T
    @fieldwise_equality T (f1, f2, ...)
    @fieldwise_equality T Twhere (f1, f2, ...)

Define `==` and `hash` for objects of type `T`, which compare fields by value with `==`
and hash them, respectively.

Types that store mutable data, like a matrix of a metric or an array of basis vectors,
or that wrap other objects which do, would otherwise be compared by identity with `===`.

The first form compares all fields of `T`, determined at runtime.
The other forms compare only the fields `f1, f2, ...` given in a tuple, which yields
type stable methods with the comparisons and the hash written out, and should be preferred
where performance matters.
The type variables in `Twhere`, either a single one or a tuple, appear in `T` and are type
parameters that are not determined by the fields: they have to agree for two objects to be
equal and they are included in the hash. Type variables that stand for types should be given
with their upper bound, like `TB <: AbstractBasis`, to avoid method ambiguities. For example
`@fieldwise_equality CachedBasis{𝔽, TB} (𝔽, TB <: AbstractBasis) (data,)` defines
`==(B::CachedBasis{𝔽, TB}, C::CachedBasis{𝔽, TB}) where {𝔽, TB <: AbstractBasis}`.
"""
macro fieldwise_equality(T)
    return esc(
        quote
            Base.:(==)(a::$T, b::$T) = ManifoldsBase._fieldwise_equal(a, b)
            Base.hash(a::$T, h::UInt) = ManifoldsBase._fieldwise_hash(a, h)
        end
    )
end
macro fieldwise_equality(T, fields::Expr)
    return _fieldwise_equality_code(T, (), fields)
end
macro fieldwise_equality(T, Twhere, fields::Expr)
    params = Twhere isa Symbol ? (Twhere,) : Tuple(Twhere.args)
    return _fieldwise_equality_code(T, params, fields)
end

function _fieldwise_equality_code(T, params, fields::Expr)
    fields.head === :tuple || throw(ArgumentError("The fields $(fields) have to be given as a tuple."))
    names = fields.args
    eq = foldl(
        (ex, f) -> :($ex && (a.$f == b.$f)), names[2:end];
        init = :(a.$(names[1]) == b.$(names[1])),
    )
    # the names of the type variables, possibly given with an upper bound
    param_names = map(p -> (p isa Expr && p.head === :<:) ? p.args[1] : p, params)
    hashed = (param_names..., (:(a.$f) for f in names)...)
    hx = foldl((ex, x) -> :(hash($x, $ex)), hashed; init = :(hash(nameof(typeof(a)), h)))
    eq_signature = :(Base.:(==)(a::$T, b::$T))
    hash_signature = :(Base.hash(a::$T, h::UInt))
    if !isempty(params)
        eq_signature = Expr(:where, eq_signature, params...)
        hash_signature = Expr(:where, hash_signature, params...)
    end
    return esc(
        Expr(
            :block,
            Expr(:function, eq_signature, Expr(:block, :(return $eq))),
            Expr(:function, hash_signature, Expr(:block, :(return $hx))),
        )
    )
end

_fieldwise_equal(a, b) = all(name -> getfield(a, name) == getfield(b, name), fieldnames(typeof(a)))

function _fieldwise_hash(a, h::UInt)
    fields = map(name -> getfield(a, name), fieldnames(typeof(a)))
    return foldr(hash, (nameof(typeof(a)), fields...); init = h)
end

# manifolds that wrap other manifolds or store data
@fieldwise_equality ConnectionManifold (manifold, connection)
@fieldwise_equality EmbeddedManifold (manifold, embedding)
@fieldwise_equality Fiber{𝔽} 𝔽 (manifold, point, fiber_type)
@fieldwise_equality MetricManifold (manifold, metric)
@fieldwise_equality PowerManifold{𝔽, <:AbstractManifold, <:Any, TPR} (
    𝔽, TPR <: AbstractPowerRepresentation,
) (manifold, size)
@fieldwise_equality ProductManifold (manifolds,)
@fieldwise_equality ValidationManifold

# bases and their data
@fieldwise_equality CachedBasis{𝔽, TB} (𝔽, TB <: AbstractBasis) (data,)
@fieldwise_equality DiagonalizingOrthonormalBasis{𝔽} 𝔽 (frame_direction,)
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
