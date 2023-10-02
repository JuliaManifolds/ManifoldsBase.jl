module ManifoldsBaseRecursiveArrayToolsExt

if isdefined(Base, :get_extension)
    using ManifoldsBase
    using RecursiveArrayTools
    using ManifoldsBase:
        AbstractBasis,
        FiberBundle,
        FiberBundleBasisData,
        FiberBundleProductVectorTransport,
        ProductBasisData,
        TangentSpaceType
    using ManifoldsBase: number_of_components
    using Random

    import ManifoldsBase:
        allocate,
        allocate_result,
        default_inverse_retraction_method,
        default_retraction_method,
        default_vector_transport_method,
        _get_dim_ranges,
        get_vector,
        get_vectors,
        inverse_retract,
        parallel_transport_direction,
        parallel_transport_to,
        project,
        _retract,
        submanifold_component,
        submanifold_components,
        vector_transport_direction,
        _vector_transport_direction,
        _vector_transport_to,
        vector_transport_to,
        ziptuples
    import Base: copyto!, getindex, setindex!, view
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..ManifoldsBase
    using ..RecursiveArrayTools
    using ..ManifoldsBase:
        AbstractBasis,
        FiberBundle,
        FiberBundleBasisData,
        FiberBundleProductVectorTransport,
        ProductBasisData,
        TangentSpaceType
    using ..ManifoldsBase: number_of_components
    using ..Random

    import ..ManifoldsBase:
        allocate,
        allocate_result,
        default_inverse_retraction_method,
        default_retraction_method,
        default_vector_transport_method,
        _get_dim_ranges,
        get_vector,
        get_vectors,
        inverse_retract,
        parallel_transport_direction,
        parallel_transport_to,
        project,
        _retract,
        submanifold_component,
        submanifold_components,
        vector_transport_direction,
        _vector_transport_direction,
        _vector_transport_to,
        vector_transport_to,
        ziptuples
    import ..Base: copyto!, getindex, setindex!, view
end

include("ProductManifoldRecursiveArrayToolsExt.jl")
include("FiberBundleRecursiveArrayToolsExt.jl")

end
