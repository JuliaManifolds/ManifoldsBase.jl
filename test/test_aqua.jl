using Aqua, ManifoldsBase, Test

@testset "Aqua.jl" begin
    Aqua.test_all(
        ManifoldsBase;
        ambiguities = (exclude = [
            allocate, # has many possible call patterns that are not supported and ambiguous
            getindex, # ambiguous call patterns are not supported
            setindex!, # ambiguous call patterns are not supported
        ], broken = false),
    )
end
