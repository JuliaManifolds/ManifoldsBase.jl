using ManifoldsBase

@testset "deprecated functions still working" begin
    x = [1.0, 2.0, 3.0]
    y = ValidationMPoint(x)
    @test (@test_deprecated ManifoldsBase.array_value(y)) == internal_value(y)
end
