using JET, ManifoldsBase, Test
# load every package ManifoldsBase has an extension for, so that JET analyses the extensions as well
using CairoMakie, Plots, Quaternions, RecursiveArrayTools, Statistics

@testset "JET.jl" begin
    # only on released Julia versions, JET follows the compiler closely
    if isempty(VERSION.prerelease)
        JET.test_package(ManifoldsBase; target_modules = (ManifoldsBase,))
    end
end
