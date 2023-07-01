#!/usr/bin/env julia
#
#

#
# (a) if docs is not the current active environment, switch to it
# (from https://github.com/JuliaIO/HDF5.jl/pull/1020/) 
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(PackageSpec(; path = (@__DIR__) * "/../"))
    Pkg.resolve()
    Pkg.instantiate()
end

# Debug Foo outside
tutorials_folder = (@__DIR__) * "/../tutorials"
# instantiate the tutorials environment if necessary
Pkg.activate(tutorials_folder)
Pkg.resolve()
Pkg.instantiate()
Pkg.build("IJulia") # build IJulia to the right version.
# Debug foo
using IJulia
@warn "$(readdir(IJulia.kerneldir(), join=true))"

Pkg.activate(@__DIR__) # but return to the docs one before
# (b) Did someone say render? Then we render!
if "--quarto" ∈ ARGS
    using CondaPkg
    CondaPkg.withenv() do
        @info "Rendering Quarto"
        run(`quarto render $(tutorials_folder)`)
        return nothing
    end
end

Pkg.activate(@__DIR__) # but return to the docs one before

using Documenter: DocMeta, HTML, MathJax3, deploydocs, makedocs
using ManifoldsBase

makedocs(
    format = HTML(;
        mathengine = MathJax3(),
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/favicon.ico"],
    ),
    modules = [ManifoldsBase],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename = "ManifoldsBase.jl",
    pages = [
        "Home" => "index.md",
        "How to define a manifold" => "tutorials/implement-a-manifold.md",
        "Design principles" => "design.md",
        "An abstract manifold" => "types.md",
        "Functions on maniolds" => [
            "Basic functions" => "functions.md",
            "Projections" => "projections.md",
            "Retractions" => "retractions.md",
            "Vector transports" => "vector_transports.md",
        ],
        "Manifolds" => "manifolds.md",
        "Decorating/Extending a Manifold" => "decorator.md",
        "Bases for tangent spaces" => "bases.md",
    ],
)
deploydocs(repo = "github.com/JuliaManifolds/ManifoldsBase.jl.git", push_preview = true)
