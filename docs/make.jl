#!/usr/bin/env julia
#
#

if "--help" ∈ ARGS
    println(
        """
docs/make.jl

Render the `ManifoldsBase.jl` documentation with optional arguments

Arguments
* `--help`              - print this help and exit without rendering the documentation
* `--prettyurls`        – toggle the prettyurls part to true (which is otherwise only true on CI)
* `--quarto`            – run the Quarto notebooks from the `tutorials/` folder before generating the documentation
  this has to be run locally at least once for the `tutorials/*.md` files to exist that are included in
  the documentation (see `--exclude-tutorials`) for the alternative.
  If they are generated ones they are cached accordingly.
  Then you can spare time in the rendering by not passing this argument.
""",
    )
    exit(0)
end

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

# (b) Did someone say render? Then we render!
if "--quarto" ∈ ARGS
    using CondaPkg
    CondaPkg.withenv() do
        @info "Rendering Quarto"
        tutorials_folder = (@__DIR__) * "/../tutorials"
        # instantiate the tutorials environment if necessary
        Pkg.activate(tutorials_folder)
        Pkg.develop(PackageSpec(; path = (@__DIR__) * "/../"))
        Pkg.resolve()
        Pkg.instantiate()
        Pkg.build("IJulia") # build IJulia to the right version.
        Pkg.activate(@__DIR__) # but return to the docs one before
        return run(`quarto render $(tutorials_folder)`)
    end
end

using Documenter
using DocumenterCitations, DocumenterInterLinks
using ManifoldsBase

# (e) ...finally! make docs
bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :alpha)
links = InterLinks("Manifolds" => ("https://juliamanifolds.github.io/Manifolds.jl/stable/"))
makedocs(;
    # for development, we disable prettyurls
    format = Documenter.HTML(;
        prettyurls = (get(ENV, "CI", nothing) == "true") || ("--prettyurls" ∈ ARGS),
        assets = ["assets/favicon.ico", "assets/citations.css", "assets/link-icons.css"],
        size_threshold_warn = 200 * 2^10, # raise slightly from 100 to 200 KiB
        size_threshold = 300 * 2^10,      # raise slightly 200 to to 300 KiB
    ),
    modules = [ManifoldsBase],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename = "ManifoldsBase.jl",
    pages = [
        "Home" => "index.md",
        "How to define a manifold" => "tutorials/implement-a-manifold.md",
        "What are manifolds?" => "tutorials/what-are-manifolds2.md",
        "Design principles" => "design.md",
        "An abstract manifold" => "types.md",
        "Functions on manifolds" => [
            "Basic functions" => "functions.md",
            "Projections" => "projections.md",
            "Retractions" => "retractions.md",
            "Vector transports" => "vector_transports.md",
        ],
        "Manifolds" => "manifolds.md",
        "Meta-Manifolds" => "metamanifolds.md",
        "Decorating/Extending a Manifold" => "decorator.md",
        "Bases for tangent spaces" => "bases.md",
        "Numerical Verification" => "numerical_verification.md",
        "References" => "references.md",
    ],
    plugins = [bib, links],
)
deploydocs(repo = "github.com/JuliaManifolds/ManifoldsBase.jl.git", push_preview = true)
