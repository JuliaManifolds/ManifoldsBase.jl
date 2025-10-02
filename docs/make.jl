#!/usr/bin/env julia
#
#

if "--help" ∈ ARGS
    println(
        """
        docs/make.jl

        Render the `ManifoldsBase.jl` documentation with optional arguments

        Arguments
        * `--exclude-tutorials` - exclude the tutorials from the menu of Documenter,
          This can be used if not all tutorials are rendered and you want to therefore exclude links
          to these, especially the corresponding menu. This option should not be set on CI.
          Locally this is also set if `--quarto` is not set and not all tutorials are rendered.
        * `--help`              - print this help and exit without rendering the documentation
        * `--prettyurls`        – toggle the pretty urls part to true, which is always set on CI
        * `--quarto`            – (re)run the Quarto notebooks from the `tutorials/` folder before
          generating the documentation. If they are generated once they are cached accordingly.
          Then you can spare time in the rendering by not passing this argument.
          If quarto is not run, some tutorials are generated as empty files, since they
          are referenced from within the documentation.
        """,
    )
    exit(0)
end

run_quarto = "--quarto" in ARGS
run_on_CI = (get(ENV, "CI", nothing) == "true")
tutorials_in_menu = !("--exclude-tutorials" ∈ ARGS)
#
#
# (a) setup the tutorials menu – check whether all files exist
tutorials_menu =
    "How to define a manifold" => "tutorials/implement-a-manifold.md"
# Check whether all tutorials are rendered, issue a warning if not (and quarto if not set)
all_tutorials_exist = true
for (name, file) in [tutorials_menu] # Adapt when we actually have a full array
    fn = joinpath(@__DIR__, "src/", file)
    if !isfile(fn) || filesize(fn) == 0 # nonexistent or empty file
        global all_tutorials_exist = false
        if !run_quarto
            @warn "Tutorial $name does not exist at $fn."
            if (!isfile(fn)) && (endswith(file, "implement-a-manifold.md"))
                @warn "Generating empty file, since this tutorial is linked to from the documentation."
                touch(fn)
            end
        end
    end
end
if !all_tutorials_exist && !run_quarto && !run_on_CI
    @warn """
    Not all tutorials exist. Run `make.jl --quarto` to generate them. For this run they are excluded from the menu.
    """
    tutorials_in_menu = false
end
if !tutorials_in_menu
    @warn """
    You are either explicitly or implicitly excluding the tutorials from the documentation.
    You will not be able to see their menu entries nor their rendered pages.
    """
    run_on_CI &&
        (@error "On CI, the tutorials have to be either rendered with Quarto or be cached.")
end
#
# (b) if docs is not the current active environment, switch to it
# (from https://github.com/JuliaIO/HDF5.jl/pull/1020/) 
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
end

# (b) If quarto is set, or we are on CI, run quarto
if run_quarto || run_on_CI
    @info "Rendering Quarto"
    tutorials_folder = (@__DIR__) * "/../tutorials"
    # instantiate the tutorials environment if necessary
    Pkg.activate(tutorials_folder)
    # For a breaking release -> also set the tutorials folder to the most recent version
    Pkg.instantiate()
    Pkg.activate(@__DIR__) # but return to the docs one before
    run(`quarto render $(tutorials_folder)`)
end

# (d) load necessary packages for the docs
using Documenter, DocumenterCitations, DocumenterInterLinks
using ManifoldsBase

function add_links(line::String, url::String = "https://github.com/JuliaManifolds/Manopt.jl")
    # replace issues (#XXXX) -> ([#XXXX](url/issue/XXXX))
    while (m = match(r"\(\#([0-9]+)\)", line)) !== nothing
        id = m.captures[1]
        line = replace(line, m.match => "([#$id]($url/issues/$id))")
    end
    # replace ## [X.Y.Z] -> with a link to the release [X.Y.Z](url/releases/tag/vX.Y.Z)
    while (m = match(r"\#\# \[([0-9]+.[0-9]+.[0-9]+)\] (.*)", line)) !== nothing
        tag = m.captures[1]
        date = m.captures[2]
        line = replace(line, m.match => "## [$tag]($url/releases/tag/v$tag) ($date)")
    end
    return line
end

# (e) add NEWS.mdto docs
generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/JuliaManifolds/Manifolds.jl/blob/master/"
isdir(generated_path) || mkdir(generated_path)
for fname in ["NEWS.md"]
    open(joinpath(generated_path, fname), "w") do io
        # Point to source license file
        println(
            io,
            """
            ```@meta
            EditURL = "$(base_url)$(fname)"
            ```
            """,
        )
        # Write the contents out below the meta block
        for line in eachline(joinpath(dirname(@__DIR__), fname))
            println(io, add_links(line))
        end
    end
end

# (f) ...finally! make docs
bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :alpha)
links = InterLinks(
    "Manifolds" => ("https://juliamanifolds.github.io/Manifolds.jl/stable/"),
)
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
        (tutorials_in_menu ? [tutorials_menu] : [])...,
        "Design principles" => "design.md",
        "An abstract manifold" => "types.md",
        "Functions on maniolds" => [
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
        "Changelog" => "NEWS.md",
        "References" => "references.md",
    ],
    plugins = [bib, links],
)
deploydocs(repo = "github.com/JuliaManifolds/ManifoldsBase.jl.git", push_preview = true)
