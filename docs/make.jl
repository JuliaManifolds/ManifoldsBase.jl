using ManifoldsBase, Documenter

makedocs(
    format = Documenter.HTML(prettyurls = false, assets = ["assets/favicon.ico"]),
    modules = [ManifoldsBase],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename = "ManifoldsBase.jl",
    pages = [
        "Home" => "index.md",
        "How to write a manifold" => "example.md",
        "Design principles" => "design.md",
        "The manifold type" => "manifold_type.md",
        "Functions on Maniolds" => [
            "Basic functions" => "functions.md",
            "Projections" => "projections.md",
            "Retractions" => "retractions.md",
            "Vector Transports" => "vector_transports.md",
        ],
        "Manifolds" => "manifolds.md",
        "Decorating/Extending a Manifold" => "decorator.md",
        "Bases for tangent spaces" => "bases.md",
    ],
)
deploydocs(repo = "github.com/JuliaManifolds/ManifoldsBase.jl.git", push_preview = true)
