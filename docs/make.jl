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
        "The manifold" => "manifold.md",
        "Basic functions" => "functions.md",
        "Meta manifolds" => "metamanifolds.md",
        "Extending Manifolds" => "decorator.md",
        "Bases for tangent spaces" => "bases.md",
    ],
)
deploydocs(repo = "github.com/JuliaManifolds/ManifoldsBase.jl.git", push_preview = true)
