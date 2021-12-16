using ManifoldsBase, Documenter

makedocs(
    format=Documenter.HTML(prettyurls=false, assets=["assets/favicon.ico"]),
    modules=[ManifoldsBase],
    authors="Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename="ManifoldsBase.jl",
    pages=[
        "Home" => "index.md",
        "How to write a Manifold" => "example.md",
        "Design Principles" => "design.md",
        "Basic functions" => "functions.md",
        "Manifolds" => "manifolds.md",
        "Extending Manifolds" => "decorator.md",
        "Bases for tangent spaces" => "bases.md",
    ],
)
deploydocs(repo="github.com/JuliaManifolds/ManifoldsBase.jl.git", push_preview=true)
