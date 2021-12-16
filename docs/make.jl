using ManifoldsBase, Documenter

makedocs(
    format=Documenter.HTML(prettyurls=false, assets=["assets/favicon.ico"]),
    modules=[ManifoldsBase],
    authors="Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename="ManifoldsBase.jl",
    pages=[
        "Home" => "index.md",
        "Basic functions" => "basicfunctions.md"
    ],
)
deploydocs(repo="github.com/JuliaManifolds/ManifoldsBase.jl.git", push_preview=true)
