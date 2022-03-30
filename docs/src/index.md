# ManifoldsBase.jl

`ManifoldsBase.jl` is a lightweight interface for manifolds.

This packages has two main purposes:
You add it as a dependency if you plan to work on manifolds (generically) or if you plan to
define own manifolds in a package.
For a package that (only) depends on `ManifoldsBase.jl`, see [Manopt.jl](https://manoptjl.org/stable/),
which implements optimization algorithms on manifolds using this interface.
These optimisation algorithms can hence be used with any manifold implemented based on `ManifoldsBase.jl`.

For a library of manifolds implemented using this interface [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/).

Your package is using `ManifoldsBase`?
We would like to add that here as well. Either [write an issue](https://github.com/JuliaManifolds/ManifoldsBase.jl/issues/new)
or add yourself by forking, editing this file and [opening a PR](https://github.com/JuliaManifolds/ManifoldsBase.jl/compare).

## Citation

If you use `ManifoldsBase.jl` in your work, please cite the following paper,
which covers both the basic interface as well as the performance for `Manifolds.jl`.

```biblatex
@online{2106.08777,
    Author = {Seth D. Axen and Mateusz Baran and Ronny Bergmann and Krzysztof Rzecki},
    Title = {Manifolds.jl: An Extensible Julia Framework for Data Analysis on Manifolds},
    Year = {2021},
    Eprint = {2106.08777},
    Eprinttype = {arXiv},
}
```

Note that the citation is in [BibLaTeX](https://ctan.org/pkg/biblatex) format.
