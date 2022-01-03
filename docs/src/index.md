# ManifoldsBase.jl

`ManifoldsBase.jl` is a lightweight interface for manifolds.

This package provides an interface, so you probably either want to add it as a dependency to your project/package to work on manifold generically or implement a new manifold.
A package that (only) depends on `ManifoldsBase.jl`, see [Manopt.jl](https://manoptjl.org/stable/), which implements optimization algorithms on manifolds using this interface, i.e. they can be used with any manifold based on `ManifoldsBase.jl`. A library of manifolds implemented using this interface is provided see [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/).

Your package is using `ManifoldsBase`? Give us a note and we add you here.

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
