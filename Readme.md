<div align="center">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="https://github.com/JuliaManifolds/ManifoldsBase.jl/raw/master/docs/src/assets/logo-text-readme-dark.png">
      <img alt="ManifoldsBase.jl logo with text on the side" src="https://github.com/JuliaManifolds/ManifoldsBase.jl/raw/master/docs/src/assets/logo-text-readme.png">
    </picture>
</div>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamanifolds.github.io/ManifoldsBase.jl/dev/)
[![CI](https://github.com/JuliaManifolds/ManifoldsBase.jl/workflows/CI/badge.svg)](https://github.com/JuliaManifolds/ManifoldsBase.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![codecov](https://codecov.io/gh/JuliaManifolds/ManifoldsBase.jl/graph/badge.svg?token=bQsBUU9knE)](https://codecov.io/gh/JuliaManifolds/ManifoldsBase.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[![ACM TOMS](https://img.shields.io/badge/ACM%20TOMS-10.1145%2F3618296-blue.svg)](http://doi.org/10.1145/3618296)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5964340.svg)](https://doi.org/10.5281/zenodo.5964340)
## Installation

In Julia you can install this package by typing

```julia
using Pkg; Pkg.add("ManifoldsBase")
```

in the Julia REPL.

Since this package provides an interface, you probably either want to add it as a dependency to your project/package to work on manifold generically or implement a new manifold.
A package that (only) depends on `ManifoldsBase.jl`, see [Manopt.jl](https://manoptjl.org/stable/), which implements optimization algorithms on manifolds using this interface, i.e. they can be used with any manifold based on `ManifoldsBase.jl`. A library of manifolds implemented using this interface is provided see [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/).

Your package is using `ManifoldsBase`?
We would be very interested to hear where you are using the interface or manifolds in general! Give us a note and we add you here.

## Citation

If you use `ManifoldsBase.jl` in your work, please cite the following open access article

```biblatex
@article{AxenBaranBergmannRzecki:2023,
    author = {Axen, Seth D. and Baran, Mateusz and Bergmann, Ronny and Rzecki, Krzysztof},
    articleno = {33},
    doi = {10.1145/3618296},
    journal = {ACM Transactions on Mathematical Software},
    month = {dec},
    number = {4},
    title = {Manifolds.Jl: An Extensible Julia Framework for Data Analysis on Manifolds},
    volume = {49},
    year = {2023},
}
```


To refer to a certain version we recommend to also cite for example

```biblatex
@software{manifoldsbasejl-zenodo-mostrecent,
    AUTHOR    = {Seth D. Axen and Mateusz Baran and Ronny Bergmann},
    TITLE     = {ManifoldsBase.jl},
    DOI       = {10.5281/ZENODO.5964340},
    URL       = {https://zenodo.org/record/5964340},
    PUBLISHER = {Zenodo},
    YEAR      = {2022},
    COPYRIGHT = {MIT License}
}
```

Note that both citations are in [BibLaTeX](https://ctan.org/pkg/biblatex) format.
