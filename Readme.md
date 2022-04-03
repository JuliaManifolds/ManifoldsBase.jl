<div align="center">
    <img src="https://github.com/JuliaManifolds/ManifoldsBase.jl/blob/master/docs/src/assets/logo-text-readme.png" alt="ManifoldsBase.jl Logo with text" width="701">
</div>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamanifolds.github.io/ManifoldsBase.jl/dev/)
[![Build Status](https://travis-ci.org/JuliaManifolds/ManifoldsBase.jl.svg?branch=master)](https://travis-ci.org/JuliaManifolds/ManifoldsBase.jl/)
[![codecov.io](http://codecov.io/github/JuliaManifolds/ManifoldsBase.jl/coverage.svg?branch=master)](https://codecov.io/gh/JuliaManifolds/ManifoldsBase.jl/)
[![arXiv](https://img.shields.io/badge/arXiv%20CS.MS-2106.08777-blue.svg)](https://arxiv.org/abs/2106.08777)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5964340.svg)](https://doi.org/10.5281/zenodo.5964340)
## Installation

In Julia you can install this package by typing

```julia
] add ManifoldsBase
```

in the Julia REPL.

Since this package provides an interface, you probably either want to add it as a dependency to your project/package to work on manifold generically or implement a new manifold.
A package that (only) depends on `ManifoldsBase.jl`, see [Manopt.jl](https://manoptjl.org/stable/), which implements optimization algorithms on manifolds using this interface, i.e. they can be used with any manifold based on `ManifoldsBase.jl`. A library of manifolds implemented using this interface is provided see [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/).

Your package is using `ManifoldsBase`? Give us a note and we add you here.

## Citation

If you use `ManifoldsBase.jl` in your work, please cite the following

```biblatex
@online{2106.08777,
Author = {Seth D. Axen and Mateusz Baran and Ronny Bergmann and Krzysztof Rzecki},
Title = {Manifolds.jl: An Extensible Julia Framework for Data Analysis on Manifolds},
Year = {2021},
Eprint = {2106.08777},
Eprinttype = {arXiv},
}
```
