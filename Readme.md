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
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

[![ACM TOMS](https://img.shields.io/badge/ACM%20TOMS-10.1145%2F3618296-blue.svg)](http://doi.org/10.1145/3618296)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5964340.svg)](https://doi.org/10.5281/zenodo.5964340)

## Installation

In Julia you can install this package by typing

```julia
using Pkg; Pkg.add("ManifoldsBase")
```

in the Julia REPL.

Since this package provides an interface, you probably either want to add it as a dependency to your project/package to work on manifold generically or implement a new manifold.
For a package that (only) depends on `ManifoldsBase.jl`, see [Manopt.jl](https://manoptjl.org/stable/), which implements optimization algorithms on manifolds using this interface, i.e. they can be used with any manifold based on `ManifoldsBase.jl`. A library of manifolds implemented using this interface is provided by [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/).

Your package is using `ManifoldsBase`?
We would be very interested to hear where you are using the interface or manifolds in general! Give us a note and we add you here.

## Citation

If you use `ManifoldsBase.jl` in your work, please cite the following open access article, which covers both the basic interface as well as the performance for `Manifolds.jl`

> _Axen, S. D., Baran, M., Bergmann, R., Rzecki, K._ (2023).
> **Manifolds.jl: An Extensible Julia Framework for Data Analysis on Manifolds**,
> ACM Transactions on Mathematical Software, Volume 49, Issue 4, Article No. 33.<br/>
> doi: [10.1145/3618296](https://doi.org/10.1145/3618296),
> arXiv: [2106.08777](https://arxiv.org/abs/2106.08777)

<details>

<summary><code>AxenBaranBergmannRzecki:2023</code> (BibLaTeX) </summary>

```biblatex
@article{AxenBaranBergmannRzecki:2023,
    AUTHOR     = {Seth D. Axen and Mateusz Baran and Ronny Bergmann and Krzysztof Rzecki},
    ARTICLENO  = {33},
    DOI        = {10.1145/3618296},
    JOURNAL    = {ACM Transactions on Mathematical Software},
    MONTH      = {dec},
    NUMBER     = {4},
    TITLE      = {Manifolds.jl: An Extensible {J}ulia Framework for Data Analysis on Manifolds},
    VOLUME     = {49},
    YEAR       = {2023},
    EPRINT     = {2106.08777},
    EPRINTTYPE = {arXiv}
}
```

</details>

To refer to a certain version or the source code in general please cite for example

> _Axen, S. D., Baran, M., Bergmann, R._ (2026). **ManifoldsBase.jl**, Zenodo.<br/>
> doi: [10.5281/ZENODO.5964340](https://doi.org/10.5281/ZENODO.5964340)

<details>

<summary><code>manifoldsbasejl-zenodo-mostrecent</code> (BibLaTeX) </summary>

```biblatex
@software{manifoldsbasejl-zenodo-mostrecent,
    AUTHOR    = {Seth D. Axen and Mateusz Baran and Ronny Bergmann},
    TITLE     = {ManifoldsBase.jl},
    DOI       = {10.5281/ZENODO.5964340},
    URL       = {https://zenodo.org/record/5964340},
    PUBLISHER = {Zenodo},
    YEAR      = {2026},
    COPYRIGHT = {MIT License}
}
```

</details>

for the most recent version.
For a corresponding version specific DOI, see [the list of all versions](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%225964340%22&sort=-version&all_versions=True).
