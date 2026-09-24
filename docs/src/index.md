```@raw html
---
layout: home

hero:
  name: ManifoldsBase.jl
  text: An interface for manifolds in Julia
  tagline: Define your own Riemannian manifolds. Use abstract manifolds in your code.
  actions:
    - theme: brand
      text: Implement a manifold
      link: tutorials/implement-a-manifold/index.html
    - theme: alt
      text: Design principles
      link: /design/index.html
    - theme: alt
      text: Functions on manifolds
      link: /functions/index.html
  image:
    src: /logo.png            # primary image (light themes)
    dark: /logo-dark.png      # variant for dark themes
    alt: ManifoldsBase.jl     # accessibility text

features:
  - icon: 🪶
    title: Lightweight
    details: Based on an `AbstractManifold` this interface allows to define own manifolds easily. Furthermore, operations and algorithms on arbitrary manifolds can be defined using this interface package.
    link: /design/index.html
  - icon: ⚡️
    title: Efficient
    details: When possible, functions are available working in-place, like `exp!` or `log!` to reduce memory allocations.
    link: /functions/index.html
  - icon: 🧩
    title: Ecosystem
    details: Several Meta-Manifolds in this package but also further packages like [ManifoldDiff.jl](https://github.com/JuliaManifolds/ManifoldDiff.jl), [ManifoldDiffEq.jl](https://github.com/JuliaManifolds/ManifoldDiffEq.jl) or [ManifoldsGPU.jl](https://github.com/JuliaManifolds/ManifoldsGPU.jl) provide a whole ecosystem to work with manifolds.
  - icon:
        light: /logo-manifolds.png
        dark: /logo-manifolds-dark.png
        alt: Manifolds.jl
        wrap: true
    title: Manifolds.jl
    details: A comprehensive library of Riemannian manifolds is available in [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) to get you started directly on the most common manifolds.
  - icon:
        src: /logo-manopt.png
        alt: Manopt.jl
        wrap: true
    title: Manopt.jl
    details: Optimisation algorithms on any manifold following this interface are available in [Manopt.jl](https://manoptjl.org/stable/). Both smooth optimization, like gradient descent, quasi-Newton and nonsmooth, like proximal-gradient or other splitting based methods are available.
  - icon:
        light: /logo-liegroups.png
        dark: /logo-liegroups-dark.png
        alt: LieGroups.jl
        wrap: true
    title: LieGroups.jl
    details: "[LieGroups.jl](https://juliamanifolds.github.io/LieGroups.jl/stable/) extends the interface to Lie groups and uses [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) to define a library of Lie groups. Further features of this package are an abstract definition of a Lie algebra and group actions as well as generic Lie groups like the semi-direct product Lie group."
---
```

```@meta
CurrentModule = ManifoldsBase
```

```@docs
ManifoldsBase.ManifoldsBase
```

This package has two main purposes.
You can add it as a dependency if you plan to work on manifolds (generically) or if you plan to
define own manifolds in a package.
For a package that (only) depends on `ManifoldsBase.jl`, see [Manopt.jl](https://manoptjl.org/stable/),
which implements optimization algorithms on manifolds using this interface.
These optimisation algorithms can hence be used with any manifold implemented based on `ManifoldsBase.jl`.

For a library of manifolds implemented using this interface see [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/).

Your package is using `ManifoldsBase`?
We would like to add that here as well. Either [write an issue](https://github.com/JuliaManifolds/ManifoldsBase.jl/issues/new)
or add yourself by forking, editing this file and [opening a PR](https://github.com/JuliaManifolds/ManifoldsBase.jl/compare).

## Citation

If you use `ManifoldsBase.jl` in your work, please cite the following open access article [AxenBaranBergmannRzecki:2023](@cite),
which covers both the basic interface as well as the performance for `Manifolds.jl`

> _Axen, S. D., Baran, M., Bergmann, R., Rzecki, K._ (2023).
> **Manifolds.jl: An Extensible Julia Framework for Data Analysis on Manifolds**,
> ACM Transactions on Mathematical Software, Volume 49, Issue 4, Article No. 33.
>
> doi: [10.1145/3618296](https://doi.org/10.1145/3618296),
> arXiv: [2106.08777](https://arxiv.org/abs/2106.08777)

```@raw html
<details><summary><code>AxenBaranBergmannRzecki:2023</code> (BibLaTeX) </summary>
```
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
```@raw html
</details><br/>
```

To refer to a certain version or the source code in general please cite for example

> _Axen, S. D., Baran, M., Bergmann, R._ (2026). **ManifoldsBase.jl**, Zenodo.
>
> doi: [10.5281/ZENODO.5964340](https://doi.org/10.5281/ZENODO.5964340)

```@raw html
<details><summary><code>manifoldsbasejl-zenodo-mostrecent</code> (BibLaTeX) </summary>
```
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
```@raw html
</details> <br/>
```

for the most recent version.
For a corresponding version specific DOI, see [the list of all versions](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%225964340%22&sort=-version&all_versions=True).
