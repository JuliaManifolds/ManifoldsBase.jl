# ManifoldsBase.jl

`ManifoldsBase.jl` is a lightweight interface for manifolds.

This packages has two main purposes.
You can add it as a dependency if you plan to work on manifolds (generically) or if you plan to
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
@article{AxenBaranBergmannRzecki:2023,
    AUTHOR    = {Axen, Seth D. and Baran, Mateusz and Bergmann, Ronny and Rzecki, Krzysztof},
    ARTICLENO = {33},
    DOI       = {10.1145/3618296},
    JOURNAL   = {ACM Transactions on Mathematical Software},
    MONTH     = {dec},
    NUMBER    = {4},
    TITLE     = {Manifolds.Jl: An Extensible Julia Framework for Data Analysis on Manifolds},
    VOLUME    = {49},
    YEAR      = {2023}
}
```

Note that the citation is in [BibLaTeX](https://ctan.org/pkg/biblatex) format.
