# ManifoldsBase.jl

`ManifoldsBase.jl` is a lightweight interface for manifolds.

You can easily implement your algorithms and even your own manifolds just using the interface.
All manifolds from the package here are also based on this interface, so any project based on the interface can benefit from all manifolds, as soon as a certain manifold provides implementations of the functions a project requires.

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
