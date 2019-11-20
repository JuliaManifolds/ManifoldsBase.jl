# ManifoldsBase.jl

[![Build Status](https://travis-ci.org/JuliaNLSolvers/ManifoldsBase.jl.svg?branch=master)](https://travis-ci.org/JuliaNLSolvers/ManifoldsBase.jl/) [![codecov.io](http://codecov.io/github/JuliaNLSolvers/ManifoldsBase.jl/coverage.svg?branch=master)](https://codecov.io/gh/JuliaNLSolvers/ManifoldsBase.jl/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://julianlsolvers.github.io/Manifolds.jl/latest/)


Basic interface for manifolds in Julia.

The project [`Manifolds.jl`](https://github.com/JuliaNLSolvers/Manifolds.jl)
is based on this interface and provides a variety of manifolds.

## Functions on a Manifold

## DefaultManifold

This interface includes a simple `DefaultManifold`, which is a reduced version
of the [`Eucliean`](https://github.com/JuliaNLSolvers/Manifolds.jl/blob/master/src/Euclidean.jl)
manifold from [`Manifolds.jl`](https://github.com/JuliaNLSolvers/Manifolds.jl)
such that the interface functions can be tested.

## Array Manifold

The `ArrayManifold` further illustrates, how one can also used types to
represent points on a manifold, tangent vectors as well as cotangent vectors,
where values are encapsulated in a certain type.

In general this might be used for manifolds, where these three types are represented
by more complex data structures or when it is neccssary to distinguish these
actually by type.

This adds a semantic layer to the interface and the default implementation of
`ArrayManifold` adds checks to all inputs and outputs of typed data.
