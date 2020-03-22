# ManifoldsBase.jl

[![Build Status](https://travis-ci.org/JuliaManifolds/ManifoldsBase.jl.svg?branch=master)](https://travis-ci.org/JuliaManifolds/ManifoldsBase.jl/) [![codecov.io](http://codecov.io/github/JuliaManifolds/ManifoldsBase.jl/coverage.svg?branch=master)](https://codecov.io/gh/JuliaManifolds/ManifoldsBase.jl/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://julianlsolvers.github.io/Manifolds.jl/latest/)

Basic interface for manifolds in Julia.

The project [`Manifolds.jl`](https://github.com/JuliaManifolds/Manifolds.jl)
is based on this interface and provides a variety of manifolds.

## `DecoratorManifold`

The decorator manifold enhances a manifold by certain, in most cases implicitly
assumed to have a standard case, properties, see for example the `EmbeddedManifold`.
The decorator acts semi transparently, i.e. `:transparent` for all functions not affected by that
decorator and `:intransparent` otherwise. Another possibility is, that the decorator just
passes to `:parent` in order to fill default values.

## `DefaultManifold`

This interface includes a simple `DefaultManifold`, which is a reduced version
of the [`Euclidean`](https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/Euclidean.jl)
manifold from [`Manifolds.jl`](https://github.com/JuliaManifolds/Manifolds.jl),
such that the interface functions can be tested.

## Embedded Manifold

The embedded manifold models the embedding of a manifold into another manifold.
This way a manifold can benefit from existing implementations.
One example is the `TransparentIsometricEmbeddingType` where a manifold uses the metric,
`inner`, from its embedding.

## `ArrayManifold`

The `ArrayManifold` further illustrates how one can also used types to
represent points on a manifold, tangent vectors, and cotangent vectors,
where values are encapsulated in a certain type.

In general, `ArrayManifold` might be used for manifolds where these three types are represented
by more complicated data structures or when it is necessary to distinguish these
by type.

This adds a semantic layer to the interface, and the default implementation of
`ArrayManifold` adds checks to all inputs and outputs of typed data.
