# ManifoldsBase.jl

[![Build Status](https://travis-ci.org/JuliaNLSolvers/ManifoldsBase.jl.svg?branch=master)](https://travis-ci.org/JuliaNLSolvers/ManifoldsBase.jl/) [![codecov.io](http://codecov.io/github/JuliaNLSolvers/ManifoldsBase.jl/coverage.svg?branch=master)](https://codecov.io/gh/JuliaNLSolvers/ManifoldsBase.jl/)

Basic interface for manifolds in Julia.

This interface includes a simple `DefaultManifold`, which is a reduced version
of the [`Eucliean`](https://github.com/JuliaNLSolvers/Manifolds.jl/blob/master/src/Euclidean.jl)
manifold from [`Manifolds.jl`](https://github.com/JuliaNLSolvers/Manifolds.jl)
such that the interface functions can be tested.

The project [`Manifolds.jl`](https://github.com/JuliaNLSolvers/Manifolds.jl)
is based on this interface and provides a variety of manifolds.
