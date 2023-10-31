# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.15.1] 30/10/2023

### Added

- `zero_vector(TpM)` to generate a zero vector in the tangent space
- a GitHub CI action that errors, when this file was not updated on a PR

### Fixed

- `is_point` and `is_vector` for the tangent space now correctly forward to
  vector checks on the corresponding manifold. The same for both `check_size`s
- add `[compat]` entries for the standard libraries.

## [0.15.0] 21/10/2023

### Added

- `ProductManifold` type was migrated from Manifolds.jl.
- `Fiber`, `VectorSpaceFiber` and `TangentSpace` types. `TangentSpace` is a generalized version of `TangentSpaceAtPoint` from Manifolds.jl.
- A keyword to `ValidationManifold` which `error=` mode to use.
  This is by default the previous `:error` mode.
- `change_representer!`, `change_metric!` and `Weingarten!` methods added to `PowerManifold`.
- `×` now also works for retractions, inverse retractions, and vector transports to create their product versions
- `retract`, `inverse_retract`, and `vector_transport_to` (and `_dir`) now also accept arbirtrary retractions on the product manifold. These act the same as the n-fold product of a retraction.

### Changed

- `retract` now behaves like `exp` in the sense that it allocates early,
  which reduces the amount of code to dispatch through levels 1-3 twice
- `inverse_retract` now behaves like `log` in the sense that it allocates early
- `Requires.jl` is added as a dependency to facilitate loading some methods related to `ProductManifolds` on Julia 1.6 to 1.8. Later versions rely on package extensions.
- `Documenter.jl` was updated to 1.0.
- `PowerManifold` can now store its size either in a field or in a type, similarly to `DefaultManifold`. By default the size is stored in a field.
- The signature of `is_point` was changed to be consistent with `isapprox.`.
  The error positional symbol (third argument) is now a keyword argument.
  We left the boolean shortcut in place.
  That means
  * `is_point(M, p, true)` works the same as before (`false` is the default anyways)
  * `is_point(M, p, :warn)` has to be changed to `is_point(M, p; error=:warn)`
- The signature of `is_vector` was changed to be consistent with `isapprox` and `is_point`.
  The error positional symbol (fourth argument) is now a keyword argument.
  The error positional boolean (fourth argument) hence moved to fifth place (after `check_base_point`)
  This means
  * `is_vector(M, p, X, true)` should now be `is_vector(M, p, X; error=:error)`
  * `is_vector(M, p, X, err, base)` for two booleans `err, base` should now be `is_vector(M, p, X, base, err)`
  * `is_vector(M, p, X, err, base)` for a symbol `err` should now be `is_vector(M, p, X, base; error=err)`

### Removed

- Julia 1.0 is no longer supported. From now on, the earliest supported Julia version is 1.6.

## [0.14.12] 23/09/2023

### Changed

- Introduce a thorough way to allocate tangent vectors for `rand`

## [0.14.11] 25/08/2023

### Added

- Make the `Weingarten` map a decorator capable function.

## [0.14.10] 17/08/2023

### Added

- introduce the `Weingarten` map and its in place variant `Weingarten!`.

## [0.14.9] 03/08/2023

### Added

- Introduce an interface that allows the static size of a manifold to be a field as well.

## [0.14.8] 07/07/2023

### Changed

- Improve `show` for cached bases and make it more robust

## [0.14.7] 07/07/2023

### Changed

- the tutorial is now written in Quarto.

## [0.14.6] 10/06/2023

### Added

- export the inplace random function `rand!`

## [0.14.5] 03/05/2023

### Added

- Allow to specify an `AbstractManifold` when converting points or tangent vector types.

## [0.14.4] 10/04/2023

### Changed

- Fix `copy` to work properly when copying `Number`s

## [0.14.3] 16/03/2023

### Changed

- Fix an allocation bug in nested power manifolds

## [0.14.2] 16/03/2023

### Added

- adds a DependaBot workflow.

### Changed

- Fix an allocation issue with `exp(M, p, X, t)` that did not respect the type of `t`.

## [0.14.1] 18/02/2023

Note that this release did not trigger a TagBot, so it appears within 0.14.2 in the tagged/created releases

### Added

- Introduce `change_representer` already in `ManifoldsBase`.

## [0.14.0] – 15/02/2023

### Added

- Type restriction for `t` in scaled retractions relaxed to `Number`.
- `embed_project(M::AbstractManifold, p)` and `embed_project(M::AbstractManifold, p, X)` that are like projections onto a manifold or a tangent space but are guaranteed to be idempotent.

### Changed

- Retractions for scaled vectors no longer dispatch to non-scaled retractions. It is now reversed for performance reasons. Please either just define `exp!(::MyManifold, q, p, X, t::Number)` or both this and `exp!(::MyManifold, q, p, X)`.
- `DefaultManifold` now stores size in a field instead of the type itself to reduce the amount of compilation needed.
- Fixed typo in `inverse_retract_caley` (now `inverse_retract_cayley`) and `retract_caley` (now `retract_cayley`).
- `retract_pade` and `retract_pade!` now receive `PadeRetraction` objects instead of just `n`.
