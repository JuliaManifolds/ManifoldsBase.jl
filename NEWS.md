# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.2.2] 02/12/2025

### Added

* `ApproximateExponentialRetraction` and `ApproximateLogarithmicInverseRetraction` as types
  to specify approximation methods for exp and log, which due to their approximation are also
  only approximate retractions and inverse retractions, respectively.

## [2.2.1] 13/11/2025

### Added

* `check_geodesic` function to numerically check whether geodesics have constant speed and are compatible with log and parallel transport.

## [2.2.0] 04/11/2025

### Fixed

* bugfixes on `get_embedding`, `get_embedding_type`, and `get_forwarding_type`, which now all accept as their second (third for the last case) argument consistently take the Type of a point (and no longer sometimes a point sometimes its type)
* `has_components` no longer propagates to the embedding

## [2.1.0] 30/10/2025

### Added

* an `StabilizedInverseRetraction` that improves numerical stability of another inverse retraction by projecting the resulting tangent vector onto its tangent space.

## [2.0.2] 23/10/2025

### Added

* Reference to `embed_project` in `project` documentation.
* Expanded docs of `inner`.

## [2.0.1] 17/10/2025

### Added

* `get_vector`, `get_coordinates`, `project` and their mutating variants now use traits for propagation.

## [2.0.0] 02/10/2025

While this release should be mostly backward compatible, especially when using defined manifolds from `Manifolds.jl`, some breaking changes were introduced.
To be precise, defining and using traits e.g. to dispatch functions to the embedding changed internally.
If you defined your own manifolds and used traits, please check the documentation of the new trait

### Added

* an interface for quotient manifolds.
  this also unified the naming a bit. Formerly `differential_canonical_project` is now `diff_canonical_project`.

### Changed

* refactor the trait system to no longer use a list of traits but single traits separately for the metric and the embedding specification
* Switch to using [Runic.jl](https://github.com/fredrikekre/Runic.jl) as code formatter

### Removed

* `ODEExponentialRetraction` was removed in favor of `solve_chart_exp_ode` implemented in `Manifolds.jl`.

## [1.2.0] 08/05/2025

### Added

* `tangent_vector_type` for converting point types to matching tangent vector types.

## [1.1.0] 29/04/2025

### Added

* `default_basis(M)` to be more flexible than a fixed `DefaultOrthonormalBasis` default.
* `StabilizedRetraction`, a retraction that improves numerical stability of another retraction by projecting the resulting point.

## [1.0.3] 08/04/2025

### Changed

* `VectorSpaceFiber` no longer requires the number system to be consistent with the wrapped manifold.

## [1.0.2] 07/04/2025

### Changed

* `Fiber` no longer requires the number system to be consistent with the wrapped manifold.

### Added

* `allocate` method that works with numeric scalars.

## [1.0.1] 05/02/2025

### Fixed

* An issue with allocation type promotion in `exp_fused`.

## [1.0] 05/02/2025

### Changed

* to avoid logical ambiguities to the forthcoming [`LieGroups.jl`](https://github.com/JuliaManifolds/LieGroups.jl),
  the “fusing” variant `exp(M, p, X, t)` has been moved to its own name `exp_fused(M, p, X, t)`
  and similarly `exp!(M, q, p, X, t)` has been moved to its own name `exp_fused!(M, q, p, X, t)`.
  Note that the new `exp_fused!` method is not exported and by default falls back to calling `exp!` with `t*X`.
  Actions to take
  * if you just implemented an own `exp(M, p, X)` or `exp!(M, q, p, X)` everything works as before.
  * if you implemented a fused variant `exp!(M, q, p, X, t)` you have to adapt two things
    1. move that implementation to `ManifoldsBase.exp_fused!(M, q, p, X, t)`
    2. Implement the default `exp!(M, q, p, X) = ManifoldBase.exp_fused!(M, q, p, one(eltype(p)), X)`,
    or an own specific implementation for the non-fused variant.
* Similar to `exp`, the “fusing” variant `retract(M, p, X, t, m)` has been moved to
  its own name `retract_fused(M, p, X, t, m)`  and similarly `retract!(M, q, p, X, t, m)`
  has been moved to its own name `retract_fused!(M, q, p, X, t, m)`.
  Note that the new `retract_fused!` method is not exported and by default falls back to calling `retract!` with `t*X`.
  Actions to take
  * if you just implemented an own `retract(M, p, X, m)` or `retract!(M, q, p, X, m)` everything works as before.
  * if you implemented a fused variant `retract!(M, q, p, X, t)` you have to adapt two things
    1. move that implementation to `ManifoldsBase.retract_fused!(M, q, p, X, t, m)`
    2. Implement the default `retract!(M, q, p, X, m) = ManifoldBase.retract_fused!(M, q, p, one(eltype(p)), X)`, or an own specific implementation for the non-fused variant.
* the `TVector` type has been renamed to `AbstractTangentVector`
* the `CoTVector` type has been renamed to `AbstractCotangentVector`

### Removed

* `parallel_transport_along(M, p, X, c)`, `vector_transport_along(M, p, X, c, m)` as well as
  their mutating variants are removed from the API for now.
  It was never specified how to actually specify a curve `c` and the method was only
  implemented for `Euclidean` in `Manifolds.jl`, where it is the identity.

## [0.15.24] 17/01/2025

### Added

* extended support for `allocate_on` with `ArrayPartition`.

## [0.15.23] 09/12/2024

### Added

* a field `point` to `ValidationFibreVector` to potentially store the point of the vector.
* a field `store_base_point` to `ValidationManifold` to indicate whether for new fibre vectors the base point should be stored.
* a keyword `ignore_contexts` to `ValidationManifold` to ignore certain contexts from validation, such as `:Input`, `:Output`, `:Point`, or `:Vector`.
* a keyword `ignore_functions` to `ValidationFibreVector` to ignore certain contexts within a single function. This is provided as a dictionary with the key being the (allocating) function and the value is a context or vector of contexts.

### Changed

* the internal function `array_value` was renamed to `internal_value` and is now exported, since it can be also used on elements that store values different from arrays,
e.g. a `ValidationMPoint` storing a subtype of a `ManifoldPoint`. `array_value` is hence deprecated.
* Minimum Julia version is now 1.10 (the LTS which replaced 1.6)

## [0.15.22] 15/11/2024

### Added

* `DefaultOrthonormalBasis()` is now the default basis for `get_vector`, `get_vector!`, `get_vectors`, `get_coordinates` and `get_coordinates!`.

## [0.15.21] 12/11/2024

### Fixed

* Coordinate allocation was improved to be more friendly with automatic differentiation.

## [0.15.20] 24/10/2024

### Changed

* `norm` function on `VectorSpaceFiber` (such as `TangentSpace`) now needs to be called without the point. The passed point was already ignored before.

## [0.15.19] 20/10/2024

### Changed

* make `has_components` introduced in the last version a decorator trait function.

## [0.15.18] 18/10/2024

### Added

* `distance(M, p, q, r)` to compute `r`-norms on manifolds that have components.
* `distance(M, p, q, m, r)` to compute (approximate) `r`-norms on manifolds that have components
  using an `AbstractInverseRetractionMethod m` within every (inner) distance call.
* `norm(M, p, X, r)` to compute `r`-norms on manifolds that have components.

## [0.15.17] 04/10/2024

### Changed

* **Mildly breaking**: the number system parameter now corresponds to the coefficients standing in front of basis vectors in a linear combination instead of components of a vector. For example, `DefaultOrthonormalBasis() == DefaultOrthonormalBasis(ℝ)` of `DefaultManifold(3, field=ℂ)` now has 6 vectors, and `DefaultOrthonormalBasis(ℂ)` of the same manifold has 3 basis vectors.

## [0.15.16] 13/09/2024

### Changed

* Adapt the traits, so that they also can be used when only `using ManifoldsBase`,
without importing internal `struct`s like `EmptyTrait` and `TraitList`

## [0.15.15] 29/08/2024

### Changed

* Refactored error message code when `ProductManifold` is used without `RecursiveArrayTools.jl`.

## [0.15.14] 27/08/2024

### Added

* A helpful error message when `ProductManifold` is used without `RecursiveArrayTools.jl`.

### Changed

* `representation_size` for `ProductManifold` now returns `nothing` instead of a one-element tuple. This change makes it easier to notice errors caused by not having `RecursiveArrayTools.jl` loaded.

## [0.15.13] 10/08/2024

### Changed

* Fixed a small bug that caused calling `get_vectors` on PowerManifolds
to sometimes cause an error, cf [#199](https://github.com/JuliaManifolds/ManifoldsBase.jl/issues/199).

## [0.15.12] 03/08/2024

### Changed

* Improved performance of power manifold creation and some cases of `get_component` on product manifold.

## [0.15.11] 28/07/2024

### Added

* Function `allocate_on` to generically allocate point and tangent vectors on a manifold without a pre-existing instance but of a particular type.
* Function `default_type` to get the default type of points and tangent vectors for a manifold.
* Package extension for the `Quaternions.jl` package that handles allocation.

### Changed

* Default allocation method was made more robust to custom promotion functions.

## [0.15.10] 19/05/2024

### Added

* Functions `fill(p, N)` and `fill!(P, p, N)` to fill values into a point on a power manifold `N`.
* introduce a `base_point(TpM)` to access the base point of a tangent space
* introduce `TpM[i]` to access tangent spaces of factors from an `AbstractPowerManifold` or a `ProductManifold`.


## [0.15.9] 02/05/2024

### Added

* Tests now also use `Aqua.jl` to spot problems in the code such as ambiguities.
* introduce a `check_inverse_retraction` function to numerically check whether an inverse retraction method is a (correct) inverse retraction.
* introduce a `check_retraction` function to numerically check whether a retraction method is a (correct) retraction.
* introduce a `check_vector_transport` function to numerically check whether a vector transport is a (correct) vector transport.

### Changed

* introduced a `ManifoldsBaseTestUtils` module to encapsulate common types and function definitions in different parts of the tests.

## [0.15.8] 13/03/2024

### Added

* `sectional_curvature` , `sectional_curvature_max` and `sectional_curvature_min` functions for obtaining information about sectional curvature of a manifold.

## [0.15.7] 24/01/2024

### Fixed

* `is_point` and `is_vector` can now more stably `:info` or `:warn` when they return false,
  since they emply `showerror` for these displays.

## [0.15.6] 15/12/2023

### Added

* An `AbstractApproximationMethod` to specify estimation methods for other more general functions,
as well as a `default_approximation_method` to specify defaults on manifolds.
* An `EmbeddedVectorTransport` to use a vector transport in the embedding and a final projection.

### Fixed

* `number_eltype` correctly returns scalar type for nested array types like `number_eltype(Vector{Vector{Float64}})`.

## [0.15.5] 13/12/2023

### Added

* Compatibility with `RecursiveArrayTools` v3.

## [0.15.4] 25/11/2023

### Fixed

* Fixed a bug reported in [Manopt#330](https://github.com/JuliaManifolds/Manopt.jl/issues/330).

## [0.15.3] 17/11/2023

### Fixed

- Pass kwargs in `rand!` for `AbstractPowerManifold` to appropriate methods on the wrapped manifold.

## [0.15.2] 8/11/2023

### Fixed

- `vee` and `hat` now use real coefficient basis for complex manifolds.

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
