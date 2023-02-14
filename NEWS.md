# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.14.0]

### Added

- Type restriction for `t` in scaled retractions relaxed to `Number`.
- `embed_project(M::AbstractManifold, p)` and `embed_project(M::AbstractManifold, p, X)` that are like projections onto a manifold or a tangent space but are guaranteed to be idempotent.

### Changed

- Retractions for scaled vectors no longer dispatch to non-scaled retractions. It is now reversed for performance reasons. Please either just define `exp!(::MyManifold, q, p, X, t::Number)` or both this and `exp!(::MyManifold, q, p, X)`.
- `DefaultManifold` now stores size in a field instead of the type itself to reduce the amount of compilation needed.
- Fixed typo in `inverse_retract_caley` (now `inverse_retract_cayley`) and `retract_caley` (now `retract_cayley`).
- `retract_pade` and `retract_pade!` now receive `PadeRetraction` objects instead of just `n`.
