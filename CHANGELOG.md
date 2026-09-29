# Changelog

## [Unreleased]

## [v3.0.0] - 2026-04-28

### Changed
- **Breaking:** now requires MimiSSPs 2, which ships updated Benveniste et al. socioeconomic
  data. SSP5 is substantially revised (global population -30% by 2100 and -52% by 2200
  relative to MimiSSPs 1), so results for `SSP_scenario = "SSP585"` change materially.
  SSP2 is only marginally revised at the country level and nets out almost exactly in the
  global aggregate, so `SSP245` results move by less than 1e-5 relative. The default
  `socioeconomics_source = :RFF` pathway is unaffected.
- Widened compat bounds for Interpolations (0.16), JSON (1) and XLSX (0.11, 0.12).
- Regenerated the regression validation data in `test/validation_data/validation_data_current`.
- Moved the validation-data regeneration script from `test/save_validation_data.jl` to
  `scripts/save_validation_data.jl`, which now has its own `Project.toml`/`Manifest.toml`
  and depends on MimiGIVE through a `..` path. The `save_*` helpers moved out of the
  `TestFunctions` test module into `scripts/validation_helpers.jl`; the test suite keeps
  only the `validate_*` functions.

### Fixed
- Fix small typo in documentation on adding a new sector

### Added
- Update compat with MimiSSPs to incorporate an update to SSP5 population data
