# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive per-gene thermodynamic parameter plots
- Enhanced diagnostics with Tm3 and Tm5 visualizations
- Excel export functionality for best mutation probes

### Changed
- Modified probe generation to retain all candidates for diagnostics
- Improved coordinate conversion and strand handling

### Fixed
- Resolved issues with missing probe data for certain mutations
- Fixed thermal filter integration in probe generation pipeline

## [2.1.0] - 2025-09-05

### Added
- Enhanced probe generation to retain all probe candidates (including failed thermal filter)
- Per-gene thermodynamic parameter plots with `gene-start-ref-alt` naming
- Advanced visualizations with Tm3 and Tm5 specific plots
- Excel export with best mutation probes for each gene

### Fixed
- Coordinate conversion issues in mutation probe design
- Strand processing problems in `_extract_padlock_arms`

### Performance
- Reduced memory usage in probe generation
- Improved diagnostic generation performance

### Technical Details
- Modified `src/mutation_probe.py` - `_generate_probe_candidates()` function to retain all candidates
- Enhanced `probe_design_BZ09_TNBC_mut.py` - Added per-gene plotting functionality
- Updated diagnostic pipeline to handle all probe candidates
- Fixed coordinate conversion and strand handling in mutation probe design

## [2.0.0] - 2025-01-10

### Added
- Updated RNA Tm calculation algorithms in search_binding.py
- Improved thermodynamic parameter calculations for RNA sequences

### Changed
- Refactored thermal filter architecture
- Enhanced probe design pipeline

## [1.9.0] - 2025-01-06

### Removed
- TC linker sequences from PRISM barcode design

### Changed
- Simplified barcode structure for better probe performance

## [1.8.0] - 2024-09-13

### Added
- Constraint that Tm difference between left and right arms should be smaller than 5°C
- Enhanced thermal filter validation for padlock probe arms

### Changed
- Improved thermal filter validation criteria
