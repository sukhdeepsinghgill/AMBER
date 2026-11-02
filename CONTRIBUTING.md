# Contributing to AMBER

Please open an issue describing the scientific or technical change before a major refactor. Include the source version, compiler, operating system, dependencies, command, parameter settings, input format and the smallest reproducer you can share.

For estimator changes, compare against an analytic model or a trusted direct triangle sum on a small grid; check normalization, units, bin edges, triangle counts and multipole conventions. A successful compilation is not a scientific validation. Keep expected results and tolerances with the example, and explain the reason for each tolerance.

Do not commit private observational data, generated binaries, FFT caches or notebook outputs. Preserve existing authorship notices. Use `make` from the repository root for the portable build; the original stage Makefiles are historical.

The publication checklist is in [docs/RELEASE_CHECKLIST.md](docs/RELEASE_CHECKLIST.md).
