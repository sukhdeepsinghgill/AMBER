# What changed in this repository draft

Scientific `.c`, `.h` and `.py` sources are preserved byte-for-byte from the supplied archive. No estimator results have been regenerated or certified by packaging.

Added a repository-root portable Makefile, documentation, paper figures with attribution, citation metadata, contribution guidance and ignore rules. BARBIE's original root makefile is preserved as `docs/legacy_makefile.txt`; the visibility packages retain their legacy stage Makefiles. The root build intentionally excludes incomplete and historical targets.

Removed bundled ELF binaries and `.o` files from the release copy. Removed ViBEr's generated `modes` file. ViBEr notebooks found only in checkpoint folders were rescued into `legacy_notebooks/`, with their stage prefix retained in the filename. Notebook code and markdown cells are retained; stored outputs and execution counts were cleared. AMBER notebook outputs were also cleared. Original uploads remain the archival copy of those outputs.

The source hash mapping is in `source_manifest.json`. Historical scripts remain for reference and must not be treated as supported quick starts.

A build-only GitHub Actions workflow is supplied. Its commands were tested locally, but the hosted workflow itself has not run; it makes no scientific-validation claim.
