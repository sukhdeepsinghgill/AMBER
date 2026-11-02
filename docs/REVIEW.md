# AMBER: release-readiness review

Source snapshot reviewed 2026-11-02. No scientific source fixes have been silently applied. Full data-driven validation remains outstanding.

## Blocking findings

1. **Triangle loops start midway through the list.** `4_fbe/amber.c:191` begins at `rp=365`; `amber_v2.c:187` begins at `rp=635`. These are partial/resume runs. Expose the start/end indices as parameters; use zero for a full run only after checking the intended partitioning and output management. A successful exit may produce no triangle files if the start exceeds `rcount`.
2. **Invalid thread counts on small machines.** Both entry points use `Nthreads=omp_get_max_threads()/2 - 2` near line 53. This is non-positive for 1–5 available threads. Accept a user thread setting, validate it, and clamp to at least one. Do not treat changing OMP_NUM_THREADS to an arbitrarily large number as a scientific fix.
3. **C11 aligned allocation sizes are not rounded.** `amber.c` allocates `ncsize*sizeof(float)` and `ncsize*sizeof(int)` with 64-byte alignment, including per-thread buffers. At default `ncsize=767²`, the float allocation is 2,353,156 bytes, not a multiple of 64. C11 aligned_alloc requires the size to be an integral multiple of the alignment; implementations may reject it. Use checked `posix_memalign` or round size up with overflow checks. Inspect all aligned_alloc helpers and both variants, not just the first allocation.
4. **Missing inputs and hard-coded configuration.** The archived UVFITS file, gridded data and `uf_chans.txt` are absent; `input_params_ess` points to a placeholder path. `uf_nc=552`, `nc10=383`, the frequency metadata and ring list are dataset-specific. Validate them against the actual gridded header and flag list before any array indexing or frequency-triplet generation.
5. **Two output formats and duplicate Python modules.** v1 writes two double headers, float32 numerator and int counts; v2 writes one double header and two float32 arrays. Use the matching reader. Root and `analysis_codes/` helpers differ and expose overlapping names; the working directory can silently choose the wrong module. Consolidate or explicitly package/version them before a stable release.

## Specific pitfalls

- `uf_chans` indices are used directly; validate bounds, uniqueness and exact count. Determine this list from the same data/flagging policy as the gridding, not an unrelated observation.
- Enabling `collapse_key` averages/sums groups of four channels in the current code. Match channel count, frequency spacing, reference channel, valid-channel list and Python lag/cosmology configuration afterward.
- v1/v2 are alternative normalization paths, not sequential stages. Numerator and normalization products should be combined consistently when merging partitions; never blindly average final normalized arrays.
- The root Python `read_bsout_float_v2` discards the stored header for its returned tuple and sets `Ntri=1`; do not interpret that as a physical count.
- The Python intensity conversion's inline unit comment needs checking: compute the direction of Q and the applied divisions explicitly before labeling final axes. Match Jy/K or Jy/mK choices at every stage.
- Fixed `rc`, `rc_prime` and beam constants are not automatically inferred from a new observing frequency or cosmology.
- Fixed per-thread buffers and full gridded channel arrays can require many GiB. Record a realistic RAM budget; allocating all fields at archived dimensions is not a smoke test.
- Most data files use native binary layouts without version/dimension metadata. Add file-size and header checks; do not assume any FITS file or ViBEr binary is compatible.

## Build and workflow defects

The historical `amber_v2` make target does not match the explicit recipe name `mabs_final_v2`; an implicit rule may build it incorrectly. The `maps` recipe refers to absent `mabs_funcs.c`. The default aggregate target also includes cleanup. The new root Makefile provides correct explicit source lists for amber/amber_v2 and omits incomplete targets, without changing the science source.

Original paths include personal cluster/library locations, and stage flags force AVX2/fast-math. The new root build removes those machine assumptions. Original scripts remain historical. Notebooks retain external datasets, path assumptions and research-specific plotting sequences; clearing outputs does not make them reproducible.

## Acceptance test before release

Create a small channel/grid fixture with a known model bispectrum. Check frequency-triplet counts independently, including flagging and channel collapse; compare both normalizations on matching supported cells; validate the lag transform with a known Fourier pair; check k-axis ordering, temperature conversion, beam factors and missing-cell treatment. Include both no-flag and flagged examples with saved reference outputs and tolerances.
