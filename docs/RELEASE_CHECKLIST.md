# Before tagging a scientific release

- [ ] Resolve the blocking findings in REVIEW.md.
- [ ] Select and document the reference implementation and supported compiler versions.
- [ ] Add a small redistributable input and saved expected outputs; run it from a clean checkout.
- [ ] Check FFT normalization, physical units, empty bins and triangle counts against a trusted reference.
- [ ] Check convergence under binning changes and document finite-grid effects.
- [ ] Record runtime, peak RAM, disk use and thread count for the example.
- [ ] Verify every notebook from a fresh kernel, with local paths replaced by user configuration.
- [ ] Agree on software and figure reuse rights, preserving third-party attribution.
- [ ] Run the root build, review compiler warnings, and archive the validation log.
- [ ] Create a versioned GitHub release and optionally archive that release on Zenodo.

Until these are completed, describe this as a research-code snapshot, not a validated turnkey reproduction package.
