# FUSPEC.jl

FUSPEC.jl is the standalone Julia physics and Monte Carlo implementation used by
the file-driven Fusion Product Energy Spectrum Calculator workflow. It loads a
self-contained HDF5 task prepared by Python and runs without Python or
JuliaCall.

The package defines particles, two-body reactions, validated run
configuration, counter-based random numbers, Maxwellian, BiMaxwellian and
tabulated distributions, cubic total and tensor-product differential cross
sections, relativistic kinematics, collective velocities, a serial DT
ReactionRate kernel, and a directional TotalSpectrum kernel. Both
kernels support bounded-workspace deterministic parallel execution with
bit-for-bit stable reduction across worker counts.

Task schema v2 supports built-in `dd-n3he` and `dt-n4he` reactions. Their
particles, Bosch-Hale total cross sections, and normalized DROSG Legendre
angular coefficients are stored in `data/*.toml` and resolved by independent
per-reaction loaders before problem construction. Explicit custom reactions
remain available and support either p3 or p4 spectra.

HDF5.jl is a direct dependency used by the validated prepared-task loader and
final-result writer. Result files contain only the normalized emission spectrum
or final reaction-rate profile; raw accumulators and continuation state are not
written. `EmissionSpectrum` keeps only raw accumulators in memory; use
`normalized_spectrum(result)` when an in-memory normalized array is explicitly
needed. That call materializes a full-size copy, so prefer the streaming
`write_result` path for production runs.

Install dependencies from the repository root:

```sh
julia --project=FUSPEC -e 'using Pkg; Pkg.instantiate()'
```

Run a prepared task with four threads and 500,000 total samples:

```sh
julia -t 4 --project=FUSPEC examples/scripts/RunTask.jl \
  --input inputs.h5 \
  --samples 500000 \
  --output output.h5
```

The task seed is read from the input bundle. The sample count is supplied only
by `--samples`, and the worker count is selected only by `julia -t N`. The CLI
does not provide a resume mode.
