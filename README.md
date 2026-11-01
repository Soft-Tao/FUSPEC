# FPESC

FPESC is a file-driven Fusion Product Energy Spectrum Calculator. Python is
used only to prepare a self-contained HDF5 task from gfile/TRANSP/NUBEAM/ENDF
inputs; the Monte Carlo calculation runs in the standalone `FPESC.jl` package.

## Layout

- `FPESC/`: Julia physics kernels, HDF5 task loader, and final-result writer.
- `scripts/prepare_fpesc_task.py`: prepare one self-contained input task.
- `scripts/run_fpesc_task.jl`: run one prepared task with Julia.
- `data/ENDF_crosssections/`: cross-section sources required by the preparer.

The removed legacy `fpesc/` Python package and JuliaCall interface are not part
of this workflow.

## Installation

Install the Python preparation dependencies:

```sh
python -m pip install -r requirements.txt
```

Install the Julia dependencies:

```sh
julia --project=FPESC -e 'using Pkg; Pkg.instantiate()'
```

## Prepare an input task

```sh
python scripts/prepare_fpesc_task.py \
  --input-dir data/inputs \
  --calculation-kind product-emission-spectrum \
  --selected-product neutron \
  --seed 20260815 \
  --output inputs.h5 \
  --resolution coarse \
  --angular-resolution medium
```

The task file contains the selected reaction, populations, cells, sampling
ranges, cross section, seed, and spectrum grid. It does not contain the Monte
Carlo sample count or Julia thread count.

Schema v2 selects the built-in `dd-n3he` or `dt-n4he` reaction with
`/reaction/name`. These tasks omit explicit particle and cross-section
definitions; Julia resolves them from the human-readable presets under
`FPESC/data/`. A task with `/reaction/name = "custom"` retains the explicit
particle and cross-section layout. The spectrum product is still selected by
`/reaction/selected_product` for emission-spectrum tasks.

The default cross-section model is the ENDF MF3 evaluation converted to the
centre-of-mass kinetic energy. Pass `--cross-section-model bosch-hale` to use
the Bosch-Hale total cross section instead; that is the parameterization
GENESIS uses and the one behind the published task-3/task-5 benchmark packages.
Either way the Legendre shape of the differential cross section comes from
ENDF MF6, and the selected model is recorded in `/provenance/task_config`.

## Run the Julia CLI

```sh
julia -t 4 --project=FPESC scripts/run_fpesc_task.jl \
  --input inputs.h5 \
  --samples 500000 \
  --output output.h5
```

`--samples` is the total sample count for the run. `julia -t N` selects the
thread and worker count. The CLI has no resume mode. The result HDF5 stores the
final normalized emission spectrum (or final reaction-rate profile) and no raw
accumulator/state datasets.

## Run the tests

```sh
julia -t 4 --project=FPESC -e 'using Pkg; Pkg.test()'
```

The suite pins the field-aligned emission convention, the `K_CM` cross-section
axis, the frozen output of both kernels and their worker-count invariance.
Start Julia with `-t N`, `N >= 2`, to exercise the parallel path.

## License and data attribution

The original FPESC software in this repository is released under the
[MIT License](LICENSE). The included ENDF-formatted nuclear data are sourced
from IAEA Nuclear Data Services and remain subject to their original terms;
see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
