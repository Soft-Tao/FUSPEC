# FUSPEC

![Julia](https://img.shields.io/badge/Julia-1.12-blue)
![Python](https://img.shields.io/badge/Python-3-blue)
![License](https://img.shields.io/github/license/Soft-Tao/FUSPEC)
![Version](https://img.shields.io/badge/version-1.0.1-blue)
![Last commit](https://img.shields.io/github/last-commit/Soft-Tao/FUSPEC?branch=main)

FUSPEC is a standalone, open-source Monte Carlo code for calculating
anisotropic fusion-product birth distributions in space, energy and pitch
angle from arbitrary reactant distributions. It applies fully relativistic
two-body reaction kinematics with angle-dependent center-of-mass differential
cross sections. Its protocol-driven architecture separates task preparation,
Monte Carlo calculation and post-processing; retaining the full phase-space
distribution allows one sampling run to provide synthetic spectra for multiple
diagnostic lines of sight and source terms for energetic-particle transport.
An efficient mode also provides reaction-rate profiles.

## Layout

- `FUSPEC/`: standalone Julia package containing the physics models and Monte
  Carlo kernels.
  - `FUSPEC/src/`: reaction kinematics, cross-section models, sampling, task
    I/O and result writing.
  - `FUSPEC/data/`: built-in D-D and D-T reaction presets, particle
    definitions and cross-section parameters.
  - `FUSPEC/test/`: package tests.
- `examples/inputs/`: example task packages, a gfile and line-of-sight data.
- `examples/outputs/`: reference spectra and example result packages.
- `examples/scripts/`: helper scripts for running tasks, post-processing and
  plotting.
- `examples/How-to-run-FUSPEC.md`: complete example workflow and expected
  outputs.
- `examples/scripts/TaskPreparation.py`: build the example task packages from a gFile
  and the thermal- and fast-ion NPZ distributions.
- `TaskPackageProtocol.md`: HDF5 task and result package specification.
- `THIRD_PARTY_NOTICES.md`: attribution for cross-section models and other
  third-party materials.

## Installation

FUSPEC requires Julia 1.12 or later. Install Julia from [here](https://julialang.org/downloads/). Install the package dependencies from the
repository root:

```sh
git clone git@github.com:Soft-Tao/FUSPEC.git
cd FUSPEC
julia --project=FUSPEC -e 'using Pkg; Pkg.instantiate()'
```

This installs the Julia dependencies declared in `FUSPEC/Project.toml`,
including HDF5.jl. Check that the package loads correctly with:

```sh
julia --project=FUSPEC -e 'using FUSPEC'
```

The Python task-preparation and plotting utilities are optional. Install
their dependencies from the requirements file:

```sh
python -m pip install -r examples/requirements.txt
```

## License and cross-section attribution

The original FUSPEC software is released under the
[MIT License](LICENSE). The built-in cross-section use the Bosch-Hale analytic
parameterization, while their angle-dependent center-of-mass differential
cross sections use the DROSG normalized Legendre expansion. The corresponding
references and attribution are provided in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
