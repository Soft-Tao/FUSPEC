# How to run FUSPEC

Run these commands from the repository root.

## Build the task packages

Build the beam-thermal and thermal-thermal D-D tasks from the gFile and NPZ profiles:

```sh
python examples/scripts/TaskPreparation.py \
  --gfile examples/inputs/gFile \
  --thermal examples/inputs/ThermalIonDistribution.npz \
  --fast examples/inputs/FastIonDistribution.npz
```

## Run the Beam-thermal task

Calculate the fast-D plus thermal-D product distribution with 1e6 samples:

```sh
julia -t 4 --project=FUSPEC examples/scripts/RunTask.jl \
  --input examples/inputs/TaskPackage1.h5 \
  --output examples/outputs/ResultPackage1.h5 \
  --samples 1000000
```

## Run the thermal-thermal task

Calculate the thermal-D plus thermal-D product distribution with 1e6 samples:

```sh
julia -t 4 --project=FUSPEC examples/scripts/RunTask.jl \
  --input examples/inputs/TaskPackage2.h5 \
  --output examples/outputs/ResultPackage2.h5 \
  --samples 1000000
```

## Generate the line-of-sight spectra

Project both product distributions onto the LOS file:

```sh
python examples/scripts/PostProcessing.py
```

## Plot the spectra

Compare the generated spectra with the reference spectra:

```sh
python examples/scripts/PlotSpectrum.py
```

The bundled `Spectrum1_ref.txt` and `Spectrum2_ref.txt` were generated with
1e8 samples. The commands above use 1e6 samples to reduce runtime, so their
total rates and spectra may show visible Monte Carlo noise. To assess
convergence, rerun them with larger `--samples` values, such as `1e7` and
`1e8`, and compare the generated spectra with the reference spectra.
