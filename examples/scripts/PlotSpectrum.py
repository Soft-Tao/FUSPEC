#!/usr/bin/env python3
"""Compare both current task spectra with their reference spectra."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT_DIRECTORY = HERE.parent / "outputs"
OUTPUT = OUTPUT_DIRECTORY / "Spectrums.png"
SERIES = (
    (
        OUTPUT_DIRECTORY / "Spectrum1.txt",
        "TaskPackage1 custom",
        "#e31a1c",
        "-",
    ),
    (
        OUTPUT_DIRECTORY / "Spectrum1_ref.txt",
        "Spectrum1_ref",
        "#9677b8",
        "--",
    ),
    (
        OUTPUT_DIRECTORY / "Spectrum2.txt",
        "TaskPackage2 custom",
        "#fdbf6f",
        ":",
    ),
    (
        OUTPUT_DIRECTORY / "Spectrum2_ref.txt",
        "Spectrum2_ref",
        "#33a02c",
        "-.",
    ),
)


def main():
    plt.style.use(HERE / "scientific.mplstyle")
    figure, axes = plt.subplots(figsize=(9, 6), constrained_layout=True)
    limits = []

    for path, label, color, linestyle in SERIES:
        data = np.loadtxt(path)
        energy, rate = data[:, 0], data[:, 2]
        nonzero = np.flatnonzero(rate > 0.0)
        limits.append((energy[nonzero[0]], energy[nonzero[-1]]))
        axes.plot(
            energy,
            rate,
            color=color,
            linestyle=linestyle,
            linewidth=3,
            label=label,
        )

    left = min(bound[0] for bound in limits)
    right = max(bound[1] for bound in limits)
    margin = 0.05 * (right - left)
    axes.set_xlim(left - margin, right + margin)
    axes.set_yscale("log")
    axes.set_ylim(bottom=1.0e-1)
    axes.set_xlabel("Neutron energy (MeV)")
    axes.set_ylabel(r"Differential rate (s$^{-1}$ MeV$^{-1}$)")
    axes.legend()
    figure.savefig(OUTPUT, bbox_inches="tight", dpi=300)
    plt.close(figure)
    print(OUTPUT)


if __name__ == "__main__":
    main()
