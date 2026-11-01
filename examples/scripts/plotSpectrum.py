#!/usr/bin/env python3
"""Plot the two FPESC line-of-sight differential spectra."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT_DIRECTORY = HERE.parent / "outputs"
SPECTRA = (OUTPUT_DIRECTORY / "Spectrum1.txt", OUTPUT_DIRECTORY / "Spectrum2.txt")
OUTPUT = OUTPUT_DIRECTORY / "Spectrums.png"
COLORS = ("#e31a1c", "#9677b8")
LINESTYLES = ("-", "--")


def main():
    plt.style.use(HERE / "scientific.mplstyle")
    figure, axes = plt.subplots()
    limits = []

    for index, (path, color, linestyle) in enumerate(
        zip(SPECTRA, COLORS, LINESTYLES), start=1
    ):
        data = np.loadtxt(path)
        energy, rate = data[:, 0], data[:, 2]
        nonzero = np.flatnonzero(rate > 0.0)
        limits.append((energy[nonzero[0]], energy[nonzero[-1]]))
        axes.plot(
            energy,
            rate,
            color=color,
            linestyle=linestyle,
            label=f"Spectrum {index}",
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
    figure.tight_layout()
    figure.savefig(OUTPUT, bbox_inches="tight", dpi=300)
    plt.close(figure)
    print(OUTPUT)


if __name__ == "__main__":
    main()
