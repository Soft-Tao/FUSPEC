#!/usr/bin/env python3
"""Standalone FPESC emission-map to line-of-sight spectrum converter."""

from pathlib import Path

import h5py
import numpy as np


# Default batch inputs. Edit this job list when processing other result files.
SCRIPT_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIRECTORY.parents[1]
GFILE = REPOSITORY_ROOT / "examples/inputs/gFile"
LOS_FILE = REPOSITORY_ROOT / "examples/inputs/LineOfSight.file"
RESULT_DIRECTORY = REPOSITORY_ROOT / "examples/outputs"
JOBS = (
    (RESULT_DIRECTORY / "ResultPackage1.h5", RESULT_DIRECTORY / "Spectrum1.txt"),
    (RESULT_DIRECTORY / "ResultPackage2.h5", RESULT_DIRECTORY / "Spectrum2.txt"),
)

MEV = 1.602176634e-13


def _fixed_values(lines, line_index, count, width):
    """Read ``count`` fixed-width Fortran values from consecutive lines."""
    values = []
    while len(values) < count:
        if line_index >= len(lines):
            raise ValueError("Unexpected end of gfile")
        line = lines[line_index]
        line_index += 1
        for start in range(0, len(line), width):
            field = line[start:start + width].strip()
            if field:
                values.append(float(field))
        if not line.strip():
            raise ValueError("Empty line in gfile numeric data")
    return values[:count], line_index


def read_gfile(path):
    """Read only the magnetic axis and LCFS boundary needed by this script."""
    lines = Path(path).read_text(encoding="ascii").splitlines()
    if not lines:
        raise ValueError(f"Empty gfile: {path}")
    header_fields = lines[0].rsplit(maxsplit=3)
    if len(header_fields) != 4:
        raise ValueError("Invalid gfile header")
    nx, ny = int(header_fields[2]), int(header_fields[3])

    first, line_index = _fixed_values(lines, 1, 20, 16)
    r_axis, z_axis = first[5], first[6]

    # Skip fpol, pressure, FF', p', psi, and qpsi.
    _, line_index = _fixed_values(lines, line_index, 4 * nx, 16)
    _, line_index = _fixed_values(lines, line_index, nx * ny, 16)
    _, line_index = _fixed_values(lines, line_index, nx, 16)

    if line_index >= len(lines):
        raise ValueError("Missing gfile boundary counts")
    count_line = lines[line_index]
    line_index += 1
    n_boundary = int(count_line[:5])
    if n_boundary <= 0:
        raise ValueError("Gfile has no LCFS boundary")

    boundary, _ = _fixed_values(lines, line_index, 2 * n_boundary, 16)
    boundary = np.asarray(boundary, dtype=np.float64)
    return r_axis, z_axis, boundary[0::2], boundary[1::2]


def read_los(path):
    """Read x, y, z, C, and the normalized ray direction."""
    rows = []
    with Path(path).open(encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line or line.startswith(("#", "/")):
                continue
            fields = line.split()
            if len(fields) != 8:
                continue
            try:
                values = [float(field) for field in fields]
            except ValueError:
                continue
            if not np.all(np.isfinite(values)):
                raise ValueError(f"Invalid LOS record: {line}")
            if values[3] < 0.0 or values[4] <= 0.0:
                raise ValueError(f"Invalid LOS coefficient: {line}")
            rows.append(values)
    if not rows:
        raise ValueError(f"No valid LOS records in {path}")

    records = np.asarray(rows, dtype=np.float64)
    direction = records[:, 5:8]
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    return (
        records[:, 0],
        records[:, 1],
        records[:, 2],
        records[:, 3],
        direction,
    )


def load_map(path):
    """Open a FPESC product-emission-spectrum HDF5 file."""
    file = h5py.File(path, "r")
    try:
        kind = file["run/calculation_kind"].asstr()[()]
        if kind != "product-emission-spectrum":
            raise ValueError(f"Expected product-emission-spectrum, got {kind}")

        shape = tuple(int(value) for value in file["result/normalized_shape"][...])
        if len(shape) != 3:
            raise ValueError(f"Unexpected result shape: {shape}")
        n_cells, n_mu, n_bins = shape

        r = file["cells/R"][...].astype(np.float64)
        z = file["cells/Z"][...].astype(np.float64)
        psi = file["cells/psi_norm"][...].astype(np.float64)
        b_hat = file["cells/B_hat_flat"][...].astype(np.float64)
        b_hat = b_hat.reshape(n_cells, 3)
        b_hat /= np.linalg.norm(b_hat, axis=1)[:, None]
        mu_grid = file["grids/mu_B"][...].astype(np.float64)
        energy = file["grids/energy"][...].astype(np.float64)
        dataset = file["result/normalized_flat"]

        if len(r) != n_cells or len(psi) != n_cells:
            raise ValueError("Cell coordinates disagree with result shape")
        if len(mu_grid) != n_mu or len(energy) != n_bins:
            raise ValueError("Pitch or energy grid disagrees with result shape")
        if dataset.shape != (n_cells * n_mu * n_bins,):
            raise ValueError("normalized_flat has an unexpected length")

        metadata = {
            "samples": int(file["run/samples_requested"][()]),
            "seed": int(file["run/seed"][()]),
            "workers": int(file["run/workers"][()]),
        }
    except Exception:
        file.close()
        raise

    return {
        "file": file,
        "dataset": dataset,
        "R": r,
        "Z": z,
        "psi": psi,
        "B_hat": b_hat,
        "mu_B": mu_grid,
        "energy": energy,
        "n_cells": n_cells,
        "n_mu": n_mu,
        "n_bins": n_bins,
        "metadata": metadata,
    }


def build_geometry(emission_map, gfile):
    """Construct the structured (rho, theta) cell index used by LOS lookup."""
    psi = emission_map["psi"]
    r = emission_map["R"]
    z = emission_map["Z"]
    levels = np.unique(psi)
    ring = np.searchsorted(levels, psi)
    counts = np.bincount(ring)
    if not np.all(counts == counts[0]):
        raise ValueError("Flux surfaces do not contain equal cell counts")
    n_theta = int(counts[0])
    if len(levels) * n_theta != len(psi):
        raise ValueError("Cells do not form a full (rho, theta) grid")

    r_axis, z_axis, boundary_r, boundary_z = read_gfile(gfile)
    angle = (np.arctan2(z - z_axis, r - r_axis) + 2.0 * np.pi) % (
        2.0 * np.pi
    )
    order = np.lexsort((angle, ring))
    cell_index = order.reshape(len(levels), n_theta)
    radius = np.hypot(
        r[cell_index] - r_axis,
        z[cell_index] - z_axis,
    )
    if np.any(np.diff(radius, axis=0) <= 0.0):
        raise ValueError("Flux-surface radii are not monotone")
    return {
        "cell_index": cell_index,
        "radius": radius,
        "r_axis": r_axis,
        "z_axis": z_axis,
        "boundary_R": boundary_r,
        "boundary_Z": boundary_z,
        "n_theta": n_theta,
    }


def point_in_polygon(r, z, boundary_r, boundary_z):
    """Ray-casting test against the gfile LCFS."""
    inside = False
    count = len(boundary_r)
    for index in range(count):
        next_index = (index + 1) % count
        r1, z1 = boundary_r[index], boundary_z[index]
        r2, z2 = boundary_r[next_index], boundary_z[next_index]
        if (z1 > z) != (z2 > z):
            crossing_r = r1 + (r2 - r1) * (z - z1) / (z2 - z1)
            if r < crossing_r:
                inside = not inside
    return inside


def pitch_weights(mu_grid, pitch):
    """Linear interpolation weights on the stored pitch grid."""
    pitch = float(np.clip(pitch, mu_grid[0], mu_grid[-1]))
    if len(mu_grid) == 1:
        return ((0, 1.0),)
    upper = int(np.searchsorted(mu_grid, pitch, side="left"))
    if upper == 0:
        return ((0, 1.0),)
    if upper == len(mu_grid):
        return ((len(mu_grid) - 1, 1.0),)
    if mu_grid[upper] == pitch:
        return ((upper, 1.0),)
    lower = upper - 1
    fraction = (pitch - mu_grid[lower]) / (mu_grid[upper] - mu_grid[lower])
    return ((lower, 1.0 - fraction), (upper, fraction))


def convolve(emission_map, geometry, los_path):
    """Convolve one full-space FPESC map with all records in a LOS file."""
    x, y, z, coefficient_c, direction = read_los(los_path)
    n_records = len(x)
    n_cells = emission_map["n_cells"]
    n_mu = emission_map["n_mu"]
    cell_index = geometry["cell_index"]
    radius = geometry["radius"]
    r_axis = geometry["r_axis"]
    z_axis = geometry["z_axis"]
    n_theta = geometry["n_theta"]
    delta_theta = 2.0 * np.pi / n_theta
    mu_grid = emission_map["mu_B"]
    b_hat = emission_map["B_hat"]

    coefficients = np.zeros((n_cells, n_mu), dtype=np.float64)
    used = outside = 0

    for record in range(n_records):
        point_r = float(np.hypot(x[record], y[record]))
        point_z = float(z[record])
        if not point_in_polygon(
            point_r,
            point_z,
            geometry["boundary_R"],
            geometry["boundary_Z"],
        ):
            outside += 1
            continue

        distance = float(np.hypot(point_r - r_axis, point_z - z_axis))
        angle = float(
            (np.arctan2(point_z - z_axis, point_r - r_axis) + 2.0 * np.pi)
            % (2.0 * np.pi)
        )
        unwrapped = (angle - 0.5 * delta_theta) % (2.0 * np.pi)
        column = int(unwrapped / delta_theta)
        j0 = column % n_theta
        j1 = (column + 1) % n_theta
        angular_fraction = unwrapped / delta_theta - column

        interpolated_radius = (
            radius[:, j0] * (1.0 - angular_fraction)
            + radius[:, j1] * angular_fraction
        )
        if distance <= interpolated_radius[0]:
            radial_cells = ((0, 1.0),)
        elif distance >= interpolated_radius[-1]:
            radial_cells = ((len(interpolated_radius) - 1, 1.0),)
        else:
            upper = int(np.searchsorted(interpolated_radius, distance))
            lower = upper - 1
            fraction = (
                (distance - interpolated_radius[lower])
                / (interpolated_radius[upper] - interpolated_radius[lower])
            )
            radial_cells = ((lower, 1.0 - fraction), (upper, fraction))

        toroidal_angle = float(np.arctan2(y[record], x[record]))
        cosine = np.cos(toroidal_angle)
        sine = np.sin(toroidal_angle)
        dx, dy, dz = direction[record]
        ray_direction = np.asarray(
            (
                dx * cosine + dy * sine,
                -dx * sine + dy * cosine,
                dz,
            ),
            dtype=np.float64,
        )
        geometry_weight = 4.0 * np.pi * coefficient_c[record]

        for radial_row, radial_weight in radial_cells:
            for column_offset, angular_weight in (
                (0, 1.0 - angular_fraction),
                (1, angular_fraction),
            ):
                if radial_weight == 0.0 or angular_weight == 0.0:
                    continue
                cell = int(
                    cell_index[
                        radial_row, (j0 + column_offset) % n_theta
                    ]
                )
                pitch = float(np.dot(ray_direction, b_hat[cell]))
                weight = geometry_weight * radial_weight * angular_weight
                for mu_index, mu_weight in pitch_weights(mu_grid, pitch):
                    coefficients[cell, mu_index] += weight * mu_weight
        used += 1

    result = np.zeros(emission_map["n_bins"], dtype=np.float64)
    block_length = n_mu * emission_map["n_bins"]
    dataset = emission_map["dataset"]
    for cell in np.flatnonzero(np.any(coefficients != 0.0, axis=1)):
        block = dataset[cell * block_length:(cell + 1) * block_length]
        block = block.reshape(n_mu, emission_map["n_bins"])
        result += coefficients[cell] @ block

    return result, {"records": n_records, "used": used, "outside": outside}


def fwhm(energy, values):
    """Interpolated full width at half maximum."""
    peak = int(np.argmax(values))
    if not values[peak] > 0.0:
        return float("nan")
    half = 0.5 * values[peak]
    left = None
    right = None
    for index in range(peak - 1, -1, -1):
        if values[index] <= half <= values[index + 1]:
            left = energy[index] + (half - values[index]) * (
                energy[index + 1] - energy[index]
            ) / (values[index + 1] - values[index])
            break
    for index in range(peak, len(values) - 1):
        if values[index] >= half >= values[index + 1]:
            right = energy[index] + (half - values[index]) * (
                energy[index + 1] - energy[index]
            ) / (values[index + 1] - values[index])
            break
    return float("nan") if left is None or right is None else right - left


def write_spectrum(
    path,
    emission_map,
    spectrum,
    stats,
    los_path,
    input_path,
    gfile,
):
    """Write the compact four-column spectrum text file."""
    energy = emission_map["energy"] / MEV
    width = float(energy[1] - energy[0])
    if not np.allclose(np.diff(energy), width, rtol=1.0e-12, atol=0.0):
        raise ValueError("Energy bins are not uniform")

    total = float(spectrum.sum())
    differential = spectrum / width
    normalized = differential / total if total > 0.0 else np.zeros_like(spectrum)
    if total > 0.0:
        mean = float(energy @ spectrum / total)
        sigma_keV = 1000.0 * float(
            np.sqrt(np.sum((energy - mean) ** 2 * spectrum) / total)
        )
    else:
        mean = sigma_keV = float("nan")
    width_keV = 1000.0 * fwhm(energy, differential)
    first_edge = float(energy[0] - 0.5 * width)
    last_edge = float(energy[-1] + 0.5 * width)
    metadata = emission_map["metadata"]

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as output:
        output.write("# FPESC line-of-sight D-D neutron spectrum\n")
        output.write(f"# Gfile = {Path(gfile).resolve()}\n")
        output.write(f"# FPESC_H5 = {Path(input_path).resolve()}\n")
        output.write(f"# LOS_file = {Path(los_path).resolve()}\n")
        output.write(f"# LOS_records = {stats['records']}\n")
        output.write(f"# LOS_cells = {stats['used']}\n")
        output.write(f"# LOS_points_outside_map = {stats['outside']}\n")
        output.write(
            f"# Neutron_energy_range_MeV = {first_edge:.8g} {last_edge:.8g}\n"
        )
        output.write(f"# Bins = {len(energy)}\n")
        output.write(f"# Bin_width_keV = {1000.0 * width:.12e}\n")
        output.write(f"# Monte_Carlo_samples = {metadata['samples']}\n")
        output.write(f"# Random_seed = {metadata['seed']}\n")
        output.write(f"# Julia_workers = {metadata['workers']}\n")
        output.write(f"# Total_reaction_rate_s^-1 = {total:.12e}\n")
        output.write(f"# Mean_neutron_energy_MeV = {mean:.12e}\n")
        output.write(f"# Standard_deviation_keV = {sigma_keV:.12e}\n")
        output.write(f"# FWHM_keV = {width_keV:.12e}\n")
        output.write(
            "# Columns: energy_center_MeV bin_rate_s^-1 "
            "differential_rate_s^-1_MeV^-1 normalized_density_MeV^-1\n"
        )
        for row in zip(energy, spectrum, differential, normalized):
            output.write("%.8f\t%.12e\t%.12e\t%.12e\n" % row)

    return {
        "output": path.resolve(),
        "total": total,
        "mean_MeV": mean,
        "std_keV": sigma_keV,
        "fwhm_keV": width_keV,
    }


def process_map(input_path, output_path, gfile=GFILE, los_path=LOS_FILE):
    """Postprocess one result HDF5 file into one text spectrum."""
    for path in (gfile, input_path, los_path):
        if not Path(path).is_file():
            raise FileNotFoundError(path)

    emission_map = load_map(input_path)
    try:
        geometry = build_geometry(emission_map, gfile)
        spectrum, stats = convolve(emission_map, geometry, los_path)
        summary = write_spectrum(
            output_path,
            emission_map,
            spectrum,
            stats,
            los_path,
            input_path,
            gfile,
        )
    finally:
        emission_map["file"].close()

    return summary, stats


def print_summary(label, summary, stats):
    """Print the main quantities for one processed result file."""
    print(f"[{label}] Spectrum: {summary['output']}")
    print(f"LOS records used: {stats['used']}/{stats['records']}")
    print(f"Total rate: {summary['total']:.6e} s^-1")
    print(f"Mean energy: {summary['mean_MeV']:.6f} MeV")
    print(f"Standard deviation: {summary['std_keV']:.3f} keV")
    print(f"FWHM: {summary['fwhm_keV']:.3f} keV")


def main():
    """Process every configured input/output pair with one invocation."""
    for job_index, (input_path, output_path) in enumerate(JOBS, start=1):
        summary, stats = process_map(input_path, output_path)
        print_summary(f"job {job_index}/{len(JOBS)}", summary, stats)


if __name__ == "__main__":
    main()
