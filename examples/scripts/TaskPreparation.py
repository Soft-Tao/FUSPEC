#!/usr/bin/env python3
"""Build the two example FUSPEC task packages from one gFile and two NPZ files."""

from __future__ import annotations

import argparse
from pathlib import Path

from freeqdsk import geqdsk
import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.spatial import cKDTree


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIRECTORY = REPOSITORY_ROOT / "examples" / "inputs"
ELEMENTARY_CHARGE = 1.602176634e-19
MEV = 1.0e6 * ELEMENTARY_CHARGE
FAST_ENERGY_J = np.linspace(500.0, 149500.0, 150) * ELEMENTARY_CHARGE
FAST_MU = np.linspace(-0.98, 0.98, 50)
THETA_COUNT = 300
DEFAULT_BINS = 128
DEFAULT_SPECTRUM_MIN_MEV = 2.0
DEFAULT_SPECTRUM_MAX_MEV = 3.0
DEFAULT_BEAM_ANGULAR_BINS = 101
DEFAULT_THERMAL_ANGULAR_BINS = 51


class FluxGeometry:
    """Small gFile helper for flux surfaces and field-aligned directions."""

    def __init__(self, path: Path):
        with path.open(encoding="ascii") as stream:
            equilibrium = geqdsk.read(stream)
        r = np.asarray(equilibrium.r_grid[:, 0], dtype=np.float64)
        z = np.asarray(equilibrium.z_grid[0, :], dtype=np.float64)
        psi_initial = np.asarray(equilibrium.psi, dtype=np.float64) - (
            equilibrium.simagx
        )
        psi_norm = psi_initial / (
            equilibrium.sibdry - equilibrium.simagx
        )
        self.r_axis = float(equilibrium.rmagx)
        self.z_axis = float(equilibrium.zmagx)
        self.boundary_r = np.asarray(equilibrium.rbdry, dtype=np.float64)
        self.boundary_z = np.asarray(equilibrium.zbdry, dtype=np.float64)
        self.psi = RegularGridInterpolator(
            (r, z), psi_norm, method="cubic",
            bounds_error=False, fill_value=np.nan,
        )
        dpsi_dr, dpsi_dz = np.gradient(psi_initial, r, z)
        self.dpsi_dr = RegularGridInterpolator(
            (r, z), dpsi_dr, method="linear",
            bounds_error=False, fill_value=np.nan,
        )
        self.dpsi_dz = RegularGridInterpolator(
            (r, z), dpsi_dz, method="linear",
            bounds_error=False, fill_value=np.nan,
        )
        self.fpol = interp1d(
            np.linspace(0.0, 1.0, len(equilibrium.fpol)), equilibrium.fpol,
            kind="linear", bounds_error=False, fill_value=np.nan,
        )

    def psi_at(self, r: np.ndarray, z: np.ndarray) -> np.ndarray:
        return np.asarray(self.psi(np.column_stack((r.ravel(), z.ravel())))).reshape(r.shape)

    def boundary_radius(self, theta: np.ndarray) -> np.ndarray:
        dx, dz = np.cos(theta)[:, None], np.sin(theta)[:, None]
        r1, z1 = self.boundary_r[:-1], self.boundary_z[:-1]
        r2, z2 = self.boundary_r[1:], self.boundary_z[1:]
        sr, sz = r2 - r1, z2 - z1
        denominator = dx * sz - dz * sr
        valid = np.abs(denominator) > 1.0e-14
        safe = np.where(valid, denominator, 1.0)
        ray = ((r1 - self.r_axis) * sz - (z1 - self.z_axis) * sr) / safe
        segment = ((r1 - self.r_axis) * dz - (z1 - self.z_axis) * dx) / safe
        hit = valid & (ray > 0.0) & (segment >= 0.0) & (segment <= 1.0)
        return np.min(np.where(hit, ray, np.inf), axis=1)

    def surface_points(self, rho: np.ndarray, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        boundary = np.broadcast_to(
            self.boundary_radius(theta), (len(rho), len(theta)),
        ).copy()
        target = np.broadcast_to(rho[:, None] ** 2, boundary.shape).copy()
        axis_psi = float(self.psi_at(
            np.asarray([self.r_axis]), np.asarray([self.z_axis]),
        )[0])
        low = np.zeros_like(boundary)
        high = boundary.copy()
        cos_theta, sin_theta = np.cos(theta)[None, :], np.sin(theta)[None, :]
        unresolved = target <= axis_psi
        high[unresolved] = 0.0
        for _ in range(50):
            radius = 0.5 * (low + high)
            psi = self.psi_at(
                self.r_axis + radius * cos_theta,
                self.z_axis + radius * sin_theta,
            )
            inside = np.nan_to_num(psi, nan=-np.inf) < target
            low = np.where(inside, radius, low)
            high = np.where(inside, high, radius)
        radius = 0.5 * (low + high)
        edge_psi = self.psi_at(
            self.r_axis + boundary * cos_theta,
            self.z_axis + boundary * sin_theta,
        )
        past_edge = ~unresolved & (
            ~np.isfinite(edge_psi) | (edge_psi < target)
        )
        radius[past_edge] = boundary[past_edge]
        for column in range(len(theta)):
            resolved = np.flatnonzero(~unresolved[:, column])
            if resolved.size:
                first = int(resolved[0])
                radius[:first, column] = (
                    radius[first, column] * rho[:first] / rho[first]
                )
            radius[:, column] = np.maximum.accumulate(radius[:, column])
            radius[:, column] += (
                np.arange(len(rho), dtype=np.float64) * 1.0e-9
            )
        return (
            self.r_axis + radius * cos_theta,
            self.z_axis + radius * sin_theta,
        )

    def field_direction(self, r: np.ndarray, z: np.ndarray) -> np.ndarray:
        psi = self.psi_at(r, z)
        points = np.column_stack((r.ravel(), z.ravel()))
        dpsi_dr = np.asarray(self.dpsi_dr(points)).reshape(r.shape)
        dpsi_dz = np.asarray(self.dpsi_dz(points)).reshape(r.shape)
        field = np.stack(
            (-dpsi_dz / r, self.fpol(np.clip(psi, 0.0, 1.0)) / r, dpsi_dr / r),
            axis=-1,
        )
        return field / np.linalg.norm(field, axis=-1, keepdims=True)


def _string(group: h5py.Group, name: str, value: str) -> None:
    group.create_dataset(name, data=value, dtype=h5py.string_dtype("utf-8"))


def _maxwellian_population(
    file: h5py.File,
    name: str,
    density: np.ndarray,
    index: np.ndarray,
    temperature_j: np.ndarray,
) -> None:
    population = file.create_group(f"populations/{name}")
    population["density"] = density
    population["distribution_index"] = index.astype(np.int64)
    distributions = population.create_group("distributions")
    distributions["kind"] = np.ones(len(temperature_j), dtype=np.uint8)
    distributions["parameters_flat"] = temperature_j
    distributions["parameter_offsets"] = np.arange(len(temperature_j) + 1, dtype=np.int64)
    for prefix in ("energy", "mu", "f"):
        distributions[f"{prefix}_values"] = np.empty(0, dtype=np.float64)
        distributions[f"{prefix}_offsets"] = np.zeros(len(temperature_j) + 1, dtype=np.int64)


def _tabulated_population(
    file: h5py.File,
    name: str,
    density: np.ndarray,
    index: np.ndarray,
    energy_j: np.ndarray,
    mu: np.ndarray,
    values: np.ndarray,
) -> None:
    population = file.create_group(f"populations/{name}")
    population["density"] = density
    population["distribution_index"] = index.astype(np.int64)
    distributions = population.create_group("distributions")
    count = values.shape[0]
    distributions["kind"] = np.full(count, 3, dtype=np.uint8)
    distributions["parameters_flat"] = np.empty(0, dtype=np.float64)
    distributions["parameter_offsets"] = np.zeros(count + 1, dtype=np.int64)
    distributions.create_dataset(
        "energy_values", data=np.tile(energy_j, count),
        compression="gzip", shuffle=True,
    )
    distributions["energy_offsets"] = np.arange(count + 1, dtype=np.int64) * len(energy_j)
    distributions.create_dataset(
        "mu_values", data=np.tile(mu, count),
        compression="gzip", shuffle=True,
    )
    distributions["mu_offsets"] = np.arange(count + 1, dtype=np.int64) * len(mu)
    distributions.create_dataset(
        "f_values", data=values.reshape(-1, order="C"),
        compression="gzip", shuffle=True,
    )
    distributions["f_offsets"] = (
        np.arange(count + 1, dtype=np.int64) * len(energy_j) * len(mu)
    )


def _write_task(
    path: Path,
    *,
    label: str,
    relation: str,
    cells: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    thermal: tuple[np.ndarray, np.ndarray, np.ndarray],
    fast: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
    sampling_energy_j: tuple[np.ndarray, np.ndarray],
    seed: int,
    bins: int,
    angular_bins: int,
    spectrum_energy_j: tuple[float, float],
) -> None:
    r, z, psi, b_hat = cells
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as file:
        file.attrs["schema_version"] = np.uint64(2)
        file.attrs["unit_system"] = "SI"
        task = file.create_group("task")
        _string(task, "calculation_kind", "product-emission-spectrum")
        _string(task, "label", label)
        file.create_group("run")["seed"] = np.uint64(seed)
        reaction = file.create_group("reaction")
        _string(reaction, "name", "dd-n3he")
        _string(reaction, "population_relation", relation)
        reaction["selected_product"] = np.int64(3)
        cell_group = file.create_group("cells")
        for name, values in (("R", r), ("Z", z), ("psi_norm", psi)):
            cell_group[name] = values
        cell_group["importance_weight"] = np.ones_like(r)
        b_dataset = cell_group.create_dataset("B_hat_flat", data=b_hat.reshape(-1))
        b_dataset.attrs["order"] = "C"
        if fast is None:
            density, index, temperature = thermal
            _maxwellian_population(file, "p0", density, index, temperature)
        else:
            fast_density, fast_index, fast_values = fast
            _tabulated_population(
                file, "p0", fast_density, fast_index, FAST_ENERGY_J, FAST_MU,
                fast_values,
            )
            density, index, temperature = thermal
            _maxwellian_population(file, "p1", density, index, temperature)
        sampling = file.create_group("sampling")
        for name, energy in zip(("reactant1", "reactant2"), sampling_energy_j):
            group = sampling.create_group(name)
            group["energy_range"] = energy
            group["mu_range"] = np.asarray([-1.0, 1.0])
        spectrum = file.create_group("spectrum")
        spectrum["energy_edges"] = np.linspace(
            spectrum_energy_j[0], spectrum_energy_j[1], bins + 1,
        )
        spectrum["mu_B"] = np.linspace(-1.0, 1.0, angular_bins)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gfile", required=True, type=Path)
    parser.add_argument("--thermal", required=True, type=Path)
    parser.add_argument("--fast", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    parser.add_argument("--seed", type=int, default=20260922)
    parser.add_argument("--bins", type=int, default=DEFAULT_BINS)
    parser.add_argument(
        "--spectrum-min-mev", type=float, default=DEFAULT_SPECTRUM_MIN_MEV,
    )
    parser.add_argument(
        "--spectrum-max-mev", type=float, default=DEFAULT_SPECTRUM_MAX_MEV,
    )
    parser.add_argument(
        "--beam-angular-bins", type=int, default=DEFAULT_BEAM_ANGULAR_BINS,
    )
    parser.add_argument(
        "--thermal-angular-bins", type=int,
        default=DEFAULT_THERMAL_ANGULAR_BINS,
    )
    args = parser.parse_args()

    with np.load(args.thermal) as data:
        thermal_density = np.asarray(data["density_m3"], dtype=np.float64)
        thermal_temperature = np.asarray(data["temperature_keV"], dtype=np.float64) * 1000.0 * ELEMENTARY_CHARGE
    with np.load(args.fast) as data:
        zones = np.asarray(data["zones_RZ_m"], dtype=np.float64)
        fast_values = np.asarray(data["distribution"], dtype=np.float64)
    if fast_values.shape != (len(zones), len(FAST_ENERGY_J), len(FAST_MU)):
        raise ValueError("FastIonDistribution distribution must have shape (n_zone, 150, 50)")

    geometry = FluxGeometry(args.gfile)
    rho = (np.arange(len(thermal_density), dtype=np.float64) + 0.5) / len(thermal_density)
    theta = (np.arange(THETA_COUNT, dtype=np.float64) + 0.5) * (2.0 * np.pi / THETA_COUNT)
    r, z = geometry.surface_points(rho, theta)
    psi = np.repeat(rho**2, THETA_COUNT)
    b_hat = geometry.field_direction(r, z)
    cells = (r.ravel(), z.ravel(), psi, b_hat.reshape(-1, 3))

    zone_rho = np.sqrt(np.clip(geometry.psi_at(zones[:, 0], zones[:, 1]), 0.0, 1.0))
    zone_theta = np.arctan2(zones[:, 1] - geometry.z_axis, zones[:, 0] - geometry.r_axis)
    zone_points = np.column_stack((zone_rho * np.cos(zone_theta), zone_rho * np.sin(zone_theta)))
    rho_flat = np.repeat(rho, THETA_COUNT)
    theta_flat = np.tile(theta, len(rho))
    cell_points = np.column_stack(
        (rho_flat * np.cos(theta_flat), rho_flat * np.sin(theta_flat)),
    )
    fast_index = cKDTree(zone_points).query(cell_points)[1]
    energy_integral = np.trapz(fast_values, x=FAST_ENERGY_J, axis=1)
    fast_density_zone = np.trapz(energy_integral, x=FAST_MU, axis=1)

    thermal_density_cells = np.repeat(thermal_density, THETA_COUNT)
    thermal_index = np.repeat(np.arange(len(thermal_density), dtype=np.int64), THETA_COUNT)
    thermal = (thermal_density_cells, thermal_index, thermal_temperature)
    thermal_range = np.asarray([0.0, 20.0 * 1000.0 * ELEMENTARY_CHARGE])
    fast_range = np.asarray([FAST_ENERGY_J[0], FAST_ENERGY_J[-1]])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_task(
        args.output_dir / "TaskPackage1.h5",
        label="EAST 13829C06 D-D beam-thermal",
        relation="distinct",
        cells=cells,
        thermal=thermal,
        fast=(fast_density_zone[fast_index], fast_index, fast_values),
        sampling_energy_j=(fast_range, thermal_range),
        seed=args.seed,
        bins=args.bins,
        angular_bins=args.beam_angular_bins,
        spectrum_energy_j=(
            args.spectrum_min_mev * MEV, args.spectrum_max_mev * MEV,
        ),
    )
    _write_task(
        args.output_dir / "TaskPackage2.h5",
        label="EAST 13829C06 D-D thermal",
        relation="same",
        cells=cells,
        thermal=thermal,
        fast=None,
        sampling_energy_j=(thermal_range, thermal_range),
        seed=args.seed,
        bins=args.bins,
        angular_bins=args.thermal_angular_bins,
        spectrum_energy_j=(
            args.spectrum_min_mev * MEV, args.spectrum_max_mev * MEV,
        ),
    )
    print(args.output_dir / "TaskPackage1.h5")
    print(args.output_dir / "TaskPackage2.h5")


if __name__ == "__main__":
    main()
