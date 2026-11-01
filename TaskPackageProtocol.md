# FPESC `input.h5` Task Package Protocol (Contract Draft)

> Status: streamlined v0.2 draft, tracked in Git, currently being implemented against.

## 1. Established Semantics

1. A single `input.h5` represents exactly one reaction and one calculation kind.
2. The calculation kind is only `reaction-rate-profile` or `product-emission-spectrum`; there is no `both`.
3. `samples` is the total sample count for one Julia CLI run, not a per-cell or per-thread count; it is not written into `input.h5`.
4. Schema v1 is fixed to relativistic kinematics; the classical branch is not retained, and no configurable kinematics field is stored in the file.
5. The physical result of a reaction-rate profile is

   ```text
   R(R,Z) = pair_factor * n1(R,Z) * n2(R,Z) * <sigma*v>(R,Z)
   ```

   in units of `m^-3 s^-1`. The result file does not store an additional `<sigma*v>` diagnostic profile; if it is needed later, the user computes it from the distributions, densities, and cross section in the task.
6. When a population reacts with itself, the pair factor is 0.5; when two distinct populations react, the pair factor is 1.0. The pair factor is not written into the task; Julia derives it uniquely from `population_relation`.
7. The pair factor applies to both the reaction-rate profile and the product emission spectrum.
8. `importance_weight` remains in the schema, but every current task package stores `1.0` for every cell (uniform). It is never multiplied into the physical results; if it were non-uniform, Julia would use it only to allocate per-cell sample counts.
9. Reaction-rate-profile tasks do not select a product; product-emission-spectrum tasks must select p3 or p4.
10. Python preprocesses cross sections into the final interpolation knots/coefficients. Julia only rebuilds and evaluates the interpolator; it does not read ENDF or re-fit.
11. The worker count is not part of the task's physics contract; it is determined by `julia -t N`.
12. Schema v1 considers only two-body fusion reactions and does not support scattering reactions, so no `link_to` parameter is stored. The differential cross section is always tabulated for the selected product relative to reactant 1 (the projectile) in the centre-of-mass frame, so Julia fixes `link_to_product1 = link_to_product2 = 1`; for a selected product 2 (`p4`) the kernel looks the table up with `mu -> -mu`, because the two products are back to back in the centre-of-mass frame.

## 2. Data Types and Unit Conventions

- Floating-point physics quantities default to IEEE `Float64`.
- Array indices use signed `Int64` and are fixed to zero-based; `-1` means the cell has no available distribution.
- Counts and seeds use `UInt64`; the Julia CLI must parse and validate `samples` as `UInt64`.
- Strings are UTF-8 variable-length strings.
- All physics quantities use SI units: energy J, mass kg, velocity m/s, density m^-3, total cross section m^2, differential cross section m^2/sr.
- Multi-dimensional arrays are stored as flat datasets with the schema-fixed `order="C"`; the logical shape is primarily derived from offsets, knots, and other array lengths, and is stored only when it cannot be derived uniquely.

## 3. Root and Task Metadata

```text
/@schema_version                 = 1
/@unit_system                    = "SI"

/task/calculation_kind           # one of the two valid enum values
/task/physics_fingerprint        # excludes CLI samples and output paths
/task/label                      # optional; for log display only
```

Schema v1 itself denotes a relativistic fusion task; Julia does not read or accept any kinematics-mode override. `task/label` does not participate in the physics fingerprint.

## 4. Run and RNG

```text
/run/seed
```

`samples`, the worker count, thread IDs, and thread shards are not written into `input.h5`. `samples` is specified by the Julia CLI; parallel thread capacity is set by `julia -t N`. Schema v1 implicitly fixes the `fpesc_counter_v1` RNG, and the collective velocity of both reactants is implicitly `[0,0,0]`, with no dataset stored. The RNG must be decoupled from the runtime choice of sample count and worker count so that deterministic acceptance with 1/2/4 threads is supported.

## 5. Reaction and Particles

```text
/reaction/population_relation        # "same" or "distinct"
/reaction/selected_product           # only for spectrum tasks; 3 or 4

/reaction/particles/r1/label         # optional, e.g. "D"
/reaction/particles/r1/mass
/reaction/particles/r1/charge_number
/reaction/particles/r1/mass_number
/reaction/particles/r2/...
/reaction/particles/p3/...
/reaction/particles/p4/...
```

Reaction-rate-profile tasks must not contain `/reaction/selected_product`. The p3/p4 information must still be present; Julia uses the four masses to compute the Q value and validate the reaction definition.

`population_relation` is the single authoritative input that decides same/distinct; it is not a provenance label. Python does not auto-infer it from array fingerprints, particle species, or labels. Julia maps populations and derives the pair factor exactly as fixed by the schema:

- `same`: r1 and r2 both use `/populations/p0`; `/populations/p1` must not exist; pair factor = 0.5.
- `distinct`: r1 uses `/populations/p0` and r2 uses `/populations/p1`; pair factor = 1.0.

For example:

- D-thermal + D-thermal: `same`, factor = 0.5.
- D-fast + D-fast, with both sides referencing the same fast population: `same`, factor = 0.5.
- D-thermal + D-fast: `distinct`, factor = 1.0.
## 6. Cells and Population References

### 6.1 Shared Spatial Cells

```text
/cells/R                         # [n_cells], m
/cells/Z                         # [n_cells], m
/cells/psi_norm                  # [n_cells], dimensionless
/cells/importance_weight         # [n_cells], > 0; currently 1.0 for every cell
/cells/B_hat_flat                # required only for spectrum tasks; fixed length 3*n_cells
```

`R/Z/psi_norm` are used for result positioning. `B_hat` never participates in the reaction-rate-profile kernel. For product-emission-spectrum tasks Julia validates it and copies it into the result bundle for the downstream line-of-sight convolution, which converts every LOS ray direction into a pitch against the local `B_hat`; the spectrum kernel itself never reads it, because it works in the field-aligned frame described in section 7.

### 6.2 Deduplicated Populations

```text
/populations/p0/label            # optional, e.g. thermal-D
/populations/p0/density          # [n_cells], m^-3
/populations/p0/distribution_index # [n_cells], zero-based; -1 means None
/populations/p0/distributions/...

/populations/p1/...              # normally needed only for distinct tasks
```

The correspondence between populations and reactant roles is uniquely determined by `population_relation` together with the p0/p1 convention above. Same-population tasks do not store two copies of the large NUBEAM arrays.

`distribution_index` maps spatial cells to the deduplicated distribution table within that population:

```text
distribution_index[i] = k
```

This means cell `i` uses the `k`-th distribution in `/populations/pN/distributions`. It must therefore have length `[n_cells]`, like `R/Z/density`; `-1` means that population has no available distribution in that cell.

For example, NUBEAM may have tens of thousands of spatial cells but only about 840 zone distribution tables. Multiple cells may point to the same `k`, avoiding a per-cell copy of a `500x50` distribution table.

Default cell sampling probabilities remain consistent with the pure-Python version:

```text
q_i = importance_weight_i * n1_i * n2_i
P(cell=i) = q_i / sum(q)
```

With the current uniform weight (`importance_weight_i = 1.0` for every cell)
this reduces to `q_i = n1_i * n2_i`, i.e. cell sampling is proportional to the
reactant density product alone.

The pair factor is a global constant within a task; it does not need to enter the cell sampling probabilities, but it must be multiplied into the final physical quantities.

## 7. Distribution Packed Schema

Each population may contain multiple deduplicated distributions:

```text
/populations/pN/distributions/kind
    # 1 = Maxwellian
    # 2 = BiMaxwellian
    # 3 = TabulatedEnergyPitch

/populations/pN/distributions/parameters_flat
/populations/pN/distributions/parameter_offsets
    # Maxwellian: [temperature]
    # BiMaxwellian: [temperature_parallel, temperature_perpendicular]

/populations/pN/distributions/energy_values
/populations/pN/distributions/energy_offsets
/populations/pN/distributions/mu_values
/populations/pN/distributions/mu_offsets
/populations/pN/distributions/f_values
/populations/pN/distributions/f_offsets
```

- Maxwellian/BiMaxwellian temperatures use J.
- Tabulated energy uses J; mu is the pitch cosine relative to the local magnetic field; the logical shape of the f table is `[n_energy,n_mu]`.
- The kernels work in a field-aligned frame: the polar axis of every sampled reactant momentum is the local magnetic field (`sample_momentum` builds `(p*sqrt(1-mu^2)*cos(phi), p*sqrt(1-mu^2)*sin(phi), p*mu)` and looks the same `mu` up in the f table). Emission directions passed into the kinematics must be expressed in that same frame. Only the scalar pitch of a line of sight, `n_los . B_hat`, couples the frame back to the laboratory geometry, and it is evaluated outside the kernel, by the post-processing convolution against the `B_hat` stored in the result bundle. Rotating the emission directions into the laboratory basis without also rotating the sampled momenta silently mixes two frames and changes the beam-thermal spectrum.
- Tabulated/NUBEAM values fully preserve the absolute `f(E,mu)` values and units from pure-Python preprocessing; they are not divided by the cell/zone density to be normalized into a shape function.
- Distribution values outside the Tabulated grid are fixed to 0.
- Distributions do not carry particle masses; Julia binds the r1/r2 masses from the reactant roles in the reaction.
## 8. Reactant Sampling Ranges

```text
/sampling/reactant1/energy_range # [Emin,Emax], J
/sampling/reactant1/mu_range     # [mu_min,mu_max], within [-1,1]
/sampling/reactant2/energy_range
/sampling/reactant2/mu_range
```

The schema stores the mu range directly rather than the polar-angle range of the Python API, avoiding the `cos(theta)` endpoint-reversal ambiguity. The azimuthal angle remains uniformly sampled in `[0,2*pi)` and is not a task parameter in v1.

## 9. Cross-Section Contract

Both tables below are functions of the centre-of-mass kinetic energy of the
reactant pair, `K_CM = sqrt(s) - (m1 + m2)*c^2`, not of the projectile energy in
the target-rest frame. Python is responsible for converting the incident
energies of the raw ENDF evaluations (which are tabulated against the
target-rest deuteron energy) before the knots are written, and Julia evaluates
every cross section at the centre-of-mass kinetic energy of the sampled pair.

### 9.1 Reaction-Rate Profile: One-Dimensional Piecewise Cubic Polynomial

```text
/cross_section/kind                   = "total-cubic-ppoly"
/cross_section/energy_knots           # [n+1], J
/cross_section/coefficients_flat      # fixed length 4*n; order defined by schema v1
```

On interval `i`:

```text
sigma(E) = c[0,i]*(E-x[i])^3 + c[1,i]*(E-x[i])^2
         + c[2,i]*(E-x[i])   + c[3,i]
```

Python must convert the MeV/barn coefficients into J/m^2 coefficients before writing. The knot count is producer-dependent: the bundled task-5 packages use 700 knots and a `[4,699]` coefficient array, while a spline on the raw 161-point ENDF grid gives 161 knots and `[4,160]`. Beyond the knot range, the piecewise-cubic extrapolation of the pure Python/SciPy version is continued; it does not return 0 or refuse to compute.

### 9.2 Product Emission Spectrum: Two-Dimensional Tensor-Product B-Spline

```text
/cross_section/kind                   = "differential-tensor-bspline"
/cross_section/energy_knots           # full B-spline knot vector, J
/cross_section/mu_knots               # full B-spline knot vector
/cross_section/degrees                # [degree_E,degree_mu]
/cross_section/coefficients_flat      # logical shape derived from knots/degrees
```

After preprocessing with the current DT SciPy cubic `RegularGridInterpolator`:

- energy knot vector: 104;
- mu knot vector: 205;
- degrees: `[3,3]`;
- coefficient logical shape: `[100,201]`;
- returns 0 outside the original energy/mu domain.

The energy axis is `K_CM` in J, as in section 9.1; `mu` is the centre-of-mass
cosine of the selected product relative to reactant 1 (see schema item 12).
The knot/degree numbers above are those of the current ENDF angular table;
other producers may tabulate a different grid, since the knots are stored
explicitly and are self-describing.

Julia does not reconstruct the differential cross section from raw ENDF or Legendre coefficients.

## 10. Required Fields for Product-Emission-Spectrum Tasks

```text
/spectrum/energy_edges           # [n_bins+1], J, strictly increasing
/spectrum/mu_B                   # [n_mu], within [-1,1], strictly increasing
```

No presets such as `medium` are stored; only the exact grid already materialized by Python is stored. A pure-Maxwellian task may use the single point `mu_B=[1]`; LOS treats it as isotropic.

Reaction-rate-profile tasks must not contain `/spectrum`.

## 11. Provenance

```text
/provenance/sources/<role>/sha256
/provenance/sources/<role>/basename   # optional
/provenance/task_config               # optional
/provenance/created_utc               # optional
```

The provenance group is entirely optional and does not participate in Julia physics branches or the physics fingerprint. Raw gfile/CDF/NUBEAM contents are not embedded in the task; only the source role, basename, and hash may be stored.

## 12. Sanity Checks Julia Must Perform

1. version/unit match strictly; schema v1 always runs the relativistic fusion algorithm and uses the fixed RNG.
2. The Julia CLI `samples` lies in `[1, 2^64-1]`.
3. All four particle masses are positive and finite; charge/mass numbers are reasonable; the Q value can be uniquely reconstructed from the masses.
4. `same` allows only p0 and derives factor = 0.5; `distinct` requires both p0 and p1 and derives factor = 1.0.
5. For same-population tasks, the r1/r2 particle definitions must agree.
6. All cell arrays have consistent lengths; densities are non-negative; importance weights are positive.
7. `distribution_index` has length `n_cells`, and each entry is `-1` or within that population's distribution range; an active cell must have both positive density and a valid distribution.
8. Energy/mu sampling ranges are finite and ordered, and mu lies in `[-1,1]`.
9. Cross-section knots are strictly increasing, and the coefficient length matches the logical shape derived from the knots/degrees; out-of-range semantics are fixed by schema v1.
10. Reaction-rate-profile tasks must not have a selected product, a spectrum grid, or a differential cross section.
11. Product-emission-spectrum tasks must have a selected product, `B_hat`, a spectrum grid, and a differential cross section.
12. Each row of `B_hat` is finite and has unit norm within the prescribed tolerance.
13. The task has at least one active cell, and `sum(importance_weight*n1*n2) > 0` (with the current unit weights this is `sum(n1*n2) > 0`).

## 13. Content Not Stored in `input.h5`

- Full gfile, LCFS boundary, and psi interpolation grids;
- raw thermal CDF time series;
- raw NUBEAM zone coordinates, KD-tree, and BMVOL;
- raw ENDF files;
- LOS files;
- Julia worker count and output file paths;
- string presets such as `coarse/medium/fine` used for physics computations.

## 14. Compatibility Decisions Frozen This Round

1. Tabulated/NUBEAM distributions preserve the absolute `f(E,mu)` values of the pure-Python version with no new normalization.
2. The total cross section beyond the knots continues the SciPy cubic extrapolation.
3. Reaction-rate results store only the physical reaction rate after multiplying in the densities and pair factor; no separate `<sigma*v>` diagnostic profile is stored.
4. same/distinct is decided entirely by `/reaction/population_relation`; Python does not auto-infer it from fingerprints.


- D-thermal + T-thermal: `distinct`, factor = 1.0.
