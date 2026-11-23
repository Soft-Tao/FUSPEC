# Third-party materials

The MIT License in [`LICENSE`](LICENSE) applies to original software authored
by the FUSPEC contributors. Third-party models and materials referenced by or
included in this repository remain subject to their original terms and should
be cited appropriately.

## Cross-section models

The built-in `dd-n3he` and `dt-n4he` presets in `FUSPEC/data/` combine two
source models. No ENDF-formatted cross-section files are bundled with the
current version.

### Bosch-Hale total cross sections

The total D-D and D-T cross sections use the analytic parameterization from:

> H.-S. Bosch and G. M. Hale, "Improved formulas for fusion cross-sections and
> thermal reactivities," *Nuclear Fusion* **32**, 611-631 (1992).
> <https://doi.org/10.1088/0029-5515/32/4/I07>

The coefficient branches stored in `FUSPEC/data/*.toml` are parameter inputs to
this model.

### DROSG differential cross sections

The center-of-mass angle-dependent differential cross sections use normalized
Legendre coefficients from:

> M. Drosg and N. Otuka, "Evaluation of the Absolute Angle-Dependent
> Differential Neutron Production Cross Sections by the Reactions
> 3H(p,n)3He, 1H(t,n)3He, 2H(d,n)3He, 3H(d,n)4He, and 2H(t,n)4He and of the
> Cross Sections of Their Time-Reversed Counterparts up to 30 MeV and Beyond"
> (2015).

The coefficients stored in `FUSPEC/data/*.toml` are used to construct the
normalized angular shape multiplying the Bosch-Hale total cross section.
