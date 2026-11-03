# vLeaf-DSSAT

vLeaf-DSSAT couples **vLeaf**, a mechanistic leaf- and canopy-scale model of
photosynthesis, stomatal conductance and leaf energy balance, into the
**CERES-Maize** crop model of the
[DSSAT Cropping System Model](https://github.com/DSSAT/dssat-csm-os) (CSM).

In stock CERES-Maize, daily carbon gain is empirical: intercepted light ×
radiation-use efficiency, reduced by stress multipliers. In vLeaf-DSSAT that
calculation is replaced by an hourly, physically based simulation of the
sunlit and shaded parts of the canopy. Water stress acts on stomatal
conductance, and nutrient stress acts on photosynthetic capacity, instead of
being multiplied onto the result.

This repository is built on the DSSAT CSM source code (v4.8.5). Everything
outside the vLeaf components is DSSAT as released upstream; other crops and
the rest of CERES-Maize behave exactly as in DSSAT.

> **Status: research prototype.** The coupling builds and runs, but vLeaf's
> default parameters are placeholders and are not yet calibrated. Calibrate
> them (see [Configuring vLeaf](#configuring-vleaf)) before using the results.

## How vLeaf works inside CERES-Maize

For every hour of every simulated day, vLeaf:

1. takes DSSAT's own hourly weather (air temperature, radiation, humidity,
   wind, solar angle, diffuse fraction) and atmospheric CO₂;
2. splits the canopy into sunlit and shaded leaf area and computes the
   radiation each absorbs and its photosynthetic capacity (Vcmax₂₅);
3. solves, for each fraction, a coupled system of:
   - C4 photosynthesis (Collatz-type),
   - Ball-Berry stomatal conductance,
   - leaf energy balance (leaf temperature, sensible and latent heat),
   - leaf boundary-layer conductance;
4. integrates hourly net assimilation (including night-time leaf
   respiration) and transpiration to daily totals.

The daily carbon gain replaces the radiation-use-efficiency calculation in
`MZ_GROSUB.for`. DSSAT's soil-water module (SPAM) still owns the water
balance: plant water stress is computed from SPAM's transpiration demand,
and vLeaf's own transpiration is reported as a diagnostic.

The full design, including units, the order of calls, and the
simplifications made, is documented in
[`Plant/CERES-Maize/VLEAF_COUPLING.md`](Plant/CERES-Maize/VLEAF_COUPLING.md).

## What differs from DSSAT

All vLeaf code lives in `Plant/CERES-Maize/`:

| File | Role |
|---|---|
| `MZ_VLEAF.for` | vLeaf driver: hourly loop, daily totals, hourly output |
| `VLEAF_UTILS.for` | Reads `MY_PARAMS.TXT`; converts DSSAT's hourly weather into vLeaf inputs |
| `VLEAF_PARAMS.f90` | Season-long vLeaf parameters |
| `can_Weather.f90`, `canopy_Photo.f90` | Sunlit/shaded radiation and photosynthetic capacity |
| `leaf_step_mod.f90` | Per-hour solver coupling the leaf modules |
| `c4_photosynth.f90`, `stomata.f90`, `enbal_mod.f90`, `boundary.f90` | Leaf physiology and physics |
| `daily_carbon_mod.f90` | Hourly → daily carbon and transpiration |

Changes to existing DSSAT files:

- `Plant/CERES-Maize/MZ_GROSUB.for`: calls vLeaf and uses its carbon gain
  in place of the radiation-use-efficiency calculation.
- `CMakeLists.txt`: adds the vLeaf source files to the build.

## Building

The build is the standard DSSAT CMake build. You need a Fortran compiler
(gfortran or Intel) and CMake. On macOS: `brew install gcc cmake`.

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=<install-dir>
cmake --build . -j8
cmake --install .
```

`-DCMAKE_BUILD_TYPE=DEBUG` (the default) adds floating-point-exception
trapping and array-bounds checks, which is useful when changing vLeaf code.

The install directory gets the `dscsm048` executable, the `DSSATPRO.L48`
profile and DSSAT's standard data (cultivar files, codes, CO₂ data).
DSSAT limits file paths to about 100 characters, so keep the install and run
directories on short paths.

## Running

vLeaf-DSSAT runs exactly like DSSAT, using a maize experiment run with the
CERES-Maize model (`MZCER048`). For example, from a directory containing
the experiment, weather and soil files:

```sh
<install-dir>/dscsm048 A UFGA8201.MZX
```

Example experiments, weather and soil files are in DSSAT's data repository,
[DSSAT/dssat-csm-data](https://github.com/DSSAT/dssat-csm-data).

In addition to DSSAT's normal outputs, each run writes:

- `VLEAF_HOURLY.OUT`: hourly sunlit and shaded leaf area, absorbed radiation,
  net assimilation, conductances, leaf temperatures, energy fluxes,
  transpiration, stress factors and solver convergence flags, with `RUN` and
  `TRNO` columns to separate treatments.
- A season summary in `WARNING.OUT` if any hourly leaf solve failed to
  converge.

## Configuring vLeaf

vLeaf parameters are read at the start of each season from an optional
`MY_PARAMS.TXT` in the run directory, one `key value` pair per line. Any key
that is left out keeps its default:

```text
vcmax25  55.0
slope    3.0
eta      1.0
```

| Key | Meaning | Default |
|---|---|---|
| `vcmax25` | Maximum carboxylation capacity at 25 °C (µmol m⁻² s⁻¹) | 55.0 |
| `rd25` | Leaf dark respiration as a fraction of Vcmax | 0.05 |
| `slope` | Ball-Berry slope | 3.0 |
| `intercept` | Ball-Berry intercept g₀ (mol m⁻² s⁻¹) | 0.008 |
| `eta` | Fraction of assimilated carbon kept as dry matter after whole-plant respiration | 1.0 |
| `switch` | `1` solves leaf temperature from the energy balance; any other value sets leaf temperature to air temperature | 1 |

The file also accepts `jmax25`, `vpmax25`, `theta`, `gbs`, `alpha`, `x`,
`sco25` and `vpr25`, but the current C4 formulation does not use them.

## Keeping up to date with DSSAT

This repository tracks the DSSAT source as a git remote:

```sh
git remote add upstream https://github.com/DSSAT/dssat-csm-os.git
git fetch upstream
git merge upstream/develop
```

Changes made upstream to `MZ_GROSUB.for` or `CMakeLists.txt` may need to be
merged by hand.

## License and acknowledgement

vLeaf-DSSAT is derived from the DSSAT Cropping System Model, Copyright (c)
DSSAT Foundation, and is distributed under the same BSD 3-Clause license;
see [`license.txt`](license.txt). Compiled binaries may include GCC runtime
libraries; see [`NOTICE.md`](NOTICE.md).

This project is independent. It is not affiliated with, maintained by or
endorsed by the DSSAT Foundation. For the official DSSAT software, see
[dssat.net](https://dssat.net) and
[DSSAT/dssat-csm-os](https://github.com/DSSAT/dssat-csm-os).

## Citing

If you use vLeaf-DSSAT, please cite DSSAT and the DSSAT Cropping System Model:

Hoogenboom, G., C.H. Porter, K.J. Boote, V. Shelia, P.W. Wilkens, U. Singh, J.W. White, S. Asseng, J.I. Lizaso, L.P. Moreno, W. Pavan, R. Ogoshi, L.A. Hunt, G.Y. Tsuji, and J.W. Jones. 2019. The DSSAT crop modeling ecosystem. In: p.173-216 [K.J. Boote, editor] Advances in Crop Modeling for a Sustainable Agriculture. Burleigh Dodds Science Publishing, Cambridge, United Kingdom (https://dx.doi.org/10.19103/AS.2019.0061.10).

Hoogenboom, G., C.H. Porter, V. Shelia, K.J. Boote, U. Singh, W. Pavan, F.A.A. Oliveira, L.P. Moreno-Cadena, T.B. Ferreira, J.W. White, J.I. Lizaso, D.N.L. Pequeno, B.A. Kimball, P.D. Alderman, K.R. Thorp, S.V. Cuadra, M.S. Vianna, F.J. Villalobos, W.D. Batchelor, S. Asseng, M.R. Jones, A. Hopf, H.B. Dias, A. Jintrawet, R. Jaikla, E. Memic, L.A. Hunt, and J.W. Jones. 2024. Decision Support System for Agrotechnology Transfer (DSSAT) Version 4.8.5 (www.DSSAT.net). DSSAT Foundation, Gainesville, Florida, USA.

Jones, J.W., G. Hoogenboom, C.H. Porter, K.J. Boote, W.D. Batchelor, L.A. Hunt, P.W. Wilkens, U. Singh, A.J. Gijsman, and J.T. Ritchie. 2003. The DSSAT cropping system model. European Journal of Agronomy 18:235-265 (https://doi.org/10.1016/S1161-0301(02)00107-7).

## Author

Antriksh Srivastava ([@GitHub-srivastava](https://github.com/GitHub-srivastava))
