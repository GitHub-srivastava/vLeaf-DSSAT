# vLeaf ↔ CERES-Maize coupling

This document describes how the vLeaf leaf/canopy photosynthesis-stomatal
conductance-energy balance submodel is wired into CERES-Maize, why it is
wired that way, and what its remaining simplifications are.

## 1. What vLeaf replaces

Stock CERES-Maize computes daily potential canopy carbon gain (`PCARB`)
empirically:

```fortran
PCARB = IPAR * RUE * PCO2        ! intercepted PAR x radiation-use efficiency
CARBO = PCARB * AMIN1(PRFT, SWFAC, NSTRES, PStres1, KSTRES) * SLPF
```

vLeaf instead computes canopy carbon gain and transpiration mechanistically,
hour by hour, from a coupled C4 photosynthesis + Ball-Berry stomatal
conductance + leaf energy balance model, run separately for sunlit and
shaded leaf fractions and integrated to a daily total:

```fortran
PCARB = CARBO_vLeaf(SWFAC=1, NSTRESS=1) / PLTPOP          ! potential
CARBO = CARBO_vLeaf(SWFAC,   NSTRESS)   / PLTPOP * SLPF   ! actual
```

`PRFT` (temperature response) and the RUE multiplication are no longer
applied for this pathway — vLeaf's own `c4_photosynth.f90` has its own
Vcmax temperature response, so temperature stress is handled inside the
mechanistic model rather than by the empirical `PRFTC` cultivar curve.
Water and nutrient stress are applied as inputs *into* vLeaf rather than as
multipliers on its output: `SWFAC` scales stomatal conductance above its
residual `g0` (`stomata.f90`), and `NSTRESS = MIN(NSTRES, PStres1, KSTRES)`
scales Vcmax25 (`canopy_Photo.f90`).

## 2. Call chain

```
MZ_CERES.for  (per-crop dispatcher)
     |
     v
MZ_GROSUB.for  (DYNAMIC = RUNINIT / SEASINIT / INTEGR / OUTPUT / SEASEND)
     |  CALL GET(WEATHER)                       [ModuleData accessor]
     |  CALL MZ_VLEAF(..., WEATHER, CO2)         one call per phase,
     |                                           1-2 calls in INTEGR (§5)
     v
MZ_VLEAF.for  (owns vLeaf's own DYNAMIC state machine)
     |  CALL MZ_VLEAF_HOURLY(WEATHER, CO2, AtmPRES, ...)   [VLEAF_UTILS.for]
     |  DO K = 1, TS   (TS = 24, ModuleDefs.for)
     |     CALL canopy_weather   [can_Weather.f90]   - sunlit/shaded PAR & NIR
     |     CALL sunshade_photo   [canopy_Photo.f90]  - sunlit/shaded Vcmax25
     |     CALL leaf_step        [leaf_step_mod.f90] - per-hour, per-fraction:
     |          energy_balance_canopy [enbal_mod.f90]
     |          call_stomata         [stomata.f90]
     |          c4_photosynthesis    [c4_photosynth.f90]
     |          call_boundary        [boundary.f90]
     |  CALL daily_carbon / daily_transpiration     [daily_carbon_mod.f90]
     v
CARBO_vLeaf (g DM/m2/d)  -->  MZ_GROSUB's PCARB / CARBO
EOPVLF      (mm/d)       -->  diagnostic only (see §5)
```

`VLEAF_PARAMS.f90` holds the season's calibration constants (Vcmax25, Jmax25,
Ball-Berry intercept/slope, `eta`, etc.), set once at `SEASINIT` from
`MZ_RPARAMS`/`MY_PARAMS.TXT` (or built-in defaults if that file is absent),
and read by the physics modules via `USE VLEAF_PARAMS`.

## 3. Weather sourcing

vLeaf reads directly from DSSAT's own `WeatherType` (declared in
`Utilities/ModuleDefs.for`), which every DSSAT run populates once per
simulated day, before any Plant module executes:

```
Weather/WEATHR.for  -->  Weather/HMET.for  -->  WEATHER%TAIRHR/RADHR/WINDHR/
                                                 RHUMHR/BETA/AZZON/FRDIFP/
                                                 FRDIFR (all DIMENSION(TS))
                     -->  CALL PUT(WEATHER)      [ModuleData]
```

`MZ_GROSUB.for` retrieves it with `CALL GET(WEATHER)` and passes it straight
through to `MZ_VLEAF` alongside `CO2` (DSSAT's daily atmospheric CO2).

`VLEAF_UTILS.for`'s `MZ_VLEAF_HOURLY` subroutine is the **single place** that
reads `WEATHER`. It converts DSSAT's fields into plain per-hour driver arrays,
and everything downstream works only from those arrays.

| vLeaf needs (hourly) | Sourced from | Notes |
|---|---|---|
| Air temperature | `WEATHER%TAIRHR(H)` | °C, direct |
| Wind speed | `WEATHER%WINDHR(H)` | m/s, floored at 0.1 |
| Vapor pressure of air | `WEATHER%RHUMHR(H)/100 * VPSAT(TAIRHR(H))` | `VPSAT` is DSSAT's saturation-vapor-pressure function (`Weather/HMET.for`, Pa) |
| Solar zenith angle | `90.0 - WEATHER%BETA(H)` | `BETA` (solar elevation, deg) computed daily by `HANG` inside `HMET` |
| Direct/diffuse PAR and NIR | `WEATHER%RADHR(H)` split 45/55 into PAR/NIR, then each split by `WEATHER%FRDIFP(H)` / `WEATHER%FRDIFR(H)` | diffuse fractions are Spitters (1986), computed by `FRACD` in `HMET.for`. At night `RADHR` is 0, so all four components go to zero |
| Atmospheric CO2 | `CO2` (already a `MZ_GROSUB` argument) | ppm, constant across the day |
| Atmospheric pressure | Computed once at `SEASINIT` from `WEATHER%XELEV` via the standard barometric formula | a missing elevation (-99) is taken as sea level |
| Downwelling longwave | Estimated each hour from air temperature and vapor pressure (Brutsaert clear-sky formula) | DSSAT has no native hourly LW estimate |

### A unit trap worth knowing about: why not `WEATHER%PARHR`?

`WEATHER%PARHR` is a **photon flux density** (µmol m⁻² s⁻¹), whereas vLeaf's
leaf physics wants PAR as an **energy flux** (W m⁻²): `c4_photosynth.f90`
applies its own `convPAR = 4.6` µmol/J conversion. Feeding `PARHR` in would
double-convert and inflate absorbed radiation ~4.6×. The adapter instead
takes `RADHR` (J m⁻² s⁻¹ = W m⁻²) and applies the conventional 45/55
PAR/NIR energy split. Using `PARHR/4.6` instead would honour a measured
PAR/SRAD ratio, but that is a calibration decision, not a bug fix.

## 4. What canopy-radiation physics vLeaf still does itself

DSSAT's `HMET`/`FRACD` compute the *atmospheric* direct/diffuse split.
vLeaf's own `can_Weather.f90` (`canopy_weather`/`sw_band`) and
`canopy_Photo.f90` (`sunshade_photo`) compute how that radiation is
partitioned *inside the canopy* between sunlit and shaded leaf area, and how
average photosynthetic capacity (Vcmax25) differs between the two fractions.
Both use the same cosine-of-zenith floor (0.05), so the sunlit/shaded split
used for Vcmax25 matches the one used for absorbed radiation at low sun.

## 5. Water stress, transpiration and the two vLeaf passes

**SPAM owns the water balance.** SPAM computes potential transpiration
`EOP` and drains the soil by `EP = MIN(EOP, TRWUP*10)` (`SPAM/SPAM.for`).
`MZ_GROSUB` computes `SWFAC`/`TURFAC` from that same `EOP` and `TRWUP`,
exactly as stock CERES-Maize does, so plant water stress and soil water
depletion are always computed from one transpiration demand.

vLeaf's own transpiration (`EOPVLF`) is **diagnostic only**: it is not fed
back into `EOP`. (An earlier version set `EOP = EOPVLF` and derived
`SWFAC = TRWUP / EOPVLF`, while SPAM still removed water using its own
`EOP` — two transpiration numbers for one crop. On UFGA8201 that version
compressed the observed 3× irrigated/rainfed biomass contrast to 1.2×.)
Making vLeaf's transpiration drive the soil water balance would require
SPAM to call vLeaf (as it does ETPHOT for `MEPHO = 'L'`); that is a larger
architectural change.

`MZ_GROSUB` calls `MZ_VLEAF` up to twice per `INTEGR` day:

1. **Potential pass** with `SWFAC = 1` and `NSTRESS = 1`. Its result sets
   `PCARB`, which drives N demand (`PDWI`/`PGRORT`). Like the RUE-based
   `PCARB` it replaces, it must be free of water *and* nutrient stress —
   otherwise an N-stressed crop would demand less N and stay stressed.
2. **Actual pass** with today's `SWFAC` and `NSTRESS`, giving the
   `CARBO_vLeaf` used for growth and the hourly values written to
   `VLEAF_HOURLY.OUT`. It is skipped when `SWFAC = 1` and `NSTRESS = 1`,
   because it would reproduce the potential pass exactly.

## 6. Daily carbon and transpiration integration

`daily_carbon_mod.f90` treats each hourly value as the mean rate over its
own hour, so all 24 hours count (the earlier trapezoid over hours 1..24
dropped one interval).

- **Leaf respiration** — including night-time dark respiration (negative
  `Anet` of the shaded canopy at night) — is kept in the daily net. The
  earlier version clipped each hour's canopy flux at zero, which discarded
  night respiration (~14% of seasonal uptake on UFGA8201 treatment 4).
- **Whole-plant respiration** of non-leaf organs, which RUE-based
  CERES-Maize accounts for implicitly, is represented by `eta`: the
  fraction of assimilated C retained as dry matter (`gDM = gC/0.45 * eta`).
  Physiologically ~0.7–0.8 for maize; the default is **1.0** because the
  other defaults are not yet calibrated for it (§8). Set it with
  `eta <value>` in `MY_PARAMS.TXT`.
- The daily total is floored at zero: CERES has no reserve pool to draw a
  negative daily carbon balance from.

## 7. Leaf-solve robustness and output

`leaf_step` iterates energy balance → stomata → photosynthesis → boundary
layer to a fixed point (max 100 iterations).

- An **unconverged** solve keeps its last iterate (far closer to the answer
  than zero) and is flagged `sunConv`/`shConv = 0` in `VLEAF_HOURLY.OUT`.
- A **NaN or out-of-range** result (|Anet| > 80 µmol m⁻² s⁻¹) is set to
  zero for that hour.
- Both are counted over the season and reported once at `SEASEND` in
  `WARNING.OUT`. On UFGA8201 about 0.5% of hourly solves hit the iteration
  limit and none were NaN.

`VLEAF_HOURLY.OUT` carries `RUN` and `TRNO` columns, so multi-treatment and
batch runs can be separated. Its header matches the data columns (an
earlier header had `shLWn` out of place, shifting every label from `sunRes`
onward by one column). File units come from DSSAT's `GETLUN`.

State persistence: values set at `SEASINIT` and read later
(`intercept`, `slope`, `switch`, `eta`) live in the saved `VLEAF_PARAMS`
module; arrays written during `OUTPUT` (`LW` and all hourly fluxes) are
`SAVE`d.

## 8. Calibration status (UFGA8201, built-in defaults)

With the built-in defaults (no `MY_PARAMS.TXT`), total biomass at maturity
(kg/ha) on the standard UFGA8201 experiment:

| TRT | Observed | Stock CERES | vLeaf |
|---|---|---|---|
| 1 rainfed, low N | 5532 | 6556 | 5526 |
| 2 rainfed, high N | 7201 | 7538 | 6406 |
| 3 irrigated, low N | 14581 | 14578 | 8422 |
| 4 irrigated, high N | 22001 | 22454 | 9931 |
| 5 veg. stress, low N | 12002 | 13458 | 7258 |
| 6 veg. stress, high N | 17146 | 19234 | 8426 |

vLeaf is close under water limitation but produces less than half the
observed biomass when water is not limiting. Diagnostics from
`VLEAF_HOURLY.OUT`:

- With the default Ball-Berry `slope = 3`, the steady state
  ci/cs ≈ 1 − 1.6/(slope·hs) goes negative at moderate humidity, so ~31% of
  sunlit daytime hours sit at the ci floor (5% of ca) in `stomata.f90`;
  mean ci/ca is 0.17 (maize is typically ~0.4).
- Raising `slope` to 6 removes that (mean ci/ca 0.59) but lifts treatment 4
  only to ~12100 kg/ha, so the main gap is photosynthetic capacity / light
  use (`vcmax25`, quantum efficiency, `rd25`), not stomata.

The defaults are therefore placeholders: calibrate `vcmax25`, `rd25`,
`slope` and `eta` (via `MY_PARAMS.TXT`) before using vLeaf results.

## 9. Known simplifications (by design, not oversights)

- **Soil surface temperature** for the canopy longwave balance is
  `tAir − 1 °C` hourly. DSSAT's `SRFTEMP` is daily and is not passed to
  the maize module.
- **Atmospheric pressure and O2** are not modeled hourly by DSSAT;
  pressure is a one-time elevation-based estimate and O2 is unused by the
  C4 formulation.
- **Downwelling longwave radiation** has no DSSAT-native hourly
  representation; vLeaf's own clear-sky estimate is the only source.
- **`PRFTC` cultivar coefficients** are still read and `PRFT` still
  computed, but not applied to `CARBO` (§1). `RUE` is likewise still read
  but unused, which is the one compiler warning (`-Wunused-dummy-argument`)
  in the vLeaf-touched code.
- **`MY_PARAMS.TXT`** is a small, self-contained override file for
  physiological constants that DSSAT's species/ecotype files don't carry
  yet. If vLeaf is promoted from research prototype to a supported option,
  these constants belong in the `.SPE`/`.CUL` file format like every other
  CERES-Maize genetic coefficient.

## 10. File manifest

| File | Role |
|---|---|
| `MZ_VLEAF.for` | vLeaf's own DYNAMIC state machine; hourly driver loop; hourly output; convergence bookkeeping |
| `VLEAF_UTILS.for` | `MZ_RPARAMS` (calibration file reader) + `MZ_VLEAF_HOURLY` (**the sole DSSAT boundary**: WEATHER → plain hourly arrays) |
| `VLEAF_PARAMS.f90` | Season-long calibration constants, set once, read by the leaf-level physics modules |
| `can_Weather.f90` | Canopy-level sunlit/shaded PAR & NIR partitioning |
| `canopy_Photo.f90` | Canopy-level sunlit/shaded Vcmax25 |
| `leaf_step_mod.f90` | Per-hour, per-leaf-fraction fixed-point iteration coupling the modules below |
| `enbal_mod.f90` | Leaf energy balance (solves leaf temperature) |
| `stomata.f90` | Ball-Berry stomatal conductance |
| `c4_photosynth.f90` | C4 (Collatz-style) photosynthesis biochemistry |
| `boundary.f90` | Leaf boundary-layer conductance |
| `daily_carbon_mod.f90` | Hourly → daily integration of carbon gain and transpiration |

All physics modules use default `REAL` throughout (no double-precision
literals), so they compile without `-Wconversion` warnings.

## 11. Verifying the build and a run

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=<prefix> ..
cmake --build . -j8 && cmake --install .
```

The Debug build traps floating-point exceptions
(`-ffpe-trap=invalid,zero,overflow`) and checks array bounds, so a full
run of a maize experiment (e.g. `UFGA8201.MZX` from DSSAT's
`dssat-csm-data`, all six treatments: `dscsm048 A UFGA8201.MZX`) is also a
robustness check of the leaf solver. Note that DSSAT limits file paths to
~100 characters, so install and run from short paths.
