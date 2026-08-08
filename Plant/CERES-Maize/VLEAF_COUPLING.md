# vLeaf ↔ CERES-Maize coupling

This document describes how the vLeaf leaf/canopy photosynthesis-stomatal
conductance-energy balance submodel is wired into CERES-Maize, why it is
wired that way, and what its remaining simplifications are.

It covers four rounds of work on the integration: repairing the build
(`CMakeLists.txt` referenced source files that do not exist), removing a
broken water-stress hook that read a file nothing ever wrote, replacing
vLeaf's ad hoc external weather files with DSSAT's own hourly weather
(§3 — the substantive coupling change), and fixing three state-persistence
bugs that made vLeaf's behaviour undefined at run time (§6).

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
CARBO = CARBO_vLeaf / PLTPOP * SLPF
```

`PRFT` (temperature response) and the RUE multiplication are no longer
applied for this pathway — vLeaf's own `c4_photosynth.f90` has its own
Vcmax temperature response, so temperature stress is now handled inside the
mechanistic model rather than by the empirical `PRFTC` cultivar curve. Water
stress (`SWFAC`) and nutrient stress (`NSTRESS`, combining N/P/K) are still
applied, but as inputs *into* vLeaf's stomatal conductance rather than as
multipliers on the output.

## 2. Call chain

```
MZ_CERES.for  (per-crop dispatcher; already holds WEATHER, CO2)
     |
     v
MZ_GROSUB.for  (DYNAMIC = RUNINIT / SEASINIT / INTEGR / OUTPUT)
     |  CALL GET(WEATHER)                       [ModuleData accessor]
     |  CALL MZ_VLEAF(..., WEATHER, CO2)         x4 call sites
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
CARBO_vLeaf (g DM/m2/d), EOPVLF (mm/d)  -->  back into MZ_GROSUB's CARBO, EOP
```

`VLEAF_PARAMS.f90` holds the season's calibration constants (Vcmax25, Jmax25,
Ball-Berry intercept/slope, etc.), set once at `SEASINIT` from
`MZ_RPARAMS`/`MY_PARAMS.TXT` (or built-in defaults if that file is absent),
and read by `c4_photosynth.f90`/`canopy_Photo.f90` via `USE VLEAF_PARAMS`.

## 3. Weather sourcing: the actual coupling fix

**Before:** `MZ_VLEAF` read its own hourly weather from two ad hoc text
files, `Daily_Points.txt` (a record count) and `MY_CLIMATE_INPUT.txt`
(hourly rows filtered by year/day), opened and parsed fresh on every call.
This meant vLeaf could not run inside a normal DSSAT simulation without
hand-prepared side files matching the simulation's actual dates and
location, and it duplicated zenith-angle/diffuse-fraction math DSSAT
already performs elsewhere. It was not a real coupling — it was two models
running side by side, coordinated only by the user manually keeping two
sets of input files in sync.

**After:** vLeaf reads directly from DSSAT's own `WeatherType` (declared in
`Utilities/ModuleDefs.for`), which every DSSAT run already populates once
per simulated day, before any Plant module executes:

```
Weather/WEATHR.for  -->  Weather/HMET.for  -->  WEATHER%TAIRHR/RADHR/WINDHR/
                                                 RHUMHR/BETA/AZZON/FRDIFP/
                                                 FRDIFR (all DIMENSION(TS))
                     -->  CALL PUT(WEATHER)      [ModuleData]
```

`MZ_GROSUB.for` retrieves it with `CALL GET(WEATHER)` (the standard DSSAT
accessor pattern, `Utilities/ModuleDefs.for`'s `ModuleData` module) rather
than adding `WEATHER` to its already-large argument list, and passes it
straight through to `MZ_VLEAF` alongside `CO2` (DSSAT's daily atmospheric
CO2, already an `MZ_GROSUB` argument).

`VLEAF_UTILS.for`'s `MZ_VLEAF_HOURLY` subroutine is the **single place** that
reads `WEATHER`. It converts DSSAT's fields into plain per-hour driver arrays,
and everything downstream (`MZ_VLEAF`'s hourly loop and all the physics
modules) works only from those arrays — no `WEATHER%` access leaks past this
boundary. That keeps the DSSAT-facing surface of vLeaf to exactly one
subroutine.

| vLeaf needs (hourly) | Sourced from | Notes |
|---|---|---|
| Air temperature | `WEATHER%TAIRHR(H)` | °C, direct |
| Wind speed | `WEATHER%WINDHR(H)` | m/s, floored at 0.1 |
| Vapor pressure of air | `WEATHER%RHUMHR(H)/100 * VPSAT(TAIRHR(H))` | `VPSAT` is DSSAT's own saturation-vapor-pressure function (`Weather/HMET.for`), reused rather than re-implemented a fourth time |
| Solar zenith angle | `90.0 - WEATHER%BETA(H)` | `BETA` (solar elevation) already computed daily by `HANG` inside `HMET`; no more re-deriving declination/hour-angle from lat/DOY in vLeaf |
| Direct/diffuse PAR and NIR | `WEATHER%RADHR(H)` split 45/55 into PAR/NIR, then each split by `WEATHER%FRDIFP(H)` / `WEATHER%FRDIFR(H)` | diffuse fractions are Spitters (1986), computed by `FRACD` in `HMET.for`. At night `FRACD` returns 1.0 for both **and** `RADHR` is 0, so all four components correctly go to zero |
| Atmospheric CO2 | `CO2` (already a `MZ_GROSUB` argument) | ppm, constant across the day |
| Atmospheric pressure | Computed once at `SEASINIT` from `WEATHER%XELEV` via the standard barometric formula | DSSAT has no hourly pressure field; station elevation does not change during a run |
| Downwelling longwave | Estimated each hour from air temperature and vapor pressure (Brutsaert clear-sky formula) | DSSAT has no native hourly LW estimate either; this was already vLeaf's fallback path, now it is the only path |

Two DSSAT-native fields that were being read from the old text files turned
out to be **dead** once traced through the call chain and were dropped
rather than replicated: `O2` (never passed into `leaf_step`/`c4_photosynth`
— the C4 model here doesn't depend on it) and `ppt`/hourly precipitation
(read, stored, never consumed).

### A unit trap worth knowing about: why not `WEATHER%PARHR`?

DSSAT does carry an hourly PAR array, `WEATHER%PARHR`, and it looks like the
obvious source for vLeaf's PAR. It is deliberately **not** used, because the
two models use different units for "PAR":

- `WEATHER%PARHR` is a **photon flux density**, µmol m⁻² s⁻¹ (see the unit
  note in `HPAR`, `Weather/HMET.for`).
- vLeaf's leaf physics wants PAR as an **energy flux**, W m⁻²:
  `c4_photosynth.f90` applies its own `convPAR = 4.6` µmol/J conversion
  internally (`Qabs = convPAR * par_leaf`).

Feeding `PARHR` straight in would therefore double-convert and inflate
absorbed radiation by roughly 4.6×. The adapter instead takes `RADHR`
(J m⁻² s⁻¹ = W m⁻², an energy flux) and applies the conventional 45/55
PAR/NIR energy split, which is what the original vLeaf code did and what its
calibration assumes.

*Possible refinement (not applied):* using `PARHR/4.6` instead of `0.45*RADHR`
would honour the actual PAR/SRAD ratio from the weather file on days when it
is measured, rather than assuming a fixed 45%. That is a change in model
behaviour, not a bug fix, so it is left as a calibration decision.

Because DSSAT's weather is fixed at `TS = 24` hourly steps (`ModuleDefs.for`),
vLeaf's hourly loop is now a fixed `DO K = 1, TS` instead of a
variable-length `DO K = 1, NREC` bounded by a record count read from
`Daily_Points.txt`. All of the file-open/IOSTAT/truncation-handling code
that existed only to cope with an external file's row count is gone.

## 4. What canopy-radiation physics vLeaf still does itself

DSSAT's `HMET`/`FRACD` compute the *atmospheric* direct/diffuse split (how
much of today's incoming radiation is beam vs. sky-diffuse). vLeaf's own
`can_Weather.f90` (`canopy_weather`/`sw_band`) and `canopy_Photo.f90`
(`sunshade_photo`) are a distinct, lower level of the model: given that
already-split incoming radiation, they compute how it is partitioned
*inside the canopy* between sunlit and shaded leaf area (a Campbell &
Norman-style two-stream/exponential-extinction treatment), and how average
photosynthetic capacity (Vcmax25) differs between the two leaf fractions.
CERES-Maize has no equivalent to this — it is the actual physiological
value vLeaf adds — so this part was left untouched by the rewire.

## 5. The two water-stress passes are intentional, not duplicated work

`MZ_GROSUB.for` calls `MZ_VLEAF` twice inside the same `DYNAMIC.EQ.INTEGR`
day:

1. First with `SWFAC = 1.0` (unstressed) to get today's **potential**
   evaporative demand, `EOPVLF`. This becomes `EOP`, which is what DSSAT's
   soil-water-uptake routine compares against `TRWUP` (actual root water
   uptake) to compute the real `SWFAC` for today.
2. Then again with the corrected `SWFAC` to get the **actual**,
   water-stressed `CARBO_vLeaf`/`EOPVLF`/hourly fluxes that get written out
   and used for the rest of today's growth calculations.

This is required because, unlike the old RUE-based `PCARB` (a pure demand
term with `SWFAC` applied afterward as a multiplier), vLeaf couples
photosynthesis and transpiration through stomatal conductance — `SWFAC` is
an *input* to the stomatal model, not a multiplier on its output, so the
potential rate has to be computed before the actual rate can be. Before
this fix, both passes also re-opened and re-parsed two text files from disk,
making the cost of this necessary two-pass design much higher than it
needed to be. Now that weather comes from an in-memory `WEATHER` struct
(no I/O), the second pass costs only the leaf-level iterative solve itself.

## 6. State-persistence bugs found and fixed in the coupling audit

DSSAT calls each Plant module repeatedly with a `DYNAMIC` phase argument
(`RUNINIT` → `SEASINIT` → `INTEGR`/`OUTPUT` each day → `SEASEND`). A value
computed in one phase and read in a later phase must therefore survive
*between subroutine calls*, which in Fortran means it needs `SAVE` (or must
live in a module). Three variables in `MZ_VLEAF` violated this:

- **`intercept`, `slope`, `switch` (critical).** These were subroutine
  locals, assigned once at `SEASINIT` from `MY_PARAMS.TXT`, then read on
  every `INTEGR` step and passed into `leaf_step`. Without `SAVE` their
  values were not guaranteed to persist, so the Ball-Berry conductance
  parameters and the energy-balance `switch` could be garbage at run time.
  A wrong `switch` is especially damaging: `enbal_mod.f90` uses it to decide
  whether to solve the leaf energy balance at all, or to short-circuit and
  force leaf temperature to air temperature. **Fix:** these values already
  had a correct, persistent home — the `VLEAF_PARAMS` module (which is
  `save`d and is what `VLEAF_SET_PARAMS` writes to). `MZ_VLEAF` now
  `USE`s them from there, and the `SEASINIT` read uses clearly-named
  `p_*` scratch copies. This removes the duplicate state rather than just
  bolting `SAVE` onto it.
- **`LW` (hourly longwave).** Computed at `INTEGR`, but written to
  `VLEAF_HOURLY.OUT` during the later `OUTPUT` phase, without `SAVE` —
  so the longwave column of the hourly output could be garbage. Every other
  array in that same `WRITE` was already `SAVE`d; `LW` was simply missed.
  **Fix:** added to the `SAVE` list.
- **`NREC`.** Read as a loop bound at `OUTPUT` but only assigned at
  `INTEGR`, so an `OUTPUT` call before the first `INTEGR` would test an
  undefined value. **Fix:** initialised to 0 at `RUNINIT`.

Also removed while auditing: `errA`/`errCi`/`errT` scalars in `MZ_VLEAF`
(assigned but never read — `leaf_step` returns its own convergence errors
into the `sunerrA`/`sherrA` arrays), and the meaningless pre-seeding of
`ci`/`eb`/`cb` (they are `INTENT(OUT)` in `leaf_step`, which re-seeds them
internally on every call — so each hour is solved from a cold start, with
no warm-start carry-over between hours).

## 7. Known simplifications (by design, not oversights)

- **The daily integration spans 23 h, not 24 h.** `daily_carbon` and
  `daily_transpiration` integrate with the trapezoid rule over hourly
  samples `HOUR = 1..24`, which covers the 23 intervals from hour 1 to
  hour 24 — the wrap-around interval (hour 24 → hour 1 of the next day) is
  not counted. Because both endpoints are night-time values (`Anet ≈ 0`),
  the effect on daily carbon is negligible; the effect on `EOPVLF` is small
  but systematically low. This behaviour is unchanged from the original
  file-driven version (which had the same gap with `HOUR = 0..23`), so it
  has been left alone rather than silently altered — changing to a
  rectangle sum over all 24 h would shift the water balance and require
  recalibration.

- **Atmospheric pressure and O2** are not modeled hourly by DSSAT at all;
  pressure is a one-time elevation-based estimate, and O2 is unused by the
  current C4 photosynthesis formulation. If a future version of vLeaf needs
  O2-sensitivity, an hourly value will have to be synthesized (DSSAT has no
  native source for it).
- **Downwelling longwave radiation** has no DSSAT-native hourly
  representation; vLeaf's own clear-sky estimate is the only source.
- **`PRFTC` cultivar coefficients** (empirical temperature response) are
  read from the species/cultivar file as before but are no longer applied
  to `CARBO` — see §1. They remain computed (`PRFT` variable) so that
  removing them entirely doesn't ripple into other output variables that
  may still reference `PRFT` in the future; they are simply not multiplied
  into vLeaf's `CARBO`.
- **`MY_PARAMS.TXT`** (vLeaf's own calibration-constant file, read once at
  `SEASINIT` via `MZ_RPARAMS`) is intentionally left as-is. It is a small,
  self-contained override file for physiological constants (Vcmax25,
  Ball-Berry slope, etc.) that DSSAT's species/ecotype files don't carry
  yet, and is not part of the weather-coupling problem this rewire fixes.
  A natural next step, if vLeaf is promoted from research prototype to a
  supported option, is to move these constants into the `.SPE`/`.CUL` file
  format like every other CERES-Maize genetic coefficient.

## 8. File manifest

| File | Role |
|---|---|
| `MZ_VLEAF.for` | vLeaf's own DYNAMIC state machine; hourly driver loop |
| `VLEAF_UTILS.for` | `MZ_RPARAMS` (calibration file reader) + `MZ_VLEAF_HOURLY` (**the sole DSSAT boundary**: WEATHER → plain hourly arrays) |
| `VLEAF_PARAMS.f90` | Season-long calibration constants, set once, read by the leaf-level physics modules |
| `can_Weather.f90` | Canopy-level sunlit/shaded PAR & NIR partitioning |
| `canopy_Photo.f90` | Canopy-level sunlit/shaded Vcmax25 |
| `leaf_step_mod.f90` | Per-hour, per-leaf-fraction fixed-point iteration coupling the three modules below |
| `enbal_mod.f90` | Leaf energy balance (solves leaf temperature) |
| `stomata.f90` | Ball-Berry stomatal conductance |
| `c4_photosynth.f90` | C4 (von Caemmerer-style) photosynthesis biochemistry |
| `boundary.f90` | Leaf boundary-layer conductance |
| `daily_carbon_mod.f90` | Hourly → daily integration of carbon gain and transpiration |

`fDiff_mod.f90` (atmospheric zenith/diffuse-fraction calculator) and
`VLEAF_DRIVER.f90`, `MZ_VLEAF_new.for`, `MZ_VLEAF_old.for`,
`MZ_GROSUB_Old.for` (dead alternate/legacy copies, never referenced by
`CMakeLists.txt`) have been removed — the table above is the complete,
live set.

## 9. Verifying the build

```
mkdir build && cd build
cmake -G "MinGW Makefiles" ..     # or your platform's generator
cmake --build . -j4
```

A successful configure + build of the `dscsm048` target (no missing source
files, no undefined-symbol link errors) is the regression check for this
coupling — it was broken (missing source files in `CMakeLists.txt`) before
this fix and is confirmed working after it.

`MZ_VLEAF.for` and `VLEAF_UTILS.for` currently compile **warning-free** under
the project's default `-Wall -Wconversion` debug flags. Worth keeping that
way: the two warnings they used to emit were both genuine type problems
(a `REAL`→`INTEGER` narrowing when parsing the `switch` key from
`MY_PARAMS.TXT`, and `REAL(8)` literals assigned to `REAL(4)` variables).

Note that the remaining vLeaf physics modules (`stomata.f90`,
`boundary.f90`, `c4_photosynth.f90`, `can_Weather.f90`, `canopy_Photo.f90`,
`enbal_mod.f90`, `leaf_step_mod.f90`) still declare `real` (single
precision) while initialising with `d0` (double) literals throughout. This
is harmless in practice — the constants are rounded to single on assignment
— but it is inconsistent, and tidying it is a worthwhile separate pass.
It was deliberately not done here: touching literals inside the physics
kernels can perturb results in the last bits, which deserves its own
before/after numerical comparison rather than being bundled into a
coupling fix.
