C-----------------------------------------------------------------------
C  MZ_VLEAF: Leaf/canopy-scale coupled photosynthesis-stomatal
C  conductance-energy balance submodel ("vLeaf") for CERES-Maize.
C
C  Sub-daily driver: runs the sunlit/shaded canopy leaf model once per
C  hour (TS = 24 hourly steps, ModuleDefs.for) using DSSAT's own hourly
C  weather (WEATHER%TAIRHR/RADHR/WINDHR/RHUMHR/BETA/FRDIFP/FRDIFR,
C  filled once per day by WEATHR/HMET before any Plant module runs -
C  see VLEAF_UTILS.for/MZ_VLEAF_HOURLY) and integrates hourly carbon
C  gain and transpiration to daily totals (CARBO_vLeaf, EOPVLF).
C-----------------------------------------------------------------------

      SUBROUTINE MZ_VLEAF(DYNAMIC,
     & LAI, LAI_eff, YR, DOY, SWFAC, NSTRESS, gDM_day, EOPVLF,
     & WEATHER, CO2)
      USE ModuleDefs
      USE ModuleData
      USE VLEAF_UTILS
C     intercept/slope/switch are taken from VLEAF_PARAMS (a saved
C     module) rather than kept as locals: they are set once at
C     SEASINIT but read every INTEGR step, and subroutine locals are
C     not guaranteed to persist between calls without SAVE.
      USE VLEAF_PARAMS, ONLY: VLEAF_SET_PARAMS,
     &                        intercept, slope, switch
      USE can_Weather, only: canopy_weather
      USE canopy_Photo, only: sunshade_photo
      USE leaf_step_mod, only: leaf_step
      USE daily_carbon_mod, only: daily_carbon, daily_transpiration

      IMPLICIT NONE
      EXTERNAL GETLUN, WARNING

      INTEGER DYNAMIC, YR, DOY, K
      REAL LAI, tSoil, SWFAC, LAI_eff
      REAL NSTRESS, CO2
      TYPE (WeatherType) WEATHER
      TYPE (ControlType) CONTROL

      REAL molCO2_day, gC_day, gDM_day
      REAL molH2O_day, EOPVLF

C Canopy radiative constants
      REAL aPARLeaf, aNIRLeaf, Omega, alpha_deg
      PARAMETER (aPARLeaf  = 0.80)   ! leaf PAR absorptivity
      PARAMETER (aNIRLeaf  = 0.23)   ! leaf NIR absorptivity
      PARAMETER (Omega     = 1.00)   ! canopy clumping index
      PARAMETER (alpha_deg = 60.00)  ! mean leaf inclination angle, deg

C     Hourly fluxes beyond this magnitude (umol m-2 s-1) are treated as
C     a failed leaf solve rather than a physical result.
      REAL ANMAX
      PARAMETER (ANMAX = 80.0)

C     Scratch copies used only to read MY_PARAMS.TXT at SEASINIT and
C     hand the values to VLEAF_SET_PARAMS; the persistent values live
C     in VLEAF_PARAMS.
      REAL p_intercept, p_slope, p_vpr25, p_vcmax25, p_jmax25
      REAL p_vpmax25, p_theta, p_gbs, p_alpha, p_x, p_rd25, p_sco25
      REAL p_eta
      INTEGER p_switch

C     RUNINIT, SEASINIT, INTEGR, OUTPUT, SEASEND come from ModuleDefs

      INTEGER NREC, LUNHR
      REAL AtmPRES, ELEV

!    Leaf-solve bookkeeping. *_D counts the latest INTEGR pass of the
!    day (the one whose results are used); *_S accumulates the season.
!    NFAIL: leaf_step hit max_iter (last iterate kept).
!    NBAD : result was NaN or out of range (hour set to zero).
      INTEGER NSOLVE_D, NFAIL_D, NBAD_D
      INTEGER NSOLVE_S, NFAIL_S, NBAD_S
      CHARACTER*6, PARAMETER :: ERRKEY = 'VLEAF '
      CHARACTER*78 MSG(3)

!    Forcing variables (TS hourly steps, matches DSSAT's WEATHER)
      REAL HOUR(TS), tAir(TS), ea(TS)
      REAL LW(TS), wind(TS)
      REAL pressure(TS), ca(TS), zenith(TS)
      LOGICAL HROUT_INIT, sunFlag, shFlag
      DATA HROUT_INIT /.FALSE./

!    Output initialize
      REAL PAR_dir(TS), PAR_dif(TS)
      REAL NIR_dir(TS), NIR_dif(TS)
      REAL ci, eb, cb

      REAL LAIsun(TS), PAR_leaf_sun(TS)
      REAL NIR_leaf_sun(TS), sun_vcmax25(TS)
      REAL sunAnet(TS), sunGs(TS), sunGb(TS)
      REAL suntLeaf(TS), sunRn(TS), sunHf(TS)
      REAL sunLWe(TS), sunTr(TS)
      REAL sunLWnet(TS), sunLEf(TS)
      REAL sunRes(TS), shRes(TS)
      REAL sunerrA(TS), sunerrCi(TS), sunerrT(TS)
      INTEGER sunConv(TS)

      REAL LAIsh(TS), PAR_leaf_sh(TS)
      REAL NIR_leaf_sh(TS), sh_vcmax25(TS)
      REAL shAnet(TS), shGs(TS), shGb(TS)
      REAL shtLeaf(TS), shRn(TS), shHf(TS)
      REAL shLWe(TS), shTr(TS)
      REAL shLWnet(TS), shLEf(TS)
      REAL sherrA(TS), sherrCi(TS), sherrT(TS)
      INTEGER shConv(TS)

!    Save Output for output
      SAVE AtmPRES
      SAVE LAIsun, PAR_leaf_sun, NIR_leaf_sun
      SAVE sun_vcmax25, sunAnet, sunGs, sunGb
      SAVE suntLeaf, sunRn, sunHf, sunLWe
      SAVE sunTr, sunLWnet, sunLEf
      SAVE sunRes, shRes
      SAVE sunerrA, sunerrCi, sunerrT, sunConv

      SAVE LAIsh, PAR_leaf_sh, NIR_leaf_sh
      SAVE sh_vcmax25, shAnet, shGs, shGb
      SAVE shtLeaf, shRn, shHf, shLWe
      SAVE shTr, shLWnet, shLEf
      SAVE sherrA, sherrCi, sherrT, shConv

C     LW is written in the OUTPUT phase, so it must persist from INTEGR
      SAVE HROUT_INIT, NREC, HOUR, LW, LUNHR
      SAVE NSOLVE_D, NFAIL_D, NBAD_D, NSOLVE_S, NFAIL_S, NBAD_S


      IF(DYNAMIC.EQ.RUNINIT) THEN

        NREC = 0

        IF (.NOT. HROUT_INIT) THEN
        CALL GETLUN('VLEAFHR', LUNHR)
        OPEN(UNIT=LUNHR, FILE='VLEAF_HOURLY.OUT',
     &     STATUS='UNKNOWN')   ! overwrite at start of run
        WRITE(LUNHR,'(A)') 'RUN TRNO YR DOY HOUR LAIsun LAIsh ' //
     &       'Vcmax25_sun Vcmax25_sh sunAnet shAnet ' //
     &       'sunLE shLE sunH shH sunRn shRn sunTleaf shTleaf ' //
     &       'sunGs shGs sunGb shGb sunTr shTr SWFAC ' //
     &       'NSTRESS PARsun PARsh NIRsun NIRsh sunLWn shLWn ' //
     &       'sunRes shRes LWd sunEA shEA sunEC shEC ' //
     &       'sunET shET sunConv shConv'
        CLOSE(LUNHR)
        HROUT_INIT = .TRUE.
        ENDIF

      ELSEIF(DYNAMIC.EQ.SEASINIT) THEN

      NREC     = 0
      NSOLVE_D = 0
      NFAIL_D  = 0
      NBAD_D   = 0
      NSOLVE_S = 0
      NFAIL_S  = 0
      NBAD_S   = 0

      CALL MZ_RPARAMS(p_intercept, p_slope, p_vpr25, p_vcmax25,
     & p_jmax25, p_vpmax25, p_theta, p_gbs, p_alpha, p_x, p_rd25,
     & p_sco25, p_eta, p_switch)

      CALL VLEAF_SET_PARAMS(p_intercept, p_slope, p_vpr25,
     &                      p_vcmax25, p_jmax25, p_vpmax25,
     &                      p_theta, p_gbs, p_alpha,
     &                      p_x, p_rd25, p_sco25, p_eta, p_switch)

C     Atmospheric pressure from station elevation (standard
C     barometric approximation); DSSAT has no hourly pressure data
C     and station elevation does not change during a run. A missing
C     elevation (-99) is taken as sea level.
      ELEV = WEATHER%XELEV
      IF (ABS(ELEV + 99.0) .LT. 0.5) ELEV = 0.0
      AtmPRES = 101325.0 * (1.0 - 2.25577E-5*ELEV)**5.25588

      ELSEIF(DYNAMIC.EQ.INTEGR) THEN
        NREC = TS
        NSOLVE_D = 0
        NFAIL_D  = 0
        NBAD_D   = 0

C       Build today's hourly forcing directly from DSSAT's WEATHER
C       (filled by WEATHR/HMET before any Plant module runs).
        CALL MZ_VLEAF_HOURLY(WEATHER, CO2, AtmPRES,
     &       HOUR, tAir, ea, LW, wind, pressure, ca, zenith,
     &       PAR_dir, PAR_dif, NIR_dir, NIR_dif)

C      ci/eb/cb are declared only to receive leaf_step's INTENT(OUT)
C      state; leaf_step re-seeds them internally every call, so each
C      hour is solved independently (no warm start between hours).

       DO 110 K = 1, NREC

C         Soil surface temperature for the canopy longwave balance.
C         DSSAT's SRFTEMP is daily and is not passed to the maize
C         module, so an hourly estimate tracking air temperature is
C         used instead.
          tSoil = tAir(K) - 1

C            Gives the absorbed leaf level radiation
          CALL canopy_weather(LAI_eff, Omega, alpha_deg, zenith(K),
     &               PAR_dir(K), PAR_dif(K), NIR_dir(K), NIR_dif(K),
     &               aPARLeaf, aNIRLeaf, LAIsun(K), LAIsh(K),
     &               PAR_leaf_sun(K), PAR_leaf_sh(K), NIR_leaf_sun(K),
     &               NIR_leaf_sh(K))

C            Computes the sun-shade canopy average photosynthetic capacity
          CALL sunshade_photo(zenith(K), LAI, Omega, NSTRESS,
     &                           sun_vcmax25(K), sh_vcmax25(K))

C            Computation of leaf scale fluxes for sunlit and then shaded
C--- SUN: if LAI is tiny, skip and set all outputs to 0
          IF (LAIsun(K) .LE. 0.001) THEN
             sunAnet(K)   = 0.0
             sunTr(K)     = 0.0
             sunLEf(K)    = 0.0
             sunHf(K)     = 0.0
             sunGs(K)     = 0.0
             sunGb(K)     = 0.0
             sunRn(K)     = 0.0
             sunLWe(K)    = 0.0
             sunLWnet(K)  = 0.0
             suntLeaf(K)  = 0.0
             sunRes(K)    = 0.0
             sunerrA(K)   = 0.0
             sunerrCi(K)  = 0.0
             sunerrT(K)   = 0.0
             sunConv(K)   = 1

          ELSE
             CALL leaf_step(
     &         PAR_leaf_sun(K), NIR_leaf_sun(K), LW(K), tAir(K),
     &         tSoil, ea(K), pressure(K), ca(K), wind(K),
     &         switch, intercept, slope, SWFAC, sun_vcmax25(K),
     &         sunAnet(K), ci, sunGs(K), sunGb(K), eb, cb,
     &         suntLeaf(K), sunRn(K), sunHf(K), sunLEf(K),
     &         sunLWe(K), sunTr(K), sunRes(K),
     &         LAI_eff, Omega, sunLWnet(K),
     &         sunerrA(K), sunerrCi(K), sunerrT(K), sunFlag)
             NSOLVE_D = NSOLVE_D + 1

C--- SUN: an unconverged solve keeps its last iterate (far closer to
C    the answer than zero); it is counted and flagged in the output.
             IF (sunFlag) THEN
                sunConv(K) = 1
             ELSE
                sunConv(K) = 0
                NFAIL_D = NFAIL_D + 1
             ENDIF

C--- SUN: discard garbage (NaN or out of range).
C    NaN check in Fortran: (x .NE. x) is TRUE only for NaN
             IF ( (sunAnet(K)  .NE. sunAnet(K))  .OR.
     &            (sunTr(K)    .NE. sunTr(K))    .OR.
     &            (suntLeaf(K) .NE. suntLeaf(K)) .OR.
     &            (ABS(sunAnet(K)) .GT. ANMAX) ) THEN
                NBAD_D = NBAD_D + 1
                sunAnet(K)   = 0.0
                sunTr(K)     = 0.0
                sunLEf(K)    = 0.0
                sunHf(K)     = 0.0
                sunGs(K)     = 0.0
                sunGb(K)     = 0.0
                sunRn(K)     = 0.0
                sunLWe(K)    = 0.0
                sunLWnet(K)  = 0.0
                suntLeaf(K)  = 0.0
                sunRes(K)    = 0.0
             ENDIF
          ENDIF

C--- SHADE: if LAI is tiny, skip and set all outputs to 0
          IF (LAIsh(K) .LE. 0.001) THEN
             shAnet(K)   = 0.0
             shTr(K)     = 0.0
             shLEf(K)    = 0.0
             shHf(K)     = 0.0
             shGs(K)     = 0.0
             shGb(K)     = 0.0
             shRn(K)     = 0.0
             shLWe(K)    = 0.0
             shLWnet(K)  = 0.0
             shtLeaf(K)  = 0.0
             shRes(K)    = 0.0
             sherrA(K)   = 0.0
             sherrCi(K)  = 0.0
             sherrT(K)   = 0.0
             shConv(K)   = 1
          ELSE
             CALL leaf_step(
     &         PAR_leaf_sh(K), NIR_leaf_sh(K), LW(K), tAir(K),
     &         tSoil, ea(K), pressure(K), ca(K), wind(K),
     &         switch, intercept, slope, SWFAC, sh_vcmax25(K),
     &         shAnet(K), ci, shGs(K), shGb(K), eb, cb,
     &         shtLeaf(K), shRn(K), shHf(K), shLEf(K),
     &         shLWe(K), shTr(K), shRes(K),
     &         LAI_eff, Omega, shLWnet(K),
     &         sherrA(K), sherrCi(K), sherrT(K), shFlag)
             NSOLVE_D = NSOLVE_D + 1

C--- SHADE: an unconverged solve keeps its last iterate.
             IF (shFlag) THEN
                shConv(K) = 1
             ELSE
                shConv(K) = 0
                NFAIL_D = NFAIL_D + 1
             ENDIF

C--- SHADE: discard garbage (NaN or out of range).
             IF ( (shAnet(K)  .NE. shAnet(K))  .OR.
     &            (shTr(K)    .NE. shTr(K))    .OR.
     &            (shtLeaf(K) .NE. shtLeaf(K)) .OR.
     &            (ABS(shAnet(K)) .GT. ANMAX) ) THEN
                NBAD_D = NBAD_D + 1
                shAnet(K)   = 0.0
                shTr(K)     = 0.0
                shLEf(K)    = 0.0
                shHf(K)     = 0.0
                shGs(K)     = 0.0
                shGb(K)     = 0.0
                shRn(K)     = 0.0
                shLWe(K)    = 0.0
                shLWnet(K)  = 0.0
                shtLeaf(K)  = 0.0
                shRes(K)    = 0.0
             ENDIF
          ENDIF

 110   CONTINUE

C     Hourly -> daily totals (per m2 ground)
      CALL daily_carbon(NREC, LAIsun, LAIsh,
     &     sunAnet, shAnet, molCO2_day, gC_day, gDM_day)

      CALL daily_transpiration(NREC, LAIsun, LAIsh, sunTr, shTr,
     &     molH2O_day, EOPVLF)


      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN

      IF (NREC .GT. 0) THEN
         CALL GET(CONTROL)
         OPEN(UNIT=LUNHR, FILE='VLEAF_HOURLY.OUT',
     &        STATUS='OLD', POSITION='APPEND')
         DO K = 1, NREC
            WRITE(LUNHR,
     &        '(I4,1X,I4,1X,I4,1X,I3,1X,F6.1,37(1X,F16.3),2(1X,I2))')
     &        CONTROL%RUN, CONTROL%TRTNUM, YR, DOY, HOUR(K),
     &        LAIsun(K), LAIsh(K),
     &        sun_vcmax25(K), sh_vcmax25(K),
     &        sunAnet(K), shAnet(K),
     &        sunLEf(K),  shLEf(K),
     &        sunHf(K),   shHf(K),
     &        sunRn(K), shRn(K),
     &        suntLeaf(K), shtLeaf(K),
     &        sunGs(K),   shGs(K),
     &        sunGb(K),   shGb(K),
     &        sunTr(K),   shTr(K), SWFAC, NSTRESS,
     &        PAR_leaf_sun(K), PAR_leaf_sh(K),
     &        NIR_leaf_sun(K), NIR_leaf_sh(K),
     &        sunLWnet(K), shLWnet(K), sunRes(K),
     &        shRes(K), LW(K), sunerrA(K), sherrA(K),
     &        sunerrCi(K), sherrCi(K), sunerrT(K), sherrT(K),
     &        sunConv(K), shConv(K)
         END DO
         CLOSE(LUNHR)
      ENDIF

      NSOLVE_S = NSOLVE_S + NSOLVE_D
      NFAIL_S  = NFAIL_S  + NFAIL_D
      NBAD_S   = NBAD_S   + NBAD_D
      NSOLVE_D = 0
      NFAIL_D  = 0
      NBAD_D   = 0
      NREC = 0

      ELSEIF (DYNAMIC .EQ. SEASEND) THEN

      IF (NFAIL_S + NBAD_S .GT. 0) THEN
        WRITE(MSG(1),'(A,I7,A,I8,A)') 'vLeaf: ', NFAIL_S, ' of ',
     &    NSOLVE_S, ' hourly leaf solves did not converge'
        WRITE(MSG(2),'(A)')
     &    '  (last iterate kept; see sunConv/shConv in VLEAF_HOURLY.OUT)'
        WRITE(MSG(3),'(A,I7,A)') 'vLeaf: ', NBAD_S,
     &    ' hourly leaf solves were NaN/out of range and set to zero'
        CALL WARNING(3, ERRKEY, MSG)
      ENDIF

      ENDIF

      RETURN
      END
