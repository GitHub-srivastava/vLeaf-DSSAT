C-----------------------------------------------------------------------
C  MZ_VLEAF: Stomatal conductance helper for CERES-Maize
C Initializes the variables at the start of season
C Compute COND using the LAI
C Daily average for the output
C-----------------------------------------------------------------------

      SUBROUTINE MZ_VLEAF(DYNAMIC,
     & LAI, LAI_eff, YR, DOY, SWFAC, NSTRESS, gDM_day, EOPVLF)
      USE VLEAF_UTILS
      USE VLEAF_PARAMS, ONLY: VLEAF_SET_PARAMS
      USE fDiff_mod, only: call_fDiff
      USE can_Weather, only: canopy_weather
      USE canopy_Photo, only: sunshade_photo
      USE leaf_step_mod, only: leaf_step
      USE daily_carbon_mod, only: daily_carbon, daily_transpiration
!     USE VLEAF_DRIVER


      IMPLICIT NONE

      INTEGER DYNAMIC, YR, DOY, K
      REAL LAI, PAR, NIR, tSoil, SWFAC, LAI_eff
      REAL NSTRESS

! Declare these somewhere in MZ_VLEAF:
      REAL molCO2_day, gC_day, gDM_day
      REAL molH2O_day, EOPVLF
!      SAVE molCO2_day, gC_day, gDM_day

C Local variables
      REAL aPARLeaf, aNIRLeaf
      REAL Omega, alpha_deg
      REAL intercept, slope, vpr25, vcmax25, jmax25
      REAL vpmax25, theta, gbs, alpha, x, rd25, sco25
      INTEGER switch

      INTEGER, PARAMETER :: RUNINIT = 1
      INTEGER, PARAMETER :: SEASINIT = 2
      INTEGER, PARAMETER :: INTEGR  = 4
      INTEGER, PARAMETER :: OUTPUT  = 5
      INTEGER, PARAMETER :: SEASEND = 6

      INTEGER NPTS, IERRN, IERRW, NREC
      INTEGER MAXPTS
      PARAMETER (MAXPTS = 288)               ! safe upper bound

!    Forcing variables
      REAL HOUR(MAXPTS), tAir(MAXPTS), ea(MAXPTS)
      REAL Rg(MAXPTS), LW(MAXPTS), ppt(MAXPTS), wind(MAXPTS)
      REAL pressure(MAXPTS), ca(MAXPTS), O2(MAXPTS)
      LOGICAL HROUT_INIT, sunFlag, shFlag
      DATA HROUT_INIT /.FALSE./


!    Output initialize
      REAL zenith, PAR_dir, PAR_dif
      REAL NIR_dir, NIR_dif
      REAL ci, eb, cb, errA, errCi, errT
!      REAL Res
      REAL LAT(MAXPTS)

      REAL LAIsun(MAXPTS), PAR_leaf_sun(MAXPTS)
      REAL NIR_leaf_sun(MAXPTS), sun_vcmax25(MAXPTS)
      REAL sunAnet(MAXPTS), sunGs(MAXPTS), sunGb(MAXPTS)
      REAL suntLeaf(MAXPTS), sunRn(MAXPTS), sunHf(MAXPTS)
      REAL sunLWe(MAXPTS), sunTr(MAXPTS)
      REAL sunLWnet(MAXPTS), sunLEf(MAXPTS)
      REAL sunRes(MAXPTS), shRes(MAXPTS)
      REAL sunerrA(MAXPTS), sunerrCi(MAXPTS), sunerrT(MAXPTS)

      REAL LAIsh(MAXPTS), PAR_leaf_sh(MAXPTS)
      REAL NIR_leaf_sh(MAXPTS), sh_vcmax25(MAXPTS)
      REAL shAnet(MAXPTS), shGs(MAXPTS), shGb(MAXPTS)
      REAL shtLeaf(MAXPTS), shRn(MAXPTS), shHf(MAXPTS)
      REAL shLWe(MAXPTS), shTr(MAXPTS)
      REAL shLWnet(MAXPTS), shLEf(MAXPTS)
      REAL sherrA(MAXPTS), sherrCi(MAXPTS), sherrT(MAXPTS)

!    Values for the LW checks
      REAL LWd_used, TaK, ea_kPa, epsA
      REAL sigmasb
      PARAMETER (sigmasb = 5.6703744E-08)

!    Save Output for output
      SAVE LAIsun, PAR_leaf_sun, NIR_leaf_sun
      SAVE sun_vcmax25, sunAnet, sunGs, sunGb
      SAVE suntLeaf, sunRn, sunHf, sunLWe
      SAVE sunTr, sunLWnet, sunLEf
      SAVE sunRes, shRes
      SAVE sunerrA, sunerrCi, sunerrT

      SAVE LAIsh, PAR_leaf_sh, NIR_leaf_sh
      SAVE sh_vcmax25, shAnet, shGs, shGb
      SAVE shtLeaf, shRn, shHf, shLWe
      SAVE shTr, shLWnet, shLEf
      SAVE sherrA, sherrCi, sherrT


      SAVE HROUT_INIT, NREC, HOUR
      INTEGER ISTAT


      IF(DYNAMIC.EQ.RUNINIT) THEN

        IF (.NOT. HROUT_INIT) THEN
        OPEN(UNIT=97, FILE='VLEAF_HOURLY.OUT',
     &     STATUS='UNKNOWN')   ! overwrite at start of run
        WRITE(97,'(A)') 'YR DOY HOUR LAIsun LAIsh ' //
     &       'Vcmax25_sun Vcmax25_sh sunAnet shAnet ' //
     &       'sunLE shLE sunH shH sunRn shRn sunTleaf shTleaf ' //
     &       'sunGs shGs sunGb shGb sunTr shTr SWFAC ' //
     &       'NSTRESS PARsun PARsh NIRsun NIRsh sunLWn ' //
     &       'sunRes shRes LWd sunEA shEA sunEC shEC ' //
     &       'shLWn sunET shET'


        CLOSE(97)
        HROUT_INIT = .TRUE.
        ENDIF

      ELSEIF(DYNAMIC.EQ.SEASINIT) THEN

      CALL MZ_RPARAMS(intercept, slope, vpr25, vcmax25, jmax25,
     & vpmax25, theta, gbs, alpha, x, rd25, sco25, switch)

      CALL VLEAF_SET_PARAMS(intercept, slope, vpr25,
     &                      vcmax25, jmax25, vpmax25,
     &                      theta, gbs, alpha,
     &                      x, rd25, sco25, switch)


      ELSEIF(DYNAMIC.EQ.INTEGR) THEN
        NREC = 0

C 1) Read NPTS from Daily_Points.txt
      CALL MZ_READ_NPOINTS(NPTS, IERRN)
        IF (IERRN .NE. 0) RETURN
        IF (NPTS .GT. MAXPTS) RETURN

      CALL MZ_READ_DIURNAL_WEATHER(YR, DOY,
     &     NPTS, NREC, HOUR, tAir, ea, Rg,
     &     LW, ppt, wind, pressure, ca, O2, LAT, IERRW)
        IF (IERRW .EQ. 1) RETURN
*        IF (IERRW .EQ. 3) then continue; you just got truncated to NPTS

       aPARLeaf = 0.80d0
       aNIRLeaf = 0.23d0
       Omega = 1.00d0
       alpha_deg = 60.00d0
       errA  = 1.0d0
       errCi = 1.0d0
       errT  = 1.0d0
       ci = 1.0d0
       eb = 1.0d0
       cb = 1.0d0

       DO 110 K = 1, NREC
C<<< NEW: compute hourly PAR (very simple) and hourly Anet >>


          PAR = 0.45d0 * MAX(Rg(K),0.0)
          NIR = 0.55d0 * MAX(Rg(K),0.0)
          tSoil = tAir(K) - 1
          wind(K) = MAX(wind(K),0.1)

C----- LW sanity check + fallback LWd from Drewry/Brutsaert-type epsA
C     Expect LWd ~ 200–550 W m-2 typically. Negative means "net LW", not incoming.
          LWd_used = LW(K)

C         If LW looks like net LW (negative) or otherwise unrealistic, overwrite it.
          IF (LW(K) .LT. 50.0 .OR. LW(K) .GT. 700.0) THEN
             TaK    = tAir(K) + 273.15
             ea_kPa = MAX(ea(K), 0.0) / 1000.0

C            Brutsaert-style clear-sky emissivity (common in canopy/LSM codes)
             epsA = 1.72 * (ea_kPa / MAX(TaK, 1.0))**(1.0/7.0)
             epsA = MIN(1.0, MAX(0.0, epsA))

             LW(K) = epsA * sigmasb * TaK**4
          ENDIF

C            Computation for the direct and diffused fraction
          call call_fDiff(HOUR(K), DOY, LAT(K), PAR, NIR,
     &         zenith, PAR_dir, PAR_dif,
     &         NIR_dir, NIR_dif)

C            Gives the absorbed leaf level radiation
          CALL canopy_weather(LAI_eff, Omega, alpha_deg, zenith,
     &               PAR_dir, PAR_dif, NIR_dir, NIR_dif,
     &               aPARLeaf, aNIRLeaf, LAIsun(K), LAIsh(K),
     &               PAR_leaf_sun(K), PAR_leaf_sh(K), NIR_leaf_sun(K),
     &               NIR_leaf_sh(K))

C            Computes the sun-shade canopy average photosynthetic capacity
          CALL sunshade_photo(zenith, LAI, Omega, NSTRESS,
     &                           sun_vcmax25(K), sh_vcmax25(K))

C            Computation of leaf scale fluxes for sunlit and then shaded
C--- SUN: if LAI is tiny, skip and set all outputs to 0
          IF (LAIsun(K) .LE. 0.001d0) THEN
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
          ENDIF
C--- SUN: check convergence flag + garbage (NaN or huge)
C    NaN check in Fortran: (x .NE. x) is TRUE only for NaN
          IF ( (.NOT. sunFlag) .OR.
     &         (sunAnet(K) .NE. sunAnet(K)) .OR.
     &         (sunTr(K)   .NE. sunTr(K))   .OR.
     &         (ABS(sunAnet(K)) .GT. 80.0) ) THEN

                sunAnet(K)   = 0.0
                sunTr(K)     = 0.0
C--- optional but recommended: keep flux bookkeeping consistent
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


          IF (LAIsh(K) .LE. 0.001d0) THEN
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
          ENDIF
C--- SUN: check convergence flag + garbage (NaN or huge)
C    NaN check in Fortran: (x .NE. x) is TRUE only for NaN
          IF ( (.NOT. shFlag) .OR.
     &         (shAnet(K) .NE. shAnet(K)) .OR.
     &         (shTr(K)   .NE. shTr(K))   .OR.
     &         (ABS(shAnet(K)) .GT. 80.0) ) THEN

                shAnet(K)   = 0.0
                shTr(K)     = 0.0
C--- optional but recommended: keep flux bookkeeping consistent
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

 110   CONTINUE

! Then at the end of OUTPUT (after hourly writes), call:
      CALL daily_carbon(NREC, HOUR, LAIsun, LAIsh,
     &     sunAnet, shAnet, molCO2_day, gC_day, gDM_day)

      CALL daily_transpiration(NREC, HOUR, LAIsun, LAIsh, sunTr, shTr,
     &     molH2O_day, EOPVLF)


      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN

      IF (NREC .GT. 0) THEN
         OPEN(UNIT=97, FILE='VLEAF_HOURLY.OUT',
     &        STATUS='OLD', POSITION='APPEND')
         DO K = 1, NREC
! OUTPUT block: write the matching values (hour + 18 reals)
            WRITE(97,'(I4,1X,I3,1X,F6.1,37(1X,F16.3))')
     &        YR, DOY, HOUR(K),
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
     &        sunerrCi(K), sherrCi(K), sunerrT(K), sherrT(K)
         END DO
         CLOSE(97)
      ENDIF

      NREC = 0
      ENDIF

      RETURN
      END

