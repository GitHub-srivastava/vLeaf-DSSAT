      MODULE VLEAF_UTILS
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: MZ_RPARAMS, MZ_VLEAF_HOURLY
      CONTAINS

C=======================================================================
C  READ_MZ_PARAMS: read stomatal params from NAMELIST file
C  File format example (MY_PARAMS.TXT):
C=======================================================================

      SUBROUTINE MZ_RPARAMS(intercept, slope, vpr25, vcmax25, jmax25,
     & vpmax25, theta, gbs, alpha, x, rd25, sco25, switch) ! Output

      IMPLICIT NONE
      REAL, INTENT(OUT) :: intercept, slope, vpr25, vcmax25, jmax25
      REAL, INTENT(OUT) :: vpmax25, theta, gbs, alpha, x
      REAL, INTENT(OUT) :: rd25, sco25
      INTEGER, INTENT(OUT) :: switch

      INTEGER :: ios
      CHARACTER(LEN=32) :: key
      REAL :: val

      ! Defaults in case file is missing or incomplete
      intercept = 0.008
      slope = 3.0
      vpr25 = 80.0
      vcmax25 = 55.0
      jmax25 = 350.0
      vpmax25 = 110.0
      theta = 0.7
      gbs = 0.003
      alpha = 0.0
      x = 0.4
      rd25 = 0.05
      sco25 = 2590.0
      switch = 1

      OPEN(UNIT=99, FILE='MY_PARAMS.TXT', STATUS='OLD', IOSTAT=ios)
      IF (ios /= 0) THEN
         RETURN   ! leave defaults
      END IF

   10 CONTINUE
      READ(99, *, IOSTAT=ios) key, val
      IF (ios /= 0) GOTO 20
      SELECT CASE (TRIM(key))
      CASE ("intercept"); intercept = val
      CASE ("slope"); slope = val
      CASE ("vpr25"); vpr25 = val
      CASE ("vcmax25"); vcmax25 = val
      CASE ("jmax25"); jmax25 = val
      CASE ("vpmax25"); vpmax25 = val
      CASE ("theta"); theta = val
      CASE ("gbs"); gbs = val
      CASE ("alpha"); alpha = val
      CASE ("x"); x = val
      CASE ("rd25"); rd25 = val
      CASE ("sco25"); sco25 = val
      CASE ("switch"); switch = NINT(val)
      END SELECT
      GOTO 10

   20 CONTINUE
      CLOSE(99)
      END SUBROUTINE MZ_RPARAMS

C=======================================================================
C  MZ_VLEAF_HOURLY
C  Builds vLeaf's hourly forcing arrays (K=1..TS) directly from DSSAT's
C  own WeatherType, instead of reading external text files. WEATHER is
C  filled once per simulation day by WEATHR/HMET (Weather/HMET.for)
C  before any Plant module runs, so TAIRHR/RADHR/WINDHR/RHUMHR/BETA/
C  FRDIFP/FRDIFR are already the real, site- and date-correct values
C  for today - vLeaf only needs to unpack and convert units.
C
C  CO2 (ppm) and ATMPRES (Pa) are treated as constant across the day:
C  CO2 is DSSAT's daily atmospheric CO2 (already available to the
C  caller); ATMPRES is estimated once per season from station elevation
C  (WEATHER%XELEV) since DSSAT does not track hourly air pressure.
C
C  Downwelling longwave (LW) has no DSSAT equivalent either, so it is
C  estimated hourly from air temperature and vapor pressure using a
C  standard clear-sky formulation (Brutsaert, 1975).
C
C  This subroutine is the ONLY place vLeaf reads WEATHER; everything
C  downstream works from the plain arrays returned here.
C=======================================================================
      SUBROUTINE MZ_VLEAF_HOURLY(WEATHER, CO2, ATMPRES,      !Input
     &  HOUR, TAIR, EA, LW, WIND, PRES, CA, ZENITH,          !Output
     &  PAR_DIR, PAR_DIF, NIR_DIR, NIR_DIF)                  !Output

      USE ModuleDefs
      IMPLICIT NONE
      TYPE (WeatherType), INTENT(IN) :: WEATHER
      REAL, INTENT(IN)  :: CO2, ATMPRES
      REAL, INTENT(OUT) :: HOUR(TS), TAIR(TS), EA(TS), LW(TS)
      REAL, INTENT(OUT) :: WIND(TS), PRES(TS), CA(TS), ZENITH(TS)
      REAL, INTENT(OUT) :: PAR_DIR(TS), PAR_DIF(TS)
      REAL, INTENT(OUT) :: NIR_DIR(TS), NIR_DIF(TS)

      INTEGER H
      REAL TAK, EA_KPA, EPSA, RG, PAR, NIR
      REAL VPSAT
      EXTERNAL VPSAT
      REAL, PARAMETER :: SIGMASB = 5.6703744E-08

C     Shortwave partitioned into PAR/NIR by the conventional 45/55
C     energy split. NOTE: WEATHER%PARHR is deliberately NOT used here -
C     it carries umol m-2 s-1, whereas vLeaf's leaf physics
C     (c4_photosynth.f90 applies its own 4.6 umol/J conversion)
C     expects PAR as an energy flux in W m-2.
      REAL, PARAMETER :: FPAR = 0.45, FNIR = 0.55

      DO H = 1, TS
        HOUR(H)   = REAL(H)
        TAIR(H)   = WEATHER % TAIRHR(H)
        WIND(H)   = MAX(WEATHER % WINDHR(H), 0.1)
        PRES(H)   = ATMPRES
        CA(H)     = CO2
        ZENITH(H) = 90.0 - WEATHER % BETA(H)

C       Vapor pressure of air [Pa] from RH and saturation vapor
C       pressure (VPSAT, Weather/HMET.for - already used by WEATHR).
        EA(H) = (WEATHER % RHUMHR(H) / 100.0) * VPSAT(TAIR(H))

C       Direct/diffuse split from DSSAT's own hourly diffuse fractions
C       (FRDIFP for PAR photons, FRDIFR for total shortwave), computed
C       by FRACD inside HMET. At night FRACD returns 1.0 for both and
C       RADHR is 0, so all four components correctly go to zero.
        RG  = MAX(WEATHER % RADHR(H), 0.0)
        PAR = FPAR * RG
        NIR = FNIR * RG

        PAR_DIR(H) = PAR * (1.0 - WEATHER % FRDIFP(H))
        PAR_DIF(H) = PAR - PAR_DIR(H)
        NIR_DIR(H) = NIR * (1.0 - WEATHER % FRDIFR(H))
        NIR_DIF(H) = NIR - NIR_DIR(H)

C       Clear-sky downwelling longwave [W m-2] (Brutsaert, 1975);
C       DSSAT has no native hourly LW estimate to draw on.
        TAK    = TAIR(H) + 273.15
        EA_KPA = MAX(EA(H), 0.0) / 1000.0
        EPSA   = 1.72 * (EA_KPA / MAX(TAK, 1.0))**(1.0/7.0)
        EPSA   = MIN(1.0, MAX(0.0, EPSA))
        LW(H)  = EPSA * SIGMASB * TAK**4
      ENDDO

      END SUBROUTINE MZ_VLEAF_HOURLY

      END MODULE VLEAF_UTILS
