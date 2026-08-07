C-----------------------------------------------------------------------
C  MZ_VLEAF: Stomatal conductance helper for CERES-Maize
C Initializes the variables at the start of season
C Compute COND using the IPAR
C Daily average for the output
C-----------------------------------------------------------------------

      SUBROUTINE MZ_VLEAF(DYNAMIC,
     & LAI,YR, DOY) ! Output
      USE VLEAF_UTILS
      USE VLEAF_DRIVER
      USE VLEAF_PARAMS, ONLY: VLEAF_SET_PARAMS

      IMPLICIT NONE

      INTEGER DYNAMIC, YR, DOY, K
      REAL LAI

C Local variables
      REAL intercept, slope, vpr25, vcmax25, jmax25
      REAL vpmax25, theta, gbs, alpha, x, rd25, sco25
      SAVE intercept, slope, vpr25, vcmax25, jmax25
      SAVE vpmax25, theta, gbs, alpha, x, rd25, sco25
      INTEGER switch
      SAVE switch

      INTEGER, PARAMETER :: RUNINIT = 1
      INTEGER, PARAMETER :: SEASINIT = 2
      INTEGER, PARAMETER :: INTEGR  = 4
      INTEGER, PARAMETER :: OUTPUT  = 5
      INTEGER, PARAMETER :: SEASEND = 6

      INTEGER NPTS, IERRN, IERRW, NREC
      INTEGER MAXPTS
      PARAMETER (MAXPTS = 288)               ! safe upper bound
      REAL HOUR(MAXPTS), tAir(MAXPTS), ea(MAXPTS)
      REAL Rg(MAXPTS), LW(MAXPTS), ppt(MAXPTS), wind(MAXPTS)
      REAL pressure(MAXPTS), ca(MAXPTS), O2(MAXPTS)
      REAL LAT(MAXPTS)
      LOGICAL HROUT_INIT
      SAVE HROUT_INIT, NREC, HOUR
      SAVE LAT
      DATA HROUT_INIT /.FALSE./

      REAL ANET(MAXPTS), TRANS(MAXPTS)
      SAVE ANET, TRANS

      REAL TLEAF_HR(MAXPTS), LE_HR(MAXPTS), H_HR(MAXPTS)
      SAVE TLEAF_HR, LE_HR, H_HR

      DOUBLE PRECISION RG_D, EA_D, TA_D, WIND_D, PRES_D, CA_D, O2_D
      DOUBLE PRECISION AN_D, TR_D, TLF_D, LE_D, H_D
      INTEGER ISTAT


      IF(DYNAMIC.EQ.RUNINIT) THEN

        IF (.NOT. HROUT_INIT) THEN
        OPEN(UNIT=97, FILE='VLEAF_HOURLY.OUT',
     &     STATUS='UNKNOWN')   ! overwrite at start of run
        WRITE(97,'(A)') 'YR DOY HOUR LAT ANET TRANS TLEAF LE H'
        CLOSE(97)
        HROUT_INIT = .TRUE.
        ENDIF

      ELSEIF(DYNAMIC.EQ.SEASINIT) THEN

      CALL MZ_RPARAMS(intercept, slope, vpr25, vcmax25, jmax25,
     & vpmax25, theta, gbs, alpha, x, rd25, sco25, switch)

      CALL VLEAF_SET_PARAMS(DBLE(intercept), DBLE(slope), DBLE(vpr25),
     &                      DBLE(vcmax25), DBLE(jmax25), DBLE(vpmax25),
     &                      DBLE(theta), DBLE(gbs), DBLE(alpha),
     &                      DBLE(x), DBLE(rd25), DBLE(sco25), switch)

      ELSEIF(DYNAMIC.EQ.INTEGR) THEN
        NREC = 0

C 1) Read NPTS from Daily_Points.txt
      CALL MZ_READ_NPOINTS(NPTS, IERRN)
        IF (IERRN .NE. 0) RETURN
        IF (NPTS .GT. MAXPTS) RETURN

      CALL MZ_READ_DIURNAL_WEATHER(YR, DOY,
     &     NPTS, NREC,
     &     HOUR, tAir, ea, Rg, LW, ppt, wind, pressure, ca, O2,
     &     LAT,
     &     IERRW)
        IF (IERRW .EQ. 1) RETURN
*        IF (IERRW .EQ. 3) then continue; you just got truncated to NPTS

       DO 110 K = 1, NREC
          RG_D   = DBLE(Rg(K))
          EA_D   = DBLE(ea(K))
          TA_D   = DBLE(tAir(K))
          WIND_D = DBLE(wind(K))
          PRES_D = DBLE(pressure(K))
          CA_D   = DBLE(ca(K))
          O2_D   = DBLE(O2(K))

          CALL VLEAF_STEP(RG_D, EA_D, TA_D, WIND_D, PRES_D, CA_D, O2_D,
     &                    AN_D, TR_D, TLF_D, LE_D, H_D, ISTAT)

          IF (ISTAT .NE. 0) THEN
             AN_D  = 0.0D0
             TR_D  = 0.0D0
             TLF_D = TA_D
             LE_D  = 0.0D0
             H_D   = 0.0D0
          ENDIF

          ANET(K)     = REAL(AN_D)
          TRANS(K)    = REAL(TR_D)
          TLEAF_HR(K) = REAL(TLF_D)
          LE_HR(K)    = REAL(LE_D)
          H_HR(K)     = REAL(H_D)

 110   CONTINUE

      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN

      IF (NREC .GT. 0) THEN
         OPEN(UNIT=97, FILE='VLEAF_HOURLY.OUT',
     &        STATUS='OLD', POSITION='APPEND')
         DO K = 1, NREC
            WRITE(97,'(I4,1X,I3,1X,F6.1,1X,F10.4,5(1X,F12.3))')
     &        YR, DOY, HOUR(K), LAT(K), ANET(K), TRANS(K),
     &        TLEAF_HR(K), LE_HR(K), H_HR(K)
         END DO
         CLOSE(97)
      ENDIF

        NREC = 0
      ENDIF

      RETURN
      END
