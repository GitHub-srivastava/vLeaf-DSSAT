C-----------------------------------------------------------------------
C  MZ_VLEAF: Stomatal conductance helper for CERES-Maize
C Initializes the variables at the start of season
C Compute COND using the LAI
C Daily average for the output
C-----------------------------------------------------------------------

      SUBROUTINE MZ_VLEAF(DYNAMIC,
     & LAI, YR, DOY) ! Output
      USE VLEAF_UTILS
      USE VLEAF_DRIVER
      USE VLEAF_PARAMS, ONLY: VLEAF_SET_PARAMS
!      USE fDiff_mod, only: call_fDiff
!      USE canopy_Photo, only: sunshade_photo
!      USE can_Weather, only: canopy_weather

      IMPLICIT NONE

      INTEGER DYNAMIC, YR, DOY, K
      REAL LAI

C Local variables
      REAL PAR, NIR, tSoil, aPARLeaf, aNIRLeaf
      REAL Omega, alpha_deg
      REAL intercept, slope, vpr25, vcmax25, jmax25
      REAL vpmax25, theta, gbs, alpha, x, rd25, sco25
      INTEGER switch

      INTEGER, PARAMETER :: RUNINIT = 1
      INTEGER, PARAMETER :: SEASINIT = 2
      INTEGER, PARAMETER :: INTEGR  = 4
      INTEGER, PARAMETER :: OUTPUT  = 5
      INTEGER, PARAMETER :: SEASEND = 6

!      SAVE intercept, slope, vpr25, vcmax25, jmax25
!      SAVE vpmax25, theta, gbs, alpha, x, rd25, sco25
!      SAVE switch

      INTEGER NPTS, IERRN, IERRW, NREC
      INTEGER MAXPTS
      PARAMETER (MAXPTS = 288)               ! safe upper bound

!    Forcing variables
      REAL HOUR(MAXPTS), tAir(MAXPTS), ea(MAXPTS)
      REAL Rg(MAXPTS), LW(MAXPTS), ppt(MAXPTS), wind(MAXPTS)
      REAL pressure(MAXPTS), ca(MAXPTS), O2(MAXPTS)
      LOGICAL HROUT_INIT
      DATA HROUT_INIT /.FALSE./


!    Output initialize
      REAL zenith(MAXPTS), PAR_dir(MAXPTS), PAR_dif(MAXPTS)
      REAL NIR_dir(MAXPTS), NIR_dif(MAXPTS)
      REAL ANET(MAXPTS), TRANS(MAXPTS)
      REAL TLEAF_HR(MAXPTS), LE_HR(MAXPTS), H_HR(MAXPTS)
      REAL LAIsun(MAXPTS), LAIsh(MAXPTS), PAR_leaf_sun(MAXPTS)
      REAL PAR_leaf_sh(MAXPTS), NIR_leaf_sun(MAXPTS)
      REAL NIR_leaf_sh(MAXPTS), sun_vcmax25(MAXPTS)
      REAL sh_vcmax25(MAXPTS)
      REAL LAT(MAXPTS)

      DOUBLE PRECISION RG_D, EA_D, TA_D, WIND_D, PRES_D, CA_D, O2_D
      DOUBLE PRECISION AN_D, TR_D, TLF_D, LE_D, H_D

!    Save Output for output
      SAVE ANET, TRANS
      SAVE TLEAF_HR, LE_HR, H_HR
      SAVE LAIsun, LAIsh, PAR_leaf_sun
      SAVE PAR_leaf_sh, NIR_leaf_sun
      SAVE NIR_leaf_sh, sun_vcmax25
      SAVE sh_vcmax25
      SAVE LAT

      SAVE HROUT_INIT, NREC, HOUR
      INTEGER ISTAT


      IF(DYNAMIC.EQ.RUNINIT) THEN

        IF (.NOT. HROUT_INIT) THEN
        OPEN(UNIT=97, FILE='VLEAF_HOURLY.OUT',
     &     STATUS='UNKNOWN')   ! overwrite at start of run
!        WRITE(97,'(A)') 'YR DOY HOUR LAIsun LAIsh PARsun PARsh NIRsun '
!        WRITE(97,'(A)') 'NIRsh ANET TRANS TLEAF LE H'
        WRITE(97,'(A)') 'YR DOY HOUR ANET TRANS TLEAF LE H'

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
     &     NPTS, NREC, HOUR, tAir, ea, Rg,
     &     LW, ppt, wind, pressure, ca, O2, LAT, IERRW)
        IF (IERRW .EQ. 1) RETURN
*        IF (IERRW .EQ. 3) then continue; you just got truncated to NPTS

       aPARLeaf = 0.80d0
       aNIRLeaf = 0.23d0
       Omega = 1.00d0
       alpha_deg = 60.00d0

       DO 110 K = 1, NREC
C<<< NEW: compute hourly PAR (very simple) and hourly Anet >>>
          RG_D   = DBLE(Rg(K))
          EA_D   = DBLE(ea(K))
          TA_D   = DBLE(tAir(K))
          WIND_D = DBLE(wind(K))
          PRES_D = DBLE(pressure(K))
          CA_D   = DBLE(ca(K))
          O2_D   = DBLE(O2(K))


          PAR = 0.45d0 * RG_D
          NIR = 0.55d0 * RG_D
          tSoil = TA_D - 1

!          call call_fDiff(HOUR(K), DOY, LAT(K), PAR, NIR,
!     &         zenith(K), PAR_dir(K), PAR_dif(K),
!     &         NIR_dir(K), NIR_dif(K))


!          CALL canopy_weather(LAI, Omega, alpha_deg, zenith(K),
!     &               PAR_dir(K), PAR_dif(K), NIR_dir(K), NIR_dif(K),
!     &               aPARLeaf, aNIRLeaf, LAIsun(K), LAIsh(K),
!     &               PAR_leaf_sun(K), PAR_leaf_sh(K), NIR_leaf_sun(K),
!     &               NIR_leaf_sh(K))

            ! Sun-shade photosynthetic capacity
!          CALL sunshade_photo(zenith(K), LAI, Omega,
!     &                           sun_vcmax25(K), sh_vcmax25(K))

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
            WRITE(97,'(I4,1X,I3,1X,F6.1,11(1X,F12.3))')
     &        YR, DOY, HOUR(K), ANET(K), TRANS(K),
     &        TLEAF_HR(K), LE_HR(K), H_HR(K)
!            WRITE(97,'(I4,1X,I3,1X,F6.1,11(1X,F12.3))')
!     &        YR, DOY, HOUR(K), LAIsun(K), LAIsh(K), PAR_leaf_sun(K),
!     &        PAR_leaf_sh(K), NIR_leaf_sun(K), NIR_leaf_sh(K),
!     &        ANET(K), TRANS(K), TLEAF_HR(K), LE_HR(K),
!     &        H_HR(K)
         END DO
         CLOSE(97)
      ENDIF

        NREC = 0
      ENDIF

      RETURN
      END

