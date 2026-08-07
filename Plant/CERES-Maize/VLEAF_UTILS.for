      MODULE VLEAF_UTILS
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: MZ_RPARAMS, MZ_READ_NPOINTS, MZ_READ_DIURNAL_WEATHER
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
      CASE ("switch"); switch = val
      END SELECT
      GOTO 10

   20 CONTINUE
      CLOSE(99)
      END SUBROUTINE MZ_RPARAMS

C=====================================================================
C  MZ_READ_NPOINTS (F77 fixed-form, fixed file path)
C  Reads a single integer from "Daily_Points.txt"
C=====================================================================
      SUBROUTINE MZ_READ_NPOINTS(NPTS, IERR)
      IMPLICIT NONE
      INTEGER NPTS, IERR

      INTEGER UNP, IOS
      CHARACTER*30 FILEPATH
      PARAMETER (UNP = 76)

      NPTS = 0
      IERR = 0
      FILEPATH = 'Daily_Points.txt'

      OPEN(UNIT=UNP, FILE=FILEPATH, STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         IERR = 1
         RETURN
      ENDIF

      READ(UNP,*,IOSTAT=IOS) NPTS
      IF (IOS .NE. 0) THEN
         IERR = 2
         CLOSE(UNP)
         RETURN
      ENDIF

      CLOSE(UNP)
      RETURN
      END SUBROUTINE MZ_READ_NPOINTS

C=====================================================================
C  MZ_READ_DIURNAL_WEATHER  (F77 fixed-form; adjustable arrays)
C
C  - Caller passes NPTS (read earlier via MZ_READ_NPOINTS in MZ_VLEAF)
C  - Opens fixed climate file 'MY_CLIMATE_INPUT.txt'
C  - Reads space-delimited climate TXT (with one header line)
C  - Filters rows by YR and DOY; fills arrays up to NPTS
C
C  IERR:
C    0 = success
C    1 = cannot open climate file or header read error
C    3 = more rows for the day than NPTS (truncated)
C
C  Expected columns:
C    YR DOY HOUR ZEN LAT LON RG LNG PPT TAIR CA O2 EA WIND PRES
C=====================================================================
      SUBROUTINE MZ_READ_DIURNAL_WEATHER(YR, DOY,
     &  NPTS, NREC,
     &  HOUR_OUT, TAIR_OUT, EA_OUT,
     &  RG_OUT, LNG_OUT, PPT_OUT, WIND_OUT,
     &  PRESS_OUT, CA_OUT, O2_OUT,
     &  LAT_OUT,
     &  IERR)

C---- Dummy args (NPTS must appear before arrays using it) ----------
      IMPLICIT NONE
      INTEGER       YR, DOY
      INTEGER       NPTS
      INTEGER       NREC
      REAL          HOUR_OUT (NPTS)
      REAL          TAIR_OUT (NPTS)
      REAL          EA_OUT   (NPTS)
      REAL          RG_OUT   (NPTS)
      REAL          LNG_OUT  (NPTS)
      REAL          PPT_OUT  (NPTS)
      REAL          WIND_OUT (NPTS)
      REAL          PRESS_OUT(NPTS)
      REAL          CA_OUT   (NPTS)
      REAL          O2_OUT   (NPTS)
      REAL          LAT_OUT  (NPTS)
      INTEGER       IERR

C---- Locals ---------------------------------------------------------
      INTEGER UWD, IOS, YY, JD, HH
      REAL ZEN, LAT, LON, RG, LNG, PPT, TAIR, CA, O2, EA, WIND, PRES
      LOGICAL DONE
      CHARACTER*40 CLIMFILE
      PARAMETER (UWD = 77)

      NREC = 0
      IERR = 0
      DONE = .FALSE.
      CLIMFILE = 'MY_CLIMATE_INPUT.txt'

C---- Open climate file
      OPEN(UNIT=UWD, FILE=CLIMFILE, STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         IERR = 1
         RETURN
      ENDIF

C---- Skip header line
      READ(UWD,'(A)',IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         CLOSE(UWD)
         IERR = 1
         RETURN
      ENDIF

C---- Read loop
 100  CONTINUE
      IF (DONE) GOTO 200

      READ(UWD,*,IOSTAT=IOS) YY, JD, HH, ZEN, LAT, LON,
     &    RG, LNG, PPT, TAIR, CA, O2, EA, WIND, PRES

      IF (IOS .NE. 0) THEN
         DONE = .TRUE.
         GOTO 100
      ENDIF

      IF (YY .EQ. YR .AND. JD .EQ. DOY) THEN
         IF (NREC .LT. NPTS) THEN
            NREC = NREC + 1
            HOUR_OUT (NREC) = FLOAT(HH)
            TAIR_OUT (NREC) = TAIR
            EA_OUT   (NREC) = EA
            RG_OUT   (NREC) = RG
            LNG_OUT  (NREC) = LNG
            PPT_OUT  (NREC) = PPT
            WIND_OUT (NREC) = WIND
            PRESS_OUT(NREC) = PRES
            CA_OUT   (NREC) = CA
            O2_OUT   (NREC) = O2
            LAT_OUT  (NREC) = LAT
         ELSE
C           More rows than NPTS -> truncate store, continue scanning
            IERR = 3
         ENDIF
      ENDIF

      GOTO 100

 200  CONTINUE
      CLOSE(UWD)
      RETURN
      END SUBROUTINE MZ_READ_DIURNAL_WEATHER

      END MODULE VLEAF_UTILS
