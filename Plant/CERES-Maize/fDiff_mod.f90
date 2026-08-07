module fDiff_mod
  implicit none
  private
  public :: call_fDiff

contains

  subroutine call_fDiff(hour, doy, lat_deg, PAR, NIR, zenith, &
                        PAR_dir, PAR_dif, NIR_dir, NIR_dif)
    implicit none

    ! ---- inputs (instant) ----
    real,    intent(in)    :: hour      ! decimal hours
    integer, intent(in)    :: doy       ! day of year
    real,    intent(in)    :: lat_deg   ! latitude (degrees)
    real,    intent(in)    :: PAR       ! W m-2
    real,    intent(in)    :: NIR       ! W m-2

    ! ---- outputs (instant) ----
    real,    intent(out)   :: zenith    ! degrees; if <= 0, it will be computed
    real,    intent(out)   :: PAR_dir, PAR_dif, NIR_dir, NIR_dif

    ! ---- constants (matching MATLAB) ----
    real, parameter :: solarConstant = 1370.0d0
    real, parameter :: pi = 3.14159265d0
    real, parameter :: deg2rad = pi/180.0
    real, parameter :: rad2deg = 180.0/pi

    ! ---- locals ----
    real :: h, tilda, elev
    real :: So, So_safe, rho, K, SW_So
    real :: sElev, lat_rad
    real :: diffFrac

    lat_rad = lat_deg * deg2rad

    ! ----------------------------
    ! 1) compute zenith if missing
    !    (treat zenith <= 0 as missing)
    ! ----------------------------
    ! hour angle (radians)
    h = 2.0d0*pi*(hour - 12.0d0)/24.0d0

    ! declination-like angle used in your code (radians)
    tilda = asin( -sin(23.45d0*deg2rad) * cos(2.0d0*pi*(doy + 10.0d0)/365.0d0) )

    ! zenith in degrees
    zenith = acos( sin(lat_rad)*sin(tilda) + cos(lat_rad)*cos(tilda)*cos(h) ) * rad2deg

    ! elevation (deg)
    elev = 90.0d0 - zenith
    sElev = sin(elev*deg2rad)

    ! default
    diffFrac = 1.0d0

    ! when sun is very low
    if (zenith > 85.0d0) then
      diffFrac = 1.0d0
    else
      ! extraterrestrial on horizontal
      So = solarConstant * (1.0d0 + 0.033d0*cos( (360.0d0*doy/365.0d0)*deg2rad )) * sElev

      ! avoid divide-by-zero / negative So at night
      So_safe = max(So, 1.0d-9)

      rho = 0.847d0 - 1.61d0*sElev + 1.04d0*(sElev*sElev)
      K   = (1.47d0 - rho)/1.66d0

      SW_So = (PAR + NIR)/So_safe

      if (SW_So <= 0.22d0) then
        diffFrac = 1.0d0
      else if (SW_So <= 0.35d0) then
        diffFrac = 1.0d0 - 6.4d0*(SW_So - 0.22d0)*(SW_So - 0.22d0)
      else if (SW_So <= K) then
        diffFrac = 1.47d0 - 1.66d0*SW_So
      else
        diffFrac = rho
      end if
    end if

    ! keep it sane (finite + clamp)
    if (diffFrac /= diffFrac) diffFrac = 1.0d0  ! NaN check
    diffFrac = max(0.0d0, min(1.0d0, diffFrac))

    ! split PAR/NIR
    PAR_dir = (1.0d0 - diffFrac) * PAR
    PAR_dif = PAR - PAR_dir

    NIR_dir = (1.0d0 - diffFrac) * NIR
    NIR_dif = NIR - NIR_dir

  end subroutine call_fDiff

end module fDiff_mod
