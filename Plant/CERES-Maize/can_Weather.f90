module can_Weather
  implicit none
  private
  public :: canopy_weather

contains

  !===========================================================
  ! can_weather
  !
  ! For ONE time step:
  !   Inputs:
  !     LAI, Omega, alpha_deg, zenith_deg
  !     PAR_dir_in, PAR_dif_in, NIR_dir_in, NIR_dif_in   (W m-2)
  !     aPARLeaf, aNIRLeaf (leaf absorptivity)
  !
  !   Outputs:
  !     LAIsun, LAIsh
  !     PAR_leaf_sun, PAR_leaf_sh, NIR_leaf_sun, NIR_leaf_sh
  !===========================================================
  subroutine canopy_weather( &
      LAI, Omega, alpha_deg, zenith, &
      PAR_dir, PAR_dif, NIR_dir, NIR_dif, &
      aPARLeaf, aNIRLeaf, &
      LAIsun, LAIsh, &
      PAR_leaf_sun, PAR_leaf_sh, NIR_leaf_sun, NIR_leaf_sh)

    implicit none

    ! ---- inputs ----
    real, intent(in) :: LAI, Omega
    real, intent(in) :: alpha_deg, zenith
    real, intent(in) :: PAR_dir, PAR_dif, NIR_dir, NIR_dif
    real, intent(in) :: aPARLeaf, aNIRLeaf

    ! ---- outputs ----
    real, intent(out) :: LAIsun, LAIsh
    real, intent(out) :: PAR_leaf_sun, PAR_leaf_sh
    real, intent(out) :: NIR_leaf_sun, NIR_leaf_sh

    ! ---- locals ----
    real, parameter :: pi = 3.141d0
    real, parameter :: deg2rad = pi/180.0d0
    real, parameter :: tiny = 1.0d-12

    real :: cosTheta, mu, cosAlpha
    real :: cosThetaBar, E

    real :: r_PAR, t_PAR, alb_PAR
    real :: r_NIR, t_NIR, alb_NIR

    ! cos(theta) from zenith
    cosTheta = cos(zenith * deg2rad)

    ! cos(thetabar) = 0.537 + 0.025*LAI, clipped
    cosThetaBar = 0.537d0 + 0.025d0*LAI
    cosThetaBar = min(1.0d0, max(1.0d-3, cosThetaBar))

    ! E = exp(-0.5*Omega*LAI/cosThetaBar)
    E = exp(-0.5d0*Omega*LAI / cosThetaBar)

    ! Sunlit/shaded LAI split
    if (cosTheta <= 0.0d0) then
      mu     = 0.05d0
      LAIsun = 0.0d0
      LAIsh  = LAI
    else
      mu = max(0.05d0, cosTheta)
      LAIsun = 2.0d0*mu * (1.0d0 - exp(-0.5d0*Omega*LAI/mu))
      LAIsun = min(LAI, max(0.0d0, LAIsun))
      LAIsh  = max(0.0d0, LAI - LAIsun)
    end if

    cosAlpha = cos(alpha_deg * deg2rad)

    ! ---- PAR optics ----
    r_PAR   = (1.0d0 - sqrt(aPARLeaf)) / (1.0d0 + sqrt(aPARLeaf))
    t_PAR   = 1.0d0 - aPARLeaf - r_PAR
    alb_PAR = r_PAR + t_PAR

    ! ---- NIR optics ----
    r_NIR   = (1.0d0 - sqrt(aNIRLeaf)) / (1.0d0 + sqrt(aNIRLeaf))
    t_NIR   = 1.0d0 - aNIRLeaf - r_NIR
    alb_NIR = r_NIR + t_NIR

    ! ---- PAR band ----
    call sw_band(PAR_dir, PAR_dif, alb_PAR, LAI, Omega, mu, cosAlpha, E, &
                 PAR_leaf_sun, PAR_leaf_sh)

    ! ---- NIR band ----
    call sw_band(NIR_dir, NIR_dif, alb_NIR, LAI, Omega, mu, cosAlpha, E, &
                 NIR_leaf_sun, NIR_leaf_sh)

  end subroutine canopy_weather


  !===========================================================
  ! sw_band: matches your MATLAB nested function
  !===========================================================
  subroutine sw_band(Sdir, Sdif, alb, LAI, Omega, mu, cosAlpha, E, Rsun, Rsh)
    implicit none

    real, intent(in)  :: Sdir, Sdif, alb, LAI, Omega, mu, cosAlpha, E
    real, intent(out) :: Rsun, Rsh

    real, parameter :: tiny = 1.0d-12
    real :: one_minus_a
    real :: Rdir, Rdif, C0, LAI_safe

    one_minus_a = 1.0d0 - alb

    Rdir = one_minus_a * Sdir * (cosAlpha / max(mu, tiny))

    C0 = 0.07d0 * Omega * Sdir * (1.1d0 - 0.1d0*LAI) * exp(-mu)
    if (C0 < 0.0d0) C0 = 0.0d0

    LAI_safe = max(LAI, tiny)

    Rdif = one_minus_a * ( Sdif * (1.0d0 - E) / LAI_safe + C0 )

    Rsun = Rdir + Rdif
    Rsh  = Rdif

  end subroutine sw_band

end module can_Weather
