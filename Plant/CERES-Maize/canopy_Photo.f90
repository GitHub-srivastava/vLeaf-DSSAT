module canopy_Photo
  use VLEAF_PARAMS, only: vcmax25
  implicit none
  private
  public :: sunshade_photo

contains

  subroutine sunshade_photo( &
      zenith, LAI, Omega, NSTRESS, &
      sun_vcmax25, sh_vcmax25)

    implicit none

    ! ---- inputs ----
    real, intent(in)  :: zenith
    real, intent(in)  :: LAI, Omega
    real, intent(in)  :: NSTRESS

    ! ---- outputs ----
    real, intent(out) :: sun_vcmax25
    real, intent(out) :: sh_vcmax25

    ! ---- local constants (same as MATLAB) ----
    real, parameter :: kn   = 0.5
    real, parameter :: G    = 0.5
    real, parameter :: pi   = 3.14159265
    real, parameter :: deg2rad = pi/180.0

    ! ---- locals ----
    real :: mu_b, kb
    real :: denom_sun, numer_sun, factor_sun
    real :: numer_sh, denom_sh
    real :: exp_kbLAI, exp_knLAI, exp_kbknLAI
    real :: tiny

    tiny = 1.0e-12

    ! Same floor as can_Weather's mu, so the sunlit/shaded split used
    ! for Vcmax25 matches the one used for absorbed radiation
    mu_b = max(0.05, cos(zenith * deg2rad))

    ! kb = (G * Omega) / mu_b
    kb  = (G * Omega) / mu_b

    ! Precompute exponentials
    exp_kbLAI   = exp(-kb * LAI)
    exp_knLAI   = exp(-kn * LAI)
    exp_kbknLAI = exp(-(kb + kn) * LAI)

    ! ---------------------------
    ! Sunlit
    ! ---------------------------
    denom_sun = 1.0 - exp_kbLAI
    if (abs(denom_sun) < tiny) denom_sun = sign(tiny, denom_sun)

    numer_sun  = 1.0 - exp_kbknLAI
    factor_sun = (kb / (kb + kn)) * (numer_sun / denom_sun)

    sun_vcmax25 = vcmax25 * NSTRESS * factor_sun

    ! ---------------------------
    ! Shaded
    ! ---------------------------
    numer_sh = ( (1.0/kn) * (1.0 - exp_knLAI) ) - &
               ( (Omega/(kb + kn)) * (1.0 - exp_kbknLAI) )

    denom_sh = LAI - 2.0 * mu_b * (1.0 - exp_kbLAI)
    if (abs(denom_sh) < tiny) denom_sh = sign(tiny, denom_sh)

    sh_vcmax25 = vcmax25 * NSTRESS * (numer_sh / denom_sh)

  end subroutine sunshade_photo

end module canopy_Photo
