module c4_photosynth
  use VLEAF_PARAMS, only: rd25
  implicit none
  private
  public :: c4_photosynthesis

  ! C4 constants
  real, parameter :: kk_C4    = 0.7d0
  real, parameter :: theta_c4 = 0.83d0
  real, parameter :: beta_c4  = 0.93d0
  real, parameter :: alpha_c4 = 0.035d0
  real, parameter :: Q10_c4   = 2.0d0
  real, parameter :: convPAR  = 4.6d0   ! PAR -> Qabs (tune if needed)

contains

  !------------------------------------------------------------
  ! C4 photosynthesis at a single time step
  !
  ! Inputs:
  !   tLeaf    - leaf temperature [degC]
  !   ci_leaf  -intercellular CO2
  !   par_leaf - incident PAR at the leaf
  !
  ! Uses from vleaf_reader:
  !   vcmax25, rd25
  !
  ! Outputs:
  !   Anet  - net assimilation
  !   M_out - light-limited rate M from the first quadratic
  !------------------------------------------------------------
  subroutine c4_photosynthesis(tLeaf, vcmax25, ci_leaf, par_leaf, Anet)
    implicit none
    real, intent(in)  :: tLeaf, ci_leaf
    real, intent(in)  :: vcmax25, par_leaf
    real, intent(out) :: Anet

    real :: Q10s, vcmax, rd_loc, kT
    real :: Qabs
    real :: aa, bb, cc, disc, sqrt_disc
    real :: M, aGross_loc

    ! Q10 temperature scaling
    Q10s = Q10_c4 ** ((tLeaf - 25.0d0) / 10.0d0)

    ! Peaked temperature response for Vcmax
    vcmax = (vcmax25 * Q10s) / &
            ((1.0d0 + exp(0.3d0 * (13.0d0 - tLeaf))) * &
             (1.0d0 + exp(0.3d0 * (tLeaf - 36.0d0))))

    ! Rd formulation: rd = rd25 * vcmax
    rd_loc = rd25 * vcmax

    ! CO2-related rate constant
    kT = kk_C4 * Q10s

    ! Absorbed radiation
    Qabs = convPAR * par_leaf

    !---------------------------
    ! Light-limited C4 rate M
    ! theta M^2 - (vcmax + alpha Qabs) M + vcmax*alpha*Qabs = 0
    !---------------------------
    aa = theta_c4
    bb = -(vcmax + alpha_c4 * Qabs)
    cc =  vcmax * alpha_c4 * Qabs

    disc = bb*bb - 4.0d0*aa*cc
    if (disc < 0.0d0) disc = 0.0d0
    sqrt_disc = sqrt(disc)

    M = (-bb - sqrt_disc) / (2.0d0 * aa)

    !---------------------------
    ! Gross assimilation A_g
    ! beta A^2 - (M + kT*ci) A + M*kT*ci = 0
    !---------------------------
    aa = beta_c4
    bb = -(M + kT * ci_leaf)
    cc =  M * kT * ci_leaf

    disc = bb*bb - 4.0d0*aa*cc
    if (disc < 0.0d0) disc = 0.0d0
    sqrt_disc = sqrt(disc)

    aGross_loc = (-bb - sqrt_disc) / (2.0d0 * aa)

    ! Net assimilation
    Anet  = aGross_loc - rd_loc

  end subroutine c4_photosynthesis

end module c4_photosynth
