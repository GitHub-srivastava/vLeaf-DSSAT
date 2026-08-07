module VLEAF_PARAMS
  implicit none
  save

  logical, public :: params_set = .false.

  ! keep DSSAT names
  real, public :: intercept, slope, vpr25
  real, public :: vcmax25, jmax25, vpmax25
  real, public :: theta, gbs, alpha, x
  real, public :: rd25, sco25
  integer, public :: switch

  contains

  subroutine VLEAF_SET_PARAMS(intercept_in, slope_in, vpr_in, vcmax_in, jmax_in, vpmax_in, &
                              theta_in, gbs_in, alpha_in, x_in, rd25_in, sco25_in, switch_in)
    implicit none
    real, intent(in) :: intercept_in, slope_in, vpr_in
    real, intent(in) :: vcmax_in, jmax_in, vpmax_in
    real, intent(in) :: theta_in, gbs_in, alpha_in, x_in
    real, intent(in) :: rd25_in, sco25_in
    integer, intent(in) :: switch_in

    intercept = intercept_in
    slope     = slope_in
    vpr25       = vpr_in
    vcmax25     = vcmax_in
    jmax25      = jmax_in
    vpmax25     = vpmax_in
    theta     = theta_in
    gbs       = gbs_in
    alpha     = alpha_in
    x         = x_in
    rd25      = rd25_in
    sco25     = sco25_in
    switch    = switch_in

    params_set = .true.
  end subroutine VLEAF_SET_PARAMS

end module VLEAF_PARAMS
