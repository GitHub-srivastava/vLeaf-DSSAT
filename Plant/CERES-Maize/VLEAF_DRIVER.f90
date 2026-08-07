!=========================================================
! File: VLEAF_DRIVER.f90
!
! Fixes vs your current version:
!  - REMOVE dependency on vleaf_reader + REMOVE vleaf_init()
!  - VLEAF_STEP interface matches your MZ_VLEAF call:
!      CALL VLEAF_STEP(RG_D, EA_D, TA_D, WIND_D, PRES_D, CA_D, O2_D,
!     &                AN_D, TR_D, TLF_D, LE_D, H_D, ISTAT)
!  - Uses parameters stored once via VLEAF_SET_PARAMS (VLEAF_PARAMS)
!  - switch is INTEGER (comes from VLEAF_PARAMS)
!
! (Also: in MZ_VLEAF you must have: IF (IERRW .EQ. 1) RETURN  ; not .NE.)
!=========================================================
module VLEAF_DRIVER
  use VLEAF_PARAMS, only: params_set, intercept, slope, switch
  use c4_photosynth, only: c4_photosynthesis
  use enbal_mod,     only: energy_balance
  use stomata,       only: call_stomata
  use boundary,      only: call_boundary
  implicit none
  private
  public :: VLEAF_STEP

  integer, parameter :: max_iter = 50
  real(8), parameter :: tol_An = 0.1d0, tol_Ci = 0.1d0, tol_T = 0.1d0

  ! energy-balance constants (keep same as your standalone)
  real(8), parameter :: Cp_air   = 29.3d0
  real(8), parameter :: sigma_sb = 5.670374419d-8
  real(8), parameter :: Hfac     = 1.0d0
  real(8), parameter :: LEfac    = 1.0d0
  real(8), parameter :: LWfac    = 2.0d0
  real(8), parameter :: epsLeaf  = 0.94d0

  ! boundary-layer constants (keep same as your standalone)
  real(8), parameter :: leafDim  = 0.08d0
  real(8), parameter :: cForced  = 4.322d-3
  real(8), parameter :: cFree    = 1.6361d-3
  real(8), parameter :: s_ratio  = 0.71d0
  real(8), parameter :: Mb       = 0.5d0 * (1.0d0 + s_ratio)**2 / (1.0d0 + s_ratio**2)

contains

  subroutine VLEAF_STEP(Rg, ea, tAir, wind, pressure, ca, O2, &
                        Anet, Tr, tLeaf_out, LE, H, istat)
    implicit none
    real(8), intent(in)  :: Rg, ea, tAir, wind, pressure, ca, O2
    real(8), intent(out) :: Anet, Tr, tLeaf_out, LE, H
    integer, intent(out) :: istat

    integer :: iter
    real(8) :: PAR, NIR, LW
    real(8) :: Rn, Hf, LEf, LWe
    real(8) :: tLeaf, ci, gs, gb, eb, cb
    real(8) :: An_prev, Ci_prev, T_prev
    real(8) :: errA, errCi, errT

    ! LW estimation helpers (if you are not passing LW in)
    real(8) :: TaK, ea_kPa, eps_sky

    if (.not. params_set) then
      istat     = 1
      Anet      = 0.0d0
      Tr        = 0.0d0
      tLeaf_out = 0.0d0
      LE        = 0.0d0
      H         = 0.0d0
      return
    end if

    ! radiation split (same as your standalone)
    PAR = 0.8d0  * 0.45d0 * Rg
    NIR = 0.23d0 * 0.55d0 * Rg

    ! LW from Ta and ea (Brutsaert-style)
    TaK = tAir
    if (TaK < 150.0d0) TaK = TaK + 273.15d0
    ea_kPa = ea
    if (ea_kPa > 20.0d0) ea_kPa = ea_kPa / 1000.0d0
    eps_sky = 1.72d0 * (ea_kPa / TaK)**(1.0d0/7.0d0)
    if (eps_sky > 1.0d0) eps_sky = 1.0d0
    LW = eps_sky * sigma_sb * TaK**4

    ! initial guesses (internal each call)
    tLeaf = tAir
    ci    = 0.7d0 * ca
    Anet  = 0.7d0 * ca
    gs    = 0.20d0
    gb    = 1.50d0
    eb    = 2.01d3
    cb    = 0.0d0

    Tr  = 0.0d0
    Hf  = 0.0d0
    LEf = 0.0d0

    do iter = 1, max_iter
      An_prev = Anet
      Ci_prev = ci
      T_prev  = tLeaf

      call energy_balance( &
           PAR, NIR, LW, tAir, ea, pressure, &
           Anet, gb, gs, &
           Cp_air, sigma_sb, Hfac, LEfac, LWfac, epsLeaf, switch, &
           tLeaf, Rn, Hf, LEf, LWe, Tr)

      call call_stomata(Anet, ca, tLeaf, eb, intercept, slope, ci, gs)

      call c4_photosynthesis(tLeaf, ci, PAR, Anet)

      call call_boundary(leafDim, cForced, cFree, Mb, &
           tAir, tLeaf, pressure, ea, ca, wind, gs, Anet, &
           eb, cb, gb)

      errA  = abs(Anet  - An_prev)
      errCi = abs(ci    - Ci_prev)
      errT  = abs(tLeaf - T_prev)

      if (errA < tol_An .and. errCi < tol_Ci .and. errT < tol_T) exit
    end do

    tLeaf_out = tLeaf
    LE        = LEf
    H         = Hf
    istat     = 0
  end subroutine VLEAF_STEP

end module VLEAF_DRIVER
