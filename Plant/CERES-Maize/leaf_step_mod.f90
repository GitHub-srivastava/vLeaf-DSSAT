module leaf_step_mod
  use c4_photosynth
  use enbal_mod
  use stomata
  use boundary
  implicit none
  private
  public :: leaf_step

contains

  subroutine leaf_step( &
      ! --- timestep drivers ---
      PAR, NIR, LW_in, tAir, tSoil, &
      ea, pressure, ca, wind, &
      ! --- model constants / options ---
      switch, intercept, slope, SWFAC, &
      ! --- Leaf photosynthetic property
      vcmax25, &
      ! --- in/out state variables (iterated) ---
      An, Ci, gs, gb, eb, cb, tLeaf, &
      ! --- outputs (fluxes + diagnostics) ---
      Rn, Hf, LEf, LWe, Tr, Res, &
      LAI_eff, Omega, LWnet, &
      errA, errCi, errT, converged )

    implicit none

    ! ===== inputs =====
    real, intent(in) :: PAR, NIR, LW_in
    real, intent(in) :: tAir, ea, SWFAC
    real, intent(in) :: tSoil, pressure, ca, wind
    integer, intent(in) :: switch
    real, intent(in) :: vcmax25
    real, intent(in) :: intercept, slope
    real, intent(in) :: LAI_eff, Omega

    integer, parameter :: max_iter = 100
    real, parameter :: tol_An   = 0.1
    real, parameter :: tol_Ci   = 0.1
    real, parameter :: tol_T    = 0.1

    ! ===== state (iterated) =====
    real, intent(out) :: An, Ci, gs, gb, eb, cb, tLeaf

    ! ===== outputs =====
    real, intent(out) :: Rn, Hf, LEf, LWe, LWnet, Tr, Res
    real, intent(out) :: errA, errCi, errT
    logical, intent(out) :: converged

    ! physical / EB constants
    real, parameter :: Cp_air   = 29.3
    real, parameter :: sigma_sb = 5.670374e-8
    real, parameter :: Hfac     = 1.0
    real, parameter :: LEfac    = 1.0
    real, parameter :: LWfac    = 2.0
    real, parameter :: epsLeaf  = 0.94
    real, parameter :: epsSoil  = 0.90

    ! Boundary layer parameters (Nikolov-style)
    real, parameter :: leafDim  = 0.08
    real, parameter :: cForced  = 4.322e-3
    real, parameter :: cFree    = 1.6361e-3
    real, parameter :: s_ratio  = 0.71
    real, parameter :: fsv      = 1.0
    real, parameter :: Mb = 0.5 * (1.0 + s_ratio)**2 / (1.0 + s_ratio**2)

    ! ===== locals =====
    integer :: iter
    real :: An_prev, Ci_prev, T_prev, gs_prev
    real :: alpha

    alpha = 0.8
    converged = .false.
    errA  = 1.0
    errCi = 1.0
    errT  = 1.0

    An    = 0.1*ca
    Ci    = 0.7*ca
    gs    = 0.2
    gb    = 1.5
    eb    = 2.01e3
    cb    = 0.7*ca
    tLeaf = tAir

    do iter = 1, max_iter

      An_prev = An
      Ci_prev = Ci
      T_prev  = tLeaf
      gs_prev = gs

      ! 1) Energy balance: updates tLeaf and fluxes
      call energy_balance_canopy( &
       PAR, NIR, LW_in, tAir, tSoil, ea, pressure, An, gb, gs, &
       Cp_air, sigma_sb, Hfac, LEfac, LWfac, epsLeaf, &
       epsSoil, LAI_eff, Omega, switch, &
       tLeaf, Rn, Hf, LEf, LWnet, LWe, Tr, Res)

      tLeaf = tLeaf*alpha + (1-alpha)*T_prev

      ! 2) Stomata: updates Ci and gs
      call call_stomata(An, cb, tLeaf, eb, &
      intercept, slope, SWFAC, Ci, gs)
      gs = gs*alpha + (1-alpha)*gs_prev
      Ci = Ci*alpha + (1-alpha)*Ci_prev

      ! 3) Photosynthesis: updates An
      call c4_photosynthesis(tLeaf, vcmax25, Ci, PAR, An)
      An = An*alpha + (1-alpha)*An_prev

      ! 4) Boundary layer: updates eb, cb, gb
      call call_boundary(leafDim, cForced, cFree, Mb, &
           tAir, tLeaf, pressure, ea, ca, wind, gs, An, &
           eb, cb, gb)



      errA  = abs(An    - An_prev)
      errCi = abs(Ci    - Ci_prev)
      errT  = abs(tLeaf - T_prev)

      if (errA < tol_An .and. errCi < tol_Ci .and. errT < tol_T) then
        converged = .true.
        exit
      end if

    end do

  end subroutine leaf_step

end module leaf_step_mod
