!========== stomata.f90 ==========
module stomata
  implicit none
  private
  public :: call_stomata

contains

  !-------------------------------------------------------------------
  ! Ball-Berry stomatal conductance model (Ball, Woodrow & Berry, 1987)
  !
  !   gs = g0 + g1 * A * hs / Ca
  !
  ! where hs is relative humidity at the leaf surface (ea/ei), derived
  ! from leaf temperature and ambient vapour pressure.
  !
  ! Inputs
  !   Anet   : net assimilation A           [umol m-2 s-1]
  !   ca_air : ambient CO2 Ca              [umol mol-1 or ppm]
  !   tleaf  : leaf temperature             [deg C]
  !   ea     : vapour pressure of air       [Pa]
  !   g0     : residual conductance         [mol m-2 s-1]
  !   g1     : Ball-Berry slope parameter   [unitless]
  !
  ! Outputs
  !   ci     : intercellular CO2            [same units as Ca]
  !   gs     : stomatal conductance (H2O)   [mol m-2 s-1]
  !-------------------------------------------------------------------
  subroutine call_stomata(Anet, ca_air, tleaf, ea, &
    g0, g1, SWFAC, ci, gs)

    implicit none
    real, intent(in)  :: Anet      ! net assimilation
    real, intent(in)  :: ca_air    ! ambient CO2
    real, intent(in)  :: tleaf     ! leaf temperature (C)
    real, intent(in)  :: ea        ! vapour pressure air (Pa)
    real, intent(in)  :: g0, g1    ! Ball-Berry parameters
    real, intent(in)  :: SWFAC
    real, intent(out) :: ci, gs

    real :: ei_pa
    real :: gs_loc, hs

    ! saturation vapour pressure at leaf T (Pa)
    ei_pa = sat_vp(tleaf)

    hs = ea / max(ei_pa, 1.0)     ! RH at leaf surface
    hs = min(1.0, max(0.05, hs))

    ! Ball-Berry conductance (for water vapour)
    gs_loc = g0 + g1 * Anet / ca_air * hs
    ! keep gs positive
    gs = max(g0, g0 + SWFAC * (gs_loc - g0))

    ! simple mass-balance for ci (assume cb ≈ ca_air)
    ci = ca_air - 1.6d0 * Anet / gs
    ci = min(0.95*ca_air, max(0.05*ca_air, ci))

  contains

    ! saturation vapour pressure (Pa) at temperature T (C)
    real function sat_vp(T) result(es)
      real, intent(in) :: T
      es = 0.611d0 * exp(17.502d0*T/(240.97d0 + T)) * 1000.0d0
    end function sat_vp

  end subroutine call_stomata

end module stomata
