!========== enbal_mod.f90 ==========
module enbal_mod
  implicit none
  private
  public :: energy_balance_canopy

contains

  !===========================================================
  ! energy_balance_canopy
  !
  ! Matches MATLAB callEnergyBalance_Canopy:
  !   Rn = PAR + NIR + LWnet_canopy(Tleaf)
  !   Resid^2 minimized over Tleaf in [2,60]
  !
  ! Inputs:
  !   PAR, NIR      [W m-2] (already per leaf area if you feed leaf-absorbed)
  !   LWd           [W m-2] downwelling LW (already epsA*sigma*Ta^4 style)
  !   tAir          [degC]
  !   ea            [Pa]
  !   pAir          [Pa]
  !   Anet          [umol m-2 s-1] (same as you used)
  !   gb, gs        [mol m-2 s-1]
  !   Cp            [J mol-1 K-1]
  !   sigma         [W m-2 K-4]
  !   Hfac, LEfac, LWfac, epsLeaf
  !   epsSoil, LAI, Omega
  !   switch        if /=1 -> force Tleaf=tAir but still compute fluxes
  !
  ! Outputs:
  !   Tleaf         [degC]
  !   Rn            [W m-2]
  !   H, LE         [W m-2]
  !   LWnet         [W m-2] (net canopy LW term per leaf area)
  !   LWemit        [W m-2] (diagnostic emission term)
  !   Tr            [umol m-2 s-1]
  !   Resid         sqrt(residual^2) [W m-2]
  !===========================================================
  subroutine energy_balance_canopy( &
       PAR, NIR, LWd, tAir, tSoil, ea, pAir, Anet, gb, gs, &
       Cp, sigma, Hfac, LEfac, LWfac, epsLeaf, &
       epsSoil, LAI, Omega, switch, &
       Tleaf, Rn, H, LE, LWnet, LWemit, Tr, Resid )

    implicit none

    ! -------- inputs --------
    real, intent(in) :: PAR, NIR, LWd
    real, intent(in) :: tAir, tSoil, ea, pAir
    real, intent(in) :: Anet
    real, intent(in) :: gb, gs
    real, intent(in) :: Cp, sigma
    real, intent(in) :: Hfac, LEfac, LWfac, epsLeaf
    real, intent(in) :: epsSoil, LAI, Omega
    integer, intent(in) :: switch

    ! -------- outputs --------
    real, intent(out) :: Tleaf, Rn, H, LE, LWnet
    real, intent(out) :: LWemit, Tr, Resid

    ! -------- locals --------
    real, parameter :: Tlo = 2.0d0, Thi = 60.0d0
    real, parameter :: tiny = 1.0d-12
    real, parameter :: aa = 273.15d0
    real, parameter :: phi = 0.5d0*(sqrt(5.0d0)-1.0d0)  ! golden ratio conjugate
    integer, parameter :: itmax = 80

    real :: gtot
    real :: cosThetaBar, E
    real :: a, b, c, d, fc, fd
    real :: t1, t2
    integer  :: it

    ! ---- total vapor conductance (LeafState.g in MATLAB) ----
    if (gb + gs > tiny) then
      gtot = (gb*gs) / (gb + gs)
    else
      gtot = 0.0d0
    end if

    ! ---- canopy LW geometric term E ----
    cosThetaBar = 0.537d0 + 0.025d0*LAI
    cosThetaBar = min(1.0d0, max(1.0d-3, cosThetaBar))
    E = exp(-0.5d0*Omega*LAI / cosThetaBar)

    !=========================================================
    ! switch: skip optimizer and force Tleaf=tAir,
    ! but compute fluxes consistently using canopy LWnet
    !=========================================================
    if (switch /= 1) then
      Tleaf = tAir
      ! LWnet  = lw_net_canopy(Tleaf, LWd, tSoil, epsLeaf, epsSoil, sigma, LAI, E, LWfac)
      call lw_net_canopy(Tleaf, LWd, tSoil, epsLeaf, epsSoil, sigma, LAI, E, LWfac, LWnet, t1, t2)
      LWemit = LWfac * epsLeaf * sigma * (aa + Tleaf)**4

      H  = sensible_heat(Tleaf, tAir, Cp, gb, Hfac)
      LE = latent_heat  (Tleaf, ea, pAir, gtot, LEfac)
      Rn = PAR + NIR + LWnet

      Tr = 0.0d0
      if (lv(Tleaf) > tiny) Tr = LE / lv(Tleaf) * 1.0d6

      Resid = abs(Rn - 0.506d0*Anet - H - LE)
      return
    end if

    !=========================================================
    ! Minimize residual^2 over [Tlo, Thi] (golden-section)
    !=========================================================
    a = Tlo
    b = Thi
    c = b - phi*(b-a)
    d = a + phi*(b-a)

    fc = residual_sq(c)
    fd = residual_sq(d)

    do it = 1, itmax
      if (fc > fd) then
        a = c
        c = d
        fc = fd
        d = a + phi*(b-a)
        fd = residual_sq(d)
      else
        b = d
        d = c
        fd = fc
        c = b - phi*(b-a)
        fc = residual_sq(c)
      end if
      if (abs(b-a) < 1.0d-6) exit
    end do

    if (fc < fd) then
      Tleaf = c
    else
      Tleaf = d
    end if

    if ( (Tleaf <= Tlo + 0.001) .or. (Tleaf >= Thi - 0.001) ) then
      Tleaf = tAir
    end if

    !=========================================================
    ! Compute final fluxes at optimal Tleaf
    !=========================================================
    !LWnet = lw_net_canopy(Tleaf, LWd, tSoil, epsLeaf, epsSoil, sigma, LAI, E, LWfac)
    call lw_net_canopy(Tleaf, LWd, tSoil, epsLeaf, epsSoil, sigma, LAI, E, LWfac, LWnet, t1, t2)
    LWemit = LWfac * epsLeaf * sigma * (aa + Tleaf)**4

    H  = sensible_heat(Tleaf, tAir, Cp, gb, Hfac)
    LE = latent_heat  (Tleaf, ea, pAir, gtot, LEfac)
    Rn = PAR + NIR + LWnet
!    write(*,'(A,1X,F7.2,4(1X,ES12.4))') 'LWin LWnet T1 T2:', LWd, LWnet, t1, t2

    Tr = 0.0d0
    if (lv(Tleaf) > tiny) Tr = LE / lv(Tleaf) * 1.0d6

    Resid = sqrt( max(residual_sq(Tleaf), 0.0d0) )

  contains

    ! es(T) in Pa (same as MATLAB): 0.611*exp(17.502*T/(240.97+T))*1000
    real function es_pa(T)
      implicit none
      real, intent(in) :: T
      es_pa = (0.611d0 * exp(17.502d0*T/(240.97d0 + T))) * 1000.0d0
    end function es_pa

    ! Lv(T) in J mol-1 (same as MATLAB)
    real function lv(T)
      implicit none
      real, intent(in) :: T
      lv = (2500.0d0 - 2.36d0*T) * 18.0d0
    end function lv

    ! Sensible heat [W m-2] (MATLAB uses 0.924*gb)
    real function sensible_heat(T, Ta, Cp_in, gb_in, Hfac_in)
      implicit none
      real, intent(in) :: T, Ta, Cp_in, gb_in, Hfac_in
      sensible_heat = Hfac_in * Cp_in * (0.924d0*gb_in) * (T - Ta)
    end function sensible_heat

    ! Latent heat [W m-2]
    real function latent_heat(T, ea_in, p_in, g_in, LEfac_in)
      implicit none
      real, intent(in) :: T, ea_in, p_in, g_in, LEfac_in
      latent_heat = 0.0d0
      if (p_in > tiny) then
        latent_heat = LEfac_in * lv(T) * g_in / p_in * (es_pa(T) - ea_in)
      end if
    end function latent_heat

    ! Residual^2 for optimizer
    real function residual_sq(T)
      implicit none
      real, intent(in) :: T
      real :: lw_local, h_local, le_local, rn_local, r, t1dum, t2dum
      call lw_net_canopy(T, LWd, tSoil, epsLeaf, epsSoil, sigma, LAI, E, LWfac, lw_local, t1dum, t2dum)
      !lw_local = lw_net_canopy(T, LWd, tSoil, epsLeaf, epsSoil, sigma, LAI, E, LWfac)
      h_local  = sensible_heat(T, tAir, Cp, gb, Hfac)
      le_local = latent_heat  (T, ea, pAir, gtot, LEfac)
      rn_local = PAR + NIR + lw_local
      r = rn_local - 0.506d0*Anet - h_local - le_local
      residual_sq = r*r
    end function residual_sq

    ! Canopy LW net per leaf area (matches your MATLAB LWnet_canopy)
!    real function lw_net_canopy(tLeaf_in, LWd_in, tSoil_in, epsC, epsG, sigma_in, LAI_in, E_in, LWfac_in)
!      implicit none
!      real, intent(in) :: tLeaf_in, LWd_in, tSoil_in
!      real, intent(in) :: epsC, epsG, sigma_in, LAI_in, E_in, LWfac_in
!      real :: TleafK, TgK, L4, Lg, term1, term2, LAI_safe!

!      TleafK = aa + tLeaf_in
!      TgK    = aa + tSoil_in

!      L4 = sigma_in * TleafK**4
!      Lg = epsG * sigma_in * TgK**4

!      term1 = ( epsC*(LWd_in + Lg) - LWfac_in*epsC*L4 ) * (1.0d0 - E_in)
!      term2 = epsC*(1.0d0 - epsG) * ( LWd_in*E_in + epsC*L4*(1.0d0 - E_in) )

!      LAI_safe = max(LAI_in, 1.0d-12)
!      lw_net_canopy = (term1 + term2) / LAI_safe
!    end function lw_net_canopy

    ! Canopy LW net per leaf area + diagnostic terms
    subroutine lw_net_canopy(tLeaf_in, LWd_in, tSoil_in, epsC, epsG, sigma_in, LAI_in, E_in, LWfac_in, &
                             LWnet_out, t1_out, t2_out)
      implicit none
      real, intent(in)  :: tLeaf_in, LWd_in, tSoil_in
      real, intent(in)  :: epsC, epsG, sigma_in, LAI_in, E_in, LWfac_in
      real, intent(out) :: LWnet_out, t1_out, t2_out
      real :: TleafK, TgK, L4, Lg, term1, term2, LAI_safe

      TleafK = aa + tLeaf_in
      TgK    = aa + tSoil_in

      L4 = sigma_in * TleafK**4
      Lg = epsG * sigma_in * TgK**4

      term1 = ( epsC*(LWd_in + Lg) - LWfac_in*epsC*L4 ) * (1.0*0 + (1.0 - E_in))
      term2 = epsC*(1.0 - epsG) * ( LWd_in*E_in + epsC*L4*(1.0 - E_in) )

      LAI_safe = max(LAI_in, 1.0e-12)

      LWnet_out = (term1 + term2) / LAI_safe
      t1_out    = term1 / LAI_safe
      t2_out    = term2 / LAI_safe
    end subroutine lw_net_canopy

  end subroutine energy_balance_canopy

end module enbal_mod
