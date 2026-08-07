module boundary
  implicit none
  private
  public :: call_boundary

contains

  subroutine call_boundary(leafDim, cForced, cFree, Mb,   &
       tAir, tLeaf, pressure, ea, ca, wind, gs, Anet,                  &
       eb, cb, gb)

    implicit none
    ! Inputs
    real, intent(in)    :: leafDim      ! leaf width [m]
    real, intent(in)    :: cForced      ! forced convection coeff
    real, intent(in)    :: cFree        ! free convection coeff
    real, intent(in)    :: Mb           ! CO2/H2O diffusivity ratio
    real, intent(in)    :: tAir         ! air temp [C]
    real, intent(in)    :: tLeaf        ! leaf temp [C]
    real, intent(in)    :: pressure     ! air pressure [Pa]
    real, intent(in)    :: ea           ! air vapour pressure [Pa]
    real, intent(in)    :: ca           ! ambient CO2 [ppm]
    real, intent(in)    :: wind         ! wind speed [m s-1]
    real, intent(in)    :: gs           ! stomatal cond. [mol m-2 s-1]
    real, intent(in)    :: Anet         ! net CO2 [µmol m-2 s-1]

    ! In/out
    real, intent(inout) :: eb           ! leaf surface vapour press [Pa]

    ! Outputs
    real, intent(out)   :: cb           ! boundary layer CO2 [ppm]
    real, intent(out)   :: gb           ! boundary layer cond. [mol m-2 s-1]

    ! Locals
    real :: tAirK, tLeafK, ei
    real :: convert
    real :: gbForced, gbFree
    real :: TDiff
    real :: gs_mps, gb_mps

    ! temperatures in Kelvin
    tAirK  = tAir  + 273.15d0
    tLeafK = tLeaf + 273.15d0

    ! conversion m s-1 -> mol m-2 s-1  (Nikolov et al. 1995)
    convert = pressure / (8.309d0 * tAirK)

    ! saturation vapour pressure at leaf temperature [Pa]
    ei = 0.611d0 * exp(17.502d0 * tLeaf / (240.97d0 + tLeaf)) * 1000.0d0

    ! ---------- forced convection (m s-1) ----------
    gbForced = cForced * tAirK**0.56d0 *                                 &
               ((tAirK + 120.d0) * (wind / leafDim / pressure))**0.5d0

    ! if eb hasn't been set yet, start with ambient vapour pressure
    if (eb <= 0.0d0) eb = ea

    ! ---------- free convection (m s-1) ----------
    TDiff = ( tLeafK / (1.0d0 - 0.378d0 * eb / pressure ) ) -            &
            ( tAirK  / (1.0d0 - 0.378d0 * ea / pressure ) )

    gbFree = cFree * tLeafK**0.56d0 *                                    &
             ((tLeafK + 120.d0)/pressure)**0.5d0 *                       &
             (abs(TDiff)/leafDim)**0.25d0

    ! choose the larger of free / forced
    gb = max(gbFree, gbForced)          ! [m s-1]

    ! convert to mol m-2 s-1
    gb = gb * convert

    ! convert gs, gb back to m s-1 for eb calculation
    gs_mps = gs / convert
    gb_mps = gb / convert

    ! leaf surface vapour pressure [Pa]
    eb = (gs_mps*ei + gb_mps*ea) / (gs_mps + gb_mps)

    ! CO2 at leaf boundary layer [ppm]
    cb = ca - 1.37d0 * Anet / (Mb * gb)

  end subroutine call_boundary

end module boundary
