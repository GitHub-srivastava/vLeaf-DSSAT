!===============================================================
! File: daily_carbon_mod.f90
! Purpose: Integrate hourly canopy CO2 uptake and transpiration to
!          daily totals
! Notes:
!   - No file writing (returns values via INTENT(OUT))
!   - Uses default REAL only (no REAL(8), no DOUBLE PRECISION)
!   - Each of the nrec hourly values is taken as the mean rate over
!     its own 86400/nrec s interval, so all 24 h of the day count
!     (the same convention DSSAT uses for its own hourly sums).
!===============================================================
module daily_carbon_mod
  use VLEAF_PARAMS, only: eta
  implicit none
  private
  public :: daily_carbon, daily_transpiration

contains

  !===============================================================
  ! daily_carbon
  !
  ! Assumptions:
  !   - asun/ash are net assimilation, umol CO2 m-2 leaf s-1
  !     (negative at night: leaf dark respiration)
  !   - laisun/laish are m2 leaf m-2 ground
  !
  ! Night-time leaf respiration is included in the daily net. The
  ! daily total is then reduced by eta (VLEAF_PARAMS) for whole-plant
  ! growth/maintenance respiration, which RUE-based CERES-Maize
  ! accounts for implicitly, and floored at zero because CERES has no
  ! reserve pool to draw a negative daily carbon balance from.
  !===============================================================
  subroutine daily_carbon(nrec, laisun, laish, asun, ash, &
                             molCO2_day, gC_day, gDM_day)
    implicit none

    integer, intent(in) :: nrec
    real,    intent(in) :: laisun(*), laish(*), asun(*), ash(*)

    real, intent(out) :: molCO2_day   ! mol CO2 m-2 ground day-1
    real, intent(out) :: gC_day       ! g C   m-2 ground day-1
    real, intent(out) :: gDM_day      ! g DM  m-2 ground day-1

    integer :: k
    real    :: umolCO2_day, dt_s
    real, parameter :: fC = 0.45      ! C fraction of dry matter

    molCO2_day = 0.0
    gC_day     = 0.0
    gDM_day    = 0.0

    if (nrec <= 0) return

    dt_s = 86400.0 / real(nrec)
    umolCO2_day = 0.0
    do k = 1, nrec
      umolCO2_day = umolCO2_day + &
                    (asun(k)*laisun(k) + ash(k)*laish(k)) * dt_s
    end do

    molCO2_day = max(0.0, umolCO2_day * 1.0e-6)
    gC_day     = molCO2_day * 12.01
    gDM_day    = (gC_day / fC) * eta

  end subroutine daily_carbon


  !===============================================================
  ! daily_transpiration
  !
  ! Integrate hourly canopy transpiration to daily totals.
  !
  ! Assumptions:
  !   - tsun/tsh are umol H2O m-2 leaf s-1
  !   - laisun/laish are m2 leaf m-2 ground
  !
  ! Outputs:
  !   - molH2O_day: mol H2O m-2 ground day-1
  !   - mm_day    : mm day-1  (DSSAT-friendly)
  !===============================================================
  subroutine daily_transpiration(nrec, laisun, laish, tsun, tsh, &
                                 molH2O_day, mm_day)
    implicit none

    integer, intent(in) :: nrec
    real,    intent(in) :: laisun(*), laish(*), tsun(*), tsh(*)

    real, intent(out) :: molH2O_day   ! mol H2O m-2 ground day-1
    real, intent(out) :: mm_day       ! mm day-1

    integer :: k
    real    :: umolH2O_day, dt_s
    real, parameter :: MW_H2O = 0.018   ! kg/mol; 1 kg/m2 = 1 mm

    molH2O_day = 0.0
    mm_day     = 0.0
    if (nrec <= 0) return

    dt_s = 86400.0 / real(nrec)
    umolH2O_day = 0.0
    do k = 1, nrec
      umolH2O_day = umolH2O_day + &
                    max(0.0, tsun(k)*laisun(k) + tsh(k)*laish(k)) * dt_s
    end do

    molH2O_day = umolH2O_day * 1.0e-6
    mm_day     = molH2O_day * MW_H2O

  end subroutine daily_transpiration

end module daily_carbon_mod
