!===============================================================
! File: mz_daily_carbon_mod.f90
! Purpose: Integrate hourly canopy CO2 uptake to daily totals
! Notes:
!   - No file writing (returns values via INTENT(OUT))
!   - Uses default REAL only (no REAL(8), no DOUBLE PRECISION)
!   - Assumes ASUN/ASH are umol CO2 m-2 leaf s-1
!   - LAISUN/LAISH are m2 leaf m-2 ground
!===============================================================
module daily_carbon_mod
  implicit none
  private
  public :: daily_carbon, daily_transpiration

contains

  subroutine daily_carbon(nrec, hour, laisun, laish, asun, ash, &
                             molCO2_day, gC_day, gDM_day)
    implicit none

    integer, intent(in) :: nrec
    real,    intent(in) :: hour(*), laisun(*), laish(*), asun(*), ash(*)

    real, intent(out) :: molCO2_day   ! mol CO2 m-2 ground day-1
    real, intent(out) :: gC_day       ! g C   m-2 ground day-1
    real, intent(out) :: gDM_day      ! g DM  m-2 ground day-1

    integer :: k
    real    :: umolCO2_day
    real    :: a1, a2, dt_s
    real, parameter :: fC  = 0.45d0
    real, parameter :: eta = 1.00d0
    logical, parameter :: uptake_only = .true.

    molCO2_day = 0.0d0
    gC_day     = 0.0d0
    gDM_day    = 0.0d0

    if (nrec <= 0) return

    umolCO2_day = 0.0d0

    if (nrec == 1) then
      dt_s = 3600.0d0
      a1 = asun(1)*laisun(1) + ash(1)*laish(1)
      if (uptake_only) a1 = max(0.0d0, a1)
      umolCO2_day = a1 * dt_s
    else
      do k = 1, nrec-1
        dt_s = (hour(k+1) - hour(k)) * 3600.0d0
        if (dt_s < 0.0d0) dt_s = dt_s + 86400.0d0   ! handle wrap if any
        if (dt_s <= 0.0d0) cycle

        a1 = asun(k)  * laisun(k)   + ash(k)  * laish(k)
        a2 = asun(k+1)* laisun(k+1) + ash(k+1)* laish(k+1)

        if (uptake_only) then
          a1 = max(0.0d0, a1)
          a2 = max(0.0d0, a2)
        end if

        umolCO2_day = umolCO2_day + 0.5*(a1 + a2)*dt_s
      end do
    end if

    molCO2_day = umolCO2_day * 1.0e-6
    gC_day     = molCO2_day * 12.01d0
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
  subroutine daily_transpiration(nrec, hour, laisun, laish, tsun, tsh, &
                                 molH2O_day, mm_day)
    implicit none

    integer, intent(in) :: nrec
    real,    intent(in) :: hour(*), laisun(*), laish(*), tsun(*), tsh(*)

    real, intent(out) :: molH2O_day   ! mol H2O m-2 ground day-1
    real, intent(out) :: mm_day       ! mm day-1

    integer :: k
    real    :: umolH2O_day
    real    :: e1, e2, dt_s
    real, parameter :: MW_H2O = 0.018d0   ! kg/mol; 1 kg/m2 = 1 mm
    logical, parameter :: nonneg_only = .true.

    molH2O_day = 0.0d0
    mm_day     = 0.0d0
    if (nrec <= 0) return

    umolH2O_day = 0.0d0

    if (nrec == 1) then
      dt_s = 3600.0d0
      e1 = tsun(1)*laisun(1) + tsh(1)*laish(1)
      if (nonneg_only) e1 = max(0.0d0, e1)
      umolH2O_day = e1 * dt_s
    else
      do k = 1, nrec-1
        dt_s = (hour(k+1) - hour(k)) * 3600.0d0
        if (dt_s < 0.0d0) dt_s = dt_s + 86400.0d0
        if (dt_s <= 0.0d0) cycle

        e1 = tsun(k)  * laisun(k)   + tsh(k)  * laish(k)
        e2 = tsun(k+1)* laisun(k+1) + tsh(k+1)* laish(k+1)

        if (nonneg_only) then
          e1 = max(0.0d0, e1)
          e2 = max(0.0d0, e2)
        end if

        umolH2O_day = umolH2O_day + 0.5*(e1 + e2)*dt_s
      end do
    end if

    molH2O_day = umolH2O_day * 1.0e-6
    mm_day     = molH2O_day * MW_H2O

  end subroutine daily_transpiration

end module daily_carbon_mod

