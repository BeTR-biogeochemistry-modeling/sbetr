module BeTR_WaterstateType
  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  implicit none
  save
  private
  !----------------------------------------------------
  ! column water state variables structure
  !----------------------------------------------------
  type, public :: betr_Waterstate_type
    real(r8), pointer :: h2osoi_liq_col(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_col(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liq_old(:,:)       !liquid water (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_old(:,:)       !ice lens (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liqvol_col(:,:)    !volumetric liquid water content
    real(r8), pointer :: h2osoi_icevol_col(:,:)    !volumetric ice water content
    real(r8), pointer :: h2osoi_vol_col(:,:)       !volumetric water content, total
    real(r8), pointer :: air_vol_col(:,:)          !volume possessed by air
    real(r8), pointer :: finundated_col(:)         ! fraction of column that is inundated, this is for bgc caclulation in betr
    real(r8), pointer :: rho_vap(:,:)              !concentration of bulk water vapor, assume in equilibrium with the liquid water
    real(r8), pointer :: rhvap_soi(:,:)            !soil relative humidity
    real(r8), pointer :: smp_l_col     (:,:)       ! col liquid phase soil matric potential, mm
    real(r8), pointer :: frac_h2osfc_col (:)       ! col fractional area with surface water greater than zero
  end type betr_Waterstate_type


end module BeTR_WaterstateType
