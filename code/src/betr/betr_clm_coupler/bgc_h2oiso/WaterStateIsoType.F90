module WaterStateIsoType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module variables for hydrology
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use clm_varctl     , only : use_vancouver, use_mexicocity, use_cn, iulog, use_luna
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno
  use clm_varcon     , only : spval
  use LandunitType   , only : lun
  use ColumnType     , only : col
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterstateiso_type
     real(r8), pointer :: h2osoi_iso_liq_col         (:,:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: h2osoi_iso_ice_col         (:,:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: snocan_iso_patch           (:,:)   ! patch canopy snow water (mm H2O)
     real(r8), pointer :: liqcan_iso_patch           (:,:)   ! patch canopy liquid water (mm H2O)
     real(r8), pointer :: h2ocan_iso_patch           (:,:)   ! patch canopy water (mm H2O)
     real(r8), pointer :: h2ocan_iso_col             (:,:)   ! col canopy water (mm H2O)
     real(r8), pointer :: h2osfc_iso_col             (:,:)   ! col surface water (mm H2O)

  end type waterstateiso_type


end module WaterStateIsoType
