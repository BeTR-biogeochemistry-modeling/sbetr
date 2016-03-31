module BeTR_WaterfluxType

  !------------------------------------------------------------------------------
  ! !USES:
  use bshr_kind_mod    , only : r8 => shr_kind_r8
  implicit none
  save
  private

  !----------------------------------------------------
  ! column water flux variables structure
  !----------------------------------------------------
  type, public :: betr_waterflux_type
    real(r8), pointer :: qflx_adv_col             (:,:) !advection velocity from one layer to another, (0:nlevgrnd), positive downward
    real(r8), pointer :: qflx_infl_col            (:)     !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf_col            (:)     !surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_rootsoi             (:,:) ! water flux between root and soil [m/s]
    real(r8), pointer :: qflx_gross_evap_soil_col (:)   ! col gross infiltration from soil, this satisfies the relationship qflx_infl_col = qflx_gross_infl_soil_col-qflx_gross_evap_soil_col
    real(r8), pointer :: qflx_gross_infl_soil_col (:)   ! col gross infiltration, before considering the evaporation, mm/s
    real(r8), pointer :: qflx_rootsoi_col         (:,:) ! col root and soil water exchange [mm H2O/s] [+ into root]
    real(r8), pointer :: qflx_drain_vr_col        (:,:) ! col liquid water losted as drainage (m /time step)
    real(r8), pointer :: qflx_totdrain_col        (:)   ! col total liquid water drainage  (m/time step), updated in betr
    real(r8), pointer :: qflx_dew_grnd_col        (:)   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
    real(r8), pointer :: qflx_dew_snow_col        (:)   ! col surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_vol_col    (:)
    real(r8), pointer :: qflx_sub_snow_col        (:)   ! col sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_h2osfc2topsoi_col   (:)   ! col liquid water coming from surface standing water top soil (mm H2O/s)
    real(r8), pointer :: qflx_snow2topsoi_col     (:)   ! col liquid water coming from residual snow to topsoil (mm H2O/s)
    real(r8), pointer :: qflx_tran_veg_patch      (:)
  end type betr_waterflux_type




end module BeTR_WaterfluxType
