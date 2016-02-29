module WaterFluxIsoType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevsoi
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use PatchType      , only : patch
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterfluxiso_type
     real(r8), pointer :: qflx_prec_grnd_iso_patch     (:,:)   ! patch water onto ground including canopy runoff [kg/(m2 s)]
     real(r8), pointer :: qflx_prec_grnd_iso_col       (:,:)   ! col water onto ground including canopy runoff [kg/(m2 s)]
     real(r8), pointer :: qflx_rain_grnd_iso_patch     (:,:)   ! patch rain on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_rain_grnd_iso_col       (:,:)   ! col rain on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_snow_grnd_iso_patch     (:,:)   ! patch snow on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_snow_grnd_iso_col       (:,:)   ! col snow on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_sub_snow_iso_patch      (:,:)   ! patch sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_sub_snow_iso_col        (:,:)   ! col sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_iso_patch      (:,:)   ! patch ground surface dew formation (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_iso_col        (:,:)   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)

     real(r8), pointer :: qflx_evap_soi_iso_patch      (:,:)   ! patch soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_soi_iso_col        (:,:)   ! col soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_iso_patch      (:,:)   ! patch vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_iso_col        (:,:)   ! col vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_iso_patch      (:,:)   ! patch evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_iso_col        (:,:)   ! col evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_tot_iso_patch      (:,:)   ! patch pft_qflx_evap_soi + pft_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_evap_tot_iso_col        (:,:)   ! col col_qflx_evap_soi + col_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_evap_grnd_iso_patch     (:,:)   ! patch ground surface evaporation rate (mm H2O/s) [+]
     real(r8), pointer :: qflx_evap_grnd_iso_col       (:,:)   ! col ground surface evaporation rate (mm H2O/s) [+]

     real(r8), pointer :: qflx_ev_snow_iso_patch       (:,:)   ! patch evaporation heat flux from snow       (W/m**2) [+ to atm]
     real(r8), pointer :: qflx_ev_snow_iso_col         (:,:)   ! col evaporation heat flux from snow         (W/m**2) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_iso_patch       (:,:)   ! patch evaporation heat flux from soil       (W/m**2) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_iso_col         (:,:)   ! col evaporation heat flux from soil         (W/m**2) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_iso_patch     (:,:)   ! patch evaporation heat flux from soil       (W/m**2) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_iso_col       (:,:)   ! col evaporation heat flux from soil         (W/m**2) [+ to atm]

     !added for betr
     real(r8), pointer :: qflx_gross_infl_soil_iso_col (:)   ! col gross infiltration, before considering the evaporation
     real(r8), pointer :: qflx_gross_evap_soil_iso_col (:)   ! col gross infiltration from soil, this satisfies the relationship qflx_infl_col = qflx_gross_infl_soil_col-qflx_gross_evap_soil_col
     real(r8), pointer :: qflx_h2osfc2topsoi_iso_col   (:)   ! col liquid water coming from surface standing water top soil (mm H2O/s)
     real(r8), pointer :: qflx_snow2topsoi_iso_col     (:)   ! col liquid water coming from residual snow to topsoil (mm H2O/s)

  end type waterfluxiso_type


end module WaterFluxIsoType
