module PlantSoilnutrientFluxType
#include "shr_assert.h"
  !!DESCRIPTION:
  ! data structure for above/below ground nutrient coupling.
  ! The vision is beyond nitrogen, which probably extends to P, S and ect.
  ! This is part of BeTRbgc
  ! Created by Jinyun Tang, Jan 11, 2015
  !
  ! !USES:
  !
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varctl             , only : iulog
  use clm_time_manager       , only : get_nstep
  use clm_varcon             , only : spval, ispval
  use decompMod              , only : bounds_type
  use ColumnType             , only : col
  use PatchType              , only : pft
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: plantsoilnutrientflux_type

    real(r8), pointer :: plant_minn_active_yield_flx_patch           (:)    !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_active_yield_flx_patch           (:)    !column level mineral phosphorus yeild from soil bgc calculation

    real(r8), pointer :: plant_minn_passive_yield_flx_patch          (:)    !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_passive_yield_flx_patch          (:)

    real(r8), pointer :: plant_minn_yield_flx_patch                  (:)
    real(r8), pointer :: plant_minp_yield_flx_patch                  (:)

    real(r8), pointer :: plant_minn_active_yield_flx_vr_patch        (:,:)  !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_active_yield_flx_vr_patch        (:,:)    !column level mineral phosphorus yeild from soil bgc calculation


    real(r8), pointer :: plant_minn_passive_yield_flx_col            (:)
    real(r8), pointer :: plant_minp_passive_yield_flx_col            (:)

    real(r8), pointer :: plant_minn_active_yield_flx_col             (:)    !column level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_active_yield_flx_col             (:)    !column level mineral phosphorus yeild from soil bgc calculation

    real(r8), pointer :: plant_minn_active_nh4_yield_flx_vr_col     (:,:)
    real(r8), pointer :: plant_minn_active_no3_yield_flx_vr_col     (:,:)
    real(r8), pointer :: plant_minn_active_yield_flx_vr_col          (:,:)  !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_active_yield_flx_vr_col        (:,:)    !column level mineral phosphorus yeild from soil bgc calculation

    real(r8), pointer :: decomp_effc_vr_patch                         (:,:)  !fine root for nutrient uptake

    real(r8), pointer :: plant_minn_nh4_uptake_vmax_vr_patch         (:,:)  !plant mineral nitrogen uptake potential for each layer
    real(r8), pointer :: plant_minn_no3_uptake_vmax_vr_patch         (:,:)  !plant mineral nitrogen uptake potential for each layer
    real(r8), pointer :: plant_effrootsc_vr_patch                    (:,:)  !fine root for nutrient uptake

    real(r8), pointer :: plant_frootsc_patch                         (:)    !fine root for nutrient uptake
    real(r8), pointer :: annavg_agnpp_patch                          (:)     ! (gC/m2/s) annual average aboveground NPP
    real(r8), pointer :: annavg_bgnpp_patch                          (:)     ! (gC/m2/s) annual average belowground NPP
    real(r8), pointer :: tempavg_agnpp_patch                         (:)     ! (gC/m2/s) temp. average aboveground NPP
    real(r8), pointer :: tempavg_bgnpp_patch                         (:)     ! (gC/m2/s) temp. average belowground NPP

    real(r8), pointer :: rr_col                                      (:)    ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: bgc_cpool_ext_loss_vr_col                 (:, :, :)  ! col-level extneral organic carbon loss gC/m3 /time step
    real(r8), pointer :: bgc_cpool_ext_inputs_vr_col               (:, :, :)  ! col-level extneral organic carbon input gC/m3 /time step
    real(r8), pointer :: bgc_ppool_ext_inputs_vr_col               (:,:,:) !col organic nitrogen input, gN/m3/time step
    real(r8), pointer :: bgc_ppool_ext_loss_vr_col                 (:,:,:) !col extneral organic nitrogen loss, gN/m3/time step

    real(r8), pointer :: bgc_ppool_inputs_col                      (:,:)   !col organic N input, gN/m2/time step
    real(r8), pointer :: bgc_npool_ext_inputs_vr_col               (:,:,:) !col organic nitrogen input, gN/m3/time step
    real(r8), pointer :: bgc_npool_ext_loss_vr_col                 (:,:,:) !col extneral organic nitrogen loss, gN/m3/time step
    real(r8), pointer :: sminn_no3_input_vr_col                    (:,:)   !col no3 input, gN/m3/time step
    real(r8), pointer :: sminn_nh4_input_vr_col                    (:,:)   !col nh4 input, gN/m3/time step
    real(r8), pointer :: sminp_input_vr_col                        (:,:)   !col minp input, gP/m3/time step
    real(r8), pointer :: biochem_pmin_vr_col                       (:,:)   ! col vertically-resolved total potential biochemical P mineralization (gP/m3/s)
    real(r8), pointer :: biochem_pmin_ppool_vr_col                 (:,:,:) ! p pool divided potential biochemical P mineralization
    real(r8), pointer :: hr_vr_col                                 (:,:)   ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)

    real(r8), pointer :: f_nit_vr_col                              (:,:)   ! col (gN/m3/s) soil nitrification flux
    real(r8), pointer :: f_denit_vr_col                            (:,:)   ! col (gN/m3/s) soil denitrification flux
    real(r8), pointer :: f_n2o_denit_vr_col                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
    real(r8), pointer :: f_n2o_nit_vr_col                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]
    real(r8), pointer :: actual_immob_no3_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
    real(r8), pointer :: actual_immob_nh4_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)
    real(r8), pointer :: smin_no3_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
    real(r8), pointer :: smin_nh4_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)

    real(r8), pointer :: km_minsurf_minnh4_vr_col                  (:,:)   !mineral NH4 adsorption affinity
    real(r8), pointer :: vmax_minsurf_minnh4_vr_col                    (:,:)
    real(r8), pointer :: r_adsorp_minp_cap_vr_col                  (:,:)   !mineral adsorption capacity, this interacts with secondary P
    real(r8), pointer :: r_adsorp_nh4_cap_vr_col                   (:,:)   !mineral adsorption capacity, this assumes equilibrium adsorption
    real(r8), pointer :: kd_desorp_minp_vr_col                     (:,:)   !desorption parameter, 1/s
    real(r8), pointer :: plant_minn_nh4_uptake_km_vr_patch         (:,:)   !
    real(r8), pointer :: plant_minn_no3_uptake_km_vr_patch         (:,:)   !
    real(r8), pointer :: plant_minp_uptake_km_vr_patch             (:,:)   !

!    real(r8), pointer :: decomp_minn_nh4_uptake_km_vr_patch         (:,:)   ! place holder
!    real(r8), pointer :: decomp_minn_no3_uptake_km_vr_patch         (:,:)   ! place holder
!    real(r8), pointer :: decomp_minp_uptake_km_vr_patch             (:,:)   ! place holder

    real(r8), pointer :: vmax_minsurf_minp_vr_col                    (:,:)
    real(r8), pointer :: km_minsurf_minp_vr_col                    (:,:)   !mineral P adsorption affinity
    real(r8), pointer :: plant_minp_uptake_vmax_vr_patch           (:,:)

    real(r8), pointer :: vmax_plant_nh4b_vr_col                      (:,:)
    real(r8), pointer :: vmax_plant_no3b_vr_col                      (:,:)
    real(r8), pointer :: vmax_plant_minpb_vr_col                     (:,:)

    real(r8), pointer :: plant_minn_nh4_uptake_km_vr_col          (:,:)
    real(r8), pointer :: plant_minn_no3_uptake_km_vr_col          (:,:)
    real(r8), pointer :: plant_minp_uptake_km_vr_col              (:,:)
    real(r8), pointer :: plant_compet_minn_vr_col                    (:,:)
    real(r8), pointer :: plant_compet_minp_vr_col                    (:,:)

    real(r8), pointer :: decomp_minn_nh4_uptake_km_vr_col          (:,:)
    real(r8), pointer :: decomp_minn_no3_uptake_km_vr_col          (:,:)
    real(r8), pointer :: decomp_minp_uptake_km_vr_col              (:,:)
    real(r8), pointer :: decomp_compet_minn_vr_col                    (:,:)
    real(r8), pointer :: decomp_compet_minp_vr_col                    (:,:)

    real(r8), pointer :: tranp_wt_patch                            (:,:)    !fraction of transpiration contributed by different pfts
    real(r8) :: decomp_km_nh4
    real(r8) :: decomp_km_no3
    real(r8) :: decomp_km_minp
   contains

     procedure , public  :: Init
     procedure , public  :: SetValues
     procedure , public  :: nutrient_flx_summary
     procedure , public  :: integrate_vr_flux_to_2D
     procedure , public  :: calc_nutrient_uptake_kinetic_pars
     procedure , public  :: init_plant_soil_feedback
     procedure , public  :: update_plant_nutrient_active_yield_patch
     procedure , public  :: do_om_phosphorus_bioextraction
     procedure , private :: InitAllocate
     procedure , private :: InitHistory
     procedure , private :: InitCold
     procedure , private :: sub_froot_prof
     procedure , private :: transp_col2patch
  end type plantsoilnutrientflux_type

 contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, lbj, ubj)
   !
   ! !DESCRIPTION:
   ! initialize data type
   !
   ! !ARGUMENTS:
    class(plantsoilnutrientflux_type) :: this
    type(bounds_type), intent(in) :: bounds

    integer          , intent(in) :: lbj, ubj

    call this%InitAllocate (bounds, lbj, ubj)
    call this%InitHistory (bounds)
    call this%InitCold (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, lbj, ubj)
    !
    ! !DESCRIPTION:
    ! Initialize pft nitrogen flux
    !
    ! USES
    use clm_varpar,  only : ndecomp_pools, nlevdecomp_full
    ! !ARGUMENTS:
    class (plantsoilnutrientflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: lbj, ubj
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc

    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%plant_minn_active_yield_flx_patch  (begp:endp          )) ; this%plant_minn_active_yield_flx_patch  (:)   = nan
    allocate(this%plant_minn_passive_yield_flx_patch (begp:endp          )) ; this%plant_minn_passive_yield_flx_patch (:)   = nan
    allocate(this%plant_minp_passive_yield_flx_patch (begp:endp          )) ; this%plant_minp_passive_yield_flx_patch (:)   = nan
    allocate(this%plant_minn_active_yield_flx_vr_patch (begp:endp,1:nlevdecomp_full )) ; this%plant_minn_active_yield_flx_vr_patch  (:,:)   = nan

    allocate(this%plant_minn_passive_yield_flx_col   (begc:endc          )) ; this%plant_minn_passive_yield_flx_col   (:)   = nan
    allocate(this%plant_minp_passive_yield_flx_col   (begc:endc          )) ; this%plant_minp_passive_yield_flx_col   (:)   = nan

    allocate(this%plant_minn_nh4_uptake_vmax_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_minn_nh4_uptake_vmax_vr_patch(:,:) = nan
    allocate(this%plant_minn_no3_uptake_vmax_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_minn_no3_uptake_vmax_vr_patch(:,:) = nan
    allocate(this%plant_minp_uptake_vmax_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_minp_uptake_vmax_vr_patch(:,:) = nan

    allocate(this%plant_minn_nh4_uptake_km_vr_patch(begp:endp,1:nlevdecomp_full)); this%plant_minn_nh4_uptake_km_vr_patch(:,:) = nan
    allocate(this%plant_minn_no3_uptake_km_vr_patch(begp:endp,1:nlevdecomp_full)); this%plant_minn_no3_uptake_km_vr_patch(:,:) = nan
    allocate(this%plant_minp_uptake_km_vr_patch(begp:endp,1:nlevdecomp_full)); this%plant_minp_uptake_km_vr_patch(:,:) = nan

!    allocate(this%decomp_minn_nh4_uptake_km_vr_patch(begp:endp,1:nlevdecomp_full)); this%decomp_minn_nh4_uptake_km_vr_patch(:,:) = nan
!    allocate(this%decomp_minn_no3_uptake_km_vr_patch(begp:endp,1:nlevdecomp_full)); this%decomp_minn_no3_uptake_km_vr_patch(:,:) = nan
!    allocate(this%decomp_minp_uptake_km_vr_patch(begp:endp,1:nlevdecomp_full)); this%decomp_minp_uptake_km_vr_patch(:,:) = nan

    allocate(this%vmax_minsurf_minnh4_vr_col (begc:endc, 1:nlevdecomp_full)); this%vmax_minsurf_minnh4_vr_col (:,:) = nan
    allocate(this%vmax_minsurf_minp_vr_col (begc:endc, 1:nlevdecomp_full));   this%vmax_minsurf_minp_vr_col (:,:) = nan

    allocate(this%plant_frootsc_patch                (begc:endp          )) ; this%plant_frootsc_patch                (:)   = nan
    allocate(this%plant_effrootsc_vr_patch             (begc:endp,lbj:ubj  )) ; this%plant_effrootsc_vr_patch            (:,:)  = nan

    allocate(this%annavg_agnpp_patch                (begp:endp))                  ; this%annavg_agnpp_patch  (:) = spval ! To detect first year
    allocate(this%annavg_bgnpp_patch                (begp:endp))                  ; this%annavg_bgnpp_patch  (:) = spval ! To detect first year
    allocate(this%tempavg_agnpp_patch               (begp:endp))                  ; this%tempavg_agnpp_patch (:) = spval
    allocate(this%tempavg_bgnpp_patch               (begp:endp))                  ; this%tempavg_bgnpp_patch (:) = spval

    allocate(this%bgc_cpool_ext_inputs_vr_col       (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_inputs_vr_col (:,:,:) = nan
    allocate(this%bgc_cpool_ext_loss_vr_col         (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_loss_vr_col   (:,:,:) = nan
    allocate(this%rr_col                            (begc:endc))                  ; this%rr_col                    (:)  =nan

    allocate(this%bgc_ppool_ext_inputs_vr_col (begc:endc,1:nlevdecomp_full,ndecomp_pools)) ;this%bgc_ppool_ext_inputs_vr_col    (:,:,:) = nan
    allocate(this%bgc_ppool_ext_loss_vr_col   (begc:endc,1:nlevdecomp_full,ndecomp_pools)) ;this%bgc_ppool_ext_loss_vr_col      (:,:,:) = nan
    allocate(this%biochem_pmin_ppool_vr_col   (begc:endc,1:nlevdecomp_full,ndecomp_pools)) ;this%biochem_pmin_ppool_vr_col      (:,:,:) = nan

    allocate(this%bgc_ppool_inputs_col        (begc:endc,ndecomp_pools))     ;this%bgc_ppool_inputs_col              (:,:) = nan

    allocate(this%bgc_npool_ext_inputs_vr_col (begc:endc,1:nlevdecomp_full,ndecomp_pools)) ;this%bgc_npool_ext_inputs_vr_col    (:,:,:) = nan
    allocate(this%bgc_npool_ext_loss_vr_col   (begc:endc,1:nlevdecomp_full,ndecomp_pools)) ;this%bgc_npool_ext_loss_vr_col      (:,:,:) = nan
    allocate(this%sminn_no3_input_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%sminn_no3_input_vr_col           (:,:) = nan
    allocate(this%sminn_nh4_input_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%sminn_nh4_input_vr_col           (:,:) = nan

    allocate(this%hr_vr_col                         (begc:endc,1:nlevdecomp_full)); this%hr_vr_col    (:,:)=nan

    allocate(this%f_n2o_denit_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_denit_vr_col               (:,:) = nan
    allocate(this%f_n2o_nit_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_nit_vr_col                 (:,:) = nan
    allocate(this%actual_immob_no3_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_no3_vr_col          (:,:) = nan
    allocate(this%actual_immob_nh4_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_nh4_vr_col          (:,:) = nan
    allocate(this%smin_no3_to_plant_vr_col    (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_to_plant_vr_col         (:,:) = nan
    allocate(this%smin_nh4_to_plant_vr_col    (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_to_plant_vr_col         (:,:) = nan
    allocate(this%f_nit_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%f_nit_vr_col                     (:,:) = nan
    allocate(this%f_denit_vr_col              (begc:endc,1:nlevdecomp_full)) ; this%f_denit_vr_col                   (:,:) = nan

    allocate(this%km_minsurf_minnh4_vr_col   (begc:endc,1:nlevdecomp_full)) ; this%km_minsurf_minnh4_vr_col(:,:) = nan
    allocate(this%km_minsurf_minp_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%km_minsurf_minp_vr_col  (:,:) = nan
    allocate(this%r_adsorp_minp_cap_vr_col       (begc:endc,1:nlevdecomp_full)) ; this%r_adsorp_minp_cap_vr_col    (:,:) = nan
    allocate(this%r_adsorp_nh4_cap_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%r_adsorp_nh4_cap_vr_col  (:,:) = nan
    allocate(this%kd_desorp_minp_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%kd_desorp_minp_vr_col      (:,:) = nan


    allocate(this%plant_compet_minp_vr_col(begc:endc, 1:nlevdecomp_full)); this%plant_compet_minp_vr_col (:,:) = nan
    allocate(this%plant_compet_minn_vr_col(begc:endc, 1:nlevdecomp_full)); this%plant_compet_minn_vr_col (:,:) = nan
    allocate(this%vmax_plant_nh4b_vr_col  (begc:endc, 1:nlevdecomp_full)); this%vmax_plant_nh4b_vr_col   (:,:) = nan
    allocate(this%vmax_plant_no3b_vr_col  (begc:endc, 1:nlevdecomp_full)); this%vmax_plant_no3b_vr_col   (:,:) = nan
    allocate(this%vmax_plant_minpb_vr_col  (begc:endc, 1:nlevdecomp_full)); this%vmax_plant_minpb_vr_col   (:,:) = nan


   allocate(this%plant_minn_nh4_uptake_km_vr_col(begc:endc,1:nlevdecomp_full));this%plant_minn_nh4_uptake_km_vr_col(:,:) = nan
   allocate(this%plant_minn_no3_uptake_km_vr_col(begc:endc,1:nlevdecomp_full));this%plant_minn_no3_uptake_km_vr_col(:,:) = nan
   allocate(this%plant_minp_uptake_km_vr_col    (begc:endc,1:nlevdecomp_full));this%plant_minp_uptake_km_vr_col    (:,:) = nan
   allocate(this%plant_minn_active_no3_yield_flx_vr_col(begc:endc,1:nlevdecomp_full));this%plant_minn_active_no3_yield_flx_vr_col(:,:) = nan
   allocate(this%plant_minn_active_nh4_yield_flx_vr_col(begc:endc,1:nlevdecomp_full));this%plant_minn_active_nh4_yield_flx_vr_col(:,:) = nan
   allocate(this%plant_minp_active_yield_flx_vr_col(begc:endc,1:nlevdecomp_full));this%plant_minp_active_yield_flx_vr_col(:,:) = nan

   allocate(this%decomp_compet_minp_vr_col(begc:endc, 1:nlevdecomp_full)); this%decomp_compet_minp_vr_col (:,:) = nan
   allocate(this%decomp_compet_minn_vr_col(begc:endc, 1:nlevdecomp_full)); this%decomp_compet_minn_vr_col (:,:) = nan

   allocate(this%decomp_minn_nh4_uptake_km_vr_col(begc:endc,1:nlevdecomp_full));this%decomp_minn_nh4_uptake_km_vr_col(:,:) = nan
   allocate(this%decomp_minn_no3_uptake_km_vr_col(begc:endc,1:nlevdecomp_full));this%decomp_minn_no3_uptake_km_vr_col(:,:) = nan
   allocate(this%decomp_minp_uptake_km_vr_col    (begc:endc,1:nlevdecomp_full));this%decomp_minp_uptake_km_vr_col    (:,:) = nan

   allocate(this%tranp_wt_patch(begc:endc,1:nlevdecomp_full)); this%tranp_wt_patch (:,:) = nan
   allocate(this%plant_minn_yield_flx_patch(begp:endp)) ; this%plant_minn_yield_flx_patch(:) = nan
   allocate(this%plant_minp_yield_flx_patch(begp:endp)) ; this%plant_minp_yield_flx_patch(:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd, crop_prog, nlevtrc_soil
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    !
    ! !ARGUMENTS:
    class(plantsoilnutrientflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer        :: k,l
    integer        :: begp, endp
    integer        :: begc, endc
    character(10)  :: active
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%plant_minn_active_yield_flx_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_MINN_ACTIVE_YIELD_FLX_PATCH', units='gN/m^2/s', &
         avgflag='A', long_name='plant nitrogen active uptake flux from soil', &
         ptr_patch=this%plant_minn_active_yield_flx_patch, default='inactive')

    this%plant_minn_passive_yield_flx_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_MINN_PASSIVE_YIELD_FLX_PATCH', units='gN/m^2/s', &
         avgflag='A', long_name='plant nitrogen passive uptake flux from soil', &
         ptr_patch=this%plant_minn_passive_yield_flx_patch, default='inactive')


    this%plant_minp_active_yield_flx_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_MINP_ACTIVE_YIELD_FLX_PATCH', units='gP/m^2/s', &
         avgflag='A', long_name='plant phosphorus active uptake flux from soil', &
         ptr_patch=this%plant_minp_active_yield_flx_patch, default='inactive')

    this%plant_minp_passive_yield_flx_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_MINP_PASSIVE_YIELD_FLX_PATCH', units='gP/m^2/s', &
         avgflag='A', long_name='plant phosphorus passive uptake flux from soil', &
         ptr_patch=this%plant_minp_passive_yield_flx_patch, default='inactive')


    this%plant_minn_passive_yield_flx_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_MINN_PASSIVE_YIELD_FLX_COL', units='gN/m^2/s', &
         avgflag='A', long_name='plant nitrogen passive uptake flux from soil', &
         ptr_col=this%plant_minn_passive_yield_flx_col)

    this%plant_minn_active_yield_flx_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_MINN_ACTIVE_YIELD_FLX_COL', units='gN/m^2/s', &
         avgflag='A', long_name='plant nitrogen active uptake flux from soil', &
         ptr_col=this%plant_minn_active_yield_flx_col)

    this%plant_minn_passive_yield_flx_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_MINP_PASSIVE_YIELD_FLX_COL', units='gP/m^2/s', &
         avgflag='A', long_name='plant phosphorus passive uptake flux from soil', &
         ptr_col=this%plant_minn_passive_yield_flx_col)

    this%plant_minp_active_yield_flx_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_MINP_ACTIVE_YIELD_FLX_COL', units='gP/m^2/s', &
         avgflag='A', long_name='plant phosphorus active uptake flux from soil', &
         ptr_col=this%plant_minn_active_yield_flx_col)


    this%f_nit_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='F_NIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='nitrification flux', &
         ptr_col=this%f_nit_vr_col)

    this%f_denit_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='F_DENIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='denitrification flux', &
         ptr_col=this%f_denit_vr_col)
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine SetValues ( this,               &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    class (plantsoilnutrientflux_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i=filter_patch(fi)
       this%plant_minn_active_yield_flx_patch(i)  = value_patch
       this%plant_minn_passive_yield_flx_patch(i) = value_patch
    enddo

    do fi = 1,num_column
       i = filter_column(fi)
       this%plant_minn_passive_yield_flx_col(i)  = value_column
    enddo

  end subroutine SetValues
  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use clm_varpar      , only : crop_prog
    use landunit_varcon , only : istsoil, istcrop
    use LandunitType   , only : lun
    !
    ! !ARGUMENTS:
    class(plantsoilnutrientflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    integer :: fp, fc                                    ! filter indices
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches


    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do


    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)
  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine nutrient_flx_summary(this, bounds, ubj,  num_soilc, filter_soilc, dz, &
     nh4_transp_vr, no3_transp_vr, minp_transp_vr)
  !
  ! !DESCRIPTION:
  ! summarize state variables from different subpools/subfluxes
  ! !USES:
  use MathfuncMod              , only : dot_sum
  use clm_time_manager         , only : get_step_size
  use clm_varcon               , only : natomw, patomw
  use clm_varpar               , only : maxpatch_pft


  ! !ARGUMENTS:
  class(plantsoilnutrientflux_type) :: this
  type(bounds_type), intent(in) :: bounds
  integer,  intent(in) :: num_soilc
  integer,  intent(in) :: filter_soilc(:)
  integer,  intent(in) :: ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,1:ubj)
  real(r8), intent(in) :: nh4_transp_vr(bounds%begc:bounds%endc, 1:ubj)
  real(r8), intent(in) :: no3_transp_vr(bounds%begc:bounds%endc, 1:ubj)
  real(r8), intent(in) :: minp_transp_vr(bounds%begc:bounds%endc,1:ubj)

  ! !LOCAL VARIABLES:
  integer :: fc, c, j, p, pi
  real(r8) :: dtime
  real(r8) :: col_minn_flx, col_minp_flx

  dtime =  get_step_size()
  !initialize to zero
  do fc = 1, num_soilc
    this%plant_minp_passive_yield_flx_col(c) = 0._r8
    this%plant_minn_passive_yield_flx_col(c) = 0._r8
    do pi = 1,maxpatch_pft
      if (pi <=  col%npfts(c)) then
        p = col%pfti(c) + pi - 1
        if (pft%active(p)) then
          this%plant_minn_passive_yield_flx_patch(p) = 0._r8
          this%plant_minp_passive_yield_flx_patch(p) = 0._r8
        endif
      endif
    enddo
  enddo
  do j = 1, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      col_minn_flx = (nh4_transp_vr(c,j) + no3_transp_vr(c,j))*natomw/dtime
      col_minp_flx = minp_transp_vr(c,j) * patomw/dtime

      this%plant_minn_passive_yield_flx_col(c) =this%plant_minn_passive_yield_flx_col(c) + col_minn_flx
      this%plant_minp_passive_yield_flx_col(c) =this%plant_minp_passive_yield_flx_col(c) + col_minp_flx

      !divide into patches based on transpiration flux of each pft
      do pi = 1,maxpatch_pft
        if (pi <=  col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p)) then
            this%plant_minn_passive_yield_flx_patch(p) =  this%plant_minn_passive_yield_flx_patch(p) + col_minn_flx*this%tranp_wt_patch(p,j)
            this%plant_minp_passive_yield_flx_patch(p) =  this%plant_minp_passive_yield_flx_patch(p) + col_minp_flx*this%tranp_wt_patch(p,j)
          endif
        endif
      enddo
    enddo
  enddo

  do fc = 1, num_soilc
    do pi = 1,maxpatch_pft
      if (pi <=  col%npfts(c)) then
        p = col%pfti(c) + pi - 1
        if (pft%active(p)) then
          this%plant_minn_yield_flx_patch(p) = this%plant_minn_passive_yield_flx_patch(p) + this%plant_minn_active_yield_flx_patch(p)
          this%plant_minp_yield_flx_patch(p) = this%plant_minp_passive_yield_flx_patch(p) + this%plant_minp_active_yield_flx_patch(p)
        endif
      endif
    enddo
  enddo
  end subroutine nutrient_flx_summary

!--------------------------------------------------------------------------------

  subroutine calc_nutrient_uptake_kinetic_pars(this, bounds, ubj, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, t_scalar, w_scalar, cnstate_vars, k_secp_to_solp, k_secp_to_occlp)
  !
  ! !DESCRIPTION:
  ! diagnose the vmax for nutrient uptake, with the vision to use ECA or something alike.
  !
  ! !USES:
  use subgridAveMod            , only : p2c
  use GridcellType             , only : grc
  use CNStateType              , only : cnstate_type
  use clm_varcon               , only : secspday
  use pftvarcon                , only : noveg
  use soilorder_varcon         , only : r_desorp, r_occlude
  use clm_time_manager         , only : get_days_per_year, get_step_size
  !
  ! !ARGUMENTS:
  class(plantsoilnutrientflux_type) :: this
  type(bounds_type) , intent(in)    :: bounds
  integer           , intent(in)    :: ubj
  integer           , intent(in)    :: num_soilc
  integer           , intent(in)    :: filter_soilc(:)
  integer           , intent(in)    :: num_soilp
  integer           , intent(in)    :: filter_soilp(:)
  real(r8)          , intent(in)    :: t_scalar(bounds%begc: , 1: )
  real(r8)          , intent(in)    :: w_scalar(bounds%begc: , 1: )
  real(r8)          , intent(inout) :: k_secp_to_solp(bounds%begc: , 1: )
  real(r8)          , intent(inout) :: k_secp_to_occlp(bounds%begc: , 1: )
  type(cnstate_type), intent(in)    :: cnstate_vars


  ! !LOCAL VARIABLES:
  integer  :: p, j, fc, c
  real(r8) :: tws
  real(r8) :: locrbc   !local root biomassc
  real(r8) :: dt, dtd, dayspyr
  real(r8) :: r_desorp_c, rr
  real(r8) :: r_occlude_c

  !calculate root nitrogen uptake potential
  SHR_ASSERT_ALL((ubound(t_scalar) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(w_scalar) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(k_secp_to_solp) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(k_secp_to_occlp) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))

  associate(                                                           &
    isoilorder                   => cnstate_vars%isoilorder          , &
    cn_scalar                    => cnstate_vars%cn_scalar           , &
    cp_scalar                    => cnstate_vars%cp_scalar             &
  )

  do j = 1, ubj
    do fc = 1, num_soilc
      c =filter_soilc(fc)
      tws = t_scalar(c,j) *w_scalar(c,j)
      do p = col%pfti(c), col%pftf(c)
        if (pft%active(p).and. (pft%itype(p) .ne. noveg)) then
          this%plant_minp_uptake_vmax_vr_patch(p, j)     = this%plant_minp_uptake_vmax_vr_patch(p, j) * cp_scalar(p) * tws
          this%plant_minn_nh4_uptake_vmax_vr_patch(p, j) = this%plant_minn_nh4_uptake_vmax_vr_patch(p, j) * cn_scalar(p) * tws
          this%plant_minn_no3_uptake_vmax_vr_patch(p, j) = this%plant_minn_no3_uptake_vmax_vr_patch(p, j)* cn_scalar(p) * tws

          locrbc=pft%wtcol(p) * this%plant_effrootsc_vr_patch(p,j)
          this%plant_minn_nh4_uptake_km_vr_col(c,j) = this%plant_minn_nh4_uptake_km_vr_col(c,j) + locrbc / this%plant_minn_nh4_uptake_km_vr_patch(p,j)
          this%plant_minn_no3_uptake_km_vr_col(c,j) = this%plant_minn_no3_uptake_km_vr_col(c,j) + locrbc / this%plant_minn_no3_uptake_km_vr_patch(p,j)
          this%plant_minp_uptake_km_vr_col(c,j)     = this%plant_minp_uptake_km_vr_col(c,j)  + locrbc / this%plant_minp_uptake_vmax_vr_patch(p,j)

          this%plant_compet_minn_vr_col(c,j) = this%plant_compet_minn_vr_col(c,j) + locrbc
          this%plant_compet_minp_vr_col(c,j) = this%plant_compet_minp_vr_col(c,j) + locrbc

          this%vmax_plant_nh4b_vr_col(c,j) = this%vmax_plant_nh4b_vr_col(c,j) + locrbc * &
            this%plant_minn_nh4_uptake_vmax_vr_patch(p,j) / this%plant_minn_nh4_uptake_km_vr_patch(p,j)
          this%vmax_plant_no3b_vr_col(c,j) = this%vmax_plant_no3b_vr_col(c,j) + locrbc * &
            this%plant_minn_no3_uptake_vmax_vr_patch(p,j) / this%plant_minn_no3_uptake_km_vr_patch(p,j)
          this%vmax_plant_minpb_vr_col(c,j)= this%vmax_plant_minpb_vr_col(c,j)+ locrbc * &
            this%plant_minp_uptake_vmax_vr_patch(p,j) / this%plant_minp_uptake_km_vr_patch(p,j)

        endif

      enddo
      this%vmax_plant_nh4b_vr_col(c,j) = this%vmax_plant_nh4b_vr_col(c,j) / this%plant_minn_nh4_uptake_km_vr_col(c,j)
      this%vmax_plant_no3b_vr_col(c,j) = this%vmax_plant_no3b_vr_col(c,j) / this%plant_minn_no3_uptake_km_vr_col(c,j)
      this%vmax_plant_minpb_vr_col(c,j)= this%vmax_plant_minpb_vr_col(c,j)/ this%plant_minp_uptake_km_vr_col(c,j)

      this%plant_minn_nh4_uptake_km_vr_col(c,j) = this%plant_minn_nh4_uptake_km_vr_col(c,j) / this%plant_compet_minn_vr_col(c,j)
      this%plant_minn_no3_uptake_km_vr_col(c,j) = this%plant_minn_no3_uptake_km_vr_col(c,j) / this%plant_compet_minn_vr_col(c,j)
      this%plant_minp_uptake_km_vr_col(c,j)     = this%plant_minp_uptake_km_vr_col(c,j) / this%plant_compet_minp_vr_col(c,j)

      this%decomp_minn_nh4_uptake_km_vr_col(c,j) = this%decomp_km_nh4
      this%decomp_minn_no3_uptake_km_vr_col(c,j) = this%decomp_km_no3
      this%decomp_minp_uptake_km_vr_col(c,j) = this%decomp_km_minp

    enddo
  enddo


  !desorption and occlusion parameters
  dayspyr = get_days_per_year()

  ! set time steps
  dt = real( get_step_size(), r8 )
  dtd = dt/(30._r8*secspday)


  do j = 1, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      !desorption rate from secondary p to dissolvable p
      r_desorp_c = r_desorp( isoilorder(c) )
      rr=-log(1._r8-r_desorp_c)
      r_desorp_c = 1._r8-exp(-rr*dtd)
      k_secp_to_solp(bounds%begc: , 1: ) = r_desorp_c/dt

      !occlusion rate parameter
      r_occlude_c = r_occlude( isoilorder(c) )
      rr=-log(1._r8-r_occlude_c)
      r_occlude_c = 1._r8-exp(-rr*dtd)
      k_secp_to_occlp(bounds%begc: , 1: ) = r_occlude_c/dt
    enddo
  enddo


  end associate
  end subroutine calc_nutrient_uptake_kinetic_pars

!--------------------------------------------------------------------------------
  subroutine sub_froot_prof(this, bounds, num_soilc, filter_soilc, frootc_patch, cnstate_vars, e_plant_scalar)

  use clm_varpar               , only : nlevtrc_soil
  use CNStateType              , only : cnstate_type
  use pftvarcon                , only : noveg
  !
  ! !ARGUMENTS:
  class(plantsoilnutrientflux_type) :: this
  integer           , intent(in)    :: num_soilc
  integer           , intent(in)    :: filter_soilc(:)
  type(bounds_type) , intent(in)    :: bounds
  real(r8)          , intent(in)    :: frootc_patch(bounds%begp: )
  real(r8)          , intent(in)    :: e_plant_scalar
  type(cnstate_type), intent(in)    :: cnstate_vars
  integer :: p, j, fc, c


  SHR_ASSERT_ALL((ubound(frootc_patch) == (/bounds%endp/)), errMsg(__FILE__,__LINE__))
  associate(                                                                        &
         froot_prof                   => cnstate_vars%froot_prof_patch              & !
  )

  do fc = 1, num_soilc
    c = filter_soilc(fc)
    do p = col%pfti(c), col%pftf(c)
      if (pft%active(p).and. (pft%itype(p) .ne. noveg)) then
        this%plant_frootsc_patch(p) = frootc_patch(p)
        !set effective nutrient uptake profile
        do j = 1, nlevtrc_soil
          this%plant_effrootsc_vr_patch(p, j) = this%plant_frootsc_patch(p) * froot_prof(p,j) * e_plant_scalar
        enddo
      endif
    enddo
  enddo
  end associate
  end subroutine sub_froot_prof

!--------------------------------------------------------------------------------
  subroutine init_plant_soil_feedback(this, bounds, num_soilc, filter_soilc, frootc_patch, cnstate_vars, &
   soilstate_vars,  waterflux_vars, ecophyscon_vars)

  use EcophysConType      , only : ecophyscon_type
  use pftvarcon           , only : noveg
  use clm_varpar          , only : nlevtrc_soil
  use CNStateType         , only : cnstate_type
  use WaterfluxType       , only : waterflux_type
  use SoilStateType       , only : soilstate_type

  class(plantsoilnutrientflux_type) :: this
  type(bounds_type)    , intent(in)    :: bounds
  integer              , intent(in)    :: num_soilc
  integer              , intent(in)    :: filter_soilc(:)
  real(r8)             , intent(in)    :: frootc_patch(bounds%begp: )
  type(ecophyscon_type), intent(in)    :: ecophyscon_vars
  type(cnstate_type)   , intent(in)    :: cnstate_vars
  type(waterflux_type) , intent(in)    :: waterflux_vars
  type(soilstate_type) , intent(in)    :: soilstate_vars

  real(r8), parameter   :: E_plant_scalar  = 0.0000125_r8
  integer :: fc, c, p, j

  SHR_ASSERT_ALL((ubound(frootc_patch) == (/bounds%endp/)), errMsg(__FILE__,__LINE__))
  associate(                                                                              &
  vmax_plant_nh4               => ecophyscon_vars%vmax_plant_nh4                        , &
  vmax_plant_no3               => ecophyscon_vars%vmax_plant_no3                        , &
  vmax_plant_p                 => ecophyscon_vars%vmax_plant_p                          , &
  vmax_minsurf_p_vr            => ecophyscon_vars%vmax_minsurf_p_vr                     , &
  ivt                          => pft%itype                                             , &
  decompmicc_patch_vr          => ecophyscon_vars%decompmicc_patch_vr                   , &
  km_decomp_nh4                => ecophyscon_vars%km_decomp_nh4                         , &
  km_decomp_no3                => ecophyscon_vars%km_decomp_no3                         , &
  km_decomp_p                  => ecophyscon_vars%km_decomp_p                             &

  )

  this%decomp_km_nh4 = km_decomp_nh4
  this%decomp_km_no3 = km_decomp_no3
  this%decomp_km_minp= km_decomp_p
  !set reference vmax
  do j = 1, nlevtrc_soil
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      this%decomp_compet_minp_vr_col(c,j) = 0._r8
      do p = col%pfti(c), col%pftf(c)
        if (pft%active(p).and. (pft%itype(p) .ne. noveg)) then
          this%plant_minp_uptake_vmax_vr_patch(p, j)     = vmax_plant_p(ivt(p))
          this%plant_minn_nh4_uptake_vmax_vr_patch(p, j) = vmax_plant_nh4(ivt(p))
          this%plant_minn_no3_uptake_vmax_vr_patch(p, j) = vmax_plant_no3(ivt(p))
          this%decomp_compet_minn_vr_col(c,j) = this%decomp_compet_minn_vr_col(c,j) + decompmicc_patch_vr(ivt(p),j)*pft%wtcol(p)

        endif
      enddo
      this%decomp_compet_minp_vr_col(c,j) = this%decomp_compet_minn_vr_col(c,j)
    enddo
  enddo

  call this%transp_col2patch(bounds, num_soilc, filter_soilc, soilstate_vars, waterflux_vars)
  !set root profile
  call this%sub_froot_prof(bounds, num_soilc, filter_soilc,  frootc_patch, cnstate_vars, e_plant_scalar)

  !set OM input profile

  !set mineral nutrient input profile

  end associate
  end subroutine init_plant_soil_feedback

  !-------------------------------------------------------------------------------
  subroutine transp_col2patch(this, bounds, num_soilc, filter_soilc, soilstate_vars, waterflux_vars)

  !
  ! USES
  use WaterfluxType         , only : waterflux_type
  use SoilStateType         , only : soilstate_type
  use clm_varpar            , only : nlevtrc_soil

  !
  !ARGUMENTS
  class(plantsoilnutrientflux_type)    :: this
  type(bounds_type)    , intent(in)    :: bounds
  integer              , intent(in)    :: num_soilc
  integer              , intent(in)    :: filter_soilc(:)
  type(waterflux_type) , intent(in)    :: waterflux_vars
  type(soilstate_type) , intent(in)    :: soilstate_vars

  real(r8) :: rootr_col
  integer  :: p, j, fc, c

  associate(                                                     &
    rootr_pft         =>    soilstate_vars%rootfr_patch         , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
    qflx_tran_veg_pft =>    waterflux_vars%qflx_tran_veg_patch   &
  )


  do j = 1,nlevtrc_soil
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      rootr_col = 0._r8
      do p = col%pfti(c), col%pftf(c)
        if (pft%active(p)) then
          rootr_col = rootr_col + rootr_pft(p,j) * qflx_tran_veg_pft(p) * pft%wtcol(p)
        end if
      enddo

      do p = col%pfti(c), col%pftf(c)
        if (pft%active(p)) then
          this%tranp_wt_patch(p,j) = rootr_pft(p,j) * qflx_tran_veg_pft(p) / rootr_col
        end if
      enddo
    end do
  end do
  end associate
  end subroutine transp_col2patch
  !-------------------------------------------------------------------------------

  subroutine update_plant_nutrient_active_yield_patch(this, bounds, num_soilc, filter_soilc)

    use clm_varpar               , only : maxpatch_pft
    use clm_varpar               , only : nlevtrc_soil
    !
    !ARGUMENTS
    class(plantsoilnutrientflux_type) :: this
    type(bounds_type)    , intent(in)    :: bounds
    integer              , intent(in)    :: num_soilc
    integer              , intent(in)    :: filter_soilc(:)

    integer :: p, pi, j, c, fc
    real(r8):: nh4_scal_denom, no3_scal_denom, minp_scal_denom
    real(r8):: fnh4b, fno3b, fminpb, nh4_scal, no3_scal, minp_scal



    do j = 1, nlevtrc_soil
      do fc = 1, num_soilc
        c =filter_soilc(fc)
        !calculate the overall flux scaling denominator
        nh4_scal_denom = this%plant_compet_minn_vr_col(c,j) * this%vmax_plant_nh4b_vr_col(c,j)/this%plant_minn_nh4_uptake_km_vr_col(c,j)
        no3_scal_denom = this%plant_compet_minn_vr_col(c,j) * this%vmax_plant_no3b_vr_col(c,j)/this%plant_minn_no3_uptake_km_vr_col(c,j)
        minp_scal_denom= this%plant_compet_minp_vr_col(c,j) * this%vmax_plant_minpb_vr_col(c,j)/this%plant_minp_uptake_km_vr_col(c,j)
        fnh4b = this%plant_minn_active_nh4_yield_flx_vr_col(c,j)
        fno3b = this%plant_minn_active_no3_yield_flx_vr_col(c,j)
        fminpb= this%plant_minp_active_yield_flx_vr_col(c,j)

        nh4_scal = fnh4b/nh4_scal_denom
        no3_scal = fno3b/no3_scal_denom
        minp_scal= fminpb/minp_scal_denom

        do pi = 1,maxpatch_pft
          if (pi <=  col%npfts(c)) then
            p = col%pfti(c) + pi - 1
            if (pft%active(p)) then
              this%plant_minn_active_yield_flx_vr_patch(p,j) = this%plant_effrootsc_vr_patch(p,j) * (nh4_scal * &
                this%plant_minn_nh4_uptake_vmax_vr_patch(p,j) / this%plant_minn_nh4_uptake_km_vr_patch(p,j) + &
                no3_scal * this%plant_minn_no3_uptake_vmax_vr_patch(p,j) / this%plant_minn_no3_uptake_km_vr_patch(p,j))
              this%plant_minp_active_yield_flx_vr_patch(p,j) = this%plant_effrootsc_vr_patch(p,j) * minp_scal * &
                this%plant_minp_uptake_vmax_vr_patch(p,j) / this%plant_minp_uptake_km_vr_patch(p,j)

            endif
          endif
        enddo
      enddo
    enddo

  end subroutine update_plant_nutrient_active_yield_patch
!-------------------------------------------------------------------------------
  subroutine do_om_phosphorus_bioextraction(this, bounds, ubj, num_soilc, filter_soilc, &
    phosphorusstate_vars, phosphorusflux_vars)

  !
  ! USES

  use clm_varpar               , only : ndecomp_pools
  use MathfuncMod              , only : safe_div
  use PhosphorusFluxType       , only : phosphorusflux_type
  use PhosphorusStateType      , only : phosphorusstate_type
  use clm_time_manager         , only : get_step_size

  !
  ! arguments
  class(plantsoilnutrientflux_type) :: this
  type(bounds_type)    , intent(in)    :: bounds
  integer              , intent(in)    :: num_soilc
  integer              , intent(in)    :: filter_soilc(:)
  integer              , intent(in)    :: ubj
  type(phosphorusstate_type)  , intent(in) :: phosphorusstate_vars
  type(phosphorusflux_type)   , intent(in) :: phosphorusflux_vars

  real(r8) :: tot_p, dtime
  integer  :: np, c, j, fc

  associate(                                                           &
    biochem_pmin_vr      => phosphorusflux_vars%biochem_pmin_vr_col  , &
    decomp_ppools_vr     => phosphorusstate_vars%decomp_ppools_vr_col  &
  )
  dtime =  get_step_size()
  do j = 1, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      this%biochem_pmin_vr_col(c,j)=biochem_pmin_vr(c,j)*dtime
      !divide biochem_pmin_vr flux to different organic p pools
      tot_p = sum(decomp_ppools_vr(c,j,1:ndecomp_pools))

      do np = 1,   ndecomp_pools
         this%biochem_pmin_ppool_vr_col(c,j,np) = &
            this%biochem_pmin_vr_col(c,j) * safe_div(decomp_ppools_vr(c,j, np),tot_p)
      enddo
    enddo
  enddo
  end associate
  end subroutine do_om_phosphorus_bioextraction


!-------------------------------------------------------------------------------
  subroutine integrate_vr_flux_to_2D(this, bounds, num_soilc, filter_soilc, carbonflux_vars, &
    nitrogenflux_vars, phosphorusflux_vars)

  use CNCarbonFluxType    , only : carbonflux_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use clm_varpar          , only : nlevdecomp
  use MathfuncMod         , only : dot_sum

  class(plantsoilnutrientflux_type) :: this
  type(bounds_type)        , intent(in)    :: bounds
  integer                  , intent(in)    :: num_soilc
  integer                  , intent(in)    :: filter_soilc(:)
  type(carbonflux_type)    , intent(inout) :: carbonflux_vars
  type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
  type(phosphorusflux_type), intent(inout) :: phosphorusflux_vars

  integer :: c, fc
!-------------------------------------------------------------------------------
  associate(                                                 &
     f_n2o_denit_col  => nitrogenflux_vars%f_n2o_denit_col , &
     f_n2o_nit_col    => nitrogenflux_vars%f_n2o_nit_col   , &
     f_nit_col        => nitrogenflux_vars%f_nit_col       , &
     f_denit_col      => nitrogenflux_vars%f_denit_col     , &
     dz               => col%dz                            , &
     hr_col           => carbonflux_vars%hr_col              &
  )


  do fc = 1, num_soilc
    c = filter_soilc(fc)
    hr_col(c)          = dot_sum(this%hr_vr_col(c,1:nlevdecomp), dz(c,1:nlevdecomp))
    f_n2o_denit_col(c) = dot_sum(this%f_n2o_denit_vr_col(c,1:nlevdecomp), dz(c,1:nlevdecomp))
    f_n2o_nit_col(c)   = dot_sum(this%f_n2o_nit_vr_col(c,1:nlevdecomp), dz(c,1:nlevdecomp))
    f_nit_col(c)       = dot_sum(this%f_nit_vr_col(c,1:nlevdecomp), dz(c,1:nlevdecomp))
    f_denit_col(c)     = dot_sum(this%f_denit_vr_col(c,1:nlevdecomp), dz(c,1:nlevdecomp))
  enddo

  end associate
  end subroutine integrate_vr_flux_to_2D
end module PlantSoilnutrientFluxType
