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

    real(r8), pointer :: plant_minn_active_yield_flx_col             (:)    !column level mineral nitrogen yield from soil bgc calculation
    real(r8), pointer :: plant_minn_passive_yield_flx_col            (:)
    real(r8), pointer :: plant_minn_active_yield_flx_patch           (:)    !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minn_passive_yield_flx_patch          (:)    !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minn_active_yield_flx_vr_col          (:, :) !layer specific active mineral nitrogen yield
    real(r8), pointer :: plant_minn_uptake_potential_patch           (:)
    real(r8), pointer :: plant_minn_uptake_potential_col             (:)
    real(r8), pointer :: plant_minn_uptake_potential_vr_patch        (:,:)  !plant mineral nitrogen uptake potential for each layer
    real(r8), pointer :: plant_minn_uptake_potential_vr_col          (:,:)  !plant mineral nitrogen uptake potential for each layer
    real(r8), pointer :: plant_totn_demand_flx_col                   (:)    !column level total nitrogen demand, g N/m2/s
    real(r8), pointer :: fppnd_col                                   (:)    !fraction of fufilled nitrogen demand
    real(r8), pointer :: plant_frootsc_vr_patch                      (:,:)  !fine root for nutrient uptake
    real(r8), pointer :: plant_frootsc_patch                         (:)    !fine root for nutrient uptake
    real(r8), pointer :: rr_col                                      (:)    ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: annavg_agnpp_patch                          (:)     ! (gC/m2/s) annual average aboveground NPP
    real(r8), pointer :: annavg_bgnpp_patch                          (:)     ! (gC/m2/s) annual average belowground NPP
    real(r8), pointer :: tempavg_agnpp_patch                         (:)     ! (gC/m2/s) temp. average aboveground NPP
    real(r8), pointer :: tempavg_bgnpp_patch                         (:)     ! (gC/m2/s) temp. average belowground NPP

    real(r8), pointer :: bgc_cpool_ext_loss_vr_col                 (:, :, :)  ! col-level extneral organic carbon loss gC/m3 /time step
    real(r8), pointer :: bgc_cpool_ext_inputs_vr_col               (:, :, :)  ! col-level extneral organic carbon input gC/m3 /time step
    real(r8), pointer :: bgc_ppool_ext_inputs_vr_col               (:,:,:) !col organic nitrogen input, gN/m3/time step
    real(r8), pointer :: bgc_ppool_ext_loss_vr_col                 (:,:,:) !col extneral organic nitrogen loss, gN/m3/time step

    real(r8), pointer :: bgc_ppool_inputs_col                      (:,:)   !col organic N input, gN/m2/time step
    real(r8), pointer :: bgc_npool_ext_inputs_vr_col               (:,:,:) !col organic nitrogen input, gN/m3/time step
    real(r8), pointer :: bgc_npool_ext_loss_vr_col                 (:,:,:) !col extneral organic nitrogen loss, gN/m3/time step
    real(r8), pointer :: sminn_no3_input_vr_col                    (:,:)   !col no3 input, gN/m3/time step
    real(r8), pointer :: sminn_nh4_input_vr_col                    (:,:)   !col nh4 input, gN/m3/time step
    real(r8), pointer :: biochem_pmin_vr_col                       (:,:)   ! col vertically-resolved total biochemical P mineralization (gP/m3/s)

    real(r8), pointer :: hr_vr_col                                 (:,:)   ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)

    real(r8), pointer :: f_nit_vr_col                              (:,:)   ! col (gN/m3/s) soil nitrification flux
    real(r8), pointer :: f_denit_vr_col                            (:,:)   ! col (gN/m3/s) soil denitrification flux
    real(r8), pointer :: f_n2o_denit_vr_col                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
    real(r8), pointer :: f_n2o_nit_vr_col                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]
    real(r8), pointer :: actual_immob_no3_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
    real(r8), pointer :: actual_immob_nh4_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)
    real(r8), pointer :: smin_no3_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
    real(r8), pointer :: smin_nh4_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)
    real(r8), pointer :: supplement_to_sminn_vr_col                (:,:)   ! col vertically-resolved supplemental N supply (gN/m3/s)

   contains

     procedure , public  :: Init
     procedure , public  :: SetValues
     procedure , public  :: summary
     procedure , public  :: calc_nutrient_uptake_potential
     procedure , private :: InitAllocate
     procedure , private :: InitHistory
     procedure , private :: InitCold
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
    allocate(this%plant_minn_uptake_potential_patch  (begp:endp          )) ; this%plant_minn_uptake_potential_patch  (:)   = nan
    allocate(this%plant_minn_uptake_potential_col    (begc:endc          )) ; this%plant_minn_uptake_potential_col    (:)   = nan
    allocate(this%fppnd_col                          (begc:endc          )) ; this%fppnd_col                          (:)   = nan


    allocate(this%plant_minn_active_yield_flx_col    (begc:endc          )) ; this%plant_minn_active_yield_flx_col    (:)   = nan
    allocate(this%plant_minn_passive_yield_flx_col   (begc:endc          )) ; this%plant_minn_passive_yield_flx_col   (:)   = nan

    allocate(this%plant_minn_active_yield_flx_vr_col (begc:endc, lbj:ubj )) ; this%plant_minn_active_yield_flx_vr_col (:,:) = nan

    allocate(this%plant_totn_demand_flx_col          (begc:endc          )) ; this%plant_totn_demand_flx_col          (:)   = nan

    allocate(this%plant_frootsc_patch                (begc:endp          )) ; this%plant_frootsc_patch                (:)   = nan
    allocate(this%plant_frootsc_vr_patch             (begc:endp,lbj:ubj  )) ; this%plant_frootsc_vr_patch            (:,:)  = nan

    allocate(this%annavg_agnpp_patch                (begp:endp))                  ; this%annavg_agnpp_patch  (:) = spval ! To detect first year
    allocate(this%annavg_bgnpp_patch                (begp:endp))                  ; this%annavg_bgnpp_patch  (:) = spval ! To detect first year
    allocate(this%tempavg_agnpp_patch               (begp:endp))                  ; this%tempavg_agnpp_patch (:) = spval
    allocate(this%tempavg_bgnpp_patch               (begp:endp))                  ; this%tempavg_bgnpp_patch (:) = spval

    allocate(this%bgc_cpool_ext_inputs_vr_col       (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_inputs_vr_col (:,:,:) = nan
    allocate(this%bgc_cpool_ext_loss_vr_col         (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_loss_vr_col   (:,:,:) = nan
    allocate(this%rr_col                            (begc:endc))                  ; this%rr_col                    (:)  =nan

    allocate(this%bgc_ppool_ext_inputs_vr_col (begc:endc,1:nlevdecomp_full,ndecomp_pools)) ;this%bgc_ppool_ext_inputs_vr_col    (:,:,:) = nan
    allocate(this%bgc_ppool_ext_loss_vr_col   (begc:endc,1:nlevdecomp_full,ndecomp_pools)) ;this%bgc_ppool_ext_loss_vr_col      (:,:,:) = nan

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
    allocate(this%supplement_to_sminn_vr_col (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminn_vr_col (:,:) = nan
    allocate(this%f_nit_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%f_nit_vr_col                     (:,:) = nan
    allocate(this%f_denit_vr_col              (begc:endc,1:nlevdecomp_full)) ; this%f_denit_vr_col                   (:,:) = nan

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

    this%fppnd_col(begc:endc) = spval
    call hist_addfld1d (fname='FPPND', units='none', &
         avgflag='A', long_name='fulfilled plant nitrogen demand from mineral nitrogen uptake', &
         ptr_col=this%fppnd_col)

    this%plant_minn_active_yield_flx_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_MINN_ACTIVE_YIELD_FLX_COL', units='gN/m^2/s', &
         avgflag='A', long_name='plant nitrogen active uptake flux from soil', &
         ptr_col=this%plant_minn_active_yield_flx_col)

    this%plant_minn_passive_yield_flx_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_MINN_PASSIVE_YIELD_FLX_COL', units='gN/m^2/s', &
         avgflag='A', long_name='plant nitrogen passive uptake flux from soil', &
         ptr_col=this%plant_minn_passive_yield_flx_col)


    this%plant_minn_active_yield_flx_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='PLANT_MINN_ACTIVE_YIELD_FLX_vr', units='gN/m^3/s', type2d='levtrc', &
            avgflag='A', long_name='plant nitrogen active_uptake flux from soil', &
            ptr_col=this%plant_minn_active_yield_flx_vr_col, default='inactive')


    this%plant_totn_demand_flx_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_TOTN_DEMAND_FLX', units='gN/m^2/s',  &
            avgflag='A', long_name='plant nitrogen demand flux', &
            ptr_col=this%plant_totn_demand_flx_col, default='inactive')

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
       this%plant_minn_active_yield_flx_col(i)   = value_column
       this%plant_minn_passive_yield_flx_col(i)  = value_column
       this%plant_totn_demand_flx_col(i)         = value_column
       this%fppnd_col(i)                         = value_column
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
  subroutine summary(this, bounds, ubj,  num_soilc, filter_soilc, dz, nh4_transp, no3_transp)
  !
  ! !DESCRIPTION:
  ! summarize state variables from different subpools/subfluxes
  ! !USES:
  use MathfuncMod              , only : dot_sum
  use clm_time_manager         , only : get_step_size
  use clm_varcon               , only : natomw

  ! !ARGUMENTS:
  class(plantsoilnutrientflux_type) :: this
  type(bounds_type), intent(in) :: bounds
  integer,  intent(in) :: num_soilc
  integer,  intent(in) :: filter_soilc(:)
  integer,  intent(in) :: ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,1:ubj)
  real(r8), intent(in) :: nh4_transp(bounds%begc:bounds%endc)
  real(r8), intent(in) :: no3_transp(bounds%begc:bounds%endc)

  ! !LOCAL VARIABLES:
  integer :: fc, c, j
  real(r8) :: dtime

  dtime =  get_step_size()

  do fc = 1, num_soilc
    c = filter_soilc(fc)
    this%plant_minn_active_yield_flx_col(c)  =dot_sum(this%plant_minn_active_yield_flx_vr_col(c,1:ubj),dz(c,1:ubj))/dtime
    this%plant_minn_passive_yield_flx_col(c) =(nh4_transp(c) + no3_transp(c))*natomw/dtime

    if (this%plant_minn_uptake_potential_col(c)>0._r8) then
       this%fppnd_col(c)   = (this%plant_minn_active_yield_flx_col(c) + &
                              this%plant_minn_passive_yield_flx_col(c))/this%plant_minn_uptake_potential_col(c)
    else
      this%fppnd_col(c) = 1._r8
    endif
  enddo

  end subroutine summary

!--------------------------------------------------------------------------------

  subroutine calc_nutrient_uptake_potential(this, bounds, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, frootc_patch)
  !
  ! !DESCRIPTION:
  ! diagnose the vmax for nutrient uptake, with the vision to use ECA or something alike.
  !
  ! !USES:
  use subgridAveMod       , only : p2c
  use GridcellType        , only : grc

  !
  ! !ARGUMENTS:
  class(plantsoilnutrientflux_type) :: this
  type(bounds_type) , intent(in)    :: bounds
  integer           , intent(in)    :: num_soilc
  integer           , intent(in)    :: filter_soilc(:)
  integer           , intent(in)    :: num_soilp
  integer           , intent(in)    :: filter_soilp(:)
  real(r8)          , intent(in)    :: frootc_patch(bounds%begp:bounds%endp)

  ! !LOCAL VARIABLES:
  real(r8) :: Vmax_minn = 1.e-6_r8  ! gN/gC/s
  integer  :: fp, p, fc, c
  real(r8) :: nscal = 1._r8

  SHR_ASSERT_ALL((ubound(frootc_patch) == (/bounds%endp/)), errMsg(__FILE__,__LINE__))

  !calculate root nitrogen uptake potential

  !default approach
  !
  do fc = 1, num_soilc
    c = filter_soilc(fc)
    this%plant_minn_uptake_potential_col(c) = this%plant_totn_demand_flx_col(c)*nscal
  enddo

  do fp = 1, num_soilp
    p = filter_soilp(fp)
    this%plant_frootsc_patch(p) = frootc_patch(p)
  enddo


end subroutine calc_nutrient_uptake_potential

end module PlantSoilnutrientFluxType
