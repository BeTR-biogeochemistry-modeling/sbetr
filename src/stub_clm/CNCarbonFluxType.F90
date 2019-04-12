module CNCarbonFluxType
  use clm_varcon             , only : spval, ispval, c14ratio
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use clm_varcon             , only : spval
  use clm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: carbonflux_type
    real(r8), pointer :: rr_col                                    (:)   => null()  ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: rr_patch                                  (:)  => null()   ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: annsum_npp_patch                          (:) => null() ! patch annual sum of NPP (gC/m2/yr)
    real(r8), pointer :: agnpp_patch                               (:)  => null()   ! (gC/m2/s) aboveground NPP
    real(r8), pointer :: bgnpp_patch                               (:)  => null()   ! (gC/m2/s) belowground NPP
    real(r8), pointer :: hr_col                                    (:)=> null()
    real(r8), pointer :: hr_vr_col                                 (:,:) => null()
    real(r8), pointer :: phenology_c_to_litr_met_c_col             (:,:)=> null()
    real(r8), pointer :: phenology_c_to_litr_cel_c_col             (:,:)=> null()
    real(r8), pointer :: phenology_c_to_litr_lig_c_col             (:,:)=> null()
    real(r8), pointer :: dwt_livecrootc_to_cwdc_col                (:,:)=> null()
    real(r8), pointer :: m_decomp_cpools_to_fire_vr_col            (:,:,:)=> null()
    real(r8), pointer :: dwt_deadcrootc_to_cwdc_col                (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_lig_c_col              (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_cel_c_col              (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_met_c_col              (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_met_c_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_cel_c_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_lig_c_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_cwdc_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_met_c_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_cel_c_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_lig_c_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_cwdc_col                     (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_met_fire_col                  (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_cel_fire_col                  (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_lig_fire_col                  (:,:)=> null()
    real(r8), pointer :: fire_mortality_c_to_cwdc_col              (:,:)=> null()
    real(r8), pointer :: fire_decomp_closs_col                     (:) => null()
    real(r8), pointer :: som_c_leached_col                         (:) => null()
    real(r8), pointer :: som_c_runoff_col                          (:) => null()
    real(r8), pointer :: cflx_input_litr_met_vr_col              (:,:) => null()
    real(r8), pointer :: cflx_input_litr_cel_vr_col              (:,:) => null()
    real(r8), pointer :: cflx_input_litr_lig_vr_col              (:,:) => null()
    real(r8), pointer :: cflx_input_litr_cwd_vr_col              (:,:) => null()
    real(r8), pointer :: tempavg_agnpp_patch                     (:) => null()
    real(r8), pointer :: tempavg_bgnpp_patch                     (:) => null()
    real(r8), pointer :: annavg_agnpp_patch                      (:) => null()
    real(r8), pointer :: annavg_bgnpp_patch                      (:) => null()
    real(r8), pointer :: co2_soi_flx_col                         (:) => null()
    real(r8), pointer :: decomp_k_col                            (:,:,:) => null()
    real(r8), pointer :: t_scalar_col(:,:) => null()
    real(r8), pointer :: w_scalar_col(:,:) => null()
    real(r8), pointer :: rr_vr_col(:,:) => null()
    real(r8), pointer :: phr_vr_col(:,:)=> null()
    real(r8), pointer :: somhr_col(:) => null()
    real(r8), pointer :: lithr_col(:)=> null()
    real(r8), pointer :: o_scalar_col(:,:)=> null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type carbonflux_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate ( bounds )

    call this%InitCold ( bounds )

  end subroutine Init
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc


    allocate(this%rr_col                  (begc:endc)) ; this%rr_col      (:)  = nan
    allocate(this%rr_patch                (begp:endp)) ; this%rr_patch    (:)  = nan
    allocate(this%tempavg_agnpp_patch(begp:endp)); this%tempavg_agnpp_patch(:) = spval
    allocate(this%tempavg_bgnpp_patch(begp:endp)); this%tempavg_bgnpp_patch(:) = spval
    allocate(this%annavg_agnpp_patch(begp:endp));  this%annavg_agnpp_patch(:) = spval
    allocate(this%annavg_bgnpp_patch(begp:endp)); this%annavg_bgnpp_patch(:)  = spval

    allocate(this%annsum_npp_patch      (begp:endp)) ; this%annsum_npp_patch      (:) = spval
    allocate(this%agnpp_patch                       (begp:endp)) ; this%agnpp_patch                               (:) = spval
    allocate(this%bgnpp_patch                       (begp:endp)) ; this%bgnpp_patch                               (:) = spval
    allocate(this%hr_col (begc:endc)); this%hr_col(:) = spval
    allocate(this%hr_vr_col (begc:endc,1:nlevdecomp_full)); this%hr_vr_col(:,:) = spval              ! added  -zlyu

    allocate(this%phenology_c_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_met_c_col(:,:) = spval
    allocate(this%phenology_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_cel_c_col(:,:) = spval
    allocate(this%dwt_livecrootc_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%dwt_livecrootc_to_cwdc_col(:,:) = spval
    allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:7)); this%m_decomp_cpools_to_fire_vr_col(:,:,:)=nan
    allocate(this%dwt_deadcrootc_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootc_to_cwdc_col(:,:) = spval
    allocate(this%dwt_frootc_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_lig_c_col(:,:) = spval
    allocate(this%dwt_frootc_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_cel_c_col(:,:) = spval
    allocate(this%dwt_frootc_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_met_c_col(:,:) = spval
    allocate(this%gap_mortality_c_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_met_c_col(:,:) = spval
    allocate(this%gap_mortality_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_cel_c_col(:,:) = spval
    allocate(this%gap_mortality_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_lig_c_col(:,:) = spval
    allocate(this%gap_mortality_c_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_cwdc_col(:,:) = spval
    allocate(this%harvest_c_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_met_c_col(:,:) = spval
    allocate(this%harvest_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_cel_c_col(:,:) = spval
    allocate(this%harvest_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_lig_c_col(:,:) = spval
    allocate(this%harvest_c_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_cwdc_col(:,:) = spval
    allocate(this%m_c_to_litr_met_fire_col(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_met_fire_col(:,:) = spval
    allocate(this%m_c_to_litr_cel_fire_col(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_cel_fire_col(:,:) = spval
    allocate(this%m_c_to_litr_lig_fire_col(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_lig_fire_col(:,:) = spval
    allocate(this%fire_mortality_c_to_cwdc_col (begc:endc,1:nlevdecomp_full)); this%fire_mortality_c_to_cwdc_col (:,:) = spval
    allocate(this%phenology_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_lig_c_col(:,:) = spval
    allocate(this%som_c_leached_col(begc:endc)); this%som_c_leached_col(:) = spval
    allocate(this%rr_vr_col(begc:endc,1:nlevdecomp_full)); this%rr_vr_col(:,:) = spval
    allocate(this%cflx_input_litr_met_vr_col(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_met_vr_col(:,:) = spval
    allocate(this%cflx_input_litr_cel_vr_col(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_cel_vr_col(:,:) = spval
    allocate(this%cflx_input_litr_lig_vr_col(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_lig_vr_col(:,:) = spval
    allocate(this%cflx_input_litr_cwd_vr_col(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_cwd_vr_col(:,:) = spval
    allocate(this%fire_decomp_closs_col(begc:endc)); this%fire_decomp_closs_col(:) = spval
    allocate(this%som_c_runoff_col(begc:endc)); this%som_c_runoff_col(:) = spval
    allocate(this%co2_soi_flx_col(begc:endc)); this%co2_soi_flx_col(:) = spval
    allocate(this%phr_vr_col(begc:endc,1:nlevdecomp_full)); this%phr_vr_col(:,:)=spval
    allocate(this%o_scalar_col(begc:endc,1:nlevdecomp_full)); this%o_scalar_col(:,:)=spval
    allocate(this%somhr_col(begc:endc)); this%somhr_col(:) = spval
    allocate(this%lithr_col(begc:endc)); this%lithr_col(:) = spval

    allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:7)); this%decomp_k_col(:,:,:)=spval                 ! start adding   -zlyu
    allocate(this%t_scalar_col(begc:endc,1:nlevdecomp_full)); this%t_scalar_col(:,:)=spval
    allocate(this%w_scalar_col(begc:endc,1:nlevdecomp_full)); this%w_scalar_col(:,:)=spval                        ! end of adding  -zlyu
    
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use clm_varctl , only : nsrest, nsrStartup
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p,n,j,m            ! indices
    real(r8) ,pointer     :: gdp (:)                  ! global gdp data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: peatf (:)                ! global peatf data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: soilorder_rdin (:)       ! global soil order data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: abm (:)                  ! global abm data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: gti (:)                  ! read in - fmax (needs to be a pointer for use in ncdio)
    integer               :: dimid                    ! dimension id
    integer               :: ier                      ! error status
    type(file_desc_t)     :: ncid                     ! netcdf id
    logical               :: readvar
    character(len=256)    :: locfn                    ! local filename
    integer               :: begc, endc
    integer               :: begg, endg


    this%phenology_c_to_litr_met_c_col(:,:) = 0._r8
    this%phenology_c_to_litr_cel_c_col(:,:) = 0._r8
    this%dwt_livecrootc_to_cwdc_col(:,:) = 0._r8
    this%m_decomp_cpools_to_fire_vr_col(:,:,:)= 0._r8
    this%dwt_deadcrootc_to_cwdc_col(:,:) = 0._r8
    this%dwt_frootc_to_litr_lig_c_col(:,:) = 0._r8
    this%dwt_frootc_to_litr_cel_c_col(:,:) = 0._r8
    this%dwt_frootc_to_litr_met_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_litr_met_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_litr_cel_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_litr_lig_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_cwdc_col(:,:) = 0._r8
    this%harvest_c_to_litr_met_c_col(:,:) = 0._r8
    this%harvest_c_to_litr_cel_c_col(:,:) = 0._r8
    this%harvest_c_to_litr_lig_c_col(:,:) = 0._r8
    this%harvest_c_to_cwdc_col(:,:) = 0._r8
    this%m_c_to_litr_met_fire_col(:,:) = 0._r8
    this%m_c_to_litr_cel_fire_col(:,:) = 0._r8
    this%m_c_to_litr_lig_fire_col(:,:) = 0._r8
    this%fire_mortality_c_to_cwdc_col (:,:) = 0._r8
    this%phenology_c_to_litr_lig_c_col(:,:) = 0._r8
    this%cflx_input_litr_met_vr_col(:,:) = 0._r8
    this%cflx_input_litr_cel_vr_col(:,:) = 0._r8
    this%cflx_input_litr_lig_vr_col(:,:) = 0._r8
    this%cflx_input_litr_cwd_vr_col(:,:) = 0._r8

  end subroutine initCold

end module CNCarbonFluxType
