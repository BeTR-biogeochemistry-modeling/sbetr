module CNCarbonFluxType
  use elm_varcon             , only : spval, ispval, c14ratio
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use elm_varcon             , only : spval
  use elm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: carbonflux_type
    real(r8), pointer :: rr                                    (:)   => null()  ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: hr                                    (:)=> null()
    real(r8), pointer :: hr_vr                                 (:,:) => null()
    real(r8), pointer :: phenology_c_to_litr_met_c             (:,:)=> null()
    real(r8), pointer :: phenology_c_to_litr_cel_c             (:,:)=> null()
    real(r8), pointer :: phenology_c_to_litr_lig_c             (:,:)=> null()
    real(r8), pointer :: dwt_livecrootc_to_cwdc                (:,:)=> null()
    real(r8), pointer :: m_decomp_cpools_to_fire_vr            (:,:,:)=> null()
    real(r8), pointer :: dwt_deadcrootc_to_cwdc                (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_lig_c              (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_cel_c              (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_met_c              (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_met_c         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_cel_c         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_lig_c         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_cwdc               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_met_c               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_cel_c               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_lig_c               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_cwdc                     (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_met_fire                  (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_cel_fire                  (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_lig_fire                  (:,:)=> null()
    real(r8), pointer :: fire_mortality_c_to_cwdc              (:,:)=> null()
    real(r8), pointer :: fire_decomp_closs                     (:) => null()
    real(r8), pointer :: som_c_leached                         (:) => null()
    real(r8), pointer :: som_c_runoff                          (:) => null()
    real(r8), pointer :: cflx_input_litr_met_vr              (:,:) => null()
    real(r8), pointer :: cflx_input_litr_cel_vr              (:,:) => null()
    real(r8), pointer :: cflx_input_litr_lig_vr              (:,:) => null()
    real(r8), pointer :: cflx_input_litr_cwd_vr              (:,:) => null()
    real(r8), pointer :: co2_soi_flx                         (:) => null()
    real(r8), pointer :: decomp_k                            (:,:,:) => null()
    real(r8), pointer :: t_scalar(:,:) => null()
    real(r8), pointer :: w_scalar(:,:) => null()
    real(r8), pointer :: rr_vr(:,:) => null()
    real(r8), pointer :: phr_vr(:,:)=> null()
    real(r8), pointer :: somhr(:) => null()
    real(r8), pointer :: lithr(:)=> null()
    real(r8), pointer :: cwdc_hr(:)=> null()
    real(r8), pointer :: o_scalar(:,:)=> null()
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


    allocate(this%rr                  (begc:endc)) ; this%rr      (:)  = nan
    allocate(this%hr (begc:endc)); this%hr(:) = spval

    if(nlevdecomp_full>0)then
      allocate(this%phenology_c_to_litr_met_c(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_met_c(:,:) = spval
      allocate(this%phenology_c_to_litr_cel_c(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_cel_c(:,:) = spval
      allocate(this%dwt_livecrootc_to_cwdc(begc:endc,1:nlevdecomp_full)); this%dwt_livecrootc_to_cwdc(:,:) = spval
      allocate(this%m_decomp_cpools_to_fire_vr(begc:endc,1:nlevdecomp_full,1:7)); this%m_decomp_cpools_to_fire_vr(:,:,:)=nan
      allocate(this%dwt_deadcrootc_to_cwdc(begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootc_to_cwdc(:,:) = spval
      allocate(this%dwt_frootc_to_litr_lig_c(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_lig_c(:,:) = spval
      allocate(this%dwt_frootc_to_litr_cel_c(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_cel_c(:,:) = spval
      allocate(this%dwt_frootc_to_litr_met_c(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_met_c(:,:) = spval
      allocate(this%gap_mortality_c_to_litr_met_c(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_met_c(:,:) = spval
      allocate(this%gap_mortality_c_to_litr_cel_c(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_cel_c(:,:) = spval
      allocate(this%gap_mortality_c_to_litr_lig_c(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_lig_c(:,:) = spval
      allocate(this%gap_mortality_c_to_cwdc(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_cwdc(:,:) = spval
      allocate(this%harvest_c_to_litr_met_c(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_met_c(:,:) = spval
      allocate(this%harvest_c_to_litr_cel_c(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_cel_c(:,:) = spval
      allocate(this%harvest_c_to_litr_lig_c(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_lig_c(:,:) = spval
      allocate(this%harvest_c_to_cwdc(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_cwdc(:,:) = spval
      allocate(this%m_c_to_litr_met_fire(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_met_fire(:,:) = spval
      allocate(this%m_c_to_litr_cel_fire(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_cel_fire(:,:) = spval
      allocate(this%m_c_to_litr_lig_fire(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_lig_fire(:,:) = spval
      allocate(this%fire_mortality_c_to_cwdc (begc:endc,1:nlevdecomp_full)); this%fire_mortality_c_to_cwdc (:,:) = spval
      allocate(this%phenology_c_to_litr_lig_c(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_lig_c(:,:) = spval
      allocate(this%rr_vr(begc:endc,1:nlevdecomp_full)); this%rr_vr(:,:) = spval
      allocate(this%cflx_input_litr_met_vr(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_met_vr(:,:) = spval
      allocate(this%cflx_input_litr_cel_vr(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_cel_vr(:,:) = spval
      allocate(this%cflx_input_litr_lig_vr(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_lig_vr(:,:) = spval
      allocate(this%cflx_input_litr_cwd_vr(begc:endc,1:nlevdecomp_full)); this%cflx_input_litr_cwd_vr(:,:) = spval
      allocate(this%phr_vr(begc:endc,1:nlevdecomp_full)); this%phr_vr(:,:)=spval
      allocate(this%o_scalar(begc:endc,1:nlevdecomp_full)); this%o_scalar(:,:)=spval
    endif
    allocate(this%som_c_leached(begc:endc)); this%som_c_leached(:) = spval
    allocate(this%fire_decomp_closs(begc:endc)); this%fire_decomp_closs(:) = spval
    allocate(this%som_c_runoff(begc:endc)); this%som_c_runoff(:) = spval
    allocate(this%co2_soi_flx(begc:endc)); this%co2_soi_flx(:) = spval
    allocate(this%somhr(begc:endc)); this%somhr(:) = spval
    allocate(this%lithr(begc:endc)); this%lithr(:) = spval
    allocate(this%cwdc_hr(begc:endc)); this%cwdc_hr(:) = spval
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use elm_varctl , only : nsrest, nsrStartup
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


    this%phenology_c_to_litr_met_c(:,:) = 0._r8
    this%phenology_c_to_litr_cel_c(:,:) = 0._r8
    this%dwt_livecrootc_to_cwdc(:,:) = 0._r8
    this%m_decomp_cpools_to_fire_vr(:,:,:)= 0._r8
    this%dwt_deadcrootc_to_cwdc(:,:) = 0._r8
    this%dwt_frootc_to_litr_lig_c(:,:) = 0._r8
    this%dwt_frootc_to_litr_cel_c(:,:) = 0._r8
    this%dwt_frootc_to_litr_met_c(:,:) = 0._r8
    this%gap_mortality_c_to_litr_met_c(:,:) = 0._r8
    this%gap_mortality_c_to_litr_cel_c(:,:) = 0._r8
    this%gap_mortality_c_to_litr_lig_c(:,:) = 0._r8
    this%gap_mortality_c_to_cwdc(:,:) = 0._r8
    this%harvest_c_to_litr_met_c(:,:) = 0._r8
    this%harvest_c_to_litr_cel_c(:,:) = 0._r8
    this%harvest_c_to_litr_lig_c(:,:) = 0._r8
    this%harvest_c_to_cwdc(:,:) = 0._r8
    this%m_c_to_litr_met_fire(:,:) = 0._r8
    this%m_c_to_litr_cel_fire(:,:) = 0._r8
    this%m_c_to_litr_lig_fire(:,:) = 0._r8
    this%fire_mortality_c_to_cwdc (:,:) = 0._r8
    this%phenology_c_to_litr_lig_c(:,:) = 0._r8
    this%cflx_input_litr_met_vr(:,:) = 0._r8
    this%cflx_input_litr_cel_vr(:,:) = 0._r8
    this%cflx_input_litr_lig_vr(:,:) = 0._r8
    this%cflx_input_litr_cwd_vr(:,:) = 0._r8

  end subroutine initCold

end module CNCarbonFluxType
