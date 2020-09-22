module PhosphorusFluxType
  use elm_varcon             , only : spval, ispval
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use elm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: phosphorusflux_type
    real(r8), pointer :: sminp_leached                         (:)    => null() !col inorganic P leaching loss, gP/m2/time step
    real(r8), pointer :: sminp_runoff                          (:)   => null()  !col inorganic P runoff loss, gP/m2/time step
    real(r8), pointer :: biochem_pmin_vr                       (:,:) => null()  ! col vertically-resolved total biochemical P mineralization (gP/m3/s)
    real(r8), pointer :: phenology_p_to_litr_met_p             (:,:)=> null()
    real(r8), pointer :: dwt_frootp_to_litr_met_p              (:,:)=> null()
    real(r8), pointer :: phenology_p_to_litr_cel_p             (:,:)=> null()
    real(r8), pointer :: dwt_livecrootp_to_cwdp                (:,:)=> null()
    real(r8), pointer :: phenology_p_to_litr_lig_p             (:,:)=> null()
    real(r8), pointer :: m_decomp_ppools_to_fire_vr            (:,:,:)=> null()
    real(r8), pointer :: dwt_deadcrootp_to_cwdp                (:,:)=> null()
    real(r8), pointer :: dwt_frootp_to_litr_lig_p              (:,:)=> null()
    real(r8), pointer :: dwt_frootp_to_litr_cel_p              (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_litr_met_p         (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_cwdp               (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_litr_lig_p         (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_litr_cel_p         (:,:)=> null()
    real(r8), pointer :: harvest_p_to_litr_met_p               (:,:)=> null()
    real(r8), pointer :: harvest_p_to_cwdp                     (:,:)=> null()
    real(r8), pointer :: harvest_p_to_litr_lig_p               (:,:)=> null()
    real(r8), pointer :: harvest_p_to_litr_cel_p               (:,:)=> null()
    real(r8), pointer :: m_p_to_litr_met_fire                  (:,:)=> null()
    real(r8), pointer :: m_p_to_litr_cel_fire                  (:,:)=> null()
    real(r8), pointer :: m_p_to_litr_lig_fire                  (:,:)=> null()
    real(r8), pointer :: fire_mortality_p_to_cwdp              (:,:)=> null()
    real(r8), pointer :: primp_to_labilep_vr                   (:,:)=> null()
    real(r8), pointer :: pdep_to_sminp                         (:)=> null()
    real(r8), pointer :: fert_p_to_sminp                       (:)=> null()
    real(r8), pointer :: sminp_to_plant_patch                      (:)=> null()
    real(r8), pointer :: supplement_to_sminp                   (:)=> null()
    real(r8), pointer :: secondp_to_occlp                      (:)=> null()
    real(r8), pointer :: fire_decomp_ploss                     (:)=> null()
    real(r8), pointer :: som_p_leached                         (:) => null()
    real(r8), pointer :: som_p_runoff                          (:) => null()
    real(r8), pointer :: pflx_input_litr_met_vr                (:,:) => null()
    real(r8), pointer :: pflx_input_litr_cel_vr                (:,:) => null()
    real(r8), pointer :: pflx_input_litr_lig_vr                (:,:) => null()
    real(r8), pointer :: pflx_input_litr_cwd_vr                (:,:) => null()
    real(r8), pointer :: pflx_minp_input_po4_vr                (:,:) => null()
    real(r8), pointer :: sminp(:)
    real(r8), pointer :: occlp(:)
    real(r8), pointer :: primp_to_labilep(:) => null()
    real(r8), pointer :: sminp_to_plant_trans_patch(:) => null()
    real(r8), pointer :: adsorb_to_labilep_vr(:,:) => null()
    real(r8), pointer :: sminp_to_plant_vr(:,:) => null()
    real(r8), pointer :: supplement_to_sminp_vr(:,:) => null()
    real(r8), pointer :: net_mineralization_p_vr(:,:) => null()
    real(r8), pointer :: col_plant_pdemand_vr                    (:,:) => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
    procedure, private :: SetValues
  end type phosphorusflux_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(phosphorusflux_type) :: this
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
    class(phosphorusflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%sminp_leached     (begc:endc              ))     ;this%sminp_leached          (:)    = spval
    allocate(this%sminp_runoff      (begc:endc              ))      ;this%sminp_runoff         (:)    = spval
    allocate(this%biochem_pmin_vr   (begc:endc,1:nlevdecomp_full))   ;this%biochem_pmin_vr      (:,:) = spval
    allocate(this%sminp_to_plant_trans_patch(begp:endp))
    allocate(this%phenology_p_to_litr_met_p(begc:endc, 1:nlevdecomp_full)); this%phenology_p_to_litr_met_p(:,:) = spval
    allocate(this%dwt_frootp_to_litr_met_p(begc:endc, 1:nlevdecomp_full)); this%dwt_frootp_to_litr_met_p(:,:) = spval
    allocate(this%phenology_p_to_litr_cel_p(begc:endc,1:nlevdecomp_full)); this%phenology_p_to_litr_cel_p(:,:) = spval
    allocate(this%dwt_livecrootp_to_cwdp(begc:endc,1:nlevdecomp_full)); this%dwt_livecrootp_to_cwdp(:,:) = spval
    allocate(this%phenology_p_to_litr_lig_p(begc:endc, 1:nlevdecomp_full)); this%phenology_p_to_litr_lig_p(:,:) = spval
    allocate(this%m_decomp_ppools_to_fire_vr(begc:endc,1:nlevdecomp_full,1:7)); this%m_decomp_ppools_to_fire_vr(:,:,:) = spval
    allocate(this%dwt_deadcrootp_to_cwdp(begc:endc, 1:nlevdecomp_full)); this%dwt_deadcrootp_to_cwdp(:,:) = spval
    allocate(this%dwt_frootp_to_litr_lig_p(begc:endc,1:nlevdecomp_full)); this%dwt_frootp_to_litr_lig_p(:,:) = spval
    allocate(this%dwt_frootp_to_litr_cel_p(begc:endc,1:nlevdecomp_full)); this%dwt_frootp_to_litr_cel_p(:,:) = spval
    allocate(this%gap_mortality_p_to_litr_met_p(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_litr_met_p(:,:) = spval
    allocate(this%gap_mortality_p_to_cwdp(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_cwdp(:,:) = spval
    allocate(this%gap_mortality_p_to_litr_lig_p(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_litr_lig_p(:,:) = spval
    allocate(this%gap_mortality_p_to_litr_cel_p(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_litr_cel_p(:,:) = spval
    allocate(this%harvest_p_to_litr_met_p(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_litr_met_p(:,:) = spval
    allocate(this%harvest_p_to_cwdp(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_cwdp(:,:) = spval
    allocate(this%harvest_p_to_litr_lig_p(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_litr_lig_p(:,:) = spval
    allocate(this%harvest_p_to_litr_cel_p(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_litr_cel_p(:,:) = spval
    allocate(this%m_p_to_litr_met_fire(begc:endc,1:nlevdecomp_full)); this%m_p_to_litr_met_fire(:,:) = spval
    allocate(this%m_p_to_litr_cel_fire(begc:endc,1:nlevdecomp_full)); this%m_p_to_litr_cel_fire(:,:) = spval
    allocate(this%m_p_to_litr_lig_fire(begc:endc,1:nlevdecomp_full)); this%m_p_to_litr_lig_fire(:,:) = spval
    allocate(this%fire_mortality_p_to_cwdp(begc:endc,1:nlevdecomp_full)); this%fire_mortality_p_to_cwdp(:,:) = spval
    allocate(this%primp_to_labilep_vr(begc:endc,1:nlevdecomp_full)); this%primp_to_labilep_vr(:,:) = spval
    allocate(this%pdep_to_sminp(begc:endc)); this%pdep_to_sminp(begc:endc) = spval
    allocate(this%pdep_to_sminp(begc:endc)); this%pdep_to_sminp(begc:endc) = spval
    allocate(this%fert_p_to_sminp(begc:endc)); this%fert_p_to_sminp(begc:endc) = spval
    allocate(this%som_p_leached(begc:endc)); this%som_p_leached(:) = spval

    allocate(this%pflx_input_litr_met_vr(begc:endc,1:nlevdecomp_full)); this%pflx_input_litr_met_vr(:,:) = spval
    allocate(this%pflx_input_litr_cel_vr(begc:endc,1:nlevdecomp_full)); this%pflx_input_litr_cel_vr(:,:) = spval
    allocate(this%pflx_input_litr_lig_vr(begc:endc,1:nlevdecomp_full)); this%pflx_input_litr_lig_vr(:,:) = spval
    allocate(this%pflx_input_litr_cwd_vr(begc:endc,1:nlevdecomp_full)); this%pflx_input_litr_cwd_vr(:,:) = spval
    allocate(this%pflx_minp_input_po4_vr(begc:endc,1:nlevdecomp_full)); this%pflx_minp_input_po4_vr(:,:) = spval
    allocate(this%supplement_to_sminp(begc:endc)); this%supplement_to_sminp(:) = spval
    allocate(this%secondp_to_occlp(begc:endc)); this%secondp_to_occlp(:) = spval
    allocate(this%fire_decomp_ploss(begc:endc)); this%fire_decomp_ploss(:) = spval
    allocate(this%occlp(begc:endc)); this%occlp(:) = spval
    allocate(this%sminp(begc:endc)); this%sminp(:) = spval
    allocate(this%som_p_runoff(begc:endc)); this%som_p_runoff(:) = spval
    allocate(this%primp_to_labilep(begc:endc)); this%primp_to_labilep(:) = spval
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
    class(phosphorusflux_type) :: this
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


    this%phenology_p_to_litr_met_p(:,:) = 0._r8
    this%dwt_frootp_to_litr_met_p(:,:) = 0._r8
    this%phenology_p_to_litr_cel_p(:,:) = 0._r8
    this%dwt_livecrootp_to_cwdp(:,:) = 0._r8
    this%phenology_p_to_litr_lig_p(:,:) = 0._r8
    this%m_decomp_ppools_to_fire_vr(:,:,:) = 0._r8
    this%dwt_deadcrootp_to_cwdp(:,:) = 0._r8
    this%dwt_frootp_to_litr_lig_p(:,:) = 0._r8
    this%dwt_frootp_to_litr_cel_p(:,:) = 0._r8
    this%gap_mortality_p_to_litr_met_p(:,:) = 0._r8
    this%gap_mortality_p_to_cwdp(:,:) = 0._r8
    this%gap_mortality_p_to_litr_lig_p(:,:) = 0._r8
    this%gap_mortality_p_to_litr_cel_p(:,:) = 0._r8
    this%harvest_p_to_litr_met_p(:,:) = 0._r8
    this%harvest_p_to_cwdp(:,:) = 0._r8
    this%harvest_p_to_litr_lig_p(:,:) = 0._r8
    this%harvest_p_to_litr_cel_p(:,:) = 0._r8
    this%m_p_to_litr_met_fire(:,:) = 0._r8
    this%m_p_to_litr_cel_fire(:,:) = 0._r8
    this%m_p_to_litr_lig_fire(:,:) = 0._r8
    this%fire_mortality_p_to_cwdp(:,:) = 0._r8
    this%primp_to_labilep_vr(:,:) = 0._r8
    this%pdep_to_sminp(:) = 0._r8
    this%pdep_to_sminp(:) = 0._r8
    this%fert_p_to_sminp(:) = 0._r8

    this%pflx_input_litr_met_vr(:,:) = 0._r8
    this%pflx_input_litr_cel_vr(:,:) = 0._r8
    this%pflx_input_litr_lig_vr(:,:) = 0._r8
    this%pflx_input_litr_cwd_vr(:,:) = 0._r8
    this%pflx_minp_input_po4_vr(:,:) = 0._r8

  end subroutine initCold


  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class (phosphorusflux_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column

    integer :: fi, i

    do fi = 1,num_column
       i = filter_column(fi)

    enddo
  end subroutine SetValues
end module PhosphorusFluxType
