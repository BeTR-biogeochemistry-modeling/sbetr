module CNNitrogenFluxType
  use elm_varcon             , only : spval, ispval
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use elm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: nitrogenflux_type
     real(r8), pointer :: smin_no3_leached                      (:)    => null() ! col soil mineral NO3 pool loss to leaching (gN/m2/s)
     real(r8), pointer :: smin_no3_runoff                       (:)   => null()  ! col soil mineral NO3 pool loss to runoff (gN/m2/s)
     real(r8), pointer :: smin_nh4_leached                      (:)    => null() ! col soil mineral NO3 pool loss to leaching (gN/m2/s)
     real(r8), pointer :: smin_nh4_runoff                       (:)   => null()  ! col soil mineral NO3 pool loss to runoff (gN/m2/s)
     real(r8), pointer :: f_n2o_denit                           (:)  => null()   ! col flux of N2o from denitrification [gN/m^2/s]
     real(r8), pointer :: f_n2o_nit                             (:)  => null()   ! col flux of N2o from nitrification [gN/m^2/s]
     real(r8), pointer :: f_nit                                 (:)=> null()
     real(r8), pointer :: f_denit                               (:)=> null()
     real(r8), pointer :: supplement_to_sminn_vr                (:,:)=> null()
     real(r8), pointer :: soyfixn_to_sminn                      (:)=> null()
     real(r8), pointer :: m_decomp_npools_to_fire_vr            (:,:,:)=> null()
     real(r8), pointer :: ndep_to_sminn                         (:)=> null()
     real(r8), pointer :: ndep_to_smin_nh3                         (:)=> null()
     real(r8), pointer :: ndep_to_smin_no3                         (:)=> null()
     real(r8), pointer :: dwt_livecrootn_to_cwdn                (:,:)=> null()
     real(r8), pointer :: phenology_n_to_litr_lig_n             (:,:)=> null()
     real(r8), pointer :: nfix_to_sminn                         (:)=> null()
     real(r8), pointer :: dwt_deadcrootn_to_cwdn                (:,:)=> null()
     real(r8), pointer :: dwt_frootn_to_litr_lig_n              (:,:)=> null()
     real(r8), pointer :: phenology_n_to_litr_cel_n             (:,:)=> null()
     real(r8), pointer :: phenology_n_to_litr_met_n             (:,:)=> null()
     real(r8), pointer :: gap_mortality_n_to_cwdn               (:,:)=> null()
     real(r8), pointer :: gap_mortality_n_to_litr_lig_n         (:,:)=> null()
     real(r8), pointer :: dwt_frootn_to_litr_cel_n              (:,:)=> null()
     real(r8), pointer :: dwt_frootn_to_litr_met_n              (:,:)=> null()
     real(r8), pointer :: harvest_n_to_cwdn                     (:,:)=> null()
     real(r8), pointer :: harvest_n_to_litr_lig_n               (:,:)=> null()
     real(r8), pointer :: gap_mortality_n_to_litr_cel_n         (:,:)=> null()
     real(r8), pointer :: gap_mortality_n_to_litr_met_n         (:,:)=> null()
     real(r8), pointer :: fire_mortality_n_to_cwdn              (:,:)=> null()
     real(r8), pointer :: m_n_to_litr_lig_fire                  (:,:)=> null()
     real(r8), pointer :: harvest_n_to_litr_cel_n               (:,:)=> null()
     real(r8), pointer :: harvest_n_to_litr_met_n               (:,:)=> null()
     real(r8), pointer :: m_n_to_litr_cel_fire                  (:,:)=> null()
     real(r8), pointer :: m_n_to_litr_met_fire                  (:,:)=> null()
     real(r8), pointer :: fert_to_sminn                         (:) => null()
     real(r8), pointer :: ndep_to_sminn_nh3                     (:) => null()
     real(r8), pointer :: ndep_to_sminn_no3                     (:) => null()

     real(r8), pointer :: fire_decomp_nloss                     (:) => null()
     real(r8), pointer :: denit                                 (:) => null()
     real(r8), pointer :: som_n_leached                         (:) => null()
     real(r8), pointer :: supplement_to_sminn                   (:) => null()
     real(r8), pointer :: som_n_runoff                          (:) => null()
     real(r8), pointer :: nflx_input_litr_met_vr              (:,:) => null()
     real(r8), pointer :: nflx_input_litr_cel_vr              (:,:) => null()
     real(r8), pointer :: nflx_input_litr_lig_vr              (:,:) => null()
     real(r8), pointer :: nflx_input_litr_cwd_vr              (:,:) => null()
     real(r8), pointer :: nflx_minn_input_nh4_vr              (:,:) => null()
     real(r8), pointer :: nflx_minn_input_no3_vr              (:,:) => null()
     real(r8), pointer :: nh3_soi_flx                           (:) => null()
     real(r8), pointer :: f_nit_vr                            (:,:) => null()
     real(r8), pointer :: f_n2o_nit_vr                        (:,:) => null()
     real(r8), pointer :: f_denit_vr                          (:,:) => null()
     real(r8), pointer :: f_n2o_denit_vr                      (:,:) => null()
     real(r8), pointer :: smin_nh4_to_plant_vr                (:,:) => null()
     real(r8), pointer :: smin_no3_to_plant_vr                (:,:) => null()
     real(r8), pointer :: actual_immob                        (:)   => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type nitrogenflux_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(nitrogenflux_type) :: this
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
    class(nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%smin_no3_runoff         (begc:endc)); this%smin_no3_runoff   (:)   = spval
    allocate(this%smin_no3_leached        (begc:endc)); this%smin_no3_leached  (:)   = spval
    allocate(this%smin_nh4_runoff         (begc:endc)); this%smin_nh4_runoff   (:)   = spval
    allocate(this%smin_nh4_leached        (begc:endc)); this%smin_nh4_leached  (:)   = spval
    allocate(this%f_n2o_denit             (begc:endc)); this%f_n2o_denit       (:)   = spval
    allocate(this%f_n2o_nit               (begc:endc)); this%f_n2o_nit         (:)   = spval
    allocate(this%f_nit                   (begc:endc)); this%f_nit             (:)   = spval
    allocate(this%f_denit                 (begc:endc)); this%f_denit           (:)   = spval
    allocate(this%fert_to_sminn          (begc:endc)); this%fert_to_sminn    (:) = spval
    allocate(this%nfix_to_sminn          (begc:endc)); this%nfix_to_sminn(:) = spval

    allocate(this%soyfixn_to_sminn(begc:endc)); this%soyfixn_to_sminn(:) = spval
    allocate(this%m_decomp_npools_to_fire_vr(begc:endc,1:nlevdecomp_full,1:7)); this%m_decomp_npools_to_fire_vr(:,:,:) =nan
    allocate(this%ndep_to_sminn(begc:endc)); this%ndep_to_sminn(:) = spval
    allocate(this%ndep_to_smin_nh3(begc:endc)); this%ndep_to_smin_nh3(:) = spval
    allocate(this%ndep_to_smin_no3(begc:endc)); this%ndep_to_smin_no3(:) = spval
    allocate(this%dwt_livecrootn_to_cwdn(begc:endc,1:nlevdecomp_full)); this%dwt_livecrootn_to_cwdn(:,:) = spval
    allocate(this%phenology_n_to_litr_lig_n(begc:endc,1:nlevdecomp_full)); this%phenology_n_to_litr_lig_n(:,:) = spval
    allocate(this%dwt_deadcrootn_to_cwdn(begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootn_to_cwdn(:,:) = spval
    allocate(this%dwt_frootn_to_litr_lig_n(begc:endc,1:nlevdecomp_full)); this%dwt_frootn_to_litr_lig_n(:,:) = spval
    allocate(this%phenology_n_to_litr_cel_n(begc:endc,1:nlevdecomp_full)); this%phenology_n_to_litr_cel_n(:,:) = spval
    allocate(this%phenology_n_to_litr_met_n(begc:endc,1:nlevdecomp_full)); this%phenology_n_to_litr_met_n(:,:) = spval
    allocate(this%gap_mortality_n_to_cwdn(begc:endc,1:nlevdecomp_full)); this%gap_mortality_n_to_cwdn(:,:) = spval
    allocate(this%gap_mortality_n_to_litr_lig_n(begc:endc,1:nlevdecomp_full)); this%gap_mortality_n_to_litr_lig_n(:,:) = spval
    allocate(this%dwt_frootn_to_litr_cel_n(begc:endc,1:nlevdecomp_full)); this%dwt_frootn_to_litr_cel_n(:,:) = spval
    allocate(this%dwt_frootn_to_litr_met_n(begc:endc,1:nlevdecomp_full)); this%dwt_frootn_to_litr_met_n(:,:) = spval
    allocate(this%harvest_n_to_cwdn(begc:endc,1:nlevdecomp_full)); this%harvest_n_to_cwdn(:,:) = spval
    allocate(this%harvest_n_to_litr_lig_n(begc:endc,1:nlevdecomp_full)); this%harvest_n_to_litr_lig_n(:,:) = spval
    allocate(this%gap_mortality_n_to_litr_cel_n(begc:endc,1:nlevdecomp_full)); this%gap_mortality_n_to_litr_cel_n(:,:) = spval
    allocate(this%gap_mortality_n_to_litr_met_n(begc:endc,1:nlevdecomp_full)); this%gap_mortality_n_to_litr_met_n(:,:) = spval
    allocate(this%fire_mortality_n_to_cwdn(begc:endc,1:nlevdecomp_full)); this%fire_mortality_n_to_cwdn(:,:) = spval
    allocate(this%m_n_to_litr_lig_fire(begc:endc,1:nlevdecomp_full)); this%m_n_to_litr_lig_fire(:,:) = spval
    allocate(this%harvest_n_to_litr_cel_n(begc:endc,1:nlevdecomp_full)); this%harvest_n_to_litr_cel_n(:,:) = spval
    allocate(this%harvest_n_to_litr_met_n(begc:endc,1:nlevdecomp_full)); this%harvest_n_to_litr_met_n(:,:) = spval
    allocate(this%m_n_to_litr_cel_fire(begc:endc,1:nlevdecomp_full)); this%m_n_to_litr_cel_fire(:,:) = spval
    allocate(this%m_n_to_litr_met_fire(begc:endc,1:nlevdecomp_full)); this%m_n_to_litr_met_fire(:,:) = spval
    allocate(this%denit(begc:endc)); this%denit(:) = spval
    allocate(this%som_n_leached(begc:endc)); this%som_n_leached(:) = spval
    allocate(this%supplement_to_sminn(begc:endc)); this%supplement_to_sminn(:) = spval
    allocate(this%f_n2o_denit_vr          (begc:endc, 1:nlevdecomp_full)); this%f_n2o_denit_vr(:,:) = spval

    allocate(this%nflx_input_litr_met_vr(begc:endc,1:nlevdecomp_full)); this%nflx_input_litr_met_vr(:,:) = spval
    allocate(this%nflx_input_litr_cel_vr(begc:endc,1:nlevdecomp_full)); this%nflx_input_litr_cel_vr(:,:) = spval
    allocate(this%nflx_input_litr_lig_vr(begc:endc,1:nlevdecomp_full)); this%nflx_input_litr_lig_vr(:,:) = spval
    allocate(this%nflx_input_litr_cwd_vr(begc:endc,1:nlevdecomp_full)); this%nflx_input_litr_cwd_vr(:,:) = spval
    allocate(this%nflx_minn_input_nh4_vr(begc:endc,1:nlevdecomp_full)); this%nflx_minn_input_nh4_vr(:,:) = spval
    allocate(this%nflx_minn_input_no3_vr(begc:endc,1:nlevdecomp_full)); this%nflx_minn_input_no3_vr(:,:) = spval
    allocate(this%fire_decomp_nloss(begc:endc)); this%fire_decomp_nloss(:) = spval
    allocate(this%som_n_runoff(begc:endc)); this%som_n_runoff(:) = spval
    allocate(this%nh3_soi_flx(begc:endc)); this%nh3_soi_flx(:) = spval
    allocate(this%actual_immob(begc:endc)); this%actual_immob(:)=spval
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
    class(nitrogenflux_type) :: this
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


    this%fert_to_sminn    (:) = 0._r8
    this%soyfixn_to_sminn(:) = 0._r8
    this%m_decomp_npools_to_fire_vr(:,:,:) = 0._r8
    this%ndep_to_sminn(:) = 0._r8
    this%ndep_to_smin_nh3(:) = 0._r8
    this%ndep_to_smin_no3(:) = 0._r8
    this%dwt_livecrootn_to_cwdn(:,:) = 0._r8
    this%phenology_n_to_litr_lig_n(:,:) = 0._r8
    this%dwt_deadcrootn_to_cwdn(:,:) = 0._r8
    this%dwt_frootn_to_litr_lig_n(:,:) = 0._r8
    this%phenology_n_to_litr_cel_n(:,:) = 0._r8
    this%phenology_n_to_litr_met_n(:,:) = 0._r8
    this%gap_mortality_n_to_cwdn(:,:) = 0._r8
    this%gap_mortality_n_to_litr_lig_n(:,:) = 0._r8
    this%dwt_frootn_to_litr_cel_n(:,:) = 0._r8
    this%dwt_frootn_to_litr_met_n(:,:) = 0._r8
    this%harvest_n_to_cwdn(:,:) = 0._r8
    this%harvest_n_to_litr_lig_n(:,:) = 0._r8
    this%gap_mortality_n_to_litr_cel_n(:,:) = 0._r8
    this%gap_mortality_n_to_litr_met_n(:,:) = 0._r8
    this%fire_mortality_n_to_cwdn(:,:) = 0._r8
    this%m_n_to_litr_lig_fire(:,:) = 0._r8
    this%harvest_n_to_litr_cel_n(:,:) = 0._r8
    this%harvest_n_to_litr_met_n(:,:) = 0._r8
    this%m_n_to_litr_cel_fire(:,:) = 0._r8
    this%m_n_to_litr_met_fire(:,:) = 0._r8
    this%nfix_to_sminn(:) = 0._r8
    this%supplement_to_sminn(:) = 0._r8
    this%nflx_input_litr_met_vr(:,:) = 0._r8
    this%nflx_input_litr_cel_vr(:,:) = 0._r8
    this%nflx_input_litr_lig_vr(:,:) = 0._r8
    this%nflx_input_litr_cwd_vr(:,:) = 0._r8
    this%nflx_minn_input_nh4_vr(:,:) = 0._r8
    this%nflx_minn_input_no3_vr(:,:) = 0._r8
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
    class (nitrogenflux_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column

    integer :: fi, i

    do fi = 1,num_column
       i = filter_column(fi)
       this%smin_no3_leached(i)       = value_column
       this%smin_no3_runoff(i)        = value_column
    enddo
  end subroutine SetValues
end module CNNitrogenFluxType
