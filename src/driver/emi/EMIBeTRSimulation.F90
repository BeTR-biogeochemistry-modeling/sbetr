module EMIBeTRSimulation
  !
  ! !DESCRIPTION:
  !  BeTR simulation class for emi.
  !  This class is independent from the shared types
  !  BeTR simulation class are API definitions, mapping data
  !  structures from a specific LSM, e.g. CLM, ALM, into BeTR data
  !  structures.
  !

  use abortutils     , only : endrun
  use clm_varctl     , only : iulog, use_cn
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use tracer_varcon  , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use BeTR_decompMod , only : betr_bounds_type
  use decompMod      , only : bounds_type
  use ColumnType     , only : column_type => column_physical_properties
  use VegetationType , only : patch_type  => vegetation_physical_properties
  use LandunitType   , only : landunit_type => landunit_physical_properties
  use VegetationDataType, only : veg_es, veg_wf
  use ColumnDataType    , only : col_es, col_ws, col_wf
  use BetrType                 , only : betr_type, create_betr_type
  use betr_ctrl                , only : max_betr_hist_type,max_betr_rest_type, biulog
  use betr_constants           , only : betr_string_length
  use betr_constants           , only : betr_filename_length
  use betr_regression_module   , only : betr_regression_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type, create_betr_biogeophys_input
  use BeTR_biogeoStateType     , only : betr_biogeo_state_type, create_betr_biogeo_state
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type, create_betr_biogeoFlux
  use BeTR_TimeMod             , only : betr_time_type
  use betr_varcon              , only : betr_maxpatch_pft
  use BetrStatusType           , only : betr_status_type, create_betr_status_type
  use BetrStatusSimType        , only : betr_status_sim_type, create_betr_status_sim_type
  use betr_columnType          , only : betr_column_type, create_betr_column_type
  use betr_patchType           , only : betr_patch_type, create_betr_patch_type
  use betr_varcon              , only : spval => bspval
  use BeTRHistVarType          , only : betr_hist_var_type
  use pftvarcon                , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
  use clm_time_manager         , only : get_curr_date
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: emi_betr_simulation_type
     type(betr_type)                     , public, pointer  :: betr(:)  => null()
     type(betr_biogeophys_input_type)    , public, pointer  :: biophys_forc(:) => null()
     type(betr_biogeo_state_type)        , public, pointer  :: biogeo_state(:) => null()
     type(betr_biogeo_flux_type)         , public, pointer  :: biogeo_flux(:) => null()
     type(betr_status_type)              , public, pointer  :: bstatus(:) => null()
     type(betr_status_sim_type)          , public, pointer  :: bsimstatus
     logical                             , public, pointer :: active_col(:) => null()
     type(betr_column_type)              , public, pointer :: betr_col(:) => null()
     type(betr_patch_type)               , public, pointer :: betr_pft(:) => null()
     type(betr_bounds_type)              , public, pointer :: betr_bounds(:) => null()
     character(len=betr_filename_length) , private :: base_filename
     character(len=betr_filename_length) , private :: hist_filename
     character(len=64)                   , private :: case_id

     type(betr_regression_type), private          :: regression
     type(betr_time_type), public                 :: betr_time
     integer, public                              :: num_soilp
     integer, public, allocatable                 :: filter_soilp(:)
     integer, public                              :: num_jtops
     integer, public, allocatable                 :: jtops(:)
     integer, public                              :: num_soilc
     integer, public, allocatable                 :: filter_soilc(:)
     integer, public :: spinup_count
     type(betr_hist_var_type), allocatable :: state_hist1d_var(:)
     type(betr_hist_var_type), allocatable :: state_hist2d_var(:)
     type(betr_hist_var_type), allocatable :: flux_hist1d_var(:)
     type(betr_hist_var_type), allocatable :: flux_hist2d_var(:)
     integer :: num_hist_state1d
     integer :: num_hist_state2d
     integer :: num_hist_flux1d
     integer :: num_hist_flux2d

     integer :: num_rest_state1d
     integer :: num_rest_state2d
     real(r8), pointer :: hist_states_2d(:,:,:)
     real(r8), pointer :: hist_states_1d(:,:)
     real(r8), pointer :: rest_states_2d(:,:,:)
     real(r8), pointer :: rest_states_1d(:,:)
     real(r8), pointer :: hist_fluxes_1d(:,:)
     real(r8), pointer :: hist_fluxes_2d(:,:,:)
     real(r8), pointer :: hist_fluxes_1d_accum(:,:)
     real(r8), pointer :: hist_fluxes_2d_accum(:,:,:)
     real(r8), pointer :: scalaravg_col(:)
     real(r8), pointer :: dom_scalar_col(:)
     logical,  private :: active_soibgc
     real(r8), private :: hist_naccum
     integer , private :: hist_record

     ! FIXME(bja, 201603) most of these types should be private!

     ! NOTE(bja, 201603) BeTR types only, no LSM specific types here!

   contains
     procedure, public :: BeTRInit
     procedure, public :: BeTRSetFilter
     procedure, public :: SetClock
     procedure, public :: InitOnline              => EMIBeTRSimulationInit
     procedure, public :: ConsistencyCheck        => EMIBeTRSimulationConsistencyCheck
     procedure, public :: PreDiagSoilColWaterFlux => EMIBeTRSimulationPreDiagSoilColWaterFlux
     procedure, public :: DiagnoseDtracerFreezeThaw=> EMIBeTRSimulationDiagnoseDtracerFreezeThaw
     procedure, public :: CalcDewSubFlux          => EMIBeTRSimulationCalcDewSubFlux
     procedure, public :: DiagAdvWaterFlux        => EMIBeTRSimulationDiagAdvWaterFlux
     procedure, public :: DiagDrainWaterFlux      => EMIBeTRSimulationDiagDrainWaterFlux
     procedure, public :: BeginSnowLayerAdjst     => EMIBeTRSimulationBeginTracerSnowLayerAdjst
     procedure, public :: EndSnowLayerAdjst       => EMIBeTRSimulationEndTracerSnowLayerAdjst
     procedure, public :: CombineSnowLayers       => EMIBeTRSimulationCombineSnowLayers
     procedure, public :: DvideSnowLayers         => EMIBeTRSimulationDvideSnowLayers
     procedure, public :: StepWithoutDrainage     => EMIBeTRSimulationStepWithoutDrainage
     procedure, public :: StepWithDrainage        => EMIBeTRSimulationStepWithDrainage
     procedure, public :: BeginMassBalanceCheck   => EMIBeTRSimulationBeginMassBalanceCheck
     procedure, public :: MassBalanceCheck        => EMIBeTRSimulationMassBalanceCheck
     procedure, public :: RetrieveBiogeoFlux      => EMIBeTRSimulationRetrieveBiogeoFlux
     procedure, public :: DiagnoseLnd2atm         => EMIBeTRSimulationDiagnoseLnd2atm

     procedure, public :: SetSpinup               => EMIBeTRSimulationSetSpinup
     procedure, public :: ReadParams              => EMIBeTRSimulationReadParams
     !the following are used to interact with lsm
     procedure, public :: do_soibgc
     procedure, public :: BeTRRestart             => EMIBeTRSimulationRestart
     procedure, public  :: HistRetrieval          => EMIBeTRSimulationHistRetrieval
     procedure, public :: BeTRSetcps              => EMIBeTRSimulationSetcps
     procedure, public :: BeTRSetBounds           => EMIBeTRSimulationSetBounds
     procedure, private :: BeTRCreateHistory      => EMIBeTRSimulationCreateHistory
     procedure, private :: BeTRRetrieveHistoryState    => EMIBeTRSimulationRetrieveHistoryState
     procedure, private :: BeTRRetrieveHistoryFlux    => EMIBeTRSimulationRetrieveHistoryFlux
     procedure, private:: hist_create_states
     procedure, private:: hist_create_fluxes
     procedure, private:: RestAlloc               => EMIBeTRSimulationRestartAlloc
     procedure, private:: HistAlloc               => EMIBeTRSimulationHistoryAlloc
     procedure, private:: set_activecol
  end type emi_betr_simulation_type
  public :: create_betr_simulation_emi
contains
!-------------------------------------------------------------------------------

  function create_betr_simulation_emi()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(emi_betr_simulation_type), pointer :: create_betr_simulation_emi
    class(emi_betr_simulation_type), pointer :: simulation

    allocate(simulation)
    create_betr_simulation_emi => simulation

  end function create_betr_simulation_emi
!-------------------------------------------------------------------------------
  subroutine SetClock(this, dtime, nelapstep)
  implicit none
  class(emi_betr_simulation_type)              , intent(inout) :: this
  real(r8), intent(in) :: dtime
  integer, intent(in) :: nelapstep


  call this%betr_time%setClock(dtime, nelapstep)

  end subroutine SetClock
  !-------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationInit(this, bounds, lun, col, pft)
    !
    ! DESCRIPTIONS
    ! Dummy routine for inheritance purposes. don't use.
    !
    !USES
    use landunit_varcon , only : istcrop, istice, istsoil
    use clm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil
    !betr types
    use betr_constants      , only : betr_filename_length
    use betr_constants      , only : betr_namelist_buffer_size
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use BeTR_landvarconType , only : betr_landvarcon
    use BeTR_decompMod      , only : betr_bounds_type
    implicit none

    class(emi_betr_simulation_type)              , intent(inout) :: this
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout) :: col
    type(patch_type)                         , intent(in) :: pft
    type(bounds_type)                        , intent(in)    :: bounds
    character(len=*), parameter :: subname = 'EMIBeTRSimulationInit'

    !grid size
    betr_nlevsoi      = nlevsoi
    betr_nlevsno      = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice
    this%spinup_count = 0

    ! now call the base simulation init to continue initialization
    call this%BeTRInit(bounds, lun, col, pft)

  end subroutine EMIBeTRSimulationInit

!-------------------------------------------------------------------------------
  subroutine BeTRSetFilter(this, maxpft_per_col, nsoilorder)
  !
  !DESCRIPTION
  ! set betr filter, only used for standalone applicaitons
  use betr_varcon                , only : betr_max_soilorder
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type), intent(inout) :: this
  integer, intent(in) :: maxpft_per_col
  integer, optional, intent(in) :: nsoilorder

    integer :: p
    !by default, surface litter layer is off
    this%num_jtops = 1
    allocate(this%jtops(this%num_jtops))
    this%jtops(:) = 1

    this%num_soilc = 1
    allocate(this%filter_soilc(this%num_soilc))
    this%filter_soilc(:) = 1

    this%num_soilp = maxpft_per_col
    allocate(this%filter_soilp(this%num_soilp))
    do p = 1, maxpft_per_col
      this%filter_soilp(p) = p
    enddo

    betr_maxpatch_pft = maxpft_per_col
    if(present(nsoilorder))then
      betr_max_soilorder= nsoilorder
    else
      betr_max_soilorder=1
    endif
  end subroutine BeTRSetFilter

!-------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationReadParams(this, bounds)

  use ncdio_pio               , only :  file_desc_t
  use ncdio_pio               , only : ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen
  use ApplicationsFactory      , only : AppLoadParameters
  use BetrStatusType           , only : betr_status_type
  use decompMod                , only : bounds_type
  use tracer_varcon            , only : bgc_param_file
  use fileutils                , only : getfil
  use spmdMod                  , only : masterproc
  implicit none
  class(emi_betr_simulation_type), intent(inout) :: this

  type(bounds_type), intent(in) :: bounds
  !temporary variables
  type(betr_status_type)   :: bstatus
  type(betr_bounds_type)   :: betr_bounds
  integer :: c
  character(len=256) :: locfn ! local file name
  type(file_desc_t)  :: ncid  ! pio netCDF file id

  !open file for parameter reading
  call getfil (bgc_param_file, locfn, 0)
  if (masterproc) then
    write(iulog,*) 'read betr bgc parameter file '//trim(locfn)
  endif
  call ncd_pio_openfile (ncid, trim(locfn), 0)

  !read in parameters
  call AppLoadParameters(ncid, bstatus)

  call ncd_pio_closefile(ncid)

  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif

  end subroutine EMIBeTRSimulationReadParams
!-------------------------------------------------------------------------------
  subroutine set_activecol(this, col)

  implicit none
    class(emi_betr_simulation_type)              , intent(inout) :: this
    type(column_type)                        , intent(inout) :: col

    logical :: do_debug
    integer :: cc
    do_debug=.false.
    if(do_debug)then
      cc=1
      col%active(:) =.false.
      col%active(cc)=.true.
      this%betr(cc)%tracers%debug=.true.
    endif
  !    col%active(c_act)=.true.

  end subroutine set_activecol
!-------------------------------------------------------------------------------
  subroutine BeTRInit(this, bounds, lun, col, pft)
    !
    ! DESCRIPTION
    ! initialize BeTR
    !
    !!USES
    use betr_constants , only : betr_filename_length
    use betr_varcon    , only : betr_maxpatch_pft
    use landunit_varcon, only : istsoil, istcrop
    implicit none
    !ARGUMENTS
    class(emi_betr_simulation_type)              , intent(inout) :: this
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout) :: col
    type(patch_type)                         , intent(in) :: pft

    !TEMPORARY VARIABLES
    character(len=*), parameter :: subname = 'BeTRInit'
    type(betr_bounds_type) :: betr_bounds
    integer :: c, l, c_l
    logical :: asoibgc
    !print*,'base_filename',trim(base_filename)

    this%hist_record=0
    this%active_soibgc=.false.
    biulog = iulog

    call this%betr_time%Init()

    !allocate memory
    allocate(this%betr(bounds%begc:bounds%endc))
    allocate(this%biophys_forc(bounds%begc:bounds%endc))
    allocate(this%biogeo_flux(bounds%begc:bounds%endc))
    allocate(this%biogeo_state(bounds%begc:bounds%endc))
    allocate(this%bstatus(bounds%begc:bounds%endc))
    allocate(this%betr_col(bounds%begc:bounds%endc))
    allocate(this%betr_pft(bounds%begc:bounds%endc))
    allocate(this%active_col(bounds%begc:bounds%endc))
    allocate(this%bsimstatus)

    call this%bsimstatus%reset()

    !grid horizontal bounds
    call this%BeTRSetBounds(betr_bounds)

    do c = bounds%begc, bounds%endc
      l = col%landunit(c)
      call this%biophys_forc(c)%Init(betr_bounds)

      call this%betr_col(c)%Init(betr_bounds)

      call this%betr_pft(c)%Init(betr_bounds)

      if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
        this%active_col(c) = .true.
      else
        this%active_col(c) = .false.
      endif
    enddo
    call this%BeTRSetcps(bounds, col, pft)

    do c = bounds%begc, bounds%endc
      call this%betr(c)%Init(betr_bounds, this%betr_col(c), &
          this%biophys_forc(c), asoibgc, this%bstatus(c))
      if(c==bounds%begc)this%active_soibgc=asoibgc
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo

    call this%set_activecol(col)

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call this%betr(c)%UpdateParas(betr_bounds, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo

    !allocate spinup factor
    if(this%do_soibgc())then
      allocate(this%scalaravg_col(bounds%begc:bounds%endc));this%scalaravg_col(:) = 0._r8
      allocate(this%dom_scalar_col(bounds%begc:bounds%endc));this%dom_scalar_col(:)=1._r8
      c_l=1
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        call this%betr(c)%set_bgc_spinup(betr_bounds, 1,  betr_nlevtrc_soil, this%biophys_forc(c))
        this%dom_scalar_col(c)=this%biophys_forc(c)%dom_scalar_col(c_l)
      enddo
    endif

    if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())

    do c = bounds%begc, bounds%endc

      call this%biogeo_state(c)%Init(betr_bounds, this%active_soibgc)

      call this%biogeo_flux(c)%Init(betr_bounds, this%active_soibgc)
    enddo
    !identify variables that are used for history output
    c = bounds%begc
    call this%betr(c)%get_hist_size(this%num_hist_state1d, this%num_hist_state2d, &
      this%num_hist_flux1d, this%num_hist_flux2d)

    call this%HistAlloc(bounds)

    call this%betr(c)%get_hist_info(this%num_hist_state1d, this%num_hist_state2d, &
      this%num_hist_flux1d, this%num_hist_flux2d, &
      this%state_hist1d_var(1:this%num_hist_state1d), this%state_hist2d_var(1:this%num_hist_state2d), &
      this%flux_hist1d_var(1:this%num_hist_flux1d), this%flux_hist2d_var(1:this%num_hist_flux2d))

    if(this%bstatus(c)%check_status())then
      call this%bsimstatus%setcol(c)
      call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
      if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())
    endif

    call this%BeTRCreateHistory(bounds, betr_nlevtrc_soil, &
         this%num_hist_state1d, this%num_hist_state2d, &
            this%num_hist_flux1d, this%num_hist_flux2d)

    !identify restart variables
    call this%betr(c)%get_restartvar_size(this%num_rest_state1d, this%num_rest_state2d)

    call this%RestAlloc(bounds)

  end subroutine BeTRInit
  !---------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationRestartAlloc(this, bounds)
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type)              , intent(inout) :: this
  type(bounds_type)                        , intent(in)    :: bounds

  integer :: begc, endc

  begc = bounds%begc; endc=bounds%endc


  allocate(this%rest_states_1d(begc:endc, 1:this%num_rest_state1d))
  this%rest_states_1d(:,:)=spval
  allocate(this%rest_states_2d(begc:endc, 1:betr_nlevtrc_soil, 1:this%num_rest_state2d))
  this%rest_states_2d(:,:,:)=spval

  end subroutine EMIBeTRSimulationRestartAlloc
  !---------------------------------------------------------------------------------

  subroutine EMIBeTRSimulationHistoryAlloc(this, bounds)
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type)              , intent(inout) :: this
  type(bounds_type)                        , intent(in)    :: bounds

  integer :: begc, endc

  begc = bounds%begc; endc=bounds%endc
  !state variables
  allocate(this%state_hist1d_var(this%num_hist_state1d))
  allocate(this%state_hist2d_var(this%num_hist_state2d))

  allocate(this%hist_states_2d(begc:endc, 1:betr_nlevtrc_soil, 1:this%num_hist_state2d))
  allocate(this%hist_states_1d(begc:endc, 1:this%num_hist_state1d))

  !flux variables
  allocate(this%flux_hist1d_var(this%num_hist_flux1d))
  allocate(this%flux_hist2d_var(this%num_hist_flux2d))

  allocate(this%hist_fluxes_2d(begc:endc, 1:betr_nlevtrc_soil, 1:this%num_hist_flux2d))
  allocate(this%hist_fluxes_1d(begc:endc, 1:this%num_hist_flux1d))

  end subroutine EMIBeTRSimulationHistoryAlloc

  !---------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationStepWithoutDrainage(this, bounds, col, pft)
  !DESCRPTION
  !interface for StepWithoutDrainage
  !
  ! USES
    use SoilStateType     , only : soilstate_type
    use WaterStateType    , only : Waterstate_Type
    use TemperatureType   , only : temperature_type
    use ChemStateType     , only : chemstate_type
    use WaterfluxType     , only : waterflux_type
    use atm2lndType       , only : atm2lnd_type
    use SoilHydrologyType , only : soilhydrology_type
    use CNCarbonFluxType  , only : carbonflux_type
    use CanopyStateType   , only : canopystate_type
    use BeTR_TimeMod      , only : betr_time_type
    use pftvarcon         , only : crop
    use clm_varctl     , only : spinup_state
    use betr_ctrl      , only : exit_spinup, enter_spinup,betr_spinup_state
    use clm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
  !ARGUMENTS
    class(emi_betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds ! bounds
    type(column_type)           , intent(in)    :: col ! column type
    type(patch_type)            , intent(in)    :: pft

    !TEMPORARY VARIABLES
    type(betr_bounds_type)     :: betr_bounds
    integer :: c, c_l, begc_l, endc_l
    integer :: year, mon, day, sec

    call get_curr_date(year, mon, day, sec)
    c_l=1
    if(this%do_soibgc())then
      if(spinup_state==1)then
        do c = bounds%begc, bounds%endc
          this%biophys_forc(c)%dom_scalar_col(c_l)=this%dom_scalar_col(c)
        enddo
      else
        betr_spinup_state=0
      endif
    endif
    call this%bsimstatus%reset()

    !pass necessary data for correct subroutine call
    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col, pft)

    c_l = 1; begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      !this%betr(c)%tracers%debug=col%debug_flag(c)
      call this%biophys_forc(c)%frac_normalize(this%betr_pft(c)%npfts, 1, betr_nlevtrc_soil)

      if(this%betr(c)%tracers%debug)call this%betr(c)%debug_info(betr_bounds, this%betr_col(c), &
         this%num_soilc, this%filter_soilc, 'bef w/o drain',this%bstatus(c))

      call this%betr(c)%step_without_drainage(this%betr_time, betr_bounds, this%betr_col(c), &
         this%betr_pft(c), this%num_soilc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc(c), this%biogeo_flux(c), this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err(),c)
        exit
      endif

      if(this%betr(c)%tracers%debug)call this%betr(c)%debug_info(betr_bounds, this%betr_col(c),&
          this%num_soilc, this%filter_soilc, 'aft w/o drain',this%bstatus(c))
    enddo
!    print*,'out without drainage'
    if(this%bsimstatus%check_status())then
      call endrun(msg=this%bsimstatus%print_msg())
    endif
    if(betr_spinup_state>0)then
      !the following needs double check for whether to keep or remove it.
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        this%dom_scalar_col(c) = this%biophys_forc(c)%dom_scalar_col(c_l)
      enddo
    endif

  end subroutine EMIBeTRSimulationStepWithoutDrainage
  !---------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationDiagnoseLnd2atm(this, bounds,  col, lnd2atm_vars)
   !
   !DESCRIPTION
   !interface for using diagnose land fluxes to atm and river copmonents
   !
   !USES
    use MathfuncMod   , only : safe_div
    use lnd2atmType    , only : lnd2atm_type
    implicit none
    !ARGUMENTS
    class(emi_betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_vars

    ! remove compiler warnings about unused dummy args
    if (this%num_soilc > 0) continue
    if (bounds%begc > 0)    continue
    if (size(col%z) > 0)    continue

  end subroutine EMIBeTRSimulationDiagnoseLnd2atm
  !---------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationStepWithDrainage(this, bounds,  col)
   !
   !DESCRIPTION
   !interface for using StepWithDrainage
   !
   !USES
    use MathfuncMod   , only : safe_div
    use clm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    !ARGUMENTS
    class(emi_betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type
    !temporary variables
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer                :: c, c_l, begc_l, endc_l

    call this%bsimstatus%reset()

    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col)

    c_l = 1; begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle

      call this%betr(c)%step_with_drainage(betr_bounds,      &
         this%betr_col(c),this%num_soilc, this%filter_soilc, this%jtops, &
         this%biogeo_flux(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif

    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=this%bsimstatus%print_msg())

  end subroutine EMIBeTRSimulationStepWithDrainage

  !------------------------------------------------------------------------

  subroutine EMIBeTRSimulationCalcDewSubFlux(this,  &
       bounds, col, num_hydrologyc, filter_soilc_hydrologyc)
   !DESCRIPTION
    ! Calculate tracer flux from dew or/and sublimation
    !External interface called by ALM

    use WaterfluxType   , only : waterflux_type
    use WaterstateType  , only : waterstate_type
    use clm_varcon      , only : denh2o,spval
    use landunit_varcon , only : istsoil, istcrop
    use betr_decompMod  , only : betr_bounds_type
    implicit none
    !ARGUMENTS
    class(emi_betr_simulation_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    type(column_type)               , intent(in)    :: col ! column type
    integer                         , intent(in)    :: num_hydrologyc ! number of column soil points in column filter_soilc
    integer                         , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points

    !temporary variables
    type(betr_bounds_type)     :: betr_bounds
    integer :: fc, c

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col)
    do fc = 1, num_hydrologyc
      c = filter_soilc_hydrologyc(fc)
      if(.not. this%active_col(c))cycle
      call this%betr(c)%calc_dew_sub_flux(this%betr_time,           &
         betr_bounds, this%betr_col(c), this%num_soilc, this%filter_soilc, &
        this%biophys_forc(c), this%betr(c)%tracers, this%betr(c)%tracerfluxes, this%betr(c)%tracerstates)
    enddo
  end subroutine EMIBeTRSimulationCalcDewSubFlux
  !---------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationBeginMassBalanceCheck(this, bounds)
    !DESCRIPTION
    !stage tracer mass balance check
    !
    !USES
    use TracerBalanceMod, only : begin_betr_tracer_massbalance
    implicit none
    !ARGUMENTS
    class(emi_betr_simulation_type), intent(inout)   :: this
    type(bounds_type), intent(in) :: bounds

    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj
    integer :: c

    !set lbj and ubj
    call this%BeTRSetBounds(betr_bounds)

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call begin_betr_tracer_massbalance(betr_bounds,             &
         this%betr_col(c), this%num_soilc, this%filter_soilc,     &
         this%num_soilp, this%filter_soilp, &
         this%betr(c)%tracers, this%betr(c)%tracerstates,         &
         this%betr(c)%tracerfluxes, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=trim(this%bsimstatus%print_msg()))
  end  subroutine EMIBeTRSimulationBeginMassBalanceCheck
  !---------------------------------------------------------------------------------

  subroutine EMIBeTRSimulationMassBalanceCheck(this, bounds)
   !DESCRIPTION
   ! do tracer mass balance check
   !
   !USES
    use TracerBalanceMod , only : betr_tracer_massbalance_check
    use BeTR_TimeMod     , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(emi_betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer :: c
    logical :: ldebug
    !set lbj and ubj
    call this%BeTRSetBounds(betr_bounds)

   do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      ldebug=(c==68) .and. .false.
      call betr_tracer_massbalance_check(this%betr_time, betr_bounds,           &
         this%betr_col(c), this%num_soilc, this%filter_soilc,           &
         this%betr(c)%tracers, this%betr(c)%tracerstates,                &
         this%betr(c)%tracerfluxes, this%bstatus(c), ldebug)
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
   enddo
   if(this%bsimstatus%check_status()) then
      print*,this%bsimstatus%cindex
      print*,trim(this%bsimstatus%print_msg())
      call endrun(msg=trim(this%bsimstatus%print_msg()))
   endif
  end subroutine EMIBeTRSimulationMassBalanceCheck

  !---------------------------------------------------------------------------------

  subroutine EMIBeTRSimulationHistRetrieval(this, bounds, numf, filter)
  use tracer_varcon  , only : betr_nlevtrc_soil
  implicit none
  !ARGUMENTS
   class(emi_betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds
   integer, intent(in) :: numf
   integer, intent(in) :: filter(:)

  call this%BeTRRetrieveHistoryState(bounds, numf, filter)

  call this%BeTRRetrieveHistoryFlux(bounds, numf, filter)

  end subroutine EMIBeTRSimulationHistRetrieval

  !---------------------------------------------------------------------------------
  subroutine EMIBeTRSimulationConsistencyCheck(this, &
     bounds, ubj, num_soilc, filter_soilc, waterstate_vars)
  ! DESCRIPTION
  ! Do consistency check, can be overwritten for varies purpose
  !
  ! USES
    use WaterStateType, only : Waterstate_Type
  implicit none
   !ARGUMENTS
    class(emi_betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc ! number of columns in column filter_soilc
    integer                     , intent(in)    :: filter_soilc(:) ! column filter_soilc
    integer                     , intent(in)    :: ubj
    type(Waterstate_Type)       , intent(in)    :: waterstate_vars ! water state variables

    ! remove compiler warnings
    if (this%num_soilc > 0)                         continue
    if (bounds%begc > 0)                            continue
    if (ubj > 0)                                    continue
    if (num_soilc > 0)                              continue
    if (size(filter_soilc) > 0)                     continue
!    if (associated(col_ws%h2osoi_liq)) continue
    if (associated(waterstate_vars%h2osoi_liq_col)) continue
  end subroutine EMIBeTRSimulationConsistencyCheck

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationRetrieveBiogeoFlux(this, bounds, lbj,ubj, carbonflux_vars,  &
    waterflux_vars)
  ! DESCRIPTIONS
  ! update and return fluxes, this eventually will be expanded to
  ! include other fluxes
  ! USES
    use WaterfluxType    , only : waterflux_type
    use CNCarbonFluxType , only : carbonflux_type
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type) , intent(inout)           :: this
  type(bounds_type)           , intent(in)              :: bounds
  integer                     , intent(in)              :: lbj, ubj
  type(carbonflux_type)       , optional, intent(inout) :: carbonflux_vars
  type(waterflux_type)        , optional, intent(inout) :: waterflux_vars

  integer :: begp, begc, endp, endc
  integer :: p, c, cc
  cc = 1
  if(present(carbonflux_vars))then
    !do nothing
  endif
  if(present(waterflux_vars))then
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
!      col_wf%qflx_infl(c)            = this%biogeo_flux(c)%qflx_infl_col(cc)
!      col_wf%qflx_adv(c,lbj-1:ubj)   = this%biogeo_flux(c)%qflx_adv_col(cc,lbj-1:ubj)
!      col_wf%qflx_totdrain(c)        = this%biogeo_flux(c)%qflx_totdrain_col(cc)
!      col_wf%qflx_gross_evap_soil(c) = this%biogeo_flux(c)%qflx_gross_evap_soil_col(cc)
!      col_wf%qflx_gross_infl_soil(c) = this%biogeo_flux(c)%qflx_gross_infl_soil_col(cc)
!      col_wf%qflx_drain_vr(c,1:ubj)  = this%biogeo_flux(c)%qflx_drain_vr_col(cc,1:ubj)
      waterflux_vars%qflx_infl_col(c)            = this%biogeo_flux(c)%qflx_infl_col(cc)
      waterflux_vars%qflx_adv_col(c,lbj-1:ubj)   = this%biogeo_flux(c)%qflx_adv_col(cc,lbj-1:ubj)
      waterflux_vars%qflx_totdrain_col(c)        = this%biogeo_flux(c)%qflx_totdrain_col(cc)
      waterflux_vars%qflx_gross_evap_soil_col(c) = this%biogeo_flux(c)%qflx_gross_evap_soil_col(cc)
      waterflux_vars%qflx_gross_infl_soil_col(c) = this%biogeo_flux(c)%qflx_gross_infl_soil_col(cc)
      waterflux_vars%qflx_drain_vr_col(c,1:ubj)  = this%biogeo_flux(c)%qflx_drain_vr_col(cc,1:ubj)
    enddo
  endif

  end subroutine EMIBeTRSimulationRetrieveBiogeoFlux

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationPreDiagSoilColWaterFlux(this, num_nolakec, filter_nolakec)
  !DESCRIPTION
  !prepare for water flux diagnosis. it is called before diagnosing the advective fluxes and applying
  !freeze-thaw tracer partition.
  !
  !USES

  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type) , intent(inout) :: this
   integer                    , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                    , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points

  !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c


   call this%BeTRSetBounds(betr_bounds)

   do fc= 1, num_nolakec
     c = filter_nolakec(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%pre_diagnose_soilcol_water_flux(betr_bounds, this%num_soilc, &
       this%filter_soilc, this%biophys_forc(c))
   enddo
  end subroutine EMIBeTRSimulationPreDiagSoilColWaterFlux

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationDiagnoseDtracerFreezeThaw(this, bounds, num_nolakec, filter_nolakec, col, lun)
  !
  ! DESCRIPTION
  ! aqueous tracer partition based on freeze-thaw
  !
  ! USES
  use WaterStateType        , only : waterstate_type
  implicit none
  !
  ! Arguments
  class(emi_betr_simulation_type), intent(inout)   :: this
  type(bounds_type)     , intent(in) :: bounds
  integer               , intent(in) :: num_nolakec                        ! number of column non-lake points in column filter
  integer               , intent(in) :: filter_nolakec(:)                  ! column filter for non-lake points
!  type(waterstate_type), intent(in) :: waterstate_vars
  type(column_type)     , intent(in) :: col                                ! column type
  type(landunit_type)   , intent(in)  :: lun

  !temporary variables
  type(betr_bounds_type)     :: betr_bounds
  integer :: fc, c

  call this%BeTRSetBounds(betr_bounds)

  call this%BeTRSetcps(bounds, col)

  do fc = 1, num_nolakec
    c = filter_nolakec(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%diagnose_dtracer_freeze_thaw(betr_bounds, this%num_soilc, this%filter_soilc,  &
      this%biophys_forc(c))
  enddo
  end subroutine EMIBeTRSimulationDiagnoseDtracerFreezeThaw

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationDiagAdvWaterFlux(this, num_hydrologyc, &
    filter_hydrologyc)

  !DESCRIPTION
  ! diagnose water fluxes for tracer advection
  !
  ! USES
  !
  implicit none
  !ARGUMENTS
   class(emi_betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c, j

   call this%BeTRSetBounds(betr_bounds)
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%diagnose_advect_water_flux(this%betr_time,              &
       betr_bounds, this%num_soilc, this%filter_soilc,                         &
       this%biophys_forc(c), this%biogeo_flux(c))
   enddo

  end subroutine EMIBeTRSimulationDiagAdvWaterFlux
  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationDiagDrainWaterFlux(this, num_hydrologyc, filter_hydrologyc)
  !DESCRIPTION
  ! diagnose water fluxes due to subsurface drainage
  !
  ! USES
  !
    use WaterfluxType     , only : waterflux_type
    use WaterStateType    , only : Waterstate_Type
    use SoilHydrologyType , only : soilhydrology_type
  implicit none
  !ARGUMENTS
   class(emi_betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%diagnose_drainage_water_flux(this%betr_time, &
       betr_bounds, this%num_soilc, this%filter_soilc,      &
       this%biophys_forc(c), this%biogeo_flux(c))
  enddo
  end subroutine EMIBeTRSimulationDiagDrainWaterFlux
  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationBeginTracerSnowLayerAdjst(this,  num_snowc, filter_snowc)
  !DESCRIPTION
  !prepare for tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(emi_betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%Enter_tracer_LayerAdjustment(betr_bounds, this%betr_col(c), &
       this%num_soilc, this%filter_soilc)
   enddo
  end subroutine EMIBeTRSimulationBeginTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationEndTracerSnowLayerAdjst(this, num_snowc, filter_snowc)
  !DESCRIPTION
  !wrap up tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(emi_betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%Exit_tracer_LayerAdjustment(betr_bounds, this%betr_col(c), &
       this%num_soilc, this%filter_soilc)
   enddo

  end subroutine EMIBeTRSimulationEndTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationDvideSnowLayers(this, bounds, num_snowc, filter_snowc, divide_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use clm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(emi_betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: divide_matrix(bounds%begc:bounds%endc , 1:nlevsno , 1:nlevsno )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%tracer_DivideSnowLayers(betr_bounds, this%betr_col(c),this%num_soilc, &
       this%filter_soilc, divide_matrix(c:c,:,:), this%bstatus(c))
     if(this%bstatus(c)%check_status())then
       call this%bsimstatus%setcol(c)
       call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
       exit
     endif
   enddo
  if(this%bsimstatus%check_status()) &
    call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine EMIBeTRSimulationDvideSnowLayers

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationCombineSnowLayers(this, bounds, num_snowc, filter_snowc, combine_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use clm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(emi_betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: combine_matrix(bounds%begc:bounds%endc,-nlevsno+1:1 ,-nlevsno+1:1 )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%tracer_CombineSnowLayers(betr_bounds, this%betr_col(c),this%num_soilc,&
       this%filter_soilc, combine_matrix(c:c,:,:),this%bstatus(c))
     if(this%bstatus(c)%check_status())then
       call this%bsimstatus%setcol(c)
       call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
       exit
     endif
   enddo
   if(this%bsimstatus%check_status()) &
    call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine EMIBeTRSimulationCombineSnowLayers

  !------------------------------------------------------------------------
  subroutine hist_create_fluxes(this, bounds, betr_nlevtrc_soil, num_flux1d, num_flux2d, ncid)
  !
  !DESCRIPTION
  !create history file for betr fluxes
  !
  use histFileMod         , only: hist_addfld1d, hist_addfld2d
  use bhistFileMod        , only : hist_def_fld2d , hist_def_fld1d
  use bncdio_pio           , only : file_desc_t, ncd_float
  implicit none
  class(emi_betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: betr_nlevtrc_soil
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d
  type(file_desc_t) ,   optional,  intent(inout)   :: ncid
  !local variables
  integer :: jj, begc, endc

  character(len=*), parameter :: subname = 'hist_create_fluxes'

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  begc=bounds%begc; endc=bounds%endc

  do jj = 1, num_flux2d
    this%hist_fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
      data2dptr => this%hist_fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj)
      call hist_addfld2d (fname=this%flux_hist2d_var(jj)%varname, &
           units=this%flux_hist1d_var(jj)%units, type2d="levtrc",  &
           avgflag=this%flux_hist2d_var(jj)%avg_flag, &
           long_name=this%flux_hist2d_var(jj)%long_name,  ptr_col=data2dptr, &
           default=this%flux_hist2d_var(jj)%use_default)
  enddo

  do jj = 1, num_flux1d
      this%hist_fluxes_1d(begc:endc,jj) = spval
      data1dptr => this%hist_fluxes_1d(begc:endc, jj)
      call hist_addfld1d (fname=this%flux_hist1d_var(jj)%varname, &
        units=this%flux_hist1d_var(jj)%units,  &
        avgflag=this%flux_hist1d_var(jj)%avg_flag, &
        long_name=this%flux_hist1d_var(jj)%long_name, &
        ptr_col=data1dptr, default=this%flux_hist1d_var(jj)%use_default)
  enddo

  end subroutine hist_create_fluxes

  !------------------------------------------------------------------------
  subroutine hist_create_states(this, bounds, betr_nlevtrc_soil, num_state1d, num_state2d, ncid)
  !
  !create history file for betr states variables
  use histFileMod   , only: hist_addfld1d, hist_addfld2d
  use bhistFileMod        , only : hist_def_fld2d , hist_def_fld1d
  use bncdio_pio , only : file_desc_t, ncd_float
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  type(file_desc_t) ,   optional,  intent(inout)   :: ncid
  !local variables
  integer :: begc, endc
  integer :: jj

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  character(len=*), parameter :: subname = 'hist_create_states'

  begc = bounds%begc; endc = bounds%endc

  do jj = 1, num_state2d
    !read namelist

      this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
      data2dptr => this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj)
      call hist_addfld2d (fname=this%state_hist2d_var(jj)%varname, &
           units=this%state_hist2d_var(jj)%units, type2d="levtrc",  &
           avgflag=this%state_hist2d_var(jj)%avg_flag, &
           long_name=this%state_hist2d_var(jj)%long_name, &
           ptr_col=data2dptr, default=this%state_hist2d_var(jj)%use_default)
  enddo

  do jj = 1, num_state1d

      this%hist_states_1d(begc:endc,jj) = spval
      data1dptr => this%hist_states_1d(begc:endc,jj)
      call hist_addfld1d (fname=this%state_hist1d_var(jj)%varname, &
          units=this%state_hist1d_var(jj)%units,      &
          avgflag=this%state_hist1d_var(jj)%avg_flag, &
          long_name=this%state_hist1d_var(jj)%long_name, &
          ptr_col=data1dptr, default=this%state_hist1d_var(jj)%use_default)
  enddo

  end subroutine hist_create_states

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationCreateHistory(this, bounds, betr_nlevtrc_soil,&
     num_state1d, num_state2d, num_flux1d, num_flux2d)
  !
  !links the variable to output
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           , intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d

  call this%hist_create_states(bounds, betr_nlevtrc_soil, num_state1d, num_state2d)

  call this%hist_create_fluxes(bounds, betr_nlevtrc_soil, num_flux1d, num_flux2d)

  end subroutine EMIBeTRSimulationCreateHistory
  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationRetrieveHistoryState(this, bounds, numf, filter)
  use tracer_varcon  , only : betr_nlevtrc_soil
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !TEMPORARY VARIABLES
  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)

  do fc = 1, numf
    c = filter(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%HistRetrieveState(betr_bounds, 1, betr_nlevtrc_soil, &
       this%num_hist_state1d, this%num_hist_state2d,&
       this%hist_states_1d(c:c,:), this%hist_states_2d(c:c,1:betr_nlevtrc_soil,:))
  enddo

  end subroutine EMIBeTRSimulationRetrieveHistoryState
  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationRetrieveHistoryFlux(this, bounds, numf, filter)
  use tracer_varcon  , only :  betr_nlevtrc_soil
  implicit none
  !ARGUMENTS
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !TEMPORARY VARIABLES
  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)
  this%hist_fluxes_1d(:,:)=spval
  this%hist_fluxes_2d(:,:,:)=spval
  do fc = 1, numf
    c = filter(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%HistRetrieveFlux(betr_bounds, 1, betr_nlevtrc_soil, &
       this%num_hist_flux1d,this%num_hist_flux2d, &
       this%hist_fluxes_1d(c:c,:),this%hist_fluxes_2d(c:c,1:betr_nlevtrc_soil,:))
  enddo

  end subroutine EMIBeTRSimulationRetrieveHistoryFlux

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationRestart(this, bounds, ncid, flag)
  !DESCRIPTION
  !create or read restart file
  use restUtilMod    , only : restartvar
  use ncdio_pio      , only : file_desc_t,ncd_double, ncd_int
  use clm_varctl     , only : spinup_state
  use clm_time_manager, only : get_nstep
  use betr_ctrl      , only : exit_spinup, enter_spinup,betr_spinup_state
  use tracer_varcon  , only : reaction_method
  implicit none
  ! !ARGUMENTS:
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(bounds_type)    , intent(in)    :: bounds
  type(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
  character(len=*)     , intent(in)    :: flag ! 'read' or 'write'

  !local variables
  integer :: nrest_1d, nrest_2d
  integer :: c, jj, fc
  character(len=255), allocatable :: rest_varname_1d(:)
  character(len=255), allocatable :: rest_varname_2d(:)
  logical :: readvar      ! determine if variable is on initial file
  real(r8), pointer :: ptr1d(:)
  real(r8), pointer :: ptr2d(:,:)
  type(betr_bounds_type)     :: betr_bounds
  integer :: recordDimID
  integer  :: idata
  integer  :: restart_file_spinup_state
  integer  :: c_l

  c_l = 1
  c = bounds%begc
  restart_file_spinup_state =0
  allocate(rest_varname_1d(this%num_rest_state1d)); rest_varname_1d=''
  allocate(rest_varname_2d(this%num_rest_state2d)); rest_varname_2d=''

  c = bounds%begc
  call this%betr(c)%get_restartvar_info(this%num_rest_state1d, &
    this%num_rest_state2d,rest_varname_1d, rest_varname_2d)

  call this%BeTRSetBounds(betr_bounds)
  if(trim(flag)=='write')then
    if(this%do_soibgc())then
      idata = spinup_state
      do c = bounds%begc, bounds%endc
        this%scalaravg_col(c) = this%biophys_forc(c)%scalaravg_col(c_l)
        this%dom_scalar_col(c) = this%biophys_forc(c)%dom_scalar_col(c_l)
      enddo
    endif
    do c = bounds%begc, bounds%endc
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d,this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo
  endif

  if(this%do_soibgc())then
    call restartvar(ncid=ncid, flag=flag, varname='spinscalar', xtype=ncd_double, &
         dim1name='column', long_name='', units='', &
         interpinic_flag = 'interp', readvar=readvar, data=this%scalaravg_col)

     call restartvar(ncid=ncid, flag=flag, varname='domspinscalar', xtype=ncd_double, &
         dim1name='column', long_name='', units='', &
         interpinic_flag = 'interp', readvar=readvar, data=this%dom_scalar_col)

    call restartvar(ncid=ncid, flag=flag, varname='betr_spinup_state', xtype=ncd_int,  &
           long_name='Spinup state of betr model that wrote this restart file: ' &
           // ' 0,1,2=not ready for spinup scalar, 3 = apply spinup scalar', units='', &
           interpinic_flag='copy', readvar=readvar,  data=betr_spinup_state)

    call restartvar(ncid=ncid, flag=flag, varname='spinup_count', xtype=ncd_int,  &
           long_name='Spinup count of the model that wrote this restart file: ' &
           // ' 0 <=2 skip mass bal check, 3 = do mass bal check', units='', &
           interpinic_flag='copy', readvar=readvar,  data=this%spinup_count)

    if(trim(flag)=='read')then
      call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
             long_name='Spinup state of the model that wrote this restart file: ' &
             // ' 0 = normal model mode, 1 = AD spinup', units='', &
             interpinic_flag='copy', readvar=readvar,  data=idata)
      if (readvar) then
        restart_file_spinup_state = idata
      else
        restart_file_spinup_state = spinup_state
      endif

    endif

  endif
  do jj = 1, this%num_rest_state1d
    ptr1d => this%rest_states_1d(:, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_1d(jj)), &
       xtype=ncd_double,  dim1name='column', long_name='',  units='', &
       interpinic_flag='interp' , readvar=readvar, data=ptr1d)
  enddo
  do jj = 1, this%num_rest_state2d
    ptr2d => this%rest_states_2d(:, :, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_2d(jj)), xtype=ncd_double,  &
      dim1name='column',dim2name='levtrc', switchdim=.true., &
      long_name='',  units='', interpinic_flag='interp',readvar=readvar, data=ptr2d)
  enddo
  if(trim(flag)=='read')then

    !assign initial conditions
    do c = bounds%begc, bounds%endc
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d,this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo

    if(this%do_soibgc())then
      exit_spinup = (spinup_state == 0 .and. restart_file_spinup_state == 1 )
      enter_spinup = (spinup_state == 1 .and. restart_file_spinup_state == 0)
      if(get_nstep() >= 2)then
        exit_spinup = .false.; enter_spinup=.false.
      endif
      do c = bounds%begc, bounds%endc
        this%biophys_forc(c)%scalaravg_col(c_l) = max(this%scalaravg_col(c),0.01_r8)
        this%biophys_forc(c)%dom_scalar_col(c_l)= this%dom_scalar_col(c)
      enddo
    endif
  endif

  deallocate(rest_varname_1d)
  deallocate(rest_varname_2d)
  end subroutine EMIBeTRSimulationRestart

  !------------------------------------------------------------------------

  subroutine EMIBeTRSimulationSetSpinup(this, bounds)
  !
  ! set spinup for betr bgc runs
  use betr_ctrl      , only : exit_spinup, enter_spinup,betr_spinup_state
  use ApplicationsFactory, only : AppSetSpinup
  implicit none
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(bounds_type), intent(in) :: bounds
  type(betr_bounds_type)     :: betr_bounds
  integer :: c

  if(exit_spinup .or. enter_spinup)then
     call AppSetSpinup()
     call this%BeTRSetBounds(betr_bounds)
     do c = bounds%begc, bounds%endc
       if(.not. this%active_col(c))cycle
       call this%betr(c)%set_bgc_spinup(betr_bounds, 1,  betr_nlevtrc_soil, this%biophys_forc(c))
     enddo
  endif
  if(exit_spinup)betr_spinup_state=0

  end subroutine EMIBeTRSimulationSetSpinup

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationSetcps(this, bounds, col, pft)
  !
  !DESCRIPTION
  ! set up columns
  !USES
  use decompMod             , only : bounds_type
  use pftvarcon             , only : noveg, crop
  use tracer_varcon         , only : betr_nlevsoi
  !ARGUMENTS
  implicit none
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(bounds_type), intent(in) :: bounds
  type(column_type), intent(in) :: col
  type(patch_type), optional, intent(in) :: pft
  integer :: c, p, pi, pp
  integer :: c_l

  c_l=1
  do c = bounds%begc, bounds%endc
    this%betr_col(c)%snl(c_l) = col%snl(c)
    if(col%snl(c)<0)this%betr_col(c)%dz_snow(c_l, col%snl(c)+1:0) = col%dz(c,col%snl(c)+1:0)
    this%betr_col(c)%zi(c_l,0:betr_nlevsoi)= col%zi(c,0:betr_nlevsoi)
    this%betr_col(c)%dz(c_l,1:betr_nlevsoi)= col%dz(c,1:betr_nlevsoi)
    this%betr_col(c)%z(c_l,1:betr_nlevsoi)= col%z(c,1:betr_nlevsoi)

    this%betr_col(c)%pfti(1)= col%pfti(c)
    this%betr_col(c)%pftf(1)= col%pftf(c)
    this%betr_col(c)%npfts(1)= col%npfts(c)

    if(present(pft))then
      this%betr_pft(c)%column(:)=1
      this%betr_pft(c)%npfts = 0
      pp = 0
      do pi = 1, betr_maxpatch_pft
        if (pi <= col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p) .and. (pft%itype(p) /= noveg)) then
            pp = pp + 1
            this%betr_pft(c)%wtcol(pp) = pft%wtcol(p)
            this%betr_pft(c)%itype(pp) = pft%itype(p)
            this%betr_pft(c)%crop(pp) = crop(pi)         !the crop looks weird here, jyt
          endif
        endif
      enddo
      this%betr_pft(c)%npfts = pp
    endif
  enddo


  end subroutine EMIBeTRSimulationSetcps

  !------------------------------------------------------------------------
  subroutine EMIBeTRSimulationSetBounds(this, betr_bounds)
  !
  !DESCRIPTION
  !set betr_bounds
  !
  use betr_varcon    , only : betr_maxpatch_pft
  use tracer_varcon  , only : betr_nlevsoi,betr_nlevsno
  implicit none
  class(emi_betr_simulation_type) , intent(inout) :: this
  type(betr_bounds_type), intent(out)  :: betr_bounds
  !the following will be adpated to simulate lake and wetlands, by allowing
  !lbj to be negative
  betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
  betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
  betr_bounds%begc = 1; betr_bounds%endc = 1
  betr_bounds%begl = 1; betr_bounds%endl = 1
  betr_bounds%begg = 1; betr_bounds%endg = 1
  betr_bounds%nlevsno=betr_nlevsno
  end subroutine EMIBeTRSimulationSetBounds

  !------------------------------------------------------------------------
  function do_soibgc(this)result(yesno)

  implicit none
  class(emi_betr_simulation_type) , intent(inout) :: this

  logical :: yesno
  yesno = this%active_soibgc
  return
  end function do_soibgc

end module EMIBeTRSimulation
