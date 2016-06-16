module BeTRSimulation
  !
  ! !DESCRIPTION:
  !  BeTR simulation base class.
  !
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
  use ColumnType     , only : col
  use PatchType      , only : pft
  ! !USES:
  use BetrType                 , only : betr_type, create_betr_type
  use betr_ctrl                , only : max_betr_hist_type, betr_offline
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
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: betr_simulation_type
     type(betr_type)                     , public, pointer  :: betr(:)
     type(betr_biogeophys_input_type)    , public, pointer  :: biophys_forc(:)
     type(betr_biogeo_state_type)        , public, pointer  :: biogeo_state(:)
     type(betr_biogeo_flux_type)         , public, pointer  :: biogeo_flux(:)
     type(betr_status_type)              , public, pointer  :: bstatus(:)
     type(betr_status_sim_type)          , public, pointer  :: bsimstatus
     character(len=betr_filename_length) , private :: base_filename
     character(len=betr_filename_length) , private :: hist_filename

     type(betr_regression_type), private          :: regression

     integer, public                              :: num_soilp
     integer, public, allocatable                 :: filter_soilp(:)
     integer, public                              :: num_jtops
     integer, public, allocatable                 :: jtops(:)
     integer, public                              :: num_soilc
     integer, public, allocatable                 :: filter_soilc(:)

     real(r8), pointer :: states_2d(:,:,:) !(col, lay,var)
     real(r8), pointer :: states_1d(:,:)   !(col, var)
     real(r8), pointer :: fluxes_2d(:,:,:) !(col, lay,var)
     real(r8), pointer :: fluxes_1d(:,:)   !(col, var)
     character(len=255) :: nmlist_hist1d_state_buffer(max_betr_hist_type)
     character(len=255) :: nmlist_hist2d_state_buffer(max_betr_hist_type)
     character(len=255) :: nmlist_hist1d_flux_buffer(max_betr_hist_type)
     character(len=255) :: nmlist_hist2d_flux_buffer(max_betr_hist_type)

     integer :: num_hist_state1d
     integer :: num_hist_state2d
     integer :: num_hist_flux1d
     integer :: num_hist_flux2d

     ! FIXME(bja, 201603) most of these types should be private!

     ! NOTE(bja, 201603) BeTR types only, no LSM specific types here!

   contains
     procedure, public :: BeTRInit
     procedure, public :: BeTRSetFilter
     procedure, public :: Init                    => BeTRSimulationInit
     procedure, public :: RestartInit             => BeTRSimulationRestartInit
     procedure, public :: ConsistencyCheck        => BeTRSimulationConsistencyCheck
     procedure, public :: PreDiagSoilColWaterFlux => BeTRSimulationPreDiagSoilColWaterFlux
     procedure, public :: DiagAdvWaterFlux        => BeTRSimulationDiagAdvWaterFlux
     procedure, public :: DiagDrainWaterFlux      => BeTRSimulationDiagDrainWaterFlux
     procedure, public :: BeginSnowLayerAdjst     => BeTRSimulationBeginTracerSnowLayerAdjst
     procedure, public :: EndSnowLayerAdjst       => BeTRSimulationEndTracerSnowLayerAdjst
     procedure, public :: CombineSnowLayers       => BeTRSimulationCombineSnowLayers
     procedure, public :: DvideSnowLayers         => BeTRSimulationDvideSnowLayers
     procedure, public :: StepWithoutDrainage     => BeTRSimulationStepWithoutDrainage
     procedure, public :: StepWithDrainage        => BeTRSimulationStepWithDrainage
     procedure, public :: BeginMassBalanceCheck   => BeTRSimulationBeginMassBalanceCheck
     procedure, public :: MassBalanceCheck        => BeTRSimulationMassBalanceCheck
     procedure, public :: BeTRSetBiophysForcing   => BeTRSimulationSetBiophysForcing
     procedure, public :: RetrieveBiogeoFlux      => BeTRSimulationRetrieveBiogeoFlux
     procedure, public :: CreateHistory           => hist_htapes_create
     procedure, public :: WriteHistory            => hist_write
     procedure, public :: WriteRegressionOutput
     !the following are used to interact with lsm
     procedure, public :: BeTRRestart             => BeTRSimulationRestart
     procedure, public :: BeTRCreateHistory       => BeTRSimulationCreateHistory
     procedure, public :: BeTRRetrieveHistory     => BeTRSimulationRetrieveHistory
     procedure, private:: hist_create_states
     procedure, private:: hist_create_fluxes
  end type betr_simulation_type

  public :: BeTRSimulationInit

contains

  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationInit(this, base_filename, namelist_buffer, &
       bounds, waterstate)
    !
    ! DESCRIPTIONS
    ! Dummy routine for inheritance purposes. don't use.
    !
    !USES
    use WaterstateType , only : waterstate_type
    use betr_constants , only : betr_namelist_buffer_size, betr_filename_length
    implicit none

    class(betr_simulation_type)              , intent(inout) :: this
    character(len=betr_filename_length)      , intent(in)    :: base_filename
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer

    type(bounds_type)                        , intent(in)    :: bounds
    type(waterstate_type)                    , intent(inout) :: waterstate

    character(len=*), parameter :: subname = 'BeTRSimulationInit'

    call endrun(msg="ERROR "//subname//" unimplemented. "//errmsg(mod_filename, __LINE__))

    if (this%num_soilc > 0)                  continue
    if (bounds%begc > 0)                     continue
    if (size(waterstate%h2osoi_liq_col) > 0) continue
    if (len(base_filename) > 0)              continue
    if (len(namelist_buffer) > 0)            continue

  end subroutine BeTRSimulationInit

!-------------------------------------------------------------------------------
  subroutine BeTRSetFilter(this)
  !
  !DESCRIPTION
  ! set betr filter, only used for standalone applicaitons

  implicit none
  !ARGUMENTS
  class(betr_simulation_type), intent(inout) :: this

    this%num_jtops = 1
    allocate(this%jtops(this%num_jtops))
    this%jtops(:) = 1

    this%num_soilc = 1
    allocate(this%filter_soilc(this%num_soilc))
    this%filter_soilc(:) = 1

    this%num_soilp = 1
    allocate(this%filter_soilp(this%num_soilp))
    this%filter_soilp(:) = 1
  end subroutine BeTRSetFilter
!-------------------------------------------------------------------------------

  subroutine BeTRInit(this, base_filename, namelist_buffer, &
       bounds, waterstate)
    !
    ! DESCRIPTION
    ! initialize BeTR
    !
    !!USES
    use WaterStateType , only : waterstate_type
    use betr_constants , only : betr_namelist_buffer_size
    use betr_constants , only : betr_filename_length
    use betr_varcon    , only : betr_maxpatch_pft
    implicit none
    !ARGUMENTS
    class(betr_simulation_type)              , intent(inout) :: this
    character(len=betr_filename_length)      , intent(in)    :: base_filename
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    type(bounds_type)                        , intent(in)    :: bounds
    type(waterstate_type)                    , intent(in)    :: waterstate

    !TEMPORARY VARIABLES
    character(len=*), parameter :: subname = 'BeTRInit'
    type(betr_bounds_type) :: betr_bounds
    integer :: c
    !print*,'base_filename',trim(base_filename)

    this%base_filename = base_filename

    !allocate memory
    allocate(this%betr(bounds%begc:bounds%endc), source=create_betr_type())
    allocate(this%biophys_forc(bounds%begc:bounds%endc), source=create_betr_biogeophys_input())
    allocate(this%biogeo_flux(bounds%begc:bounds%endc), source=create_betr_biogeoFlux())
    allocate(this%biogeo_state(bounds%begc:bounds%endc), source=create_betr_biogeo_state())
    allocate(this%bstatus(bounds%begc:bounds%endc), source=create_betr_status_type())
    allocate(this%bsimstatus, source = create_betr_status_sim_type())
    call this%bsimstatus%reset()
    !grid horizontal bounds
    betr_bounds%lbj  = 1 ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begc = 1 ; betr_bounds%endc = 1
    betr_bounds%begp = 1 ; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begl = 1 ; betr_bounds%endl = 1
    betr_bounds%begg = 1 ; betr_bounds%endg = 1

    do c = bounds%begc, bounds%endc

      call this%biophys_forc(c)%Init(betr_bounds)

      call this%biogeo_state(c)%Init(betr_bounds)

      call this%biogeo_flux(c)%Init(betr_bounds)

    enddo

    call this%BeTRSetBiophysForcing(bounds, betr_bounds%lbj, betr_bounds%ubj, &
        waterstate_vars = waterstate)

    do c = bounds%begc, bounds%endc
      call this%betr(c)%Init(namelist_buffer, betr_bounds, this%biophys_forc(c), this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())
    !print*,'create hist files',betr_offline
    if(betr_offline)then
      call this%CreateHistory(betr_nlevtrc_soil, this%num_soilc)
    else
      c = bounds%begc
      call this%betr(c)%get_hist_size(this%num_hist_state1d, this%num_hist_state2d, &
         this%num_hist_flux1d, this%num_hist_flux2d, &
         this%nmlist_hist1d_state_buffer, this%nmlist_hist2d_state_buffer, &
         this%nmlist_hist1d_flux_buffer, this%nmlist_hist2d_flux_buffer)

      call this%BeTRCreateHistory(bounds, betr_nlevtrc_soil, &
         this%num_hist_state1d, this%num_hist_state2d, &
            this%num_hist_flux1d, this%num_hist_flux2d)
    endif
    call this%regression%Init(base_filename, namelist_buffer, this%bsimstatus)
    if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())
  end subroutine BeTRInit

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartInit(this, bounds, ncid, flag)
    !
    !! DESCRIPTION
    ! initialize for restart run
    ! !USES:
    use ncdio_pio, only : file_desc_t
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    class(file_desc_t)          , intent(inout) :: ncid ! netcdf id
    character(len=*)            , intent(in)    :: flag ! 'read' or 'write'

    !TEMPORARY VARIABLES
    type(betr_bounds_type)     :: betr_bounds
    integer :: lbj, ubj

    !set lbj and ubj
    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp =  betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1

   !  call this%betr%tracerstates%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betr%tracers)

   !  call this%betr%tracerfluxes%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betr%tracers)

   !  call this%betr%tracercoeffs%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betr%tracers)
  end subroutine BeTRSimulationRestartInit


  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithoutDrainage(this, betr_time, bounds, col)
  !DESCRPTION
  !interface for StepWithoutDrainage
  !
  ! USES
    use SoilStateType     , only : soilstate_type
    use WaterStateType    , only : Waterstate_Type
    use TemperatureType   , only : temperature_type
    use ChemStateType     , only : chemstate_type
    use WaterfluxType     , only : waterflux_type
    use ColumnType        , only : column_type
    use atm2lndType       , only : atm2lnd_type
    use SoilHydrologyType , only : soilhydrology_type
    use CNCarbonFluxType  , only : carbonflux_type
    use CanopyStateType   , only : canopystate_type
    use BeTR_PatchType    , only : betr_pft
    use BeTR_TimeMod      , only : betr_time_type
    use PatchType         , only : pft
    use pftvarcon         , only : crop
    implicit none
  !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    class(betr_time_type)       , intent(in)    :: betr_time
    type(bounds_type)           , intent(in)    :: bounds ! bounds
    type(column_type)           , intent(in)    :: col ! column type

    ! remove compiler warnings about unused dummy args
    if (this%num_soilc > 0)                           continue
    if (betr_time%tstep > 0)                          continue
    if (bounds%begc > 0)                              continue
    if (size(col%z) > 0)                              continue

  end subroutine BeTRSimulationStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithDrainage(this, bounds,  col)
   !
   !DESCRIPTION
   !interface for using StepWithDrainage
   !
   !USES
    use ColumnType    , only : column_type
    use MathfuncMod   , only : safe_div
    use WaterFluxType , only : waterflux_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type

    ! remove compiler warnings about unused dummy args
    if (this%num_soilc > 0) continue
    if (bounds%begc > 0)    continue
    if (size(col%z) > 0)    continue

  end subroutine BeTRSimulationStepWithDrainage

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationBeginMassBalanceCheck(this, bounds)
    !DESCRIPTION
    !stage tracer mass balance check
    !
    !USES
    use TracerBalanceMod, only : begin_betr_tracer_massbalance
    implicit none
    !ARGUMENTS
    class(betr_simulation_type), intent(inout)   :: this
    type(bounds_type), intent(in) :: bounds

    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj
    integer :: c

    !set lbj and ubj
    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp =  betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1

    do c = bounds%begc, bounds%endc
      call begin_betr_tracer_massbalance(betr_bounds,         &
         this%num_soilc, this%filter_soilc,                   &
         this%betr(c)%tracers, this%betr(c)%tracerstates,     &
         this%betr(c)%tracerfluxes, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo

  end  subroutine BeTRSimulationBeginMassBalanceCheck
  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationMassBalanceCheck(this, betr_time, bounds)
   !DESCRIPTION
   ! do tracer mass balance check
   !
   !USES
    use TracerBalanceMod , only : betr_tracer_massbalance_check
    use BeTR_TimeMod     , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    class(betr_time_type)       , intent(in)    :: betr_time
    type(bounds_type)           , intent(in)    :: bounds
    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer :: c

    !set lbj and ubj
    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp =  betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1

   do c = bounds%begc, bounds%endc
      call betr_tracer_massbalance_check(betr_time, betr_bounds,           &
         this%num_soilc, this%filter_soilc,                              &
         this%betr(c)%tracers, this%betr(c)%tracerstates,                &
         this%betr(c)%tracerfluxes, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
   enddo
  end subroutine BeTRSimulationMassBalanceCheck

!-------------------------------------------------------------------------------

  subroutine hist_htapes_create(this, nlevtrc_soil, ncol)
  !
  ! DESCRIPTIONS
  ! create history file and define output variables, only for standalone applicaitons
  !
  ! USES
    use netcdf          , only : nf90_float
    use ncdio_pio       , only : file_desc_t
    use ncdio_pio       , only : ncd_pio_createfile
    use ncdio_pio       , only : ncd_pio_closefile
    use ncdio_pio       , only : ncd_enddef
    use ncdio_pio       , only : ncd_defvar
    use ncdio_pio       , only : ncd_putvar
    use histFileMod     , only : hist_file_create, hist_def_fld1d, hist_def_fld2d
    use betr_varcon     , only : bspval
    use BeTR_ColumnType , only : betr_col
    !
    !ARGUMENTS
    implicit none

    class(betr_simulation_type) , intent(inout) :: this
    integer                     , intent(in)    :: nlevtrc_soil, ncol
  !TEMPORARY VARIABLES
    integer                     :: jj, kk, c
    type(file_desc_t)           :: ncid
    character(len=*), parameter :: subname = 'hist_htapes_create'
   c = 1
    associate(                                                     &
         ntracers          => this%betr(c)%tracers%ntracers,          &
         ngwmobile_tracers => this%betr(c)%tracers%ngwmobile_tracers, &
         is_volatile       => this%betr(c)%tracers%is_volatile,       &
         is_h2o            => this%betr(c)%tracers%is_h2o,            &
         is_isotope        => this%betr(c)%tracers%is_isotope,        &
         volatileid        => this%betr(c)%tracers%volatileid,        &
         tracernames       => this%betr(c)%tracers%tracernames        &
         )

      this%hist_filename = trim(this%base_filename) // '.output.nc'
      call ncd_pio_createfile(ncid, this%hist_filename)

      call hist_file_create(ncid,nlevtrc_soil, ncol)

     call ncd_defvar(ncid, "ZSOI", nf90_float, &
        dim1name="ncol",dim2name="levgrnd",    &
        long_name="grid center"           ,    &
        units="m", missing_value=bspval, fill_value=bspval)

      call  hist_def_fld2d(ncid, varname="TRACER_P_GAS", nf90_type=nf90_float,dim1name="ncol", &
           dim2name="levgrnd",long_name="total gas pressure", units="Pa")

      call  hist_def_fld2d(ncid, varname="QFLX_ADV", nf90_type=nf90_float,dim1name="ncol",     &
           dim2name="levgrnd",long_name="advective flux / velocity", units="m/s")

      do jj = 1, ntracers
         if(jj<= ngwmobile_tracers)then
            call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_CONC_MOIBLE',    &
                 nf90_type=nf90_float, dim1name="ncol",                                        &
                 dim2name="levgrnd", long_name=trim(tracernames(jj))//"tracer concentrations", &
                 units="mol m-3")

            if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then
               call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_P_GAS_FRAC',                     &
                    nf90_type=nf90_float, dim1name="ncol",                                                        &
                    dim2name="levgrnd", long_name='fraction of gas phase contributed by '//trim(tracernames(jj)), &
                    units="none")
            endif

            if(is_volatile(jj))then
               call hist_def_fld1d (ncid, varname=trim(tracernames(jj))//'_FLX_SURFEMI', units='mol/m2/s', &
                    nf90_type=nf90_float,  dim1name="ncol",                                                &
                    long_name='loss from surface emission for '//trim(tracernames(jj)))
            endif
         else
            call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_CONC_SOLID_PASSIVE', &
                 nf90_type=nf90_float, dim1name="ncol",                                            &
                 dim2name="levgrnd", long_name=trim(tracernames(jj))//"tracer concentrations",     &
                 units="mol m-3")
         endif
      enddo

      call ncd_enddef(ncid)
      call ncd_putvar(ncid,"ZSOI",1,betr_col%z(1:1,1:nlevtrc_soil))
      call ncd_pio_closefile(ncid)

    end associate
  end subroutine hist_htapes_create

  !-------------------------------------------------------------------------------
  subroutine hist_write(me, record, lbj, ubj, time_vars, velocity)
    !
    ! DESCRIPTION
    ! output hist file, only for standalone applications
    !
    ! USES
    use ncdio_pio    , only : file_desc_t
    use ncdio_pio    , only : ncd_pio_openfile_for_write
    use ncdio_pio    , only : ncd_putvar
    use ncdio_pio    , only : ncd_pio_closefile
    use BeTR_TimeMod , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: me
    integer                     , intent(in)    :: record
    integer                     , intent(in)    :: lbj,ubj
    type(betr_time_type)        , intent(in)    :: time_vars
    real(r8)                    , intent(in)    :: velocity(:, :)
    !TEMPORARY VARIABLES
    type(file_desc_t)           :: ncid
    integer                     :: jj
    integer                     :: c
    character(len=*), parameter :: subname='hist_write'
      c = 1
    associate(                                                   &
         ntracers          => me%betr(c)%tracers%ntracers,          &
         ngwmobile_tracers => me%betr(c)%tracers%ngwmobile_tracers, &
         is_volatile       => me%betr(c)%tracers%is_volatile,       &
         is_h2o            => me%betr(c)%tracers%is_h2o,            &
         is_isotope        => me%betr(c)%tracers%is_isotope,        &
         volatileid        => me%betr(c)%tracers%volatileid,        &
         tracernames       => me%betr(c)%tracers%tracernames        &
         )

      call ncd_pio_openfile_for_write(ncid, me%hist_filename)

      if (mod(time_vars%time, 86400._r8)==0) then
         print*,'day', time_vars%time/86400._r8
      end if
      call ncd_putvar(ncid, "time", record, time_vars%time)

      do jj = 1, ntracers
         if(jj<= ngwmobile_tracers)then
            call ncd_putvar(ncid,trim(tracernames(jj))//'_TRACER_CONC_MOIBLE',         &
                 record,me%betr(c)%tracerstates%tracer_conc_mobile_col(1:1,lbj:ubj,jj))

            if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then
               call ncd_putvar(ncid, trim(tracernames(jj))//'_TRACER_P_GAS_FRAC',      &
                    record,me%betr(c)%tracerstates%tracer_P_gas_frac_col(1:1,lbj:ubj,volatileid(jj)))
            endif
            if(is_volatile(jj))then
               call ncd_putvar(ncid, trim(tracernames(jj))//'_FLX_SURFEMI',            &
                    record, me%betr(c)%tracerfluxes%tracer_flx_surfemi_col(1:1, volatileid(jj)))
            endif
         else
            call ncd_putvar(ncid, trim(tracernames(jj))//'_TRACER_CONC_SOLID_PASSIVE', &
                 record, me%betr(c)%tracerstates%tracer_conc_solid_passive_col(1:1,lbj:ubj,jj-ngwmobile_tracers))
         endif
      enddo

      call ncd_putvar(ncid, 'TRACER_P_GAS', record, me%betr(c)%tracerstates%tracer_P_gas_col(1:1,lbj:ubj))

      call ncd_putvar(ncid, 'QFLX_ADV', record, velocity(1:1, lbj:ubj))

      call ncd_pio_closefile(ncid)
    end associate
  end subroutine hist_write

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationConsistencyCheck(this, &
     bounds, ubj, num_soilc, filter_soilc, waterstate_vars)
  ! DESCRIPTION
  ! Do consistency check, can be overwritten for varies purpose
  !
  ! USES
    use WaterStateType, only : Waterstate_Type
  implicit none
   !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
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
    if (associated(waterstate_vars%h2osoi_liq_col)) continue

  end subroutine BeTRSimulationConsistencyCheck

  !---------------------------------------------------------------------------------

  subroutine WriteRegressionOutput(this, velocity)
  !DESCRIPTION
  ! output for regression test
  !
  ! USES
    use bshr_kind_mod  , only : r8 => shr_kind_r8
    use betr_constants , only : betr_string_length
    implicit none
   !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    real(r8)                    , intent(in)    :: velocity(:, :)
  !TEMPORARY VARIABLES
    integer                           :: jj, tt, begc, endc, c
    character(len=betr_string_length) :: category
    character(len=betr_string_length) :: name

    ! FIXME(bja, 201603) should we output units as well...?
    begc = 1
    endc = 1
    c= 1
    if (this%regression%write_regression_output) then
       call this%regression%OpenOutput()
       ! NOTE(bja, 201603) currently we are allocating all tracer
       ! state vars all the time.
       do tt = 1, this%betr(c)%tracers%ntracers
          if (tt <= this%betr(c)%tracers%ngwmobile_tracers) then
             category = 'concentration'
             name = trim(this%betr(c)%tracers%tracernames(tt)) // '_total_aqueous_conc'
             call this%regression%WriteData(category, name, &
                  this%betr(c)%tracerstates%tracer_conc_mobile_col(begc, :, tt))
          end if
          if (tt <= this%betr(c)%tracers%nvolatile_tracers) then
             category = 'pressure'
             name = trim(this%betr(c)%tracers%tracernames(tt)) // '_gas_partial_pressure'
             call this%regression%WriteData(category, name, &
                  this%betr(c)%tracerstates%tracer_P_gas_frac_col(begc, :, tt))
          end if
       end do

       name = 'total_gas_pressure'
       category = 'pressure'
       call this%regression%WriteData(category, name, &
            this%betr(c)%tracerstates%tracer_P_gas_col(begc, :))

       name = 'advective flux'
       category = 'velocity'
       call this%regression%WriteData(category, name, &
            velocity(begc, :))

       call this%regression%CloseOutput()
    end if
  end subroutine WriteRegressionOutput

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcing(this, bounds,  lbj, ubj, carbonflux_vars, waterstate_vars, &
    waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars)
  !DESCRIPTION
  !pass in biogeophysical variables for running betr
  !USES
    use SoilStateType     , only : soilstate_type
    use WaterStateType    , only : Waterstate_Type
    use TemperatureType   , only : temperature_type
    use ChemStateType     , only : chemstate_type
    use WaterfluxType     , only : waterflux_type
    use atm2lndType       , only : atm2lnd_type
    use SoilHydrologyType , only : soilhydrology_type
    use CNCarbonFluxType  , only : carbonflux_type
    use CanopyStateType   , only : canopystate_type
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout)        :: this
  type(bounds_type)           , intent(in)           :: bounds
  integer                     , intent(in)           :: lbj, ubj
  type(carbonflux_type)       , optional, intent(in) :: carbonflux_vars
  type(Waterstate_Type)       , optional, intent(in) :: Waterstate_vars
  type(waterflux_type)        , optional, intent(in) :: waterflux_vars
  type(temperature_type)      , optional, intent(in) :: temperature_vars
  type(soilhydrology_type)    , optional, intent(in) :: soilhydrology_vars
  type(atm2lnd_type)          , optional, intent(in) :: atm2lnd_vars
  type(canopystate_type)      , optional, intent(in) :: canopystate_vars
  type(chemstate_type)        , optional, intent(in) :: chemstate_vars
  type(soilstate_type)        , optional, intent(in) :: soilstate_vars

  !TEMPORARY VARIABLES
  integer :: p, pi, cc, c

  cc = 1

  do c = bounds%begc, bounds%endc
  if(present(carbonflux_vars))then
    do pi = 1, betr_maxpatch_pft
       this%biophys_forc(c)%annsum_npp_patch(pi) = 0._r8
       this%biophys_forc(c)%agnpp_patch(pi) = 0._r8
       this%biophys_forc(c)%bgnpp_patch(pi)  = 0._r8
       if (pi <= col%npfts(c)) then
         p = col%pfti(c) + pi - 1
         if (pft%active(p)) then
           this%biophys_forc(c)%annsum_npp_patch(pi) = carbonflux_vars%annsum_npp_patch(p)
           this%biophys_forc(c)%agnpp_patch(pi)      = carbonflux_vars%agnpp_patch(p)
           this%biophys_forc(c)%bgnpp_patch(pi)      = carbonflux_vars%bgnpp_patch(p)
         endif
       endif
    enddo
  endif
  !assign waterstate
  if(present(waterstate_vars))then
      this%biophys_forc(c)%finundated_col(cc)            = waterstate_vars%finundated_col(c)
      this%biophys_forc(c)%frac_h2osfc_col(cc)           = waterstate_vars%frac_h2osfc_col(c)
      this%biophys_forc(c)%h2osoi_liq_col(cc,lbj:ubj)    = waterstate_vars%h2osoi_liq_col(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_ice_col(cc,lbj:ubj)    = waterstate_vars%h2osoi_ice_col(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_liq_old(cc,lbj:ubj)    = waterstate_vars%h2osoi_liq_old(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_ice_old(cc,lbj:ubj)    = waterstate_vars%h2osoi_ice_old(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_liqvol_col(cc,lbj:ubj) = waterstate_vars%h2osoi_liqvol_col(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_icevol_col(cc,lbj:ubj) = waterstate_vars%h2osoi_icevol_col(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_vol_col(cc,lbj:ubj)    = waterstate_vars%h2osoi_vol_col(c,lbj:ubj)
      this%biophys_forc(c)%air_vol_col(cc,lbj:ubj)       = waterstate_vars%air_vol_col(c,lbj:ubj)
      this%biophys_forc(c)%rho_vap(cc,lbj:ubj)           = waterstate_vars%rho_vap(c,lbj:ubj)
      this%biophys_forc(c)%rhvap_soi(cc,lbj:ubj)         = waterstate_vars%rhvap_soi(c,lbj:ubj)
      this%biophys_forc(c)%smp_l_col(cc,lbj:ubj)         = waterstate_vars%smp_l_col(c,lbj:ubj)
  endif
  if(present(waterflux_vars))then
      this%biogeo_flux(c)%qflx_infl_col(cc)             = waterflux_vars%qflx_infl_col(c)
      this%biogeo_flux(c)%qflx_totdrain_col(cc)         = waterflux_vars%qflx_totdrain_col(c)
      this%biogeo_flux(c)%qflx_gross_evap_soil_col(cc)  = waterflux_vars%qflx_gross_evap_soil_col(c)
      this%biogeo_flux(c)%qflx_gross_infl_soil_col(cc)  = waterflux_vars%qflx_gross_infl_soil_col(c)
      this%biophys_forc(c)%qflx_surf_col(cc)            = waterflux_vars%qflx_surf_col(c)
      this%biophys_forc(c)%qflx_dew_grnd_col(cc)        = waterflux_vars%qflx_dew_grnd_col(c)
      this%biophys_forc(c)%qflx_dew_snow_col(cc)        = waterflux_vars%qflx_dew_snow_col(c)
      this%biophys_forc(c)%qflx_sub_snow_vol_col(cc)    = waterflux_vars%qflx_sub_snow_vol_col(c)
      this%biophys_forc(c)%qflx_sub_snow_col(cc)        = waterflux_vars%qflx_sub_snow_col(c)
      this%biophys_forc(c)%qflx_h2osfc2topsoi_col(cc)   = waterflux_vars%qflx_h2osfc2topsoi_col(c)
      this%biophys_forc(c)%qflx_snow2topsoi_col(cc)     = waterflux_vars%qflx_snow2topsoi_col(c)
      this%biophys_forc(c)%qflx_rootsoi_col(cc,lbj:ubj) = waterflux_vars%qflx_rootsoi_col(c,lbj:ubj)
      this%biogeo_flux(c)%qflx_adv_col(cc,lbj-1:ubj)    = waterflux_vars%qflx_adv_col(c,lbj-1:ubj)
      this%biogeo_flux(c)%qflx_drain_vr_col(cc,lbj:ubj) = waterflux_vars%qflx_drain_vr_col(c,lbj:ubj)

    do pi = 1, betr_maxpatch_pft
       this%biophys_forc(c)%qflx_tran_veg_patch(pi) = 0._r8
       if (pi <= col%npfts(c)) then
         p = col%pfti(c) + pi - 1
         if (pft%active(p)) then
           this%biophys_forc(c)%qflx_tran_veg_patch(pi)     = waterflux_vars%qflx_tran_veg_patch(p)
         endif
       endif
    enddo
  endif
  if(present(temperature_vars))then
      this%biophys_forc(c)%t_soi_10cm(cc)           = temperature_vars%t_soi_10cm(c)
      this%biophys_forc(c)%t_soisno_col(cc,lbj:ubj) = temperature_vars%t_soisno_col(c,lbj:ubj)
    do pi = 1, betr_maxpatch_pft
       this%biophys_forc(c)%t_veg_patch(pi) = 0._r8
       if (pi <= col%npfts(c)) then
         p = col%pfti(c) + pi - 1
         if (pft%active(p)) then
           this%biophys_forc(c)%t_veg_patch(pi)         = temperature_vars%t_veg_patch(p)
         endif
       endif
    enddo
  endif
  if(present(soilhydrology_vars))then
      this%biophys_forc(c)%qflx_bot_col(cc)        = soilhydrology_vars%qflx_bot_col(c)
      this%biophys_forc(c)%fracice_col(cc,lbj:ubj) = soilhydrology_vars%fracice_col(c,lbj:ubj)
  endif

  if(present(atm2lnd_vars))then
      this%biophys_forc(c)%forc_pbot_downscaled_col(cc) = atm2lnd_vars%forc_pbot_downscaled_col(c)
      this%biophys_forc(c)%forc_t_downscaled_col(cc)    = atm2lnd_vars%forc_t_downscaled_col(c)
  endif

  if(present(canopystate_vars))then
    this%biophys_forc(c)%altmax_col(cc)          = canopystate_vars%altmax_col(c)
    this%biophys_forc(c)%altmax_lastyear_col(cc) = canopystate_vars%altmax_lastyear_col(c)

    do pi = 1, betr_maxpatch_pft
       this%biophys_forc(c)%lbl_rsc_h2o_patch(pi) = 0._r8
       this%biophys_forc(c)%elai_patch(pi)        = 0._r8
       if (pi <= col%npfts(c)) then
         p = col%pfti(c) + pi - 1
         if (pft%active(p)) then
           this%biophys_forc(c)%lbl_rsc_h2o_patch(pi) = canopystate_vars%lbl_rsc_h2o_patch(p)
           this%biophys_forc(c)%elai_patch(pi)        = canopystate_vars%elai_patch(p)
         endif
       endif
    enddo
  endif
  if(present(chemstate_vars))then
    this%biophys_forc(c)%soil_pH(cc,lbj:ubj) = chemstate_vars%soil_pH(c,lbj:ubj)
  endif
  if(present(soilstate_vars))then

      this%biophys_forc(c)%bsw_col(cc,lbj:ubj)          = soilstate_vars%bsw_col(c,lbj:ubj)
      this%biophys_forc(c)%watsat_col(cc,lbj:ubj)       = soilstate_vars%watsat_col(c,lbj:ubj)
      this%biophys_forc(c)%eff_porosity_col(cc,lbj:ubj) = soilstate_vars%eff_porosity_col(c,lbj:ubj)
      this%biophys_forc(c)%soilpsi_col(cc,lbj:ubj)      = soilstate_vars%soilpsi_col(c,lbj:ubj)
      this%biophys_forc(c)%cellorg_col(cc,lbj:ubj)      = soilstate_vars%cellorg_col(c,lbj:ubj)
      this%biophys_forc(c)%cellclay_col(cc,lbj:ubj)     = soilstate_vars%cellclay_col(c,lbj:ubj)
      this%biophys_forc(c)%cellsand_col(cc,lbj:ubj)     = soilstate_vars%cellsand_col(c,lbj:ubj)
      this%biophys_forc(c)%bd_col(cc,lbj:ubj)           = soilstate_vars%bd_col(c,lbj:ubj)
      this%biophys_forc(c)%watfc_col(cc,lbj:ubj)        = soilstate_vars%watfc_col(c,lbj:ubj)
      this%biophys_forc(c)%sucsat_col(cc,lbj:ubj)       = soilstate_vars%sucsat_col(c,lbj:ubj)

    do pi = 1, betr_maxpatch_pft
       this%biophys_forc(c)%rootfr_patch(pi,lbj:ubj) = 0._r8
       if (pi <= col%npfts(c)) then
         p = col%pfti(c) + pi - 1
         if (pft%active(p)) then
           this%biophys_forc(c)%rootfr_patch(pi,lbj:ubj) = soilstate_vars%rootfr_patch(p,lbj:ubj)
         endif
       endif
    enddo
   endif
  enddo
  end subroutine BeTRSimulationSetBiophysForcing

  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveBiogeoFlux(this, bounds, lbj,ubj, carbonflux_vars,  &
    waterflux_vars)
  ! DESCRIPTIONS
  ! update and return fluxes, this eventually will be expanded to
  ! include other fluxes
  ! USES
    use WaterfluxType    , only : waterflux_type
    use CNCarbonFluxType , only : carbonflux_type
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout)           :: this
  type(bounds_type)           , intent(in)              :: bounds
  integer                     , intent(in)              :: lbj, ubj
  type(carbonflux_type)       , optional, intent(inout) :: carbonflux_vars
  type(waterflux_type)        , optional, intent(inout) :: waterflux_vars

  integer :: begp, begc, endp, endc
  integer :: p, c

  if(present(carbonflux_vars))then
    !do nothing
  endif
  if(present(waterflux_vars))then
    do c = bounds%begc, bounds%endc
      waterflux_vars%qflx_infl_col(c)            = this%biogeo_flux(c)%qflx_infl_col(c)
      waterflux_vars%qflx_adv_col(c,lbj-1:ubj)   = this%biogeo_flux(c)%qflx_adv_col(c,lbj-1:ubj)
      waterflux_vars%qflx_totdrain_col(c)        = this%biogeo_flux(c)%qflx_totdrain_col(c)
      waterflux_vars%qflx_gross_evap_soil_col(c) = this%biogeo_flux(c)%qflx_gross_evap_soil_col(c)
      waterflux_vars%qflx_gross_infl_soil_col(c) = this%biogeo_flux(c)%qflx_gross_infl_soil_col(c)
    enddo
  endif

  end subroutine BeTRSimulationRetrieveBiogeoFlux

  !------------------------------------------------------------------------
  subroutine BeTRSimulationPreDiagSoilColWaterFlux(this, num_nolakec, filter_nolakec)
  !DESCRIPTION
  !prepare for water flux diagnosis. it is called before diagnosing the advective fluxes and applying
  !freeze-thaw tracer partition.
  !
  !USES

  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
   integer                    , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                    , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points

  !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp =  betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1

   do fc= 1, num_nolakec
     c = filter_nolakec(fc)
     call this%betr(c)%pre_diagnose_soilcol_water_flux(betr_bounds, this%num_soilc, &
       this%filter_soilc, this%biophys_forc(c))
   enddo
  end subroutine BeTRSimulationPreDiagSoilColWaterFlux

  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagAdvWaterFlux(this, betr_time, num_hydrologyc, &
    filter_hydrologyc)

  !DESCRIPTION
  ! diagnose water fluxes for tracer advection
  !
  ! USES
  !

  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   class(betr_time_type)       , intent(in)    :: betr_time
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1

   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     call this%betr(c)%diagnose_advect_water_flux(betr_time,                      &
       betr_bounds, this%num_soilc, this%filter_soilc,                         &
       this%biophys_forc(c), this%biogeo_flux(c))
   enddo

  end subroutine BeTRSimulationDiagAdvWaterFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagDrainWaterFlux(this, betr_time, &
       num_hydrologyc, filter_hydrologyc)
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
   class(betr_simulation_type) , intent(inout) :: this
   class(betr_time_type)       , intent(in)    :: betr_time
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1

   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     call this%betr(c)%diagnose_drainage_water_flux(betr_time, &
       betr_bounds, this%num_soilc, this%filter_soilc,      &
       this%biophys_forc(c), this%biogeo_flux(c))
  enddo
  end subroutine BeTRSimulationDiagDrainWaterFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationBeginTracerSnowLayerAdjst(this,  num_snowc, filter_snowc)
  !DESCRIPTION
  !prepare for tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1
   do fc = 1, num_snowc
     c = filter_snowc(fc)
     call this%betr(c)%Enter_tracer_LayerAdjustment(betr_bounds, this%num_soilc, this%filter_soilc)
   enddo
  end subroutine BeTRSimulationBeginTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine BeTRSimulationEndTracerSnowLayerAdjst(this, num_snowc, filter_snowc)
  !DESCRIPTION
  !wrap up tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1
  do fc = 1, num_snowc
    c = filter_snowc(fc)
    call this%betr(c)%Exit_tracer_LayerAdjustment(betr_bounds, this%num_soilc, this%filter_soilc)
  enddo

  end subroutine BeTRSimulationEndTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine BeTRSimulationDvideSnowLayers(this, bounds, num_snowc, filter_snowc, divide_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use clm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: divide_matrix(bounds%begc:bounds%endc , 1:nlevsno , 1:nlevsno )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1
  do fc = 1, num_snowc
    c = filter_snowc(fc)
    call this%betr(c)%tracer_DivideSnowLayers(betr_bounds, this%num_soilc, &
      this%filter_soilc, divide_matrix(c:c,:,:), this%bstatus(c))
    if(this%bstatus(c)%check_status())then
      call this%bsimstatus%setcol(c)
      call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
      exit
    endif
  enddo

  end subroutine BeTRSimulationDvideSnowLayers

  !------------------------------------------------------------------------
  subroutine BeTRSimulationCombineSnowLayers(this, bounds, num_snowc, filter_snowc, combine_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use clm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: combine_matrix(bounds%begc:bounds%endc,-nlevsno+1:1 ,-nlevsno+1:1 )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c
    betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begc = 1; betr_bounds%endc = 1
    betr_bounds%begl = 1; betr_bounds%endl = 1
    betr_bounds%begg = 1; betr_bounds%endg = 1

  do fc = 1, num_snowc
    c = filter_snowc(fc)
    call this%betr(c)%tracer_CombineSnowLayers(betr_bounds, this%num_soilc,&
      this%filter_soilc, combine_matrix(c:c,:,:),this%bstatus(c))
    if(this%bstatus(c)%check_status())then
      call this%bsimstatus%setcol(c)
      call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
      exit
    endif
  enddo

  end subroutine BeTRSimulationCombineSnowLayers
  !------------------------------------------------------------------------
  subroutine hist_create_fluxes(this, bounds, betr_nlevtrc_soil, num_flux1d, num_flux2d)
  !
  !DESCRIPTION
  !create history file for betr fluxes
  !
  use betr_varcon         , only : spval => bspval
  use histFileMod   , only: hist_addfld1d, hist_addfld2d
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: betr_nlevtrc_soil
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d

  !local variables
  integer :: jj, begc, endc
  character(len=100) :: fname
  character(len=30) :: units
  character(len=20) :: avgflag
  character(len=20) :: type2d
  character(len=200) :: long_name
  character(len=20) :: default
  integer :: nml_error
  character(len=200):: ioerror_msg

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  namelist /hist2d_fmt/    &
  fname, units, avgflag,type2d,long_name, default

  namelist /hist1d_fmt/    &
  fname, units, avgflag, long_name, default

  begc=bounds%begc; endc=bounds%endc
  if(num_flux2d >0)then
    allocate(this%fluxes_2d(begc:endc, 1:betr_nlevtrc_soil, 1:num_flux2d))
  endif
  if(num_flux1d>0)then
    allocate(this%fluxes_1d(begc:endc, 1:num_flux1d))
  endif

  do jj = 1, num_flux2d
    !read name list
    read(this%nmlist_hist2d_flux_buffer(jj), nml=hist2d_fmt, iostat=nml_error, iomsg=ioerror_msg)
    if(nml_error/=0)then
      write(*,*)'reading ',jj,'-th namelist failed'//ioerror_msg
    endif
    this%fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
    data2dptr => this%states_2d(begc:endc,1:betr_nlevtrc_soil, jj)
    call hist_addfld2d (fname=fname, units=units, type2d=type2d,  &
           avgflag=avgflag, long_name=long_name,  ptr_col=data2dptr, default=default)
  enddo

  do jj = 1, num_flux1d
    !read name list
    read(this%nmlist_hist1d_flux_buffer(jj), nml=hist1d_fmt, iostat=nml_error, iomsg=ioerror_msg)
    if(nml_error/=0)then
      write(*,*)'reading ',jj,'-th namelist failed'//ioerror_msg
    endif
    this%fluxes_1d(begc:endc,jj) = spval
    data1dptr => this%fluxes_1d(begc:endc, jj)
    call hist_addfld1d (fname=fname, units=units,  avgflag=avgflag, long_name=long_name, &
      ptr_col=data1dptr, default=default)
  enddo
  end subroutine hist_create_fluxes
  !------------------------------------------------------------------------
  subroutine hist_create_states(this, bounds, betr_nlevtrc_soil, num_state1d, num_state2d)
  !
  !create history file for betr states variables
  use histFileMod   , only: hist_addfld1d, hist_addfld2d
  use betr_varcon         , only : spval => bspval
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d

  !local variables
  integer :: begc, endc
  integer :: jj
  character(len=100) :: fname
  character(len=30) :: units
  character(len=20) :: avgflag
  character(len=20) :: type2d
  character(len=200) :: long_name
  character(len=20) :: default
  integer :: nml_error
  character(len=200):: ioerror_msg
  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  namelist /hist2d_fmt/    &
  fname, units, avgflag,type2d,long_name, default

  namelist /hist1d_fmt/    &
  fname, units, avgflag, long_name, default

  begc = bounds%begc; endc = bounds%endc

  if(num_state2d>0)then
    allocate(this%states_2d(bounds%begc:bounds%endc, 1:betr_nlevtrc_soil, 1:num_state2d))
  endif
  if(num_state1d>0)then
    allocate(this%states_1d(bounds%begc:bounds%endc, 1:num_state1d))
  endif

  do jj = 1, num_state2d
    !read namelist
    read(this%nmlist_hist2d_state_buffer(jj), nml=hist2d_fmt, iostat=nml_error, iomsg=ioerror_msg)
    if(nml_error/=0)then
      write(*,*)'reading ',jj,'-th namelist failed'//ioerror_msg
    endif
    this%states_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
    data2dptr => this%states_2d(begc:endc,1:betr_nlevtrc_soil, jj)
    call hist_addfld2d (fname=fname, units=units, type2d=type2d,  &
           avgflag=avgflag, long_name=long_name,  ptr_col=data2dptr, default=default)
  enddo

  do jj = 1, num_state1d
    !read namelist
    read(this%nmlist_hist1d_state_buffer(jj), nml=hist1d_fmt, iostat=nml_error, iomsg=ioerror_msg)
    if(nml_error/=0)then
      write(*,*)'reading ',jj,'-th namelist failed'//ioerror_msg
    endif
    this%states_1d(begc:endc,jj) = spval
    data1dptr => this%states_1d(begc:endc,jj)
    call hist_addfld1d (fname=fname, units=units, avgflag=avgflag, &
      long_name=long_name, ptr_col=data1dptr, default=default)
  enddo
  end subroutine hist_create_states
  !------------------------------------------------------------------------
  subroutine BeTRSimulationCreateHistory(this, bounds, betr_nlevtrc_soil,&
     num_state1d, num_state2d, num_flux1d, num_flux2d)
  use betr_varcon         , only : spval => bspval
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           , intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d

  call this%hist_create_states(bounds, betr_nlevtrc_soil, num_state1d, num_state2d)

  call this%hist_create_fluxes(bounds, betr_nlevtrc_soil, num_flux1d, num_flux2d)

  end subroutine BeTRSimulationCreateHistory
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveHistory(this, numf, filter)
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !TEMPORARY VARIABLES
  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds


  betr_bounds%lbj  = 1 ; betr_bounds%ubj  = betr_nlevsoi
  betr_bounds%begc = 1 ; betr_bounds%endc = 1
  betr_bounds%begp = 1 ; betr_bounds%endp = betr_maxpatch_pft
  betr_bounds%begl = 1 ; betr_bounds%endl = 1
  betr_bounds%begg = 1 ; betr_bounds%endg = 1

  do fc = 1, numf
    c = filter(fc)
    call this%betr(c)%HistRetrieve(betr_bounds, 1, betr_nlevtrc_soil, &
       this%num_hist_state1d, this%num_hist_state2d, this%num_hist_flux1d,&
       this%num_hist_flux2d, this%states_1d(c:c,:), &
       this%states_2d(c:c,1:betr_nlevtrc_soil,:), this%fluxes_1d(c:c,:),&
       this%fluxes_2d(c:c,1:betr_nlevtrc_soil,:))
  enddo

  end subroutine BeTRSimulationRetrieveHistory
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRestart(this, bounds, ncid, flag, numf, filter)
  !DESCRIPTION
  !create or read restart file
  use ncdio_pio      , only : file_desc_t
  use betr_varcon         , only : spval => bspval
  use restUtilMod    , only : restartvar
  use ncdio_pio      , only : file_desc_t
  use ncdio_pio      , only : ncd_double
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)    , intent(in)    :: bounds
  class(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
  character(len=*)     , intent(in)    :: flag                                         ! 'read' or 'write'
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !local variables
  real(r8), pointer :: states_1d(:,:)
  real(r8), pointer :: states_2d(:,:,:)
  integer :: nrest_1d, nrest_2d
  integer :: c, jj, fc
  character(len=255) :: rest_varname_1d(200)
  character(len=255) :: rest_varname_2d(200)
  logical :: readvar      ! determine if variable is on initial file
  real(r8), pointer :: ptr1d(:)
  real(r8), pointer :: ptr2d(:,:)
  type(betr_bounds_type)     :: betr_bounds

  c = bounds%begc
  call this%betr(c)%get_restartvar(nrest_1d, nrest_2d,rest_varname_1d, &
     rest_varname_2d)

  allocate(states_1d(bounds%begc:bounds%endc, 1:nrest_1d))
  allocate(states_2d(bounds%begc:bounds%endc, 1:betr_nlevtrc_soil, 1:nrest_2d))

  if(trim(flag)/='define')then
    !assign initial conditions
    betr_bounds%lbj  = 1 ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begc = 1 ; betr_bounds%endc = 1
    betr_bounds%begp = 1 ; betr_bounds%endp = betr_maxpatch_pft
    betr_bounds%begl = 1 ; betr_bounds%endl = 1
    betr_bounds%begg = 1 ; betr_bounds%endg = 1
    do fc = 1, numf
      c = filter(fc)
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, nrest_1d,&
        nrest_2d, states_1d(c:c,:), states_2d(c:c,:,:), flag)
    enddo
  endif

  do jj = 1, nrest_1d
    ptr1d => states_1d(:, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_1d(jj)), &
       xtype=ncd_double,  dim1name='column', long_name='',  units='', &
       interpinic_flag='interp' , readvar=readvar, data=ptr1d)
  enddo
  do jj = 1, nrest_2d
    ptr2d => states_2d(:, :, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_2d(jj)), xtype=ncd_double,  &
      dim1name='column',dim2name='levtrc', switchdim=.true., &
      long_name='',  units='', fill_value=spval, &
      interpinic_flag='interp', readvar=readvar, data=ptr2d)
  enddo




  deallocate(states_1d)
  deallocate(states_2d)
  end subroutine BeTRSimulationRestart

end module BeTRSimulation
