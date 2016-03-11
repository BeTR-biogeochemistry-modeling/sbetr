module BeTRSimulation
  !
  ! !DESCRIPTION:
  !  factory to load the specific betr simulator
  !
  ! !USES:
  !
  use abortutils                  , only : endrun
  use clm_varctl                  , only : iulog
  use shr_log_mod                 , only : errMsg => shr_log_errMsg
  use tracer_varcon               , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use BeTR_decompMod              , only : betr_bounds_type
  use decompMod                   , only : bounds_type

  ! !USES:
  use BeTRTracerType            , only : BeTRtracer_type
  use TracerCoeffType           , only : TracerCoeff_type
  use TracerFluxType            , only : TracerFlux_type
  use TracerStateType           , only : TracerState_type
  use tracerboundarycondType    , only : tracerboundarycond_type
  use BGCReactionsMod           , only : bgc_reaction_type
  use PlantSoilBGCMod           , only : plant_soilbgc_type
  use BeTR_aerocondType         , only : betr_aerecond_type

  implicit none
  private
  character(len=*), private, parameter :: mod_filename = __FILE__

  type, public :: betr_simulation_type
     character(len=128) :: reaction_method

     type(BeTRtracer_type), public :: betrtracer_vars
     type(TracerCoeff_type), public :: tracercoeff_vars
     type(TracerFlux_type), public :: tracerflux_vars
     type(TracerState_type), public :: tracerState_vars
     type(tracerboundarycond_type), public :: tracerboundarycond_vars
     class(plant_soilbgc_type), allocatable,public :: plant_soilbgc
     class(bgc_reaction_type), allocatable,public :: bgc_reaction
     type(betr_aerecond_type), public :: betr_aerecond_vars
   contains
     procedure, public :: Init => BeTRSimulationInit
     procedure, public :: ReadNameList => BeTRSimulationReadNameList
     procedure, public :: RestartInit => BeTRSimulationRestartInit
     procedure, public :: StepWithoutDrainage => BeTRSimulationStepWithoutDrainage
     procedure, public :: StepWithDrainage => BeTRSimulationStepWithDrainage
     procedure, public :: BeginMassBalanceCheck => BeTRSimulationBeginMassBalanceCheck
     procedure, public :: MassBalanceCheck      => BeTRSimulationMassBalanceCheck
  end type betr_simulation_type

  public :: BeTRSimulationInit

contains

!-------------------------------------------------------------------------------

  subroutine BeTRSimulationInit(this, reaction_method, bounds, &
       waterstate)
    !
    use TransportMod          , only : init_transportmod
    use TracerParamsMod       , only : tracer_param_init
    use WaterstateType        , only : waterstate_type
    use BeTR_WaterstateType   , only : betr_waterstate_type
    implicit none

    class(betr_simulation_type) :: this
    character(len=*), intent(in) :: reaction_method
    type(bounds_type)    , intent(in) :: bounds

    type(waterstate_type), intent(inout) :: waterstate

    character(len=32) :: subname='BeTRSimulationInit'
    type(betr_waterstate_type) :: betr_waterstate
    type(betr_bounds_type)     :: betr_bounds
    integer :: lbj, ubj

    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj
    betr_waterstate%h2osoi_liq_col => waterstate%h2osoi_liq_col
    betr_waterstate%h2osoi_ice_col    => waterstate%h2osoi_ice_col


    call this%betrtracer_vars%init_scalars()

    call this%bgc_reaction%Init_betrbgc(betr_bounds, lbj, ubj, this%betrtracer_vars)

    call this%betr_aerecond_vars%Init(betr_bounds)

    call init_transportmod()

    call this%tracerState_vars%Init(betr_bounds, lbj, ubj, this%betrtracer_vars)

    call this%tracerflux_vars%Init(betr_bounds,  lbj, ubj, this%betrtracer_vars)

    call this%tracercoeff_vars%Init(betr_bounds, lbj, ubj, this%betrtracer_vars)

    call this%tracerboundarycond_vars%Init(betr_bounds, this%betrtracer_vars)

    !inside Init_plant_soilbgc, specific plant soil bgc coupler data type will be created
    call this%plant_soilbgc%Init_plant_soilbgc(betr_bounds, lbj, ubj)

    !initialize state variable
    call this%bgc_reaction%initCold(betr_bounds,  this%betrtracer_vars, betr_waterstate, this%tracerstate_vars)

    !initialize boundary condition type
    call this%bgc_reaction%init_boundary_condition_type(betr_bounds, this%betrtracer_vars, this%tracerboundarycond_vars)

    !initialize the betr parameterization module
    call tracer_param_init(betr_bounds)

    !initialize the betrBGC module
    !X!call betrbgc_init(bounds) - NOTE(bja, 2016-03) empty subroutine...

  end subroutine BeTRSimulationInit

  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationReadNameList(this, filename)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_type) :: this
    character(len=*), intent(IN) :: filename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                    ! error code
    integer :: unitn                   ! unit for namelist file
    character(len=32) :: subname = 'BeTRSimulationReadNameList'
    character(len=128) :: reaction_method

    !-----------------------------------------------------------------------

    namelist / betr_inparm / reaction_method

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in betr_inparm  namelist'
       call opnfil (filename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'betr_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, betr_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading betr_inparm namelist"//errmsg(filename, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if
    ! Broadcast namelist variables read in
    this%reaction_method = reaction_method
    call shr_mpi_bcast(reaction_method, mpicom)

  end subroutine BeTRSimulationReadNameList

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartInit(this, bounds, ncid, flag)
    !
    !! DESCRIPTION
    ! initialize for restart run
    ! !USES:
    use ncdio_pio, only : file_desc_t
    implicit none

    class(betr_simulation_type) :: this
    type(bounds_type), intent(in) :: bounds
    class(file_desc_t), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in)    :: flag ! 'read' or 'write'

    type(betr_bounds_type)     :: betr_bounds
    integer :: lbj, ubj

    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    call this%tracerstate_vars%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betrtracer_vars)

    call this%tracerflux_vars%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betrtracer_vars)

    call this%tracercoeff_vars%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betrtracer_vars)
  end subroutine BeTRSimulationRestartInit


  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithoutDrainage(this, bounds, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)

    !this will be overridden in actual application
    use BetrBGCMod         , only : run_betr_one_step_without_drainage
    use SoilStateType      , only : soilstate_type
    use WaterStateType     , only : Waterstate_Type
    use TemperatureType    , only : temperature_type
    use ChemStateType      , only : chemstate_type
    use WaterfluxType      , only : waterflux_type
    use ColumnType         , only : column_type
    use BGCReactionsMod    , only : bgc_reaction_type
    use atm2lndType        , only : atm2lnd_type
    use SoilHydrologyType  , only : soilhydrology_type
    use BeTR_CNStateType   , only : betr_cnstate_type
    use BeTR_CarbonFluxType, only : betr_carbonflux_type
    use CNStateType        , only : cnstate_type
    use CNCarbonFluxType   , only : carbonflux_type
    use CanopyStateType    , only : canopystate_type
    use BeTR_PatchType     , only : betr_pft
    use PatchType          , only : pft
    use pftvarcon          , only : crop

    implicit none
    class(betr_simulation_type) :: this
    type(bounds_type), intent(in) :: bounds ! bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: num_soilp
    integer, intent(in) :: filter_soilp(:) ! pft filter


    type(column_type), intent(in) :: col ! column type
    type(Waterstate_Type), intent(in) :: waterstate_vars ! water state variables
    type(soilstate_type), intent(in) :: soilstate_vars ! column physics variable
    type(temperature_type), intent(in) :: temperature_vars ! energy state variable
    type(chemstate_type), intent(in) :: chemstate_vars
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    type(cnstate_type), intent(inout) :: cnstate_vars
    type(canopystate_type), intent(in) :: canopystate_vars
    type(carbonflux_type), intent(in) :: carbonflux_vars
    type(waterflux_type), intent(inout) :: waterflux_vars !temporary variables
    type(betr_cnstate_type) :: betr_cnstate_vars
    type(betr_carbonflux_type) :: betr_carbonflux_vars !pass necessary data for correct subroutine call


  end subroutine BeTRSimulationStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithDrainage(this, bounds, num_soilc, &
       filter_soilc, jtops, waterflux_vars, col)
    use BetrBGCMod, only : run_betr_one_step_with_drainage
    use ColumnType, only : column_type
    use MathfuncMod, only : safe_div
    use WaterFluxType, only : waterflux_type

    implicit none
    class(betr_simulation_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: jtops(bounds%begc: )
    type(waterflux_type)    , intent(in) :: waterflux_vars
    type(column_type), intent(in) :: col ! column type

  end subroutine BeTRSimulationStepWithDrainage


  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationBeginMassBalanceCheck(this, bounds, num_soilc, &
       filter_soilc)
  use TracerBalanceMod, only : begin_betr_tracer_massbalance
  implicit none
    class(betr_simulation_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc

    type(betr_bounds_type)     :: betr_bounds
    integer  :: lbj, ubj
    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    call begin_betr_tracer_massbalance(betr_bounds, lbj, ubj, num_soilc, filter_soilc, &
         this%betrtracer_vars, this%tracerstate_vars, &
         this%tracerflux_vars)

  end  subroutine BeTRSimulationBeginMassBalanceCheck
  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationMassBalanceCheck(this, bounds, num_soilc, &
       filter_soilc)
  use TracerBalanceMod, only : betr_tracer_massbalance_check
  implicit none
    class(betr_simulation_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc

    integer :: lbj, ubj
    type(betr_bounds_type)     :: betr_bounds

    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    call betr_tracer_massbalance_check(betr_bounds, lbj, ubj, num_soilc, filter_soilc, &
         this%betrtracer_vars, this%tracerstate_vars, &
         this%tracerflux_vars)
  end subroutine BeTRSimulationMassBalanceCheck
end module BeTRSimulation
