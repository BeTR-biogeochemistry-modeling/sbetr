module betr_initializeMod
  !
  ! !DESCRIPTION:
  !  subroutines to initialize betr based on namelist

  ! !USES:
  use BeTRTracerType            , only : BeTRtracer_type
  use TracerCoeffType           , only : TracerCoeff_type
  use TracerFluxType            , only : TracerFlux_type
  use TracerStateType           , only : TracerState_type
  use tracerboundarycondType    , only : tracerboundarycond_type
  use BGCReactionsMod           , only : bgc_reaction_type
  use PlantSoilBGCMod           , only : plant_soilbgc_type
  use BeTR_aerocondType         , only : betr_aerecond_type
  use clm_varctl                , only : iulog
  use abortutils                , only : endrun
  use shr_log_mod               , only : errMsg => shr_log_errMsg
  implicit none
  save
  private   ! By default everything is public

  public :: betr_initialize
  public :: betr_readNL
  public :: betr_rest
  character(len=32) :: bgc_method='mock_run'

  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  type(BeTRtracer_type)                , public :: betrtracer_vars
  type(TracerCoeff_type)               , public :: tracercoeff_vars
  type(TracerFlux_type)                , public :: tracerflux_vars
  type(TracerState_type)               , public :: tracerState_vars
  type(tracerboundarycond_type)        , public :: tracerboundarycond_vars
  class(plant_soilbgc_type), allocatable,public :: plant_soilbgc
  class(bgc_reaction_type) , allocatable,public :: bgc_reaction
  type(betr_aerecond_type)              ,public :: betr_aerecond_vars
contains

  !-------------------------------------------------------------------------------
  subroutine betr_readNL(NLFilename)
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
    character(len=*), intent(IN) :: NLFilename              ! Namelist filename
                                                            !
                                                            ! !LOCAL VARIABLES:
    integer                      :: ierr                    ! error code
    integer                      :: unitn                   ! unit for namelist file
    character(len=32)            :: subname = 'betr_readNL' ! subroutine name
    !-----------------------------------------------------------------------

    namelist / betr_inparm / bgc_method

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in betr_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_CanopyHydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, betr_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading betr_inparm namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if
    ! Broadcast namelist variables read in
    call shr_mpi_bcast(bgc_method, mpicom)

  end subroutine betr_readNL

  !-------------------------------------------------------------------------------
  subroutine betr_initialize(bounds, lbj, ubj, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Initialize BeTR
    !
    ! !USES:
    use decompMod             , only : bounds_type
    use BGCReactionsFactoryMod, only : create_betr_application
    use BetrBGCMod            , only : betrbgc_init
    use TransportMod          , only : init_transportmod
    use TracerParamsMod       , only : tracer_param_init
    use WaterstateType        , only : waterstate_type

    implicit none
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in) :: bounds
    integer              , intent(in) :: lbj, ubj
    type(waterstate_type), intent(in) :: waterstate_vars

    character(len=32) :: subname='betr_initialize'

    call betrtracer_vars%init_scalars()


    call create_betr_application(bgc_reaction, plant_soilbgc, bgc_method)

    call bgc_reaction%Init_betrbgc(bounds, lbj, ubj, betrtracer_vars)

    call betr_aerecond_vars%Init(bounds)

    call init_transportmod

    call tracerState_vars%Init(bounds, lbj, ubj, betrtracer_vars)

    call tracerflux_vars%Init(bounds,  lbj, ubj, betrtracer_vars)

    call tracercoeff_vars%Init(bounds, lbj, ubj, betrtracer_vars)

    call tracerboundarycond_vars%Init(bounds, betrtracer_vars)

    !inside Init_plant_soilbgc, specific plant soil bgc coupler data type will be created
    call plant_soilbgc%Init_plant_soilbgc(bounds, lbj, ubj)

    !initialize state variable
    call bgc_reaction%initCold(bounds,  betrtracer_vars, waterstate_vars, tracerstate_vars)

    !initialize boundary condition type
    call bgc_reaction%init_boundary_condition_type(bounds, betrtracer_vars, tracerboundarycond_vars)

    !initialize the betr parameterization module
    call tracer_param_init(bounds)

    !initialize the betrBGC module
    call betrbgc_init(bounds, betrtracer_vars)

  end subroutine betr_initialize
  !---------------------------------------------------------------------------------


  subroutine betr_rest(bounds, ncid, flag)
  !
  !! DESCRIPTION
  ! initialize for restart run
  ! !USES:
  use ncdio_pio             , only : file_desc_t
  use decompMod             , only : bounds_type
  implicit none
  type(bounds_type)    , intent(in) :: bounds
  class(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
  character(len=*)     , intent(in)    :: flag                                         ! 'read' or 'write'


  call tracerstate_vars%Restart(bounds, ncid, flag=flag, betrtracer_vars=betrtracer_vars)
  call tracerflux_vars%Restart( bounds, ncid, flag=flag, betrtracer_vars=betrtracer_vars)
  call tracercoeff_vars%Restart(bounds, ncid, flag=flag, betrtracer_vars=betrtracer_vars)
  end subroutine betr_rest
end module betr_initializeMod
