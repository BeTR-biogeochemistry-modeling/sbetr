module BeTRSimulationALM
  !
  ! !DESCRIPTION:
  !  factory to load the specific betr simulator
  !
  ! !USES:
  !
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use shr_log_mod , only : errMsg => shr_log_errMsg

  use BeTRSimulation  , only : betr_simulator_type
  use decompMod       , only : bounds_type
  use BeTRSimulation  , only : betr_simulation_type
  use BeTR_CNStateType, only : betr_cnstate_type
  use EcophysConType  , only : ecophyscon_type

  implicit none

  private
  character(len=*), private, parameter :: mod_filename = __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_alm_type

   contains
     procedure :: Init => ALMInit
     procedure, public :: StepWithoutDrainage => ALMStepWithoutDrainage
     procedure, public :: StepWithDrainage    => ALMStepWithDrainage

     procedure, public :: BetrReadParams            => ALMBetrReadParams
     procedure, public :: DiagnoseDtracerFreezeThaw => ALMDiagnoseDtracerFreezeThaw
     procedure, public :: CalcDewSubFlux            => ALMCalcDewSubFlux
     procedure, public :: BetrFluxStateReceive      => ALMBetrFluxStateReceive
     procedure, public :: CalcSmpL                  => ALMCalcSmpL
     procedure, public :: BetrPlantSoilBGCSend      => ALMBetrPlantSoilBGCSend
  end type betr_simulation_alm_type

  public :: create_betr_simulation_alm

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_alm()
    implicit none
    class(betr_simulation_alm_type), pointer :: create_betr_simulation_alm
    class(betr_simulation_alm_type), pointer :: simulation

    allocate(simulation)
    call simulation%Init()
    create_betr_simulation_alm => simulation

  end function create_betr_simulation_alm

!-------------------------------------------------------------------------------

  subroutine ALMInit(this, reaction_method, bounds, lbj, ubj, waterstate)

    implicit none
    class(betr_simulation_alm_type)   :: this
    character(len=*)     , intent(in) :: reaction_method
    type(bounds_type)    , intent(in) :: bounds
    integer              , intent(in) :: lbj, ubj
    type(waterstate_type), intent(inout) :: waterstate

  end subroutine ALMInit
!-------------------------------------------------------------------------------
  subroutine ALMStepWithoutDrainage(this, bounds, lbj, ubj, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)

    implicit none

    class(betr_simulation_alm_type) :: this

    ! !ARGUMENTS :
    type(bounds_type)       , intent(in) :: bounds ! bounds
    integer                 , intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer                 , intent(in) :: filter_soilc(:) ! column filter_soilc
    integer                 , intent(in) :: num_soilp
    integer                 , intent(in) :: filter_soilp(:) ! pft filter
    integer                 , intent(in) :: lbj, ubj ! lower and upper bounds, make sure they are > 0

    type(column_type)       , intent(in) :: col ! column type
    type(Waterstate_Type)   , intent(in) :: waterstate_vars ! water state variables
    type(soilstate_type)    , intent(in) :: soilstate_vars ! column physics variable
    type(temperature_type)  , intent(in) :: temperature_vars ! energy state variable
    type(chemstate_type)    , intent(in) :: chemstate_vars
    type(atm2lnd_type)      , intent(in) :: atm2lnd_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    type(cnstate_type)      , intent(inout) :: cnstate_vars
    type(canopystate_type)  , intent(in) :: canopystate_vars
    type(carbonflux_type)   , intent(in) :: carbonflux_vars
    type(waterflux_type)    , intent(inout) :: waterflux_vars

    !temporary variables
    type(betr_carbonflux_type) :: betr_carbonflux_vars

  end subroutine ALMStepWithoutDrainage


  !---------------------------------------------------------------------------------
  subroutine ALMStepWithDrainage(this, bounds, lbj, ubj, &
       num_soilc, filter_soilc, jtops, waterflux_vars, col)

    implicit none

    ! !ARGUMENTS:
    class(betr_simulation_standalone_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: lbj, ubj
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: jtops(bounds%begc: )
    type(waterflux_type)    , intent(in) :: waterflux_vars
    type(column_type), intent(in) :: col ! column type


  end subroutine ALMStepWithDrainage


  !------------------------------------------------------------------------
  subroutine ALMBetrPlantSoilBGCSend(bounds, num_soilc, &
           filter_soilc)

  implicit none
  type(bounds_type)         , intent(in)   :: bounds
  integer                   , intent(in)   :: num_soilc
  integer                   , intent(in)   :: filter_soilc(:)


  call plant_soilbgc%init_plant_soil_feedback(bounds, num_soilc, &
               filter_soilc,  carbonstate_vars%frootc_patch, cnstate_vars, &
               soilstate_vars,  waterflux_vars, ecophyscon_vars)

  end subroutine ALMBetrPlantSoilBGCSend

  !------------------------------------------------------------------------

  subroutine ALMDiagnoseDtracerFreezeThaw(bounds, num_nolakec, filter_nolakec, col, lun, &
    waterstate_vars)
  !
  ! DESCRIPTION
  ! aqueous tracer partition based on freeze-thaw
  !
  ! USES
  !
  use ColumnType            , only : column_type
  use LandunitType          , only : landunit_type
  use WaterStateType        , only : waterstate_type
  use landunit_varcon       , only : istsoil
  use TracerParamsMod       , only : diagnose_dtracer_freeze_thaw
  !
  ! Arguments
  implicit none
  type(bounds_type)      , intent(in)    :: bounds
  integer                , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
  integer                , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points
  type(landunit_type)    , intent(in)    :: lun
  type(waterstate_type)  , intent(in)    :: waterstate_vars
  type(column_type)      , intent(in)    :: col                                ! column type


  call diagnose_dtracer_freeze_thaw(bounds, num_nolakec, filter_nolakec, lun, &
    waterstate_vars, betrtracer_vars, tracerstate_vars)

  end subroutine ALMDiagnoseDtracerFreezeThaw


  !------------------------------------------------------------------------
  subroutine ALMCalcDewSubFlux(bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars)

  use ColumnType            , only : col
  use LandunitType          , only : lun
  use WaterfluxType         , only : waterflux_type
  use BetrBGCMod            , only : calc_dew_sub_flux, check_mass_err
  use WaterstateType        , only : waterstate_type
  use clm_varcon            , only : denh2o,spval
  use landunit_varcon       , only : istsoil, istcrop

  ! !ARGUMENTS:
  type(bounds_type)         , intent(in)    :: bounds
  integer                   , intent(in)    :: num_hydrologyc             ! number of column soil points in column filter_soilc
  integer                   , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
  type(waterstate_type)     , intent(in)    :: waterstate_vars
  type(waterflux_type)      , intent(in)    :: waterflux_vars

  call calc_dew_sub_flux(bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars, betrtracer_vars, tracerflux_vars, tracerstate_vars)

  end subroutine ALMCalcDewSubFlux

  !------------------------------------------------------------------------
  subroutine ALMBetrFluxStateReceive(bounds, num_soilc, filter_soilc)

  implicit none
  ! !ARGUMENTS:
  type(bounds_type)         , intent(in)   :: bounds
  integer                   , intent(in)   :: num_soilc
  integer                   , intent(in)   :: filter_soilc(:)

  call bgc_reaction%betr_lsm_flux_state_sendback(bounds, &
       num_soilc, filter_soilc,                    &
       tracerstate_vars, tracerflux_vars,  betrtracer_vars)

  end subroutine ALMBetrFluxStateReceive

  !------------------------------------------------------------------------
  subroutine ALMCalcSmpL(bounds, lbj, ubj, numf, filter, t_soisno, soilstate_vars, waterstate_vars, soil_water_retention_curve)

  use SoilStateType              , only : soilstate_type
  use WaterStateType             , only : waterstate_type
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
  use clm_varcon                 , only : grav,hfus,tfrz

  implicit none
  type(bounds_type)         , intent(in)    :: bounds  ! bounds
  integer                   , intent(in)    :: lbj, ubj                                          ! lower and upper bounds, make sure they are > 0
  integer                   , intent(in)    :: numf                                              ! number of columns in column filter
  integer                   , intent(in)    :: filter(:)                                         ! column filter
  real(r8)                  , intent(in)    :: t_soisno(bounds%begc: , lbj: )                    ! soil temperature
  type(soilstate_type)      , intent(in)    :: soilstate_vars
  type(waterstate_type)     , intent(inout) :: waterstate_vars
  class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

  !local variables
  real(r8) :: s_node
  integer  :: fc, c, j

  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(__FILE__,__LINE__))


  associate(                                                     & !
    h2osoi_vol        =>    waterstate_vars%h2osoi_vol_col     , & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
    smp_l             =>    waterstate_vars%smp_l_col          , & ! Output: [real(r8) (:,:) ]  soil suction (mm)
    bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
    watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
    sucsat            =>    soilstate_vars%sucsat_col            & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
  )


  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(j>=1)then

        if(t_soisno(c,j)<tfrz)then
          smp_l(c,j)= hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
        else
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_l(c,j))
        endif

      endif
    enddo
  enddo
  end associate
  end subroutine ALMCalcSmpL


  !------------------------------------------------------------------------
  subroutine ALMBetrReadParams(ncid)
  !
  ! DESCRIPTIONS
  ! read in parameters from namelist
  !
  ! USES
  use ncdio_pio                , only : file_desc_t

  !
  !ARGUMENTS
  !
  type(file_desc_t)                 , intent(inout) :: ncid  ! pio netCDF file id

  call bgc_reaction%readParams(ncid, betrtracer_vars)
  end subroutine ALMBetrReadParams

end module BeTRSimulationALM
