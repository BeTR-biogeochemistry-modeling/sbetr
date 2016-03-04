module BeTRSimulationStandalone
  !
  ! !DESCRIPTION:
  !  factory to load the specific betr simulator
  !
  ! !USES:
  !
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_log_mod, only : errMsg => shr_log_errMsg

  use decompMod, only : bounds_type
  use BeTRSimulation, only : betr_simulation_type
  use BeTR_CNStateType, only : betr_cnstate_type

  use EcophysConType, only : ecophyscon_type

  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_standalone_type
     type(betr_cnstate_type) :: betr_cnstate_vars
     type(ecophyscon_type) :: ecophyscon


   contains
     procedure :: Init => StandaloneInit
     procedure, public :: StepWithoutDrainage => StandaloneStepWithoutDrainage
     procedure, public :: StepWithDrainage => StandaloneStepWithDrainage
  end type betr_simulation_standalone_type

  public :: create_betr_simulation_standalone
  
contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_standalone()

    implicit none

    class(betr_simulation_standalone_type), pointer :: create_betr_simulation_standalone
    class(betr_simulation_standalone_type), pointer :: simulation

    allocate(simulation)

    create_betr_simulation_standalone => simulation

  end function create_betr_simulation_standalone

  !-------------------------------------------------------------------------------

  subroutine StandaloneInit(this, reaction_method, bounds, lbj, ubj, waterstate)

    use BeTRSimulation, only : BeTRSimulationInit
    use ReactionsFactoryStandalone, only : create_standalone_bgc_reaction_type, &
         create_standalone_plant_soilbgc_type

    use BeTR_PatchType, only : betr_pft
    use BeTR_pftvarconType, only : betr_pftvarcon
    use PatchType, only : pft
    use pftvarcon, only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use clm_instMod, only : cnstate_vars
    use WaterStateType, only : waterstate_type
    
    implicit none
    class(betr_simulation_standalone_type) :: this
    character(len=*), intent(in) :: reaction_method
    type(bounds_type)    , intent(in) :: bounds
    integer              , intent(in) :: lbj, ubj
    type(waterstate_type), intent(inout) :: waterstate
    
    
    betr_pft%wtcol                        => pft%wtcol
    betr_pft%column                       => pft%column
    betr_pft%itype                        => pft%itype
    betr_pft%landunit                     => pft%landunit
    
    ! allocate the reaction types that may only be known to this
    ! simulation type.
    allocate(this%bgc_reaction, source=create_standalone_bgc_reaction_type(reaction_method))
    allocate(this%plant_soilbgc, source=create_standalone_plant_soilbgc_type(reaction_method))

    ! now call the base simulation init to continue initialization
    ! FIXME(bja, 2016-03) missing water state vars!
    call BeTRSimulationInit(this, reaction_method, bounds, lbj, ubj, waterstate)

    !pass necessary data
    this%betr_cnstate_vars%isoilorder  => cnstate_vars%isoilorder
    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    call this%bgc_reaction%init_betr_lsm_bgc_coupler(bounds, this%plant_soilbgc, &
         this%betrtracer_vars, this%tracerstate_vars, this%betr_cnstate_vars, &
         this%ecophyscon)

  end subroutine StandaloneInit

  
  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithoutDrainage(this, bounds, lbj, ubj, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)

    use BetrBGCMod, only : run_betr_one_step_without_drainage
    use SoilStateType, only : soilstate_type
    use WaterStateType, only : Waterstate_Type
    use TemperatureType, only : temperature_type
    use ChemStateType, only : chemstate_type
    use WaterfluxType, only : waterflux_type
    use ColumnType, only : column_type
    use BGCReactionsMod, only : bgc_reaction_type
    use atm2lndType, only : atm2lnd_type
    use SoilHydrologyType, only : soilhydrology_type
    use BeTR_CarbonFluxType, only : betr_carbonflux_type
    use CNStateType, only : cnstate_type
    use CNCarbonFluxType, only : carbonflux_type
    use CanopyStateType, only : canopystate_type
    use BeTR_PatchType, only : betr_pft
    use PatchType, only : pft
    use pftvarcon, only : crop

    implicit none

    class(betr_simulation_standalone_type) :: this

    ! !ARGUMENTS :
    type(bounds_type), intent(in) :: bounds ! bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: num_soilp
    integer, intent(in) :: filter_soilp(:) ! pft filter
    integer, intent(in) :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    
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
    type(waterflux_type), intent(inout) :: waterflux_vars

    !temporary variables
    type(betr_carbonflux_type) :: betr_carbonflux_vars

    !pass necessary data for correct subroutine call
    
    this%betr_cnstate_vars%isoilorder          => cnstate_vars%isoilorder
    
    betr_carbonflux_vars%annsum_npp_patch => carbonflux_vars%annsum_npp_patch
    betr_carbonflux_vars%agnpp_patch      => carbonflux_vars%agnpp_patch
    betr_carbonflux_vars%bgnpp_patch      => carbonflux_vars%bgnpp_patch
    
    betr_pft%wtcol                        => pft%wtcol
    betr_pft%column                       => pft%column
    betr_pft%itype                        => pft%itype
    betr_pft%landunit                     => pft%landunit
    betr_pft%crop                         => crop
    
    call run_betr_one_step_without_drainage(bounds, lbj, ubj, &
         num_soilc, filter_soilc, num_soilp, filter_soilp, col,   &
         atm2lnd_vars, soilhydrology_vars, soilstate_vars, &
         waterstate_vars, temperature_vars, waterflux_vars, &
         chemstate_vars, this%betr_cnstate_vars, canopystate_vars, &
         betr_carbonflux_vars, this%betrtracer_vars, this%bgc_reaction, &
         this%betr_aerecond_vars, this%tracerboundarycond_vars, this%tracercoeff_vars, &
         this%tracerstate_vars, this%tracerflux_vars, this%plant_soilbgc)


    
  end subroutine StandaloneStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithDrainage(this, bounds, lbj, ubj, &
       num_soilc, filter_soilc, jtops, waterflux_vars, col)

    use BetrBGCMod, only : run_betr_one_step_with_drainage
    use ColumnType, only : column_type
    use MathfuncMod, only : safe_div
    use WaterFluxType, only : waterflux_type

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
    
    call run_betr_one_step_with_drainage(bounds, lbj, ubj, &
         num_soilc, filter_soilc, &
         jtops, waterflux_vars, col, this%betrtracer_vars, this%tracercoeff_vars, &
         this%tracerstate_vars,  this%tracerflux_vars)

  
  end subroutine StandaloneStepWithDrainage


end module BeTRSimulationStandalone
