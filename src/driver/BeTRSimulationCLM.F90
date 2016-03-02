module BeTRSimulationCLM
  !
  ! !DESCRIPTION:
  !  factory to load the specific betr simulator
  !
  ! !USES:
  !
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_log_mod, only : errMsg => shr_log_errMsg
  
  use BeTRSimulation, only : betr_simulator_type
  
  implicit none

  private
  character(len=*), private, parameter :: mod_filename = __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_clm_type

   contains
     
  end type betr_simulation_clm_type

  public :: create_betr_simulation_clm
  
contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_clm()
    implicit none
    class(betr_simulation_clm_type), pointer :: create_betr_simulation_clm
    class(betr_simulation_clm_type), pointer :: simulation
    allocate(simulation)
    call simulation%Init()
    create_betr_simulation_clm => simulation
  end function create_betr_simulation_clm

  
  
end module BeTRSimulationCLM
