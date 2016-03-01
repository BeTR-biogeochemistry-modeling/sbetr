module BeTRSimulationALM
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
  character(len=*), private, parameter :: filename = "__FILE__"

  type, public, extends(betr_simulation_type) :: betr_simulation_alm_type

   contains
     
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

  
  
end module BeTRSimulationALM
