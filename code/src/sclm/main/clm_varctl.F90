module clm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  ! Unit Numbers
  !
  use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CL
  implicit none
  private
  save
  integer , parameter, public ::  iundef = -9999999
  real(r8), parameter, public ::  rundef = -9999999._r8
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6
  logical, public :: single_column = .true. !do single column run
  logical, public :: use_cn              = .false.  

  ! Type of run
  integer, public :: nsrest             = iundef

  ! Startup from initial conditions
  integer, public, parameter :: nsrStartup  = 0

  ! Continue from restart files
  integer, public, parameter :: nsrContinue = 1

  integer, public :: spinup_state = 0

  logical :: carbon_only  = .false.

  public :: CNAllocate_Carbon_only
  public :: cnallocate_carbon_only_set
 contains

  ! Set module carbon_only flag
  subroutine cnallocate_carbon_only_set(carbon_only_in)
    logical, optional, intent(in) :: carbon_only_in

    if(present(carbon_only_in))then
      carbon_only = carbon_only_in
    else
      carbon_only = .false.
    endif
  end subroutine cnallocate_carbon_only_set

  ! Get module carbon_only flag
  logical function CNAllocate_Carbon_only()
    cnallocate_carbon_only = carbon_only
  end function CNAllocate_Carbon_only
end module clm_varctl
