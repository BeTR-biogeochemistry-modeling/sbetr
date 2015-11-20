module clm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  ! Unit Numbers
  !
  use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CL
  implicit none
  integer , parameter, public ::  iundef = -9999999
  real(r8), parameter, public ::  rundef = -9999999._r8
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6
  logical, public :: single_column = .true. !do single column run

  ! Type of run
  integer, public :: nsrest             = iundef

  ! Startup from initial conditions
  integer, public, parameter :: nsrStartup  = 0

  ! Continue from restart files
  integer, public, parameter :: nsrContinue = 1

  integer, public :: spinup_state = 0
end module clm_varctl
