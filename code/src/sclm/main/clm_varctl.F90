module clm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  ! Unit Numbers
  !
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6
  logical, public :: single_column = .true. !do single column run  
  
end module clm_varctl  
