module betr_ctrl
!
!logical switches for global control of betr functions
implicit none

  logical, public :: use_pH_data = .false.
  logical, public :: betr_use_cn =.false.
  integer, public :: biulog = 6        ! "stdout" log file unit number, default is 6
  logical, public :: do_betr_otuput = .true.
  integer, parameter, public :: max_betr_hist_type=300
  integer, parameter, public :: max_betr_rest_type=200
  logical, public :: betr_offline = .true.
end module betr_ctrl
