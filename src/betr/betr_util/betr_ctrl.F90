module betr_ctrl
!
!logical switches for global control of betr functions
implicit none


  logical, public :: use_pH_data = .false.
  logical, public :: is_active_betr_bgc = .false.
  logical, public :: betr_use_cn =.false.
  integer, public :: biulog = 6        ! "stdout" log file unit number, default is 6

end module betr_ctrl
