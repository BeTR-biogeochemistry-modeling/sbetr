module TracerBaseType
use betr_ctrl, only : max_betr_hist_type
implicit none

  type, public :: tracerbase_type
     integer :: num_hist1d
     integer :: num_hist2d
     character(len=255) :: nmlist_hist1d_buffer(max_betr_hist_type)
     character(len=255) :: nmlist_hist2d_buffer(max_betr_hist_type)
  contains
     procedure, public :: tracer_base_init
     procedure, public :: add_hist_var2d
     procedure, public :: add_hist_var1d
  end type tracerbase_type
contains

  !-----------------------------------------------------------------------
  subroutine tracer_base_init(this)
  implicit none
  class(tracerbase_type) :: this

  !set number of history files to zero
  this%num_hist1d = 0
  this%num_hist2d = 0
  end subroutine tracer_base_init
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine add_hist_var2d(this, fname, units, type2d,  &
     avgflag, long_name, default)
  !
  !DESCRIPTION
  !build the namelist for an output 2d variable
  implicit none
  class(tracerbase_type) :: this
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: type2d
  character(len=*), intent(in) :: avgflag
  character(len=*), intent(in) :: long_name
  character(len=*),optional, intent(in) :: default

  character(len=20) :: default_loc = "active"

  character(len=1)  :: quote = ''''

  !use the following format to read
  !x namelist /hist2d_fmt/    &
  !x fname, units, avgflag,type2d,long_name, default

  if(present(default))then
    default_loc=trim(default)
  endif
  this%num_hist2d = this%num_hist2d  + 1

  this%nmlist_hist2d_buffer(this%num_hist2d)='&hist2d_fmt'//new_line('A') &
    //' fname='//quote//trim(fname)//quote//new_line('A') &
    //' units='//quote//trim(units)//quote//new_line('A') &
    //' avgflag='//quote//trim(avgflag)//quote//new_line('A') &
    //' type2d='//quote//trim(type2d)//quote//new_line('A') &
    //' long_name='//quote//trim(long_name)//quote//new_line('A') &
    //' default='//quote//trim(default_loc)//quote//new_line('A') &
    //'/'
  end subroutine add_hist_var2d

  !-----------------------------------------------------------------------
  subroutine add_hist_var1d(this, fname, units,  &
     avgflag, long_name, default)
  !
  !DESCRIPTION
  !build the namelist for an output 2d variable
  implicit none
  class(tracerbase_type) :: this
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: avgflag
  character(len=*), intent(in) :: long_name
  character(len=*),optional, intent(in) :: default

  character(len=20) :: default_loc = "active"

  character(len=1)  :: quote = ''''

  !use the following format to read
  !x namelist /hist1d_fmt/    &
  !x fname, units, avgflag,type2d,long_name, default

  if(present(default))then
    default_loc=trim(default)
  endif
  this%num_hist1d = this%num_hist1d  + 1

  this%nmlist_hist1d_buffer(this%num_hist1d)='&hist1d_fmt'//new_line('A') &
    //' fname='//quote//trim(fname)//quote//new_line('A') &
    //' units='//quote//trim(units)//quote//new_line('A') &
    //' avgflag='//quote//trim(avgflag)//quote//new_line('A') &
    //' long_name='//quote//trim(long_name)//quote//new_line('A') &
    //' default='//quote//trim(default_loc)//quote//new_line('A') &
    //'/'
  end subroutine add_hist_var1d

end module TracerBaseType
