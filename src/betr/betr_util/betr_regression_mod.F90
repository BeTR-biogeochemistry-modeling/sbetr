module betr_regression_module

  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  
  type, public :: betr_regression_type
     integer, private :: num_cells
   contains
     procedure, public :: Init
     procedure, public :: WriteOutput 
     procedure, private :: ReadNamelist 
  end type betr_regression_type

contains

  subroutine Init(this, namelist_buffer)

    use betr_constants, only : betr_namelist_buffer_size

    implicit none

    class(betr_regression_type), intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    call  this%ReadNamelist(namelist_buffer)
  end subroutine Init
    
  !--------------------------------------------------------------------
  
  subroutine ReadNamelist(this, namelist_buffer)
    use babortutils, only : endrun
    use bshr_log_mod, only : errMsg => shr_log_errMsg
    
    use betr_constants, only : stdout, betr_string_length_long, betr_namelist_buffer_size

    implicit none

    class(betr_regression_type), intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer

    character(len=*), parameter :: subname = 'betr_regression:ReadNamelist'
    ! !LOCAL VARIABLES:
    integer :: nml_error
    integer :: cells
    character(len=betr_string_length_long) :: ioerror_msg
    logical :: write_regression_output
    !-----------------------------------------------------------------------

    namelist / regression_test / cells, write_regression_output

    cells = 0
    write_regression_output = .false.

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=regression_test, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading betr_regression_test namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr regression test type :'
       write(stdout, *)
       write(stdout, *) ' regression_test namelist settings :'
       write(stdout, *)
       write(stdout, regression_test)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    this%num_cells = cells
  end subroutine ReadNamelist

  !---------------------------------------------------------------------------------

  subroutine WriteOutput(this)

    use betr_constants, only : stdout
    
    implicit none

    class(betr_regression_type), intent(inout) :: this

    write(stdout, *) 'Writing regression output.'
    
  end subroutine WriteOutput

  
end module betr_regression_module
