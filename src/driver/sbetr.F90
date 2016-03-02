program main

  use decompMod            , only : bounds_type
  use clmgridMod           , only : init_clm_vertgrid
  use clm_varpar           , only : nlevgrnd
  use clm_initializeMod    , only : initialize
  use sbetrDriverMod, only : time_type
  use sbetrDriverMod, only : sbetrBGC_driver

  use betr_constants, only : betr_filename_length, betr_namelist_buffer_size
  
  
  implicit none

  type(bounds_type) :: bounds
  type(time_type) :: ttime
  integer :: numf
  integer, pointer :: filter(:)

  integer :: arg_count
  integer :: args
  character(len=betr_filename_length) :: namelist_filename
  character(len=betr_namelist_buffer_size) :: namelist_buffer
  

  arg_count = command_argument_count()
  if (arg_count /= 1) then
     write(*, *) 'ERROR: must pass exactly one arguement one the command line, received ', arg_count
     call usage()
     call abort()
  end if

  call get_command_argument(1, namelist_filename)
  write(*, *) 'Reading namelist filename : ', trim(namelist_filename)
  namelist_buffer = ''
  call namelist_to_buffer(namelist_filename, namelist_buffer)
  
  !set up mask
  bounds%begc = 1
  bounds%endc = 1
  bounds%begp = 1
  bounds%endp = 1
  bounds%begl = 1
  bounds%endl = 1
  bounds%lbj  = 1
  bounds%ubj  = nlevgrnd

  numf = 1
  ttime%time_end = 1800*48*10
  ttime%restart_dtime=1800*2
  allocate(filter(numf))
  filter(:) = 1

  !set up grid
  call init_clm_vertgrid(nlevgrnd)

  call initialize(bounds)

  call sbetrBGC_driver(namelist_buffer, bounds, numf, filter, ttime)

end program main


subroutine usage()
  use betr_constants, only : stdout
  
  write(stdout, *) 'sbetr - standalone driver for BeTR reactive transport library.'
  write(stdout, *) 'usage: sbetr namelist_filename'
end subroutine usage

subroutine namelist_to_buffer(namelist_filename, namelist_buffer)

  use betr_constants, only : betr_string_length_long, betr_namelist_buffer_size, stdout
  
  character(len=*), intent(in) :: namelist_filename
  character(len=betr_namelist_buffer_size), intent(out) :: namelist_buffer

  character(len=*), parameter :: subname = 'namelist_to_buffer'
  character(len=betr_string_length_long) :: ioerror_msg
  integer :: nml_unit, nml_error

  nml_unit = 16
  
  ! read the namelist file into a buffer.
  open(unit=nml_unit, file=trim(namelist_filename), action='read', access='stream', &
       form='unformatted', iostat=nml_error)
  if (nml_error == 0) then
     read(unit=nml_unit, iostat=nml_error, iomsg=ioerror_msg) namelist_buffer

     ! we should always reach the EOF to capture the entire file...
     if (.not. is_iostat_end(nml_error)) then
        write(stdout, '(a, a, i8)') subname, &
             ": IO ERROR reading namelist file into buffer: ", nml_error
        write(stdout, '(a)') ioerror_msg
        call abort()
     else
        write(stdout, '(a, a, a)') "Read '", trim(namelist_filename), "' until EOF."
     end if

     write(stdout, '(a, a, i7, a)') subname, ": Read buffer of ", &
          len_trim(namelist_buffer), " characters."

     write(stdout, '(a)') "  If it looks like part of the namelist is missing, "
     write(stdout, '(a)') "  compare the number of characters read to the actual "
     write(stdout, '(a,a,a)') "  size of your file ($ wc -c ", trim(namelist_filename), ") and increase "
     write(stdout, '(a)') "  the buffer size if necessary."
     write(stdout, '(a)') "------------------------------"
     write(stdout, '(a)') trim(namelist_buffer)
     write(stdout, '(a)') "------------------------------"
  else
     write(stdout, '(a, a, i8, a, a)') subname, ": IO ERROR ", nml_error, &
          " opening namelist file : ", trim(namelist_filename)
     call abort()
  end if
  close(nml_unit)


end subroutine namelist_to_buffer
