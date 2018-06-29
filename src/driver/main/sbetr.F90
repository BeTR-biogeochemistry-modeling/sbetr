program main
  !DESCRIPTION
  !main driver file for standalone betr
  !
  !USES
  use sbetrDriverMod , only : sbetrBGC_driver
  use betr_constants , only : stdout, betr_filename_length, betr_namelist_buffer_size
  use betr_utils     , only : remove_filename_extension, namelist_to_buffer
  implicit none
  integer :: arg_count
  integer :: args
  character(len=betr_filename_length)      :: namelist_filename
  character(len=betr_namelist_buffer_size) :: namelist_buffer
  character(len=betr_filename_length)      :: base_filename


  arg_count = command_argument_count()
  if (arg_count /= 1) then
     write(stdout, '(a, i3)') 'ERROR: must pass exactly one arguement one the command line, received ', arg_count
     call usage()
     call abort()
  end if

  call get_command_argument(1, namelist_filename)
  base_filename = remove_filename_extension(namelist_filename)
  write(stdout, '(a, a)') 'Using base filename for output : ', trim(base_filename)

  write(stdout, '(a, a)') 'Reading namelist filename : ', trim(namelist_filename)
  namelist_buffer = ''
  call namelist_to_buffer(namelist_filename, namelist_buffer)

  call sbetrBGC_driver(base_filename, namelist_buffer)

  ! return correct error code to caller
  call exit(0)

end program main


! ----------------------------------------------------------------------
subroutine usage()
  !DESCRIPTION
  !display something
  use betr_constants, only : stdout
 implicit none
  write(stdout, *) 'sbetr - standalone driver for BeTR reactive transport library.'
  write(stdout, *) 'usage: sbetr namelist_filename'
end subroutine usage
