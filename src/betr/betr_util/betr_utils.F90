module betr_utils

  implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
contains
  ! ----------------------------------------------------------------------
  function remove_filename_extension(filename) result(basename)
    !
    ! Remove the extension from a file name to get a base filename.
    !
    ! We start at the end of the filename and assume that the extension
    ! is marked by a period.
    !

    use betr_constants, only : betr_filename_length

    implicit none

    character(len=betr_filename_length), intent(in) :: filename
    character(len=betr_filename_length) :: basename
    integer :: ext_index

    ext_index = scan(filename, '.', .true.)
    if (ext_index == 0) then
       ! no period marking an extension...
       ext_index = len(trim(filename)) + 1

    end if
    basename = filename(1:ext_index-1)
  end function remove_filename_extension

  ! ----------------------------------------------------------------------
  function num2str(num, fmt)result(ans)
  use betr_constants, only : betr_string_length_long

  implicit none
  integer :: num
  character(len=*), intent(in) :: fmt
  character(len=betr_string_length_long) :: ans

  write(ans,fmt)num
  return
  end function num2str
  ! ----------------------------------------------------------------------
  function log2str(logval)result(str)
  implicit none
  logical, intent(in) :: logval
  character(len=8) :: str

  if(logval)then
    str='.true.'
  else
    str='.false.'
  endif
  end function log2str


  ! ----------------------------------------------------------------------

  subroutine namelist_to_buffer(namelist_filename, namelist_buffer)
  !DESCRIPTION
  !read in namelist
  !USES
  use betr_constants, only : betr_string_length_long, betr_namelist_buffer_size, stdout
  implicit none
  character(len=*)                         , intent(in)  :: namelist_filename
  character(len=betr_namelist_buffer_size) , intent(out) :: namelist_buffer

  character(len=*), parameter                            :: subname = 'namelist_to_buffer'
  character(len=betr_string_length_long)                 :: ioerror_msg
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
end module betr_utils
