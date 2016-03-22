module betr_utils

  implicit none

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
     
end module betr_utils
   
