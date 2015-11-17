MODULE shr_sys_mod

   use shr_kind_mod  ! defines real & integer kinds
   use shr_log_mod, only: s_logunit => shr_log_Unit   
   implicit none
! PUBLIC: Public interfaces

   private
   
   public :: shr_sys_abort   ! abort a program
   public :: shr_sys_flush   !flush
   
   contains
   
SUBROUTINE shr_sys_abort(string,rc)

   IMPLICIT none

   character(*)        ,optional :: string  ! error message string
   integer(SHR_KIND_IN),optional :: rc      ! error code

   !----- local -----
   integer(SHR_KIND_IN) :: ierr
   logical              :: flag

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_abort) '
   character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   if (present(string)) then
      if (len_trim(string) > 0) write(s_logunit,F00) 'ERROR: '//trim(string)
   end if
   write(s_logunit,F00) 'WARNING: calling shr_mpi_abort() and stopping'

   stop

END SUBROUTINE shr_sys_abort

SUBROUTINE shr_sys_flush(iulog)
  implicit none
  integer, intent(in) :: iulog
 
  flush(iulog)
  
END SUBROUTINE shr_sys_flush

end module shr_sys_mod    
