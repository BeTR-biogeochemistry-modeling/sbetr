module ListMod
!
! list types
implicit none
  private

  integer, parameter :: loc_str_len=16
  type, public :: list_s
     character(len=loc_str_len) :: value
     integer :: id
     type(list_s), pointer :: next
   end type list_s

   public :: copy_name, list_init, list_insert, list_free
contains


  !-----------------------------------------------------------------------
  subroutine list_init(self, name, id)
  implicit none
  type(list_s), pointer :: self
  character(len=*), intent(in) :: name
  integer, intent(inout) :: id

  allocate(self)
  nullify(self%next)
  id=id+1;
  write(self%value,'(A)')trim(name)
  self%id=id
  end subroutine list_init
  !-----------------------------------------------------------------------
  subroutine list_insert(self, name, id)

  implicit none
  type(list_s), pointer :: self
  character(len=*), intent(in) :: name
  integer, intent(inout) :: id
  type(list_s), pointer :: next

  allocate(next)
  id=id+1
  write(next%value,'(A)')trim(name)
  next%id=id
  next%next=> self
  self => next

  end subroutine list_insert
  !-----------------------------------------------------------------------
  subroutine list_free(self)
  implicit none
  type(list_s), pointer :: self
  type(list_s), pointer :: current
  type(list_s), pointer :: elem

  elem => self
  do while(associated(elem))
    current => elem
    elem => current%next
    deallocate(current)
  enddo
  end subroutine list_free
  !-----------------------------------------------------------------------
  function list_next(self)result(next)

  implicit none
  type(list_s), pointer :: self
  type(list_s), pointer :: next

  next => self%next
  end function list_next

  !-------------------------------------------------------------------------------
  subroutine copy_name(num_names, list_name, outnames)

  implicit none
  integer, intent(in) :: num_names
  type(list_s), pointer :: list_name
  character(len=*), intent(out) :: outnames(num_names)

  type(list_s), pointer :: next
  integer :: jj
  next => list_name
  do jj = num_names, 1, -1
    write(outnames(jj),'(A)')trim(next%value)
    next=>list_next(next)
  enddo
  end subroutine copy_name

end module ListMod
