module ListMod
!
! list types
implicit none
  private

  integer, parameter :: loc_str_len=64
  type, public :: list_s
     character(len=loc_str_len) :: value
     integer :: id
     integer :: itype
     type(list_s), pointer :: next
   end type list_s

   public :: copy_name, list_init, list_insert, list_free, copy_name_type
contains


  !-----------------------------------------------------------------------
  subroutine list_init(self, name, id, itype)
  implicit none
  type(list_s), pointer :: self
  character(len=*), intent(in) :: name
  integer, intent(inout) :: id
  integer, optional, intent(in) :: itype
  allocate(self)
  nullify(self%next)
  id=id+1;
  write(self%value,'(A)')trim(name)
  self%id=id
  if(present(itype))then
    self%itype=itype
  else
    self%itype=0
  endif
  end subroutine list_init
  !-----------------------------------------------------------------------
  subroutine list_insert(self, name, id, itype)

  implicit none
  type(list_s), pointer :: self
  character(len=*), intent(in) :: name
  integer, intent(inout) :: id
  integer, optional, intent(in) :: itype
  type(list_s), pointer :: next

  allocate(next)
  id=id+1
  write(next%value,'(A)')trim(name)
  next%id=id
  if(present(itype))then
    next%itype=itype
  else
    next%itype=0
  endif
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

  !-------------------------------------------------------------------------------
  subroutine copy_name_type(num_names, list_name, vartypes)

  implicit none
  integer, intent(in) :: num_names
  type(list_s), pointer :: list_name
  integer, intent(out) :: vartypes(num_names)

  type(list_s), pointer :: next
  integer :: jj
  next => list_name
  do jj = num_names, 1, -1
    vartypes(jj) = next%itype
    next=>list_next(next)
  enddo
  end subroutine copy_name_type

  !-------------------------------------------------------------------------------
  subroutine list_disp(list)
  implicit none
  type(list_s), pointer :: list
  type(list_s), pointer :: next

  next => list
  do while(associated(next))
    write(*,'(A30,A,I0)')trim(next%value),' =',next%id
    next=>list_next(next)
  enddo

  end subroutine list_disp

end module ListMod
