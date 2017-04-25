module TracerBaseType
use betr_ctrl, only : max_betr_hist_type
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  integer, parameter :: loc_str_len=255
  type, private :: list_t
    character(len=loc_str_len)  :: name
    type(list_t), pointer :: next => null()
  end type list_t

  type, public :: tracerbase_type
     integer :: num_hist1d
     integer :: num_hist2d
     character(len=loc_str_len), allocatable :: nmlist_hist1d_buffer(:)
     character(len=loc_str_len), allocatable :: nmlist_hist2d_buffer(:)
  contains
     procedure, public :: tracer_base_init
     procedure, public :: add_hist_var2d
     procedure, public :: add_hist_var1d
     procedure, public :: sort_hist_list
  end type tracerbase_type

  type(list_t), pointer :: head1d => null()
  type(list_t), pointer :: tail1d => null()
  type(list_t), pointer :: head2d => null()
  type(list_t), pointer :: tail2d => null()


contains
  !-----------------------------------------------------------------------
  subroutine list_init(self, name)
  implicit none
  type(list_t), pointer :: self
  character(len=loc_str_len), intent(in) :: name

  allocate(self)
  nullify(self%next)
  self%name=name

  end subroutine list_init
  !-----------------------------------------------------------------------
  subroutine list_insert(self, name)

  implicit none
  type(list_t), pointer :: self
  character(len=loc_str_len), intent(in) :: name
  type(list_t), pointer :: next

  allocate(next)
  next%name=name
  next%next=> self
  self => next

  end subroutine list_insert
  !-----------------------------------------------------------------------
  subroutine list_free(self)
  implicit none
  type(list_t), pointer :: self
  type(list_t), pointer :: current
  type(list_t), pointer :: elem

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
  type(list_t), pointer :: self
  type(list_t), pointer :: next

  next => self%next
  end function list_next

  !-----------------------------------------------------------------------
  subroutine tracer_base_init(this)
  implicit none
  class(tracerbase_type), intent(inout) :: this

  !set number of history files to zero
  this%num_hist1d = 0
  this%num_hist2d = 0
  end subroutine tracer_base_init

  !-----------------------------------------------------------------------
  subroutine add_hist_var2d(this, fname, units, type2d,  &
     avgflag, long_name, default)
  !
  !DESCRIPTION
  !build the namelist for an output 2d variable
  implicit none
  class(tracerbase_type), intent(inout) :: this
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: type2d
  character(len=*), intent(in) :: avgflag
  character(len=*), intent(in) :: long_name
  character(len=*),optional, intent(in) :: default

  character(len=20) :: default_loc = "active"
  character(len=loc_str_len):: tempstr
  character(len=1)  :: quote = ''''

  !use the following format to read
  !x namelist /hist2d_fmt/    &
  !x fname, units, avgflag,type2d,long_name, default

  if(present(default))then
    default_loc=trim(default)
  endif
  this%num_hist2d = this%num_hist2d  + 1

  tempstr='&hist2d_fmt'//new_line('A') &
    //' fname='//quote//trim(fname)//quote//new_line('A') &
    //' units='//quote//trim(units)//quote//new_line('A') &
    //' avgflag='//quote//trim(avgflag)//quote//new_line('A') &
    //' type2d='//quote//trim(type2d)//quote//new_line('A') &
    //' long_name='//quote//trim(long_name)//quote//new_line('A') &
    //' default='//quote//trim(default_loc)//quote//new_line('A') &
    //'/'
  if(this%num_hist2d==1)then
    call list_init(head2d, tempstr)
  else
    call list_insert(head2d, tempstr)
  endif

  end subroutine add_hist_var2d

  !-----------------------------------------------------------------------
  subroutine add_hist_var1d(this, fname, units,  &
     avgflag, long_name, default)
  !
  !DESCRIPTION
  !build the namelist for an output 2d variable
  implicit none
  class(tracerbase_type), intent(inout) :: this
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: avgflag
  character(len=*), intent(in) :: long_name
  character(len=*),optional, intent(in) :: default

  character(len=20) :: default_loc = "active"

  character(len=1)  :: quote = ''''
  character(len=loc_str_len):: tempstr
  !use the following format to read
  !x namelist /hist1d_fmt/    &
  !x fname, units, avgflag,type2d,long_name, default

  if(present(default))then
    default_loc=trim(default)
  endif
  this%num_hist1d = this%num_hist1d  + 1

  tempstr='&hist1d_fmt'//new_line('A') &
    //' fname='//quote//trim(fname)//quote//new_line('A') &
    //' units='//quote//trim(units)//quote//new_line('A') &
    //' avgflag='//quote//trim(avgflag)//quote//new_line('A') &
    //' long_name='//quote//trim(long_name)//quote//new_line('A') &
    //' default='//quote//trim(default_loc)//quote//new_line('A') &
    //'/'

  if(this%num_hist1d==1)then
    call list_init(head1d, tempstr)
  else
    call list_insert(head1d, tempstr)
  endif

  end subroutine add_hist_var1d
  !-----------------------------------------------------------------------
  subroutine sort_hist_list(this)

  implicit none
  class(tracerbase_type), intent(inout) :: this
  type(list_t), pointer :: temp_node, temp1_node, next
  integer :: jj

  !allocate memories
  allocate(this%nmlist_hist1d_buffer(this%num_hist1d))
  allocate(this%nmlist_hist2d_buffer(this%num_hist2d))

  !flip 1d history fields
  call list_init(tail1d, head1d%name)

  temp_node=>list_next(head1d)
  do while(associated(temp_node))
    call list_insert(tail1d, temp_node%name)
    if(associated(temp_node%next))then
      temp_node=list_next(temp_node)
    else
      exit
    endif
  enddo

  jj = 0
  temp1_node=>tail1d
  do while(associated(temp1_node))
    jj = jj + 1
    this%nmlist_hist1d_buffer(jj)=temp1_node%name
    if(associated(temp1_node%next))then
      temp1_node=list_next(temp1_node)
    else
      exit
    endif
  enddo

  !flip 2d history fields
  call list_init(tail2d, head2d%name)

  temp_node=>head2d%next
  do while(associated(temp_node))
    call list_insert(tail2d, temp_node%name)
    if(associated(temp_node%next))then
      temp_node=> list_next(temp_node)
    else
      exit
    endif
  enddo

  jj = 0
  temp1_node=>tail2d
  do while(associated(temp1_node))
    jj = jj + 1

    this%nmlist_hist2d_buffer(jj)=temp1_node%name
    if(associated(temp1_node%next))then
      temp1_node=> list_next(temp1_node)
    else
      exit
    endif
  enddo

!x  if(jj/=this%num_hist2d)then
!x    print*,'bang'
!x  endif
!x  print*,'free lists'
  call list_free(tail1d)
  call list_free(tail2d)
  call list_free(head1d)
  call list_free(head2d)
  end subroutine sort_hist_list

end module TracerBaseType
