module TracerBaseType
use betr_ctrl, only : max_betr_hist_type
implicit none



  type, private :: node
    character(len=255)  :: name
    type(node), pointer :: next
  end type node

  type, public :: tracerbase_type
     integer :: num_hist1d
     integer :: num_hist2d
     character(len=255), allocatable :: nmlist_hist1d_buffer(:)
     character(len=255), allocatable :: nmlist_hist2d_buffer(:)
  contains
     procedure, public :: tracer_base_init
     procedure, public :: add_hist_var2d
     procedure, public :: add_hist_var1d
     procedure, public :: sort_hist_list
  end type tracerbase_type

  type(node), private, pointer :: head1d, tail1d
  type(node), private, pointer :: head2d, tail2d


contains

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
  character(len=255):: tempstr
  type(node), pointer :: temp_node
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
    allocate(head2d)
    head2d%name=tempstr
    nullify(head2d%next)
  else
    allocate(temp_node)
    temp_node%name=tempstr
    temp_node%next=>head2d
    head2d=>temp_node
    nullify(temp_node)
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
  character(len=255):: tempstr
  type(node), pointer :: temp_node
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
    allocate(head1d)
    head1d%name=tempstr
    nullify(head1d%next)
  else
    allocate(temp_node)
    temp_node%name=tempstr
    temp_node%next=>head1d
    head1d=>temp_node
    nullify(temp_node)
  endif

  end subroutine add_hist_var1d
  !-----------------------------------------------------------------------
  subroutine sort_hist_list(this)

  implicit none
  class(tracerbase_type), intent(inout) :: this
  type(node), pointer :: temp_node, temp1_node
  integer :: jj

  !flip 1d history fields
  allocate(tail1d)
  tail1d%name=head1d%name
  nullify(tail1d%next)

  temp_node=>head1d%next
  do while(associated(temp_node))
    allocate(temp1_node)
    temp1_node%name=temp_node%name
    temp1_node%next=>tail1d
    tail1d => temp1_node
    temp_node=> temp_node%next
  enddo

  allocate(this%nmlist_hist1d_buffer(this%num_hist1d))
  jj = 0
  nullify(temp1_node)
  temp1_node=>tail1d
  do while(associated(temp1_node))
    jj = jj + 1
    this%nmlist_hist1d_buffer(jj)=temp1_node%name
    temp1_node=> temp1_node%next
  enddo
!x  if(jj/=this%num_hist1d)then
!x    print*,'bang'
!x  endif

  !flip 2d history fields
  allocate(tail2d)
  tail2d%name=head2d%name
  nullify(tail2d%next)

  temp_node=>head2d%next
  do while(associated(temp_node))
    allocate(temp1_node)
    temp1_node%name=temp_node%name
    temp1_node%next=>tail2d
    tail2d => temp1_node
    temp_node=> temp_node%next
  enddo

  allocate(this%nmlist_hist2d_buffer(this%num_hist2d))
  jj = 0
  nullify(temp1_node)
  temp1_node=>tail2d
  do while(associated(temp1_node))
    jj = jj + 1
    this%nmlist_hist2d_buffer(jj)=temp1_node%name
    temp1_node=> temp1_node%next
  enddo
  nullify(temp1_node)

!x  if(jj/=this%num_hist2d)then
!x    print*,'bang'
!x  endif
  nullify(tail1d)
  nullify(tail2d)
  nullify(head1d)
  nullify(head2d)
  end subroutine sort_hist_list

end module TracerBaseType
