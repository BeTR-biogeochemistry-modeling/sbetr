module SparseMatMod


  ! !USES:
  ! sparse matrix capability
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
implicit none
  private
  type, public :: spm_list_type
    real(r8) :: val
    integer  :: icol
    integer  :: irow
    type(spm_list_type), pointer :: next => null()
  end type spm_list_type


  type, public :: sparseMat_type
    real(r8), pointer :: val(:)      !nonzero val
    integer , pointer :: icol(:)     !col id for each val
    integer , pointer :: pB(:)       !first nonzero column id (as in val) in each row
    integer , pointer :: ncol(:)     !number of nonzero columns in eacho row
    integer :: szcol                 ! number of columns of the matrix
    integer :: szrow
  contains
    procedure, public :: init
  end type sparseMat_type
  real(r8), parameter :: tiny_val=1.e-14_r8
  public :: spm_list_init, spm_list_insert, spm_list_to_mat
  public :: spm_axpy
contains

!---------------------------------------------------------------
  function create_spm_type()
    !
    ! ! DESCRIPTION:
    !
    ! create an object of type sparseMat_type.
    ! Right now it is purposely empty
   class(sparseMat_type), pointer :: spm
   class(sparseMat_type), pointer ::  create_spm_type
   allocate(spm)
   create_spm_type => spm
  end function create_spm_type

!---------------------------------------------------------------
  subroutine init(this,nelms,nrows,szcol)
  !
  !DESCRIPTION
  ! initialize the sparse matrix
  implicit none
  class(sparseMat_type), intent(inout) :: this
  integer, intent(in) :: nelms
  integer, intent(in) :: nrows
  integer, intent(in) :: szcol

  allocate(this%val(nelms)) ;  this%val(:)  = 0._r8
  allocate(this%icol(nelms));  this%icol(:) = 0
  allocate(this%pB(nrows));    this%pB(:) = 0
  allocate(this%ncol(nrows));  this%ncol(:) = 0
  this%szcol = szcol
  this%szrow= nrows
  end subroutine init

!---------------------------------------------------------------
  subroutine spm_pack(Amat, nelms)

  implicit none
  real(r8), dimension(:,:),intent(in) :: Amat
  integer, intent(in) :: nelms

  class(sparseMat_type), pointer :: spm
  integer :: nrows, szcol
  integer :: ii,jj,kk
  logical :: first
  nrows=size(Amat,1)
  szcol=size(Amat,2)
  allocate(spm, source=create_spm_type())

  call spm%init(nelms, nrows, szcol)

  kk = 0
  do ii = 1, nrows
    first =.true.
    do jj = 1, szcol
      if(abs(amat(ii,jj))>tiny_val)then
        kk = kk + 1
        spm%val(kk) = amat(ii,jj)
        spm%icol(kk)= jj
        if(first)then
          spm%pB(ii) = jj
        endif
        spm%ncol(ii)=spm%ncol(ii)+1
        first=.false.
      endif
    enddo
  enddo

  end subroutine spm_pack

!---------------------------------------------------------------

  subroutine spm_unpack(spm, bmat)

  implicit none
  class(sparseMat_type), intent(in) :: spm
  real(r8), pointer :: bmat(:,:)

  integer :: nrs, ncs
  integer :: i, j, id

  nrs = size(spm%ncol)
  ncs = spm%szcol
  allocate(bmat(nrs,ncs))

  do i = 1, nrs
    do j = 1, spm%ncol(i)
      id=j+spm%pB(i)
      bmat(i,spm%icol(id))=spm%val(id)
    enddo
  enddo

  end subroutine spm_unpack

!---------------------------------------------------------------

  subroutine spm_axpy(nx, ny, a, x, spm, y, errinfo)
  implicit none
  real(r8), intent(in) :: a
  integer , intent(in) :: nx
  integer , intent(in) :: ny
  real(r8), intent(in) :: x(nx)
  type(sparseMat_type), intent(in) :: spm
  real(r8), intent(inout) :: y(ny)
  integer, intent(out) :: errinfo

  integer :: ii, jj, id
  errinfo=0
  if(nx /=spm%szcol .or. ny /= spm%szrow)then
    errinfo=-1
    return
  endif

  do ii = 1, spm%szrow
    do jj = 1, spm%ncol(ii)
      id = spm%pB(ii)+jj-1
      y(ii)=y(ii)+a*spm%val(id)*x(spm%icol(id))
    enddo
  enddo
  end subroutine spm_axpy
!---------------------------------------------------------------


  subroutine spm_list_init(self, val, icol, irow, nelms)
  implicit none
  type(spm_list_type), pointer :: self
  real(r8), intent(in) :: val
  integer, intent(in) :: icol
  integer, intent(in) :: irow
  integer, intent(out):: nelms
  allocate(self)
  nullify(self%next)

  self%val=val
  self%icol=icol
  self%irow=irow
  nelms=1
  end subroutine spm_list_init
!---------------------------------------------------------------

  subroutine spm_list_insert(self, val, icol, irow, nelms)
  implicit none
  type(spm_list_type), pointer :: self
  real(r8), intent(in) :: val
  integer, intent(in) :: icol
  integer, intent(in) :: irow
  integer, intent(inout):: nelms
  type(spm_list_type), pointer :: next

  allocate(next)

  next%val=val
  next%icol=icol
  next%irow=irow
  next%next=> self
  self => next
  nelms=nelms+1
  end subroutine spm_list_insert
!---------------------------------------------------------------

  function spm_list_next(self)result(next)

  implicit none
  type(spm_list_type), pointer :: self
  type(spm_list_type), pointer :: next

  next => self%next
  end function spm_list_next
!---------------------------------------------------------------

  subroutine spm_list_free(self)
  implicit none
  type(spm_list_type), pointer :: self
  type(spm_list_type), pointer :: current
  type(spm_list_type), pointer :: elem

  elem => self
  do while(associated(elem))
    current => elem
    elem => current%next
    deallocate(current)
  enddo
  end subroutine spm_list_free
!---------------------------------------------------------------

  subroutine spm_list_to_mat(spm_list, spm, nelms, szcol)
  !
  !copy the list to sparse matrix, and free the list
  implicit none
  type(spm_list_type), pointer :: spm_list
  class(sparseMat_type), pointer :: spm
  integer, intent(in) :: nelms    !total number of nonzero elements
  integer, intent(in) :: szcol    !total columns of the matrix

  type(spm_list_type), pointer :: next
  integer :: nrows
  integer :: irow, jj

  nrows=spm_list%irow

  allocate(spm, source=create_spm_type())

  call spm%init(nelms, nrows, szcol)

  jj = nelms
  irow=nrows
  next => spm_list
  do while(associated(next))
    if (irow /= next%irow)then
       spm%pB(irow)=jj+1
       irow=next%irow
    endif
    spm%val(jj) = next%val
    spm%icol(jj)= next%icol
    spm%ncol(irow)=spm%ncol(irow)+1
    next=>spm_list_next(next)
    jj=jj-1
  enddo
  !fill the last row
  spm%pB(irow)=jj+1
  call spm_list_free(spm_list)
  end subroutine spm_list_to_mat
end module SparseMatMod
