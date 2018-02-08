module summsecaOutType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: summseca_out_type
    real(r8), pointer :: ystates0(:)  => null()
    real(r8), pointer :: ystatesf(:)  => null()
  contains
    procedure, public :: Init
    procedure, private :: InitAllocate

  end type summseca_out_type

  public :: create_summs_out_type
contains
  !-------------------------------------------------------------------------------
  function create_summs_out_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(summseca_out_type), pointer :: create_summs_out_type
    class(summseca_out_type), pointer :: summsout

    allocate(summsout)
    create_summs_out_type => summsout

  end function create_summs_out_type
  !-------------------------------------------------------------------------------

  subroutine Init(this, nstates)
  implicit none
  class(summseca_out_type), intent(inout) :: this
  integer, intent(in) :: nstates


  call this%InitAllocate(nstates)
  end subroutine Init
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, nstates)
  implicit none
  class(summseca_out_type), intent(inout) :: this
  integer, intent(in) :: nstates


  allocate(this%ystates0(nstates))
  allocate(this%ystatesf(nstates))
  end subroutine InitAllocate


end module summsecaOutType
