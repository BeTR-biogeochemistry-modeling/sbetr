module BeTR_phosphorusstateRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  type, public :: betr_phosphorusstate_recv_type
    real(r8), pointer :: cwdp_col(:) => null()
    real(r8), pointer :: totlitp_col(:) => null()
    real(r8), pointer :: totsomp_col(:) => null()
    real(r8), pointer :: totlitp_1m_col(:) => null()
    real(r8), pointer :: totsomp_1m_col(:) => null()
    real(r8), pointer :: solutionp_col(:) => null()
    real(r8), pointer :: labilep_col(:) => null()
    real(r8), pointer :: secondp_col(:) => null()
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
  end type betr_phosphorusstate_recv_type

 contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_phosphorusstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_phosphorusstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%cwdp_col(begc:endc))
  allocate(this%totlitp_col(begc:endc))
  allocate(this%totsomp_col(begc:endc))
  allocate(this%totlitp_1m_col(begc:endc))
  allocate(this%totsomp_1m_col(begc:endc))
  allocate(this%solutionp_col(begc:endc))
  allocate(this%labilep_col(begc:endc))
  allocate(this%secondp_col(begc:endc))
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_phosphorusstate_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%cwdp_col(:) = value_column
  this%totlitp_col(:) = value_column
  this%totsomp_col(:) = value_column
  this%totlitp_1m_col(:) = value_column
  this%totsomp_1m_col(:) = value_column
  this%solutionp_col(:) = value_column
  this%labilep_col(:) = value_column
  this%secondp_col(:) = value_column
  end subroutine reset
end module BeTR_phosphorusstateRecvType
