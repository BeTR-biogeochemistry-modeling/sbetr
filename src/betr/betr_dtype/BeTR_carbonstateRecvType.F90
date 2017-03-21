module BeTR_carbonstateRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  type, public :: betr_carbonstate_recv_type
     real(r8), pointer :: cwdc_col(:) => null()
     real(r8), pointer :: totlitc_col(:)  => null()
     real(r8), pointer :: totsomc_col(:) => null()
     real(r8), pointer :: totlitc_1m_col(:) => null()
     real(r8), pointer :: totsomc_1m_col(:) => null()
 contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
  end type betr_carbonstate_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_carbonstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_carbonstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%cwdc_col(begc:endc))
  allocate(this%totlitc_col(begc:endc))
  allocate(this%totsomc_col(begc:endc))
  allocate(this%totlitc_1m_col(begc:endc))
  allocate(this%totsomc_1m_col(begc:endc))
  end subroutine InitAllocate



  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_carbonstate_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%cwdc_col(:) = value_column
  this%totlitc_col(:) = value_column
  this%totsomc_col(:) = value_column
  this%totlitc_1m_col(:) = value_column
  this%totsomc_1m_col(:) = value_column
  end subroutine reset
end module BeTR_carbonstateRecvType
