module BeTR_carbonfluxRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  type, public :: betr_carbonflux_recv_type
    real(r8), pointer :: hr_col(:) => null()
    real(r8), pointer :: hr_vr_col(:,:) => null()
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
  end type betr_carbonflux_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%hr_col(begc:endc))
  allocate(this%hr_vr_col(begc:endc,lbj:ubj))

  end subroutine InitAllocate



  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_carbonflux_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%hr_col(:) = value_column
  this%hr_vr_col(:,:) = value_column

  end subroutine reset
end module BeTR_carbonfluxRecvType
