module BeTR_nitrogenstateRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_nitrogenstate_recv_type
    real(r8), pointer :: cwdn_col(:) => null()
    real(r8), pointer :: totlitn_col(:) => null()
    real(r8), pointer :: totsomn_col(:) => null()
    real(r8), pointer :: sminn_col(:) => null()
    real(r8), pointer :: totlitn_1m_col(:) => null()
    real(r8), pointer :: totsomn_1m_col(:) => null()
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
  end type betr_nitrogenstate_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_nitrogenstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_nitrogenstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%cwdn_col(begc:endc))
  allocate(this%totlitn_col(begc:endc))
  allocate(this%totsomn_col(begc:endc))
  allocate(this%sminn_col(begc:endc))
  allocate(this%totlitn_1m_col(begc:endc))
  allocate(this%totsomn_1m_col(begc:endc))
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_nitrogenstate_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%cwdn_col(:) = value_column
  this%totlitn_col(:) = value_column
  this%totsomn_col(:) = value_column
  this%sminn_col(:) = value_column
  this%totlitn_1m_col(:) = value_column
  this%totsomn_1m_col(:) = value_column
  end subroutine reset
end module BeTR_nitrogenstateRecvType
