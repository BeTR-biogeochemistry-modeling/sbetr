module BeTR_biogeoStateType

  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use betr_decompMod  , only : betr_bounds_type
implicit none

  type betr_biogeo_state_type
    real(r8), pointer :: zwts_col           (:)   => null() ! the shallower between zwt_perch and zwt
    contains
      procedure, public  :: Init
      procedure, private :: InitAllocate
  end type betr_biogeo_state_type

contains

  subroutine Init(this, bounds)

  class(betr_biogeo_state_type) :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  class(betr_biogeo_state_type) :: this
  type(betr_bounds_type), intent(in) :: bounds

  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp; endp=bounds%endp
  begc = bounds%begc; endc=bounds%endc
  lbj = bounds%lbj; ubj=bounds%ubj

  !soilhydrology
  allocate(this%zwts_col           (begc:endc) ) ! the shallower between zwt_perch and zwt

  end subroutine InitAllocate


end module BeTR_biogeoStateType
