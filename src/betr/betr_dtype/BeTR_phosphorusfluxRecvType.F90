module BeTR_phosphorusfluxRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  type, public :: betr_phosphorusflux_recv_type
    real(r8), pointer :: sminp_leached_col(:) => null()
    real(r8), pointer :: sminp_to_plant_patch(:) => null() !integrated phosphate goes to plant at patch (gN/m2/s)
    real(r8), pointer :: fire_decomp_ploss_col(:) => null()
    real(r8), pointer :: supplement_to_sminp_col(:) => null()
    real(r8), pointer :: secondp_to_occlp_col(:) => null()
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
  end type betr_phosphorusflux_recv_type

 contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_phosphorusflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_phosphorusflux_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%sminp_leached_col(begc:endc))
  allocate(this%sminp_to_plant_patch(begp:endp))
  allocate(this%fire_decomp_ploss_col(begc:endc))
  allocate(this%supplement_to_sminp_col(begc:endc))
  allocate(this%secondp_to_occlp_col(begc:endc))
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_phosphorusflux_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%sminp_leached_col(:) = value_column
  this%sminp_to_plant_patch(:) = value_column
  this%fire_decomp_ploss_col(:) = value_column
  this%supplement_to_sminp_col(:) = value_column
  this%secondp_to_occlp_col(:) = value_column
  end subroutine reset
end module BeTR_phosphorusfluxRecvType
