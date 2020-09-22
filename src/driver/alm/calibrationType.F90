module calibrationType
  use shr_kind_mod        , only : r8 => shr_kind_r8
implicit none

  type, public :: calibration_type

    real(r8), pointer :: plant_nh4_vmax_scalar(:,:) => null()
    real(r8), pointer :: plant_no3_vmax_scalar(:,:) => null()
    real(r8), pointer :: plant_p_vmax_scalar  (:,:) => null()
    real(r8), pointer :: plant_nh4_km_scalar  (:,:) => null()
    real(r8), pointer :: plant_no3_km_scalar  (:,:) => null()
    real(r8), pointer :: plant_p_km_scalar  (:,:) => null()
  contains
    procedure, public :: Init
  end type calibration_type

  contains

  subroutine Init(this, bounds, npfts)
  use decompMod       , only : bounds_type
  implicit none
  class(calibration_type), intent(inout) :: this
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: npfts

  integer :: begg, endg
  begg = bounds%begg; endg=bounds%endg
  allocate(this%plant_nh4_vmax_scalar(begg:endg, 0:npfts))
  allocate(this%plant_no3_vmax_scalar(begg:endg, 0:npfts))
  allocate(this%plant_p_vmax_scalar(begg:endg, 0:npfts))
  allocate(this%plant_nh4_km_scalar(begg:endg, 0:npfts))
  allocate(this%plant_no3_km_scalar(begg:endg, 0:npfts))
  allocate(this%plant_p_km_scalar(begg:endg, 0:npfts))

  end subroutine Init

end module calibrationType
