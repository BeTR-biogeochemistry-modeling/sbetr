module BeTR_aerocondType

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use decompMod              , only : bounds_type
  use clm_varcon             , only : spval, ispval
implicit none
save
private


!
! !PUBLIC DATA:
!
  type, public :: betr_aerecond_type


   real(r8), pointer :: tempavg_agnpp_patch  (:) !  temporary average above-ground NPP (gC/m2/s)
   real(r8), pointer :: annavg_agnpp_patch   (:) !  annual average above-ground NPP (gC/m2/s)
   real(r8), pointer :: tempavg_bgnpp_patch  (:) !  temporary average below-ground NPP (gC/m2/s)
   real(r8), pointer :: annavg_bgnpp_patch   (:) !  annual average below-ground NPP (gC/m2/s)
   real(r8), pointer :: plant_frootsc_patch  (:) !
  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type betr_aerecond_type

 contains

  subroutine Init(this, bounds)

  class(betr_aerecond_type) :: this
  type(bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)


  class(betr_aerecond_type) :: this
  type(bounds_type), intent(in) :: bounds

  integer :: begp, endp

  begp = bounds%begp; endp=bounds%endp

  allocate(this%plant_frootsc_patch (begp:endp)); this%plant_frootsc_patch (:)   = nan
  allocate(this%annavg_agnpp_patch  (begp:endp)); this%annavg_agnpp_patch  (:) = spval ! To detect first year
  allocate(this%annavg_bgnpp_patch  (begp:endp)); this%annavg_bgnpp_patch  (:) = spval ! To detect first year
  allocate(this%tempavg_agnpp_patch (begp:endp)); this%tempavg_agnpp_patch (:) = spval
  allocate(this%tempavg_bgnpp_patch (begp:endp)); this%tempavg_bgnpp_patch (:) = spval

  end subroutine InitAllocate
end module BeTR_aerocondType
