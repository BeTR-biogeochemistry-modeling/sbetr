module BeTR_CarbonFluxType
  !carbon flux data type to pass in data
  use bshr_kind_mod, only: r8 => shr_kind_r8
implicit none
save
private

!
! !PUBLIC DATA:
!
  type, public :: betr_carbonflux_type

   real(r8), pointer :: annsum_npp_patch(:)
   real(r8), pointer :: agnpp_patch(:)
   real(r8), pointer :: bgnpp_patch(:)

  end type betr_carbonflux_type

end module BeTR_CarbonFluxType
