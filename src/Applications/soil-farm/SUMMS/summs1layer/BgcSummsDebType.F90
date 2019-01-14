module BgcSummsDebType
  !DESCRIPTION
!module defines the data type for microbial growth
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  
  implicit none

  ! Kitchen sink container that can be used to pass data into call
  ! back functions. To keep the size down:
  !
  ! * all large arrays should be dynamically allocated.
  !
  ! * variables in the struct should be kept generic, then use
  ! associate blocks in the calling routine and call back functions to
  ! give them meaningful names.
  !
  type, public :: debs
    real(r8), pointer :: gmax_mic  => null()                 ! Maximum microbial growth rate (1/day)
    real(r8), pointer :: yld_mic   => null()               ! Yield rate of structural biomass from reserve metabolistes (g microbep/g micc)
    real(r8), pointer :: yld_enz   => null()                  ! Yield rate of enzymes from reserve metabolistes (g enzymes/g micc)
    real(r8), pointer :: mr_mic   => null()
    real(r8), pointer :: pmax_enz    => null()                     
    real(r8), pointer :: je         => null()                   
    real(r8), pointer :: ec     => null()
    real(r8), pointer :: gB      => null()
    real(r8), pointer :: pE    => null()
  end type debs

end module BgcSummsDebType


