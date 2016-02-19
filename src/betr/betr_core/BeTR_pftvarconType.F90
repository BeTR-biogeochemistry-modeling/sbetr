module BeTR_pftvarconType

implicit none

  type, public :: betr_pftvarcon_type
   integer :: nc3_arctic_grass       !value for C3 arctic grass
   integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
   integer :: nc4_grass              !value for C4 grass
   integer :: noveg                  !value for not vegetated
  end type betr_pftvarcon_type

  type(betr_pftvarcon_type), public ::betr_pftvarcon

end module BeTR_pftvarconType
