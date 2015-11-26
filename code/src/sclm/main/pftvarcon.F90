module pftvarcon
  use shr_kind_mod, only : r8 => shr_kind_r8
 use clm_varpar  , only : mxpft
implicit none


  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  real(r8), allocatable :: crop(:)        !crop pft: 0. = not crop, 1. = crop pft
  integer :: noveg                  !value for not vegetated 
contains



!-----------------------------------------------------------------------
subroutine pftconrd

! !ARGUMENTS:
implicit none



    allocate( crop          (0:mxpft) )
end subroutine pftconrd
end module pftvarcon
