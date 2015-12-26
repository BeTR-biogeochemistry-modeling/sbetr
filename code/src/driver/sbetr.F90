program main

use sbetrDriverMod

use decompMod            , only: bounds_type
use clmgridMod           , only : init_clm_vertgrid
use clm_varpar           , only : nlevgrnd
use clm_initializeMod    , only : initialize1, initialize2

implicit none

  type(bounds_type) :: bounds
  type(time_type) :: ttime
  integer :: numf
  integer, pointer :: filter(:)

  !set up mask
  bounds%begc = 1
  bounds%endc = 1
  bounds%begp = 1
  bounds%endp = 1
  bounds%lbj  = 1
  bounds%ubj  = nlevgrnd

  numf = 1
  ttime%time_end = 1800*48*10
  ttime%restart_dtime=1800*2
  allocate(filter(numf))
  filter(:) = 1

  !set up grid
  call init_clm_vertgrid(nlevgrnd)

  call initialize1(bounds)

  call initialize2(bounds)

  call sbetrBGC_driver(bounds, numf, filter, ttime)

end program main
