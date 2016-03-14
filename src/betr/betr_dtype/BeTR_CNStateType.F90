module BeTR_CNStateType
  !data type to pass in cn information
implicit none
save
private

!
! !PUBLIC DATA:
!
type, public :: betr_cnstate_type
  integer, pointer :: isoilorder(:)
end type betr_cnstate_type

end module BeTR_CNStateType
