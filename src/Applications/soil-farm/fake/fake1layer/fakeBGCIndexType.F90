module fakeBGCIndexType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_varcon    , only : var_flux_type, var_state_type
implicit none
private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

integer, parameter :: loc_name_len=64

  type, public :: fake_bgc_index_type
    integer           :: nelms
    integer           :: nstvars          !number of equations for the state variabile vector
    integer           :: nreactions
  contains
    procedure, public  :: Init
  end type fake_bgc_index_type
contains
 !-------------------------------------------------------------------------------
  subroutine Init(this)
  !
  ! DESCRIPTION:
  ! Initialize fake_bgc_index_type
  ! !USES:
  implicit none
  ! !ARGUMENTS:
  class(fake_bgc_index_type), intent(inout) :: this
  
  this%nelms=1
  this%nstvars=1
  this%nreactions=1
  end subroutine Init

end module fakeBGCIndexType

