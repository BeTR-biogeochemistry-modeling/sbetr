module PfTemperatureType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use elm_varcon      , only : spval
  implicit none
  save
  private

!----------------------------------------------------
! column energy state variables structure
!----------------------------------------------------
  type, public :: pf_temperature_type
    real(r8), pointer :: t_veg              (:) => null()  ! patch vegetation temperature (Kelvin)
  contains
     procedure, public  :: Init
     procedure, private :: InitAllocate

  end type pf_temperature_type

  contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(pf_temperature_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(pf_temperature_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    allocate(this%t_veg  (begp:endp))        ; this%t_veg              (:)   = spval
  end subroutine InitAllocate
end module PfTemperatureType
