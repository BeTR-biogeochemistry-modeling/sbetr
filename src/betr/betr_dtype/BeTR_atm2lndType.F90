module betr_atm2lndType

! DESCRIPTION
! derived data for atmospheric/land exchange variables

  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use betr_decompMod  , only : betr_bounds_type

implicit none

  type, public :: betr_atm2lnd_type

  real(r8), pointer :: forc_pbot_downscaled_col      (:)   => null() ! downscaled atm pressure (Pa)
  real(r8), pointer :: forc_t_downscaled_col         (:)   => null() ! downscaled atm temperature (Kelvin)
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
  end type betr_atm2lnd_type
  !----------------------------------------------------

  contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(betr_atm2lnd_type) :: this
    type(betr_bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    ! !ARGUMENTS:
    class(betr_atm2lnd_type) :: this
    type(betr_bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg

    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc


    allocate(this%forc_pbot_downscaled_col(begc:endc));    this%forc_pbot_downscaled_col(:) = nan

  end subroutine InitAllocate
end module betr_atm2lndType
