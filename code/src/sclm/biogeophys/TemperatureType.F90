module TemperatureType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun

  implicit none
  save
  private

!----------------------------------------------------
! column energy state variables structure
!----------------------------------------------------
  type, public :: temperature_type
    real(r8), pointer :: t_soisno(:,:)         !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: t_soi_10cm(:)         !soil temperature in top 10cm of soil (Kelvin)
  contains
     procedure, public  :: Init         
     procedure, private :: InitAllocate
     
  end type temperature_type

  contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(temperature_type) :: this
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
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    lbj  = bounds%lbj;  ubj = bounds%ubj
    
    allocate(this%t_soisno(begc:endc, lbj:ubj))
    allocate(this%t_soi_10cm(begc:endc))
    
  end subroutine InitAllocate  
end module TemperatureType  
