module SoilStateType
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
  ! column physical state variables structure
  !----------------------------------------------------
  type, public :: soilstate_type
    real(r8), pointer :: bsw(:,:)              !Clapp and Hornberger "b" (nlevgrnd)  
    real(r8), pointer :: watsat(:,:)           !volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: eff_porosity(:,:)     !effective porosity = porosity - vol_ice (nlevgrnd)
    
  contains
    procedure, public  :: Init         
    procedure, private :: InitAllocate 
  end type soilstate_type
  
  contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilstate_type) :: this
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
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    lbj  = bounds%lbj; ubj = bounds%ubj

    allocate(this%bsw(begc:endc, lbj:ubj))
    allocate(this%watsat(begc:endc, lbj:ubj))
    allocate(this%eff_porosity(begc:endc, lbj:ubj))
    
  end subroutine InitAllocate    
end module SoilStateType
