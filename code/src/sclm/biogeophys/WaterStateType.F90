module WaterstateType
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
  ! column water state variables structure
  !----------------------------------------------------
  type, public :: Waterstate_type

    real(r8), pointer :: h2osoi_liq(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
    real(r8), pointer :: h2osoi_ice(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liq_old(:,:)   !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
    real(r8), pointer :: h2osoi_ice_old(:,:)   !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)   
    real(r8), pointer :: h2osoi_liqvol(:,:)    !volumetric liquid water content
    real(r8), pointer :: h2osoi_icevol(:,:)    !volumetric ice water content
    real(r8), pointer :: h2osoi_vol(:,:)       !volumetric water content, total
    real(r8), pointer :: air_vol(:,:)          !volume possessed by air
    real(r8), pointer :: rho_vap(:,:)                !concentration of bulk water vapor, assume in equilibrium with the liquid water
    real(r8), pointer :: rhvap_soi(:,:)                 !soil relative humidity
  contains
    procedure          :: Init
    procedure, private :: InitAllocate      
  end type Waterstate_type
  
  contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(Waterstate_type) :: this
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
    class(Waterstate_type) :: this
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
  
    allocate(this%h2osoi_liq(begc:endc, lbj:ubj))
    allocate(this%h2osoi_ice(begc:endc, lbj:ubj))
    allocate(this%h2osoi_liq_old(begc:endc, lbj:ubj))
    allocate(this%h2osoi_ice_old(begc:endc, lbj:ubj))
    allocate(this%h2osoi_liqvol(begc:endc, lbj:ubj))
    allocate(this%h2osoi_icevol(begc:endc, lbj:ubj))
    allocate(this%air_vol(begc:endc, lbj:ubj))
    allocate(this%rho_vap(begc:endc, lbj:ubj))
    allocate(this%rhvap_soi(begc:endc, lbj:ubj))
    
  end subroutine InitAllocate  

end module WaterstateType
