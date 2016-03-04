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

    real(r8), pointer :: h2osoi_liq_col(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_col(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liq_old(:,:)       !liquid water (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_old(:,:)       !ice lens (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liqvol_col(:,:)    !volumetric liquid water content
    real(r8), pointer :: h2osoi_icevol_col(:,:)    !volumetric ice water content
    real(r8), pointer :: h2osoi_vol_col(:,:)       !volumetric water content, total
    real(r8), pointer :: air_vol_col(:,:)          !volume possessed by air
    real(r8), pointer :: finundated_col         (:)   ! fraction of column that is inundated, this is for bgc caclulation in betr
    real(r8), pointer :: rho_vap(:,:)                !concentration of bulk water vapor, assume in equilibrium with the liquid water
    real(r8), pointer :: rhvap_soi(:,:)                 !soil relative humidity
    real(r8), pointer :: smp_l_col              (:,:) ! col liquid phase soil matric potential, mm
    real(r8), pointer :: frac_h2osfc_col        (:)   ! col fractional area with surface water greater than zero
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

    allocate(this%h2osoi_liq_col(begc:endc, lbj:ubj))     ; this%h2osoi_liq_col(:,:) = nan
    allocate(this%h2osoi_ice_col(begc:endc, lbj:ubj))     ; this%h2osoi_ice_col(:,:) = nan
    allocate(this%h2osoi_liq_old(begc:endc, lbj:ubj))     ; this%h2osoi_liq_old(:,:) = nan
    allocate(this%h2osoi_ice_old(begc:endc, lbj:ubj))     ; this%h2osoi_ice_old(:,:) = nan
    allocate(this%h2osoi_liqvol_col(begc:endc, lbj:ubj))  ; this%h2osoi_liqvol_col(:,:) = nan
    allocate(this%h2osoi_icevol_col(begc:endc, lbj:ubj))  ; this%h2osoi_icevol_col(:,:) = nan
    allocate(this%h2osoi_vol_col(begc:endc, lbj:ubj))     ; this%h2osoi_vol_col(:,:) = nan
    allocate(this%air_vol_col(begc:endc, lbj:ubj))        ; this%air_vol_col(:,:) = nan
    allocate(this%rho_vap(begc:endc, lbj:ubj))            ; this%air_vol_col(:,:) = nan
    allocate(this%rhvap_soi(begc:endc, lbj:ubj))          ; this%rhvap_soi(:,:) = nan
    allocate(this%smp_l_col  (begc:endc,lbj:ubj))         ; this%smp_l_col              (:,:) = nan
    allocate(this%finundated_col         (begc:endc))                     ; this%finundated_col         (:)   = nan
    allocate(this%frac_h2osfc_col        (begc:endc))                     ; this%frac_h2osfc_col        (:)   = nan
  end subroutine InitAllocate

end module WaterstateType
