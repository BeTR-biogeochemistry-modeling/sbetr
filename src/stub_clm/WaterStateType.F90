module WaterstateType
  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varcon      , only : spval
  implicit none
  save
  private
  !----------------------------------------------------
  ! column water state variables structure
  !----------------------------------------------------
  type, public :: Waterstate_type

    real(r8), pointer :: h2osoi_liq_col(:,:)    => null()   !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_col(:,:)    => null()   !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liq_old_col(:,:)    => null()   !liquid water (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_old_col(:,:)    => null()   !ice lens (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liqvol_col(:,:)  => null()  !volumetric liquid water content
    real(r8), pointer :: h2osoi_icevol_col(:,:)  => null()  !volumetric ice water content
    real(r8), pointer :: h2osoi_vol_col(:,:)     => null()  !volumetric water content, total
    real(r8), pointer :: air_vol_col(:,:)       => null()   !volume possessed by air
    real(r8), pointer :: finundated_col         (:)  => null() ! fraction of column that is inundated, this is for bgc caclulation in betr
    real(r8), pointer :: rho_vap_col(:,:)              => null()  !concentration of bulk water vapor, assume in equilibrium with the liquid water
    real(r8), pointer :: rhvap_soi_col(:,:)             => null()    !soil relative humidity
    real(r8), pointer :: smp_l_col              (:,:) => null()! col liquid phase soil matric potential, mm
    real(r8), pointer :: frac_h2osfc_col        (:)  => null() ! col fractional area with surface water greater than zero
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

    allocate(this%h2osoi_liq_col(begc:endc, lbj:ubj))     ; this%h2osoi_liq_col(:,:) = spval
    allocate(this%h2osoi_ice_col(begc:endc, lbj:ubj))     ; this%h2osoi_ice_col(:,:) = spval
    allocate(this%h2osoi_liq_old_col(begc:endc, lbj:ubj))     ; this%h2osoi_liq_old_col(:,:) = spval
    allocate(this%h2osoi_ice_old_col(begc:endc, lbj:ubj))     ; this%h2osoi_ice_old_col(:,:) = spval
    allocate(this%h2osoi_liqvol_col(begc:endc, lbj:ubj))  ; this%h2osoi_liqvol_col(:,:) = spval
    allocate(this%h2osoi_icevol_col(begc:endc, lbj:ubj))  ; this%h2osoi_icevol_col(:,:) = spval
    allocate(this%h2osoi_vol_col(begc:endc, lbj:ubj))     ; this%h2osoi_vol_col(:,:) = spval
    allocate(this%air_vol_col(begc:endc, lbj:ubj))        ; this%air_vol_col(:,:) = spval
    allocate(this%rho_vap_col(begc:endc, lbj:ubj))            ; this%air_vol_col(:,:) = spval
    allocate(this%rhvap_soi_col(begc:endc, lbj:ubj))          ; this%rhvap_soi_col(:,:) = spval
    allocate(this%smp_l_col  (begc:endc,lbj:ubj))         ; this%smp_l_col              (:,:) = spval
    allocate(this%finundated_col         (begc:endc))                     ; this%finundated_col         (:)   = spval
    allocate(this%frac_h2osfc_col        (begc:endc))                     ; this%frac_h2osfc_col        (:)   = spval
  end subroutine InitAllocate

end module WaterstateType
