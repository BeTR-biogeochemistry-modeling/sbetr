module WaterstateType
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
  ! column water state variables structure
  !----------------------------------------------------
  type, public :: Waterstate_type

    real(r8), pointer :: h2osoi_liq(:,:)    => null()   !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice(:,:)    => null()   !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liq_old(:,:)    => null()   !liquid water (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_old(:,:)    => null()   !ice lens (kg/m2) (old) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liqvol(:,:)  => null()  !volumetric liquid water content
    real(r8), pointer :: h2osoi_icevol(:,:)  => null()  !volumetric ice water content
    real(r8), pointer :: h2osoi_vol(:,:)     => null()  !volumetric water content, total
    real(r8), pointer :: air_vol(:,:)       => null()   !volume possessed by air
    real(r8), pointer :: finundated         (:)  => null() ! fraction of column that is inundated, this is for bgc caclulation in betr
    real(r8), pointer :: rho_vap(:,:)              => null()  !concentration of bulk water vapor, assume in equilibrium with the liquid water
    real(r8), pointer :: rhvap_soi(:,:)             => null()    !soil relative humidity
    real(r8), pointer :: smp_l              (:,:) => null()! col liquid phase soil matric potential, mm
    real(r8), pointer :: frac_h2osfc        (:)  => null() ! col fractional area with surface water greater than zero
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

    allocate(this%h2osoi_liq(begc:endc, lbj:ubj))     ; this%h2osoi_liq(:,:) = spval
    allocate(this%h2osoi_ice(begc:endc, lbj:ubj))     ; this%h2osoi_ice(:,:) = spval
    allocate(this%h2osoi_liq_old(begc:endc, lbj:ubj))     ; this%h2osoi_liq_old(:,:) = spval
    allocate(this%h2osoi_ice_old(begc:endc, lbj:ubj))     ; this%h2osoi_ice_old(:,:) = spval
    allocate(this%h2osoi_liqvol(begc:endc, lbj:ubj))  ; this%h2osoi_liqvol(:,:) = spval
    allocate(this%h2osoi_icevol(begc:endc, lbj:ubj))  ; this%h2osoi_icevol(:,:) = spval
    allocate(this%h2osoi_vol(begc:endc, lbj:ubj))     ; this%h2osoi_vol(:,:) = spval
    allocate(this%air_vol(begc:endc, lbj:ubj))        ; this%air_vol(:,:) = spval
    allocate(this%rho_vap(begc:endc, lbj:ubj))            ; this%air_vol(:,:) = spval
    allocate(this%rhvap_soi(begc:endc, lbj:ubj))          ; this%rhvap_soi(:,:) = spval
    allocate(this%smp_l  (begc:endc,lbj:ubj))         ; this%smp_l              (:,:) = spval
    allocate(this%finundated         (begc:endc))                     ; this%finundated         (:)   = spval
    allocate(this%frac_h2osfc        (begc:endc))                     ; this%frac_h2osfc        (:)   = spval
  end subroutine InitAllocate

end module WaterstateType
