module PfWaterfluxType

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
  ! column water flux variables structure
  !----------------------------------------------------
  type, public :: pf_waterflux_type
    real(r8), pointer :: qflx_tran_veg      (:)=> null()
    real(r8), pointer :: qflx_rootsoi       (:,:) => null()! pft root and soil water exchange [mm H2O/s] [+ into atmosphere]
    real(r8), pointer :: qflx_rootsoi_frac  (:,:) => null()
  contains
    procedure          :: Init
    procedure, private :: InitAllocate
  end type pf_waterflux_type

  contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(pf_waterflux_type) :: this
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
    class(pf_waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------


    begp = bounds%begp; endp= bounds%endp
    lbj  = bounds%lbj; ubj = bounds%ubj

    allocate(this%qflx_tran_veg      (begp:endp))              ; this%qflx_tran_veg      (:)   = nan
    allocate(this%qflx_rootsoi       (begp:endp,lbj:ubj))      ; this%qflx_rootsoi       (:,:) = nan
    allocate(this%qflx_rootsoi_frac  (begp:endp,lbj:ubj))      ; this%qflx_rootsoi_frac(:,:) = nan
  end subroutine InitAllocate


end module PfWaterfluxType
