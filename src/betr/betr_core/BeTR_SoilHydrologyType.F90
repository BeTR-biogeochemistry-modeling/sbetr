Module BeTR_SoilHydrologyType

! DESCRIPTION
! derived data for soilhydrology

 use bshr_kind_mod   , only : r8 => shr_kind_r8
 use betr_decompMod      , only : betr_bounds_type
implicit none

  type, public :: betr_soilhydrology_type

  real(r8), pointer :: fracice_col       (:,:)   ! col fractional impermeability (-)
  real(r8), pointer :: zwts_col           (:)    ! the shallower between zwt_perch and zwt
  real(r8), pointer :: qflx_bot_col      (:)     ! bottom of soil col flux, (mm/s)
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
  end type betr_soilhydrology_type

  contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(betr_soilhydrology_type) :: this
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
    class(betr_soilhydrology_type) :: this
    type(betr_bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    lbj  = bounds%lbj; ubj = bounds%ubj

    allocate(this%fracice_col       (begc:endc,lbj:ubj))        ; this%fracice_col       (:,:)   = nan
    allocate(this%zwts_col           (begc:endc))                ; this%zwts_col         (:)     = nan
    allocate(this%qflx_bot_col      (begc:endc))                 ; this%qflx_bot_col      (:)     = nan
  end subroutine InitAllocate
end Module BeTR_SoilHydrologyType
