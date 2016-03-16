module BeTR_CanopyStateType
  use betr_varcon         , only : spval => bspval, ispval => bispval, c14ratio => bc14ratio
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use betr_decompMod      , only : betr_bounds_type
implicit none

  type, public :: betr_canopystate_type
     real(r8) , pointer :: altmax_col               (:)   ! col maximum annual depth of thaw
     real(r8) , pointer :: altmax_lastyear_col      (:)   ! col prior year maximum annual depth of thaw
     real(r8),  pointer :: lbl_rsc_h2o_patch        (:)   ! laminar boundary layer resistance for water over dry leaf (s/m)
     real(r8) , pointer :: elai_patch               (:)   ! patch canopy one-sided leaf area index with burying by snow
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type betr_canopystate_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(betr_canopystate_type) :: this
    type(betr_bounds_type), intent(in) :: bounds

    call this%InitAllocate ( bounds )

    call this%InitCold ( bounds )

  end subroutine Init
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(betr_canopystate_type) :: this
    type(betr_bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    allocate(this%altmax_col               (begc:endc))           ; this%altmax_col               (:)   = spval
    allocate(this%altmax_lastyear_col      (begc:endc))           ; this%altmax_lastyear_col      (:)   = spval
    allocate(this%lbl_rsc_h2o_patch        (begp:endp))           ; this%lbl_rsc_h2o_patch        (:)   = nan
    allocate(this%elai_patch               (begp:endp))           ; this%elai_patch               (:)   = nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use clm_varctl , only : nsrest, nsrStartup
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(betr_canopystate_type) :: this
    type(betr_bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p,n,j,m            ! indices
    real(r8) ,pointer     :: gdp (:)                  ! global gdp data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: peatf (:)                ! global peatf data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: soilorder_rdin (:)       ! global soil order data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: abm (:)                  ! global abm data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: gti (:)                  ! read in - fmax (needs to be a pointer for use in ncdio)
    integer               :: dimid                    ! dimension id
    integer               :: ier                      ! error status
    type(file_desc_t)     :: ncid                     ! netcdf id
    logical               :: readvar
    character(len=256)    :: locfn                    ! local filename
    integer               :: begc, endc
    integer               :: begg, endg



  end subroutine initCold

end module BeTR_CanopyStateType
