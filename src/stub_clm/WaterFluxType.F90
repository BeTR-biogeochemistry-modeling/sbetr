module WaterfluxType

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
  type, public :: waterflux_type
    real(r8), pointer :: qflx_adv(:,:)    => null()   !advection velocity from one layer to another, (0:nlevgrnd), positive downward
    real(r8), pointer :: qflx_infl(:)	    => null()   !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)	    => null()   !surface runoff (mm H2O /s)
    real(r8), pointer :: h2oliq_vol_tendency(:,:)   => null()      !temporal change of water during the solution of soil water movement
    real(r8), pointer :: qflx_gross_evap_soil (:) => null()  ! col gross infiltration from soil, this satisfies the relationship qflx_infl = qflx_gross_infl_soil-qflx_gross_evap_soil
    real(r8), pointer :: qflx_gross_infl_soil (:) => null()  ! col gross infiltration, before considering the evaporation, mm/s
    real(r8), pointer :: qflx_rootsoi         (:,:)=> null() ! col root and soil water exchange [mm H2O/s] [+ into root]
    real(r8), pointer :: qflx_drain_vr        (:,:)=> null() ! col liquid water losted as drainage (m /time step)
    real(r8), pointer :: qflx_totdrain        (:) => null()  ! col total liquid water drainage  (m/time step), updated in betr
    real(r8), pointer :: qflx_dew_grnd        (:)=> null()   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
    real(r8), pointer :: qflx_dew_snow        (:) => null()  ! col surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_vol    (:)=> null()
    real(r8), pointer :: qflx_sub_snow        (:)=> null()   ! col sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_h2osfc2topsoi   (:) => null()  ! col liquid water coming from surface standing water top soil (mm H2O/s)
    real(r8), pointer :: qflx_snow2topsoi     (:) => null()  ! col liquid water coming from residual snow to topsoil (mm H2O/s)
    real(r8), pointer :: qflx_runoff          (:)   ! col total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)s
  contains
    procedure          :: Init
    procedure, private :: InitAllocate
  end type waterflux_type

  contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(waterflux_type) :: this
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
    class(waterflux_type) :: this
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

    allocate(this%qflx_adv(begc:endc, lbj-1:ubj))
    allocate(this%qflx_infl(begc:endc))
    allocate(this%qflx_surf(begc:endc))
    allocate(this%qflx_gross_evap_soil (begc:endc))              ; this%qflx_gross_evap_soil (:)   = nan
    allocate(this%qflx_gross_infl_soil (begc:endc))              ; this%qflx_gross_infl_soil (:)   = nan
    allocate(this%qflx_rootsoi         (begc:endc,lbj:ubj))      ; this%qflx_rootsoi         (:,:) = nan
    allocate(this%qflx_drain_vr        (begc:endc,lbj:ubj))      ; this%qflx_drain_vr        (:,:) = nan
    allocate(this%qflx_dew_grnd        (begc:endc))              ; this%qflx_dew_grnd        (:)   = nan
    allocate(this%qflx_dew_snow        (begc:endc))              ; this%qflx_dew_snow        (:)   = nan
    allocate(this%qflx_sub_snow_vol    (begc:endc))              ; this%qflx_sub_snow_vol    (:)   = 0._r8
    allocate(this%qflx_sub_snow        (begc:endc))              ; this%qflx_sub_snow        (:)   = 0.0_r8
    allocate(this%qflx_snow2topsoi     (begc:endc))              ; this%qflx_snow2topsoi     (:)   = nan
    allocate(this%qflx_h2osfc2topsoi   (begc:endc))              ; this%qflx_h2osfc2topsoi   (:)   = nan
    allocate( this%qflx_totdrain       (begc:endc))              ; this%qflx_totdrain        (:)   = nan
    allocate(this%qflx_runoff     (begc:endc))              ; this%qflx_runoff     (:)   = nan
  end subroutine InitAllocate




end module WaterfluxType
