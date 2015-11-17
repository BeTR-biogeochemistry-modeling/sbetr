module WaterfluxType

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
  ! column water flux variables structure
  !----------------------------------------------------
  type, public :: waterflux_type
    real(r8), pointer :: qflx_adv(:,:)         !advection velocity from one layer to another, (0:nlevgrnd)
    real(r8), pointer :: qflx_infl(:)	       !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)	       !surface runoff (mm H2O /s)
    real(r8), pointer :: h2oliq_vol_tendency(:,:)         !temporal change of water during the solution of soil water movement
    real(r8), pointer :: qflx_rootsoi(:,:)     ! water flux between root and soil [m/s]
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
    allocate(this%qflx_rootsoi(begc:endc,lbj:ubj)) 
  end subroutine InitAllocate


  

end module WaterfluxType
