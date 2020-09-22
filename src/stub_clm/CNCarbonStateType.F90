module CNCarbonStateType
  use elm_varcon     , only : spval, ispval, c14ratio
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use elm_varpar     , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: carbonstate_type
    real(r8), pointer :: decomp_cpools_vr    (:,:,:) => null() ! col (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    real(r8), pointer :: frootc_patch             (:)    => null() ! (gC/m2) fine root C
    real(r8), pointer :: totlitc_1m          (:) => null()
    real(r8), pointer :: totlitc             (:) => null()
    real(r8), pointer :: totsomc             (:) => null()
    real(r8), pointer :: cwdc                (:) => null()
    real(r8), pointer :: totsomc_1m          (:) => null()
    real(r8), pointer :: decomp_som2c_vr     (:,:)=> null()
    real(r8), pointer :: som1c               (:) => null()
    real(r8), pointer :: som2c               (:) => null()
    real(r8), pointer :: som3c               (:) => null()
    real(r8), pointer :: domc                (:) => null()
    real(r8), pointer :: beg_totomc          (:) => null()
    real(r8), pointer :: totomc              (:) => null()
    real(r8), pointer :: cmass_residual       (:) => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type carbonstate_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds

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
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    if(nlevdecomp_full>0 .and. ndecomp_pools>0)then
      allocate(this%decomp_cpools_vr(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%decomp_cpools_vr(:,:,:)= spval
    endif
    allocate(this%frootc_patch (begp :endp));     this%frootc_patch(:)   = spval
    allocate(this%cwdc(begc:endc)); this%cwdc(:) = spval
    allocate(this%totlitc(begc:endc)); this%totlitc(:) = spval
    allocate(this%totsomc(begc:endc)); this%totsomc(:) = spval
    allocate(this%totlitc_1m(begc:endc)); this%totlitc_1m(:) = spval
    allocate(this%totsomc_1m(begc:endc)); this%totsomc_1m(:) = spval
    allocate(this%som1c(begc:endc)); this%som1c(:) = spval
    allocate(this%som2c(begc:endc)); this%som2c(:) = spval
    allocate(this%som3c(begc:endc)); this%som3c(:) = spval
    allocate(this%domc(begc:endc)); this%domc(:) = spval
    allocate(this%beg_totomc(begc:endc)); this%beg_totomc(:) = spval
    allocate(this%totomc(begc:endc)); this%totomc(:) = spval
    allocate(this%cmass_residual(begc:endc)); this%cmass_residual(:) = spval
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use elm_varctl , only : nsrest, nsrStartup
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds
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

    begc=bounds%begc; endc=bounds%endc

    this%totomc(begc:endc) = 0._r8

  end subroutine initCold

end module CNCarbonStateType
