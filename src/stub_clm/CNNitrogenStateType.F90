module CNNitrogenStateType
  use clm_varcon     , only : spval, ispval, c14ratio
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varpar     , only : ndecomp_pools, nlevdecomp_full
implicit none

  type, public :: nitrogenstate_type
    real(r8), pointer :: decomp_npools_vr         (:,:,:) => null()! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
    real(r8), pointer :: smin_no3_vr              (:,:)  => null() ! col (gN/m3) vertically-resolved soil mineral NO3
    real(r8), pointer :: smin_no3                 (:)   => null()  ! col (gN/m2) soil mineral NO3 pool
    real(r8), pointer :: smin_nh4_vr              (:,:)=> null()   ! col (gN/m3) vertically-resolved soil mineral NH4
    real(r8), pointer :: smin_nh4                 (:)   => null()  ! col (gN/m2) soil mineral NH4 pool
    real(r8), pointer :: sminn_vr                 (:,:)=> null()   ! col (gN/m3) vertically-resolved soil mineral N
    real(r8), pointer :: sminn                    (:) => null()
    real(r8), pointer :: cwdn                     (:) => null()
    real(r8), pointer :: totsomn                  (:) => null()
    real(r8), pointer :: totsomn_1m               (:) => null()
    real(r8), pointer :: totlitn_1m               (:) => null()
    real(r8), pointer :: totlitn                  (:) => null()
    real(r8), pointer :: som1n                    (:) => null()
    real(r8), pointer :: som2n                    (:) => null()
    real(r8), pointer :: som3n                    (:) => null()
    real(r8), pointer :: domn                     (:) => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type nitrogenstate_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(nitrogenstate_type) :: this
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
    class(nitrogenstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    allocate(this%decomp_npools_vr(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%decomp_npools_vr(:,:,:)= spval
    allocate(this%smin_no3_vr          (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_vr          (:,:) = spval
    allocate(this%smin_nh4_vr          (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_vr          (:,:) = spval
    allocate(this%smin_no3             (begc:endc))                   ; this%smin_no3             (:)   = spval
    allocate(this%smin_nh4             (begc:endc))                   ; this%smin_nh4             (:)   = spval
    allocate(this%sminn_vr             (begc:endc,1:nlevdecomp_full)) ; this%sminn_vr             (:,:) = spval
    allocate(this%cwdn(begc:endc)); this%cwdn(:) = spval
    allocate(this%totlitn(begc:endc)); this%totlitn(:) = spval
    allocate(this%totsomn(begc:endc)); this%totsomn(:) = spval
    allocate(this%totlitn_1m(begc:endc)); this%totlitn_1m(:) = spval
    allocate(this%totsomn_1m(begc:endc)); this%totsomn_1m(:) = spval
    allocate(this%som1n(begc:endc)); this%som1n(:) = spval
    allocate(this%som2n(begc:endc)); this%som2n(:) = spval
    allocate(this%som3n(begc:endc)); this%som3n(:) = spval
    allocate(this%sminn(begc:endc)); this%sminn(:) = spval
    allocate(this%domn(begc:endc)); this%domn(:) = spval

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
    class(nitrogenstate_type) :: this
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


  end subroutine initCold

end module CNNitrogenStateType
