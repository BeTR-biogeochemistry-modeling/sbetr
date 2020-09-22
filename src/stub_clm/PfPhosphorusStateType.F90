module PhosphorusStateType
  use elm_varcon     , only : spval, ispval, c14ratio
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use elm_varpar     , only : ndecomp_pools, nlevdecomp_full
implicit none

  type, public :: phosphorusstate_type
    real(r8), pointer :: decomp_ppools            (:,:)  => null() ! col (gP/m2)  decomposing (litter, cwd, soil) P pools
    real(r8), pointer :: sminp                    (:)   => null()  ! col (gP/m2) soil mineral P
    real(r8), pointer :: solutionp                (:)  => null()       ! col (gP/m2) soil solution P
    real(r8), pointer :: labilep                  (:)   => null()      ! col (gP/m2) soil labile mineral P
    real(r8), pointer :: secondp                  (:)   => null()      ! col (gP/m2) soil secondary mineralP
    real(r8), pointer :: occlp                    (:)  => null()       ! col (gP/m2) soil occluded mineral P
    real(r8), pointer :: primp                    (:)  => null()       ! col (gP/m2) soil parimary mineral P
    real(r8), pointer :: totlitp                  (:)  => null()   ! col (gP/m2) total litter phosphorus
    real(r8), pointer :: totsomp                  (:) => null()    ! col (gP/m2) total soil organic matter phosphorus
    real(r8), pointer :: totecosysp               (:) => null()    ! col (gP/m2) total ecosystem phosphorus, incl veg
    real(r8), pointer :: totcolp                  (:) => null()    ! col (gP/m2) total column phosphorus, incl veg
    real(r8), pointer :: cwdp                     (:) => null()    ! col (gP/m2) Diagnostic: coarse woody debris P
    real(r8), pointer :: totlitp_1m               (:) => null()
    real(r8), pointer :: totsomp_1m               (:) => null()
    ! patch averaged to column variables
    real(r8), pointer :: totvegp                  (:)  => null()   ! col (gP/m2) total vegetation phosphorus (p2c)

    real(r8), pointer :: decomp_ppools_vr         (:,:,:) => null()    ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
    real(r8), pointer :: solutionp_vr             (:,:)  => null()     ! col (gP/m3) vertically-resolved soil solution P
    real(r8), pointer :: labilep_vr               (:,:) => null()      ! col (gP/m3) vertically-resolved soil labile mineral P
    real(r8), pointer :: secondp_vr               (:,:)  => null()     ! col (gP/m3) vertically-resolved soil secondary mineralP
    real(r8), pointer :: occlp_vr                 (:,:) => null()      ! col (gP/m3) vertically-resolved soil occluded mineral P
    real(r8), pointer :: primp_vr                 (:,:)  => null()     ! col (gP/m3) vertically-resolved soil parimary mineral P
    real(r8), pointer :: sminp_vr                 (:,:)  => null()     ! col (gP/m3) vertically-resolved soil total mineral P, diagnostic
    real(r8), pointer :: som1p                    (:) => null()
    real(r8), pointer :: som2p                    (:) => null()
    real(r8), pointer :: som3p                    (:) => null()
    real(r8), pointer :: domp                     (:) => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type phosphorusstate_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(phosphorusstate_type) :: this
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
    class(phosphorusstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    allocate(this%decomp_ppools        (begc:endc,1:ndecomp_pools))   ; this%decomp_ppools        (:,:) = spval
    allocate(this%solutionp            (begc:endc))                   ; this%solutionp            (:)   = spval
    allocate(this%solutionp            (begc:endc))                   ; this%solutionp            (:)   = spval
    allocate(this%labilep              (begc:endc))                   ; this%labilep              (:)   = spval
    allocate(this%secondp              (begc:endc))                   ; this%secondp              (:)   = spval
    allocate(this%occlp                (begc:endc))                   ; this%occlp                (:)   = spval
    allocate(this%primp                (begc:endc))                   ; this%primp                (:)   = spval
    allocate(this%cwdp                 (begc:endc))                   ; this%cwdp                 (:)   = spval
    allocate(this%sminp                (begc:endc))                   ; this%sminp                (:)   = spval
    allocate(this%totecosysp           (begc:endc))                   ; this%totecosysp           (:)   = spval
    allocate(this%totcolp              (begc:endc))                   ; this%totcolp              (:)   = spval

    allocate(this%decomp_ppools_vr(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%decomp_ppools_vr(:,:,:)= spval
    allocate(this%solutionp_vr         (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr         (:,:) = spval
    allocate(this%labilep_vr           (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr           (:,:) = spval
    allocate(this%secondp_vr           (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr           (:,:) = spval
    allocate(this%occlp_vr             (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr             (:,:) = spval
    allocate(this%primp_vr             (begc:endc,1:nlevdecomp_full)) ; this%primp_vr             (:,:) = spval
    allocate(this%sminp_vr             (begc:endc,1:nlevdecomp_full)) ; this%sminp_vr             (:,:) = spval
    allocate(this%totlitp (begc:endc)); this%totlitp(:) = spval
    allocate(this%totsomp(begc:endc)); this%totsomp(:) = spval
    allocate(this%totlitp_1m(begc:endc)); this%totlitp_1m(:) = spval
    allocate(this%totsomp_1m(begc:endc)); this%totsomp_1m(:) = spval
    allocate(this%som1p(begc:endc)); this%som1p(:) = spval
    allocate(this%som2p(begc:endc)); this%som2p(:) = spval
    allocate(this%som3p(begc:endc)); this%som3p(:) = spval
    allocate(this%domp(begc:endc)); this%domp(:) = spval
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
    class(phosphorusstate_type) :: this
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

end module PhosphorusStateType
