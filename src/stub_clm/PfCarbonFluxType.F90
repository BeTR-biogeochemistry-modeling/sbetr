module PfCarbonFluxType
  use elm_varcon             , only : spval, ispval, c14ratio
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use elm_varcon             , only : spval
  use elm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: pf_carbonflux_type
    real(r8), pointer :: rr                                  (:)  => null()   ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: annsum_npp                          (:) => null() ! patch annual sum of NPP (gC/m2/yr)
    real(r8), pointer :: agnpp                               (:)  => null()   ! (gC/m2/s) aboveground NPP
    real(r8), pointer :: bgnpp                               (:)  => null()   ! (gC/m2/s) belowground NPP
    real(r8), pointer :: tempavg_agnpp                     (:) => null()
    real(r8), pointer :: tempavg_bgnpp                     (:) => null()
    real(r8), pointer :: annavg_agnpp                      (:) => null()
    real(r8), pointer :: annavg_bgnpp                      (:) => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type pf_carbonflux_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(pf_carbonflux_type) :: this
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
    class(pf_carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    allocate(this%rr                (begp:endp)) ; this%rr    (:)  = nan
    allocate(this%tempavg_agnpp(begp:endp)); this%tempavg_agnpp(:) = spval
    allocate(this%tempavg_bgnpp(begp:endp)); this%tempavg_bgnpp(:) = spval
    allocate(this%annavg_agnpp(begp:endp));  this%annavg_agnpp(:) = spval
    allocate(this%annavg_bgnpp(begp:endp)); this%annavg_bgnpp(:)  = spval

    allocate(this%annsum_npp      (begp:endp)) ; this%annsum_npp      (:) = spval
    allocate(this%agnpp                       (begp:endp)) ; this%agnpp                               (:) = spval
    allocate(this%bgnpp                       (begp:endp)) ; this%bgnpp                               (:) = spval
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
    class(pf_carbonflux_type) :: this
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

end module PfCarbonFluxType
