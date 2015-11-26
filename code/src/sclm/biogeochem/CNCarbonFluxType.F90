module CNCarbonFluxType
  use clm_varcon     , only : spval, ispval, c14ratio
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varcon             , only : spval
  use clm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: carbonflux_type
    real(r8), pointer :: o_scalar_col                              (:,:)   ! fraction by which decomposition is limited by anoxia
    real(r8), pointer :: w_scalar_col                              (:,:)   ! fraction by which decomposition is limited by moisture availability
    real(r8), pointer :: t_scalar_col                              (:,:)   ! fraction by which decomposition is limited by temperature
    real(r8), pointer :: hr_vr_col                                 (:,:)   ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(r8), pointer :: bgc_cpool_ext_loss_vr_col                 (:, :, :)  ! col-level extneral organic carbon loss gC/m3 /time step
    real(r8), pointer :: bgc_cpool_ext_inputs_vr_col               (:, :, :)  ! col-level extneral organic carbon input gC/m3 /time step
    real(r8), pointer :: rr_col                                    (:)     ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: annsum_npp_patch            (:) ! patch annual sum of NPP (gC/m2/yr)
    real(r8), pointer :: annavg_agnpp_patch                        (:)     ! (gC/m2/s) annual average aboveground NPP
    real(r8), pointer :: annavg_bgnpp_patch                        (:)     ! (gC/m2/s) annual average belowground NPP
    real(r8), pointer :: agnpp_patch                               (:)     ! (gC/m2/s) aboveground NPP
    real(r8), pointer :: bgnpp_patch                               (:)     ! (gC/m2/s) belowground NPP
    real(r8), pointer :: tempavg_agnpp_patch                       (:)     ! (gC/m2/s) temp. average aboveground NPP
    real(r8), pointer :: tempavg_bgnpp_patch                       (:)     ! (gC/m2/s) temp. average belowground NPP
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type carbonflux_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(carbonflux_type) :: this
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
    class(carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%t_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%t_scalar_col (:,:)=spval
    allocate(this%w_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%w_scalar_col (:,:)=spval
    allocate(this%o_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%o_scalar_col (:,:)=spval
    allocate(this%hr_vr_col                         (begc:endc,1:nlevdecomp_full)); this%hr_vr_col    (:,:)=nan
    allocate(this%bgc_cpool_ext_inputs_vr_col       (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_inputs_vr_col (:,:,:) = nan
    allocate(this%bgc_cpool_ext_loss_vr_col         (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_loss_vr_col   (:,:,:) = nan
    allocate(this%rr_col                            (begc:endc))                  ; this%rr_col                    (:)  =nan
    allocate(this%annsum_npp_patch      (begp:endp)) ; this%annsum_npp_patch      (:) = nan
    allocate(this%annavg_agnpp_patch                (begp:endp))                  ; this%annavg_agnpp_patch  (:) = spval ! To detect first year
    allocate(this%annavg_bgnpp_patch                (begp:endp))                  ; this%annavg_bgnpp_patch  (:) = spval ! To detect first year
    allocate(this%agnpp_patch                               (begp:endp)) ; this%agnpp_patch                               (:) = nan
    allocate(this%bgnpp_patch                               (begp:endp)) ; this%bgnpp_patch                               (:) = nan
    allocate(this%tempavg_agnpp_patch               (begp:endp))                  ; this%tempavg_agnpp_patch (:) = spval
    allocate(this%tempavg_bgnpp_patch               (begp:endp))                  ; this%tempavg_bgnpp_patch (:) = spval
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
    class(carbonflux_type) :: this
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

end module CNCarbonFluxType
