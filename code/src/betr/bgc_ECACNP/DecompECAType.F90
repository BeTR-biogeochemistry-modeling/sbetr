module DecompECAType

!
! DESCRIPTIONS
! contains data structure for ECA-century decomposition.

! USES

  use shr_kind_mod , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type
implicit none


  type, public :: DecompECA_Type
  real(r8), pointer :: o_scalar_col                              (:,:)   ! fraction by which decomposition is limited by anoxia
  real(r8), pointer :: w_scalar_col                              (:,:)   ! fraction by which decomposition is limited by moisture availability
  real(r8), pointer :: t_scalar_col                              (:,:)   ! fraction by which decomposition is limited by temperature


contains
  procedure, public  :: Init
  procedure, private :: InitCold
  procedure, private :: InitAllocate

 end type DecompECA_Type


 contains

 !------------------------------------------------------------------------
 subroutine Init(this, bounds)

   class(DecompECA_Type) :: this
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
   use clm_varpar     , only : nlevdecomp_full
   use clm_varcon             , only : spval   
   !
   ! !ARGUMENTS:
   class(DecompECA_Type) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   integer :: begp, endp
   integer :: begc, endc
   allocate(this%t_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%t_scalar_col (:,:)=spval
   allocate(this%w_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%w_scalar_col (:,:)=spval
   allocate(this%o_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%o_scalar_col (:,:)=spval
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
    class(DecompECA_Type) :: this
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

end module DecompECAType
