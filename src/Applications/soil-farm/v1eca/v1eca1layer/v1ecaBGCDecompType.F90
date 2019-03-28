module v1ecaBGCDecompType

!
! DESCRIPTIONS
! contains data structure for ECA-century decomposition.

! USES

  use bshr_kind_mod  , only : r8 => shr_kind_r8
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: DecompCent_type
  real(r8) :: o_scalar       ! fraction by which decomposition is limited by anoxia
  real(r8) :: decomp_k(7)
  !parameters
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding
  real(r8) :: minpsi
  real(r8) :: k_m_o2
  real(r8) :: vmax_decomp_n
  real(r8) :: vmax_decomp_p
  contains
    procedure, public  :: Init
    procedure, public  :: set_decompk_scalar
    procedure, private :: InitCold
    procedure, private :: InitAllocate
    procedure, public  :: UpdateParas
  end type DecompCent_type

 contains

 !------------------------------------------------------------------------
 subroutine Init(this, biogeo_con)

  use v1ecaParaType, only : v1eca_para_type
  implicit none
  class(DecompCent_type), intent(inout) :: this
  type(v1eca_para_type),intent(in) :: biogeo_con

  call this%InitAllocate ()

  call this%InitCold ()

 end subroutine Init
 !------------------------------------------------------------------------
 subroutine InitAllocate(this)
   !
   ! !DESCRIPTION:
   ! Initialize module data structure
   !
   ! !USES:
   use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
   use betr_varcon    , only : spval  => bspval
   !
   ! !ARGUMENTS:
   class(DecompCent_type), intent(inout) :: this

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this)
    !
    ! !USES:

    implicit none
    !
    ! !ARGUMENTS:
    class(DecompCent_type), intent(inout) :: this

    this%o_scalar = 1._r8

  end subroutine initCold

  !-----------------------------------------------------------------------
  subroutine UpdateParas(this, biogeo_con)

  use v1ecaParaType, only : v1eca_para_type
  implicit none
  class(DecompCent_type) , intent(inout) :: this
  type(v1eca_para_type)   , intent(in) :: biogeo_con

  ! set "Q10" parameter
  this%Q10 = biogeo_con%Q10
  this%froz_q10 = biogeo_con%froz_q10
  this%decomp_depth_efolding = biogeo_con%decomp_depth_efolding
  this%minpsi= biogeo_con%minpsi_bgc
  this%k_m_o2=biogeo_con%k_m_o2_bgc
  this%vmax_decomp_n=biogeo_con%vmax_decomp_n
  this%vmax_decomp_p=biogeo_con%vmax_decomp_p
  end subroutine UpdateParas
  !-----------------------------------------------------------------------
  subroutine set_decompk_scalar(this, o2b, bgc_forc)

  use JarBgcForcType , only : JarBGC_forc_type
  use bshr_const_mod     , only : SHR_CONST_TKFRZ
  implicit none
  ! !ARGUMENTS:
  class(DecompCent_type)     , intent(inout) :: this
  real(r8)                   , intent(in) :: o2b
  type(JarBGC_forc_type) , intent(in) :: bgc_forc

  ! !LOCAL VARIABLES:
  real(r8)   :: o2w
  integer    :: jj
  associate(                                      &
    o2_w2b        => bgc_forc%o2_w2b            , &
    k_m_o2        => this%k_m_o2                , &
    decomp_k      => bgc_forc%decomp_k            &
  )

  !oxygen scalar, this is different from what CLM4.5bgc does, I use a M-M formulation to indicate O2 stress
  !and the O2 budget is done on the fly
  o2w = o2b / o2_w2b
  this%o_scalar = max(o2w/(o2w+k_m_o2),1.e-20_r8)   !the value 0.22 mol O3/m3 is from Arah and Kirk, 2000
  do jj = 1, 7
    this%decomp_k(jj) = decomp_k(jj)
  enddo
  end associate
  end subroutine set_decompk_scalar


end module v1ecaBGCDecompType
