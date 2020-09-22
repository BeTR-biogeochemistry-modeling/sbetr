module BeTR_phosphorusstateRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
  use tracer_varcon  , only : reaction_method
  use betr_varcon     , only : spval => bspval
  use betr_ctrl, only : bgc_type
implicit none
#include "bshr_alloc.h"
  private
  character(len=*), parameter :: mod_filename = &
       __FILE__
  type, public :: betr_phosphorusstate_recv_type
    real(r8), pointer :: cwdp_col(:) => null()
    real(r8), pointer :: totlitp_col(:) => null()
    real(r8), pointer :: totsomp_col(:) => null()
    real(r8), pointer :: totlitp_1m_col(:) => null()
    real(r8), pointer :: totsomp_1m_col(:) => null()
    real(r8), pointer :: sminp_col(:) => null()
    real(r8), pointer :: occlp_col(:) => null()
    real(r8), pointer :: domp_col(:) => null()

    real(r8), pointer :: cwdp_vr_col(:,:) => null()
    real(r8), pointer :: totlitp_vr_col(:,:) => null()
    real(r8), pointer :: totsomp_vr_col(:,:) => null()
    real(r8), pointer :: sminp_vr_col(:,:) => null()
    real(r8), pointer :: occlp_vr_col(:,:) => null()
    real(r8), pointer :: som1p_col(:) => null()
    real(r8), pointer :: som2p_col(:) => null()
    real(r8), pointer :: som3p_col(:) => null()
    real(r8), pointer :: som1p_vr_col(:,:) => null()
    real(r8), pointer :: som2p_vr_col(:,:) => null()
    real(r8), pointer :: som3p_vr_col(:,:) => null()
    real(r8), pointer :: domp_vr_col(:,:) => null()
    real(r8), pointer :: decomp_ppools_vr(:,:,:) => null()
    real(r8), pointer :: solutionp_vr_col(:,:) => null()
    real(r8), pointer :: totsoip(:) => null()
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_phosphorusstate_recv_type

 contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_phosphorusstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  implicit none
  class(betr_phosphorusstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  SPVAL_ALLOC(this%cwdp_col(begc:endc))
  SPVAL_ALLOC(this%totlitp_col(begc:endc))
  SPVAL_ALLOC(this%totsomp_col(begc:endc))
  SPVAL_ALLOC(this%totlitp_1m_col(begc:endc))
  SPVAL_ALLOC(this%totsomp_1m_col(begc:endc))
  SPVAL_ALLOC(this%sminp_col(begc:endc))
  SPVAL_ALLOC(this%occlp_col(begc:endc))
  SPVAL_ALLOC(this%domp_col(begc:endc))
  SPVAL_ALLOC(this%som1p_col(begc:endc))
  SPVAL_ALLOC(this%som2p_col(begc:endc))
  SPVAL_ALLOC(this%som3p_col(begc:endc))
  SPVAL_ALLOC(this%totsoip(begc:endc))
!  if(index(bgc_type,'type1_bgc')/=0)then
    SPVAL_ALLOC(this%decomp_ppools_vr(begc:endc, lbj:ubj, 1:7))
    SPVAL_ALLOC(this%solutionp_vr_col(begc:endc, lbj:ubj))
!  else
    SPVAL_ALLOC(this%som1p_vr_col(begc:endc, lbj:ubj))
    SPVAL_ALLOC(this%som2p_vr_col(begc:endc, lbj:ubj))
    SPVAL_ALLOC(this%som3p_vr_col(begc:endc, lbj:ubj))

    SPVAL_ALLOC(this%cwdp_vr_col(begc:endc,lbj:ubj))
    SPVAL_ALLOC(this%totlitp_vr_col(begc:endc,lbj:ubj))
    SPVAL_ALLOC(this%totsomp_vr_col(begc:endc,lbj:ubj))
    SPVAL_ALLOC(this%domp_vr_col(begc:endc,lbj:ubj))
!  endif
  SPVAL_ALLOC(this%sminp_vr_col(begc:endc,lbj:ubj))
  SPVAL_ALLOC(this%occlp_vr_col(begc:endc,lbj:ubj))

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_phosphorusstate_recv_type)  :: this
  real(r8), intent(in) :: value_column


  if(index(bgc_type,'type1_bgc')/=0)then
    this%decomp_ppools_vr(:,:,:) = value_column
    this%solutionp_vr_col(:,:) = value_column
  else
    this%cwdp_vr_col(:,:) = value_column
    this%totlitp_vr_col(:,:) = value_column
    this%totsomp_vr_col(:,:) = value_column

    this%som1p_vr_col(:,:) = value_column
    this%som2p_vr_col(:,:) = value_column
    this%som3p_vr_col(:,:) = value_column

    this%domp_vr_col(:,:) = value_column
  endif
  this%sminp_vr_col(:,:) = value_column
  this%occlp_vr_col(:,:) = value_column

  end subroutine reset

  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz,zs)

  implicit none
  class(betr_phosphorusstate_recv_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)
  real(r8), intent(in) :: zs(bounds%begc:bounds%endc,lbj:ubj)


  integer :: c, j

  if(index(bgc_type,'type1_bgc')/=0)return
  this%cwdp_col(:) = 0._r8
  this%totlitp_col(:) = 0._r8
  this%totsomp_col(:) = 0._r8
  this%totlitp_1m_col(:) = 0._r8
  this%totsomp_1m_col(:) = 0._r8
  this%sminp_col(:) = 0._r8
  this%occlp_col(:) = 0._r8
  this%domp_col(:) = 0._r8
  this%som1p_col(:) = 0._r8
  this%som2p_col(:) = 0._r8
  this%som3p_col(:) = 0._r8
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      this%som1p_col(c) =   this%som1p_col(c) + dz(c,j)*this%som1p_vr_col(c,j)
      this%som2p_col(c) =   this%som2p_col(c) + dz(c,j)*this%som2p_vr_col(c,j)
      this%som3p_col(c) =   this%som3p_col(c) + dz(c,j)*this%som3p_vr_col(c,j)
      this%domp_col(c)  =   this%domp_col(c) + dz(c,j)*this%domp_vr_col(c,j)
    enddo
  enddo

  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      this%cwdp_col(c) = this%cwdp_col(c) + dz(c,j) * this%cwdp_vr_col(c,j)
      this%totlitp_col(c) = this%totlitp_col(c) + dz(c,j)*this%totlitp_vr_col(c,j)
      this%totsomp_vr_col(c,j) = this%som1p_vr_col(c,j) + this%som2p_vr_col(c,j) + &
         this%som3p_vr_col(c,j) + this%domp_vr_col(c,j)
      this%totsomp_col(c) = this%totsomp_col(c) + dz(c,j)*this%totsomp_vr_col(c,j)
      this%sminp_col(c) = this%sminp_col(c) + dz(c,j)*this%sminp_vr_col(c,j)
      this%occlp_col(c) = this%occlp_col(c) + dz(c,j)*this%occlp_vr_col(c,j)

      this%som1p_col(c) =   this%som1p_col(c) + dz(c,j)*this%som1p_vr_col(c,j)
      this%som2p_col(c) =   this%som2p_col(c) + dz(c,j)*this%som2p_vr_col(c,j)
      this%som3p_col(c) =   this%som3p_col(c) + dz(c,j)*this%som3p_vr_col(c,j)
      if(zs(c,j)<1._r8)then
        if(zs(c,j+1)>1._r8)then
          this%totlitp_1m_col(c) = this%totlitp_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totlitp_vr_col(c,j)
          this%totsomp_1m_col(c) = this%totsomp_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totsomp_vr_col(c,j)
        else
          this%totlitp_1m_col(c) = this%totlitp_1m_col(c) + dz(c,j)*this%totlitp_vr_col(c,j)
          this%totsomp_1m_col(c) = this%totsomp_1m_col(c) + dz(c,j)*this%totsomp_vr_col(c,j)
        endif
      endif
    enddo
  enddo
  do c = bounds%begc, bounds%endc
     this%totsoip(c) = this%totlitp_col(c) + this%totsomp_col(c) + &
        this%cwdp_col(c) + this%sminp_col(c) + this%occlp_col(c)
  enddo
  end subroutine summary
end module BeTR_phosphorusstateRecvType
