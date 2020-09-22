module BeTR_carbonfluxRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
  use betr_varcon     , only : spval => bspval
  use betr_ctrl, only : bgc_type
implicit none
#include "bshr_alloc.h"
  private
  character(len=*),  parameter :: mod_filename = &
       __FILE__
  type, public :: betr_carbonflux_recv_type
    real(r8), pointer :: hr_col(:) => null()
    real(r8), pointer :: som_c_leached_col(:) => null()
    real(r8), pointer :: som_c_runoff_col(:) => null()
    real(r8), pointer :: som_c_qdrain_col(:) => null()
    real(r8), pointer :: hr_vr_col(:,:) => null()
    real(r8), pointer :: fire_decomp_closs_vr_col(:,:) => null()  !will be summarized from the specific bgc model
    real(r8), pointer :: fire_decomp_closs_col(:) => null()  !will be summarized from the specific bgc model
    real(r8), pointer :: tempavg_agnpp_patch(:) => null()
    real(r8), pointer :: tempavg_bgnpp_patch(:) => null()
    real(r8), pointer :: annavg_agnpp_patch(:) => null()
    real(r8), pointer :: annavg_bgnpp_patch(:) => null()
    real(r8), pointer :: co2_soi_flx_col(:) => null()
    real(r8), pointer :: decomp_k(:,:,:) => null()
    real(r8), pointer :: phr_vr_col(:,:) => null()
    real(r8), pointer :: somhr_vr_col(:,:) => null()
    real(r8), pointer :: lithr_vr_col(:,:) => null()
    real(r8), pointer :: o_scalar_col(:,:) => null()
    real(r8), pointer :: cwdhr_vr_col(:,:) => null()
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_carbonflux_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  SPVAL_ALLOC(this%hr_col(begc:endc))
  SPVAL_ALLOC(this%fire_decomp_closs_col(begc:endc))
  SPVAL_ALLOC(this%som_c_leached_col(begc:endc))
  SPVAL_ALLOC(this%som_c_runoff_col(begc:endc))
  SPVAL_ALLOC(this%som_c_qdrain_col(begc:endc))
  SPVAL_ALLOC(this%hr_vr_col(begc:endc,lbj:ubj))
  SPVAL_ALLOC(this%fire_decomp_closs_vr_col(begc:endc,lbj:ubj))

  SPVAL_ALLOC(this%tempavg_agnpp_patch(begp:endp))
  SPVAL_ALLOC(this%tempavg_bgnpp_patch(begp:endp))
  SPVAL_ALLOC(this%annavg_agnpp_patch(begp:endp))
  SPVAL_ALLOC(this%annavg_bgnpp_patch(begp:endp))
  SPVAL_ALLOC(this%co2_soi_flx_col(begc:endc))
  SPVAL_ALLOC(this%phr_vr_col(begc:endc,lbj:ubj))
  SPVAL_ALLOC(this%o_scalar_col(begc:endc,lbj:ubj))
  SPVAL_ALLOC(this%somhr_vr_col(begc:endc,lbj:ubj))
  SPVAL_ALLOC(this%lithr_vr_col(begc:endc,lbj:ubj))
  SPVAL_ALLOC(this%cwdhr_vr_col(begc:endc,lbj:ubj))
!  if(index(bgc_type,'type1_bgc')/=0)then
    SPVAL_ALLOC(this%decomp_k(begc:endc, lbj:ubj,7)) !decomposition k parameter
!  endif
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_carbonflux_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%hr_vr_col(:,:) = value_column
  this%fire_decomp_closs_vr_col(:,:) = value_column
  this%som_c_leached_col(:) = value_column
  this%som_c_qdrain_col(:) = value_column
  this%som_c_runoff_col(:) = value_column
  this%co2_soi_flx_col(:)=value_column
  this%phr_vr_col(:,:)= value_column
  this%o_scalar_col(:,:)=value_column
  this%somhr_vr_col(:,:)=value_column
  this%lithr_vr_col(:,:)=value_column
  this%cwdhr_vr_col(:,:)=value_column
  end subroutine reset


  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)

  integer :: c, j

  this%hr_col(:) = 0._r8
  this%fire_decomp_closs_col(:) = 0._r8
  do j = lbj, ubj
    do c=bounds%begc, bounds%endc
      this%hr_col(c) = this%hr_col(c) + dz(c,j) * this%hr_vr_col(c,j)
      this%fire_decomp_closs_col(c) = this%fire_decomp_closs_col(c) + dz(c,j) * &
           this%fire_decomp_closs_vr_col(c,j)
    enddo
  enddo
  end subroutine summary
end module BeTR_carbonfluxRecvType
