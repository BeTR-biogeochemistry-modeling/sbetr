module BeTR_carbonstateRecvType
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
  type, public :: betr_carbonstate_recv_type
     real(r8), pointer :: cwdc_col(:) => null()
     real(r8), pointer :: totlitc_col(:)  => null()
     real(r8), pointer :: totsomc_col(:) => null()
     real(r8), pointer :: totlitc_1m_col(:) => null()
     real(r8), pointer :: totsomc_1m_col(:) => null()
     real(r8), pointer :: cwdc_vr_col(:,:) => null()
     real(r8), pointer :: totlitc_vr_col(:,:) => null()
     real(r8), pointer :: totsomc_vr_col(:,:) => null()
     real(r8), pointer :: som1c_col(:) => null()
     real(r8), pointer :: som2c_col(:) => null()
     real(r8), pointer :: som3c_col(:) => null()
     real(r8), pointer :: domc_col(:) => null()
     real(r8), pointer :: som1c_vr_col(:,:) => null()
     real(r8), pointer :: som2c_vr_col(:,:) => null()
     real(r8), pointer :: som3c_vr_col(:,:) => null()
     real(r8), pointer :: domc_vr_col(:,:) => null()
     real(r8), pointer :: decomp_cpools_vr(:,:,:) => null()
     real(r8), pointer :: totomc(:) => null()
 contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_carbonstate_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_carbonstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)

  implicit none
  class(betr_carbonstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  SPVAL_ALLOC(this%cwdc_col(begc:endc))
  SPVAL_ALLOC(this%totlitc_col(begc:endc))
  SPVAL_ALLOC(this%totsomc_col(begc:endc))
  SPVAL_ALLOC(this%totlitc_1m_col(begc:endc))
  SPVAL_ALLOC(this%totsomc_1m_col(begc:endc))

  SPVAL_ALLOC(this%som1c_col(begc:endc))
  SPVAL_ALLOC(this%som2c_col(begc:endc))
  SPVAL_ALLOC(this%som3c_col(begc:endc))
  SPVAL_ALLOC(this%domc_col(begc:endc))
  SPVAL_ALLOC(this%totomc(begc:endc))

    SPVAL_ALLOC(this%decomp_cpools_vr(begc:endc,lbj:ubj,7))
    SPVAL_ALLOC(this%som1c_vr_col(begc:endc, lbj:ubj))
    SPVAL_ALLOC(this%som2c_vr_col(begc:endc, lbj:ubj))
    SPVAL_ALLOC(this%som3c_vr_col(begc:endc, lbj:ubj))

    SPVAL_ALLOC(this%cwdc_vr_col(begc:endc,lbj:ubj))
    SPVAL_ALLOC(this%totlitc_vr_col(begc:endc,lbj:ubj))
    SPVAL_ALLOC(this%totsomc_vr_col(begc:endc,lbj:ubj))
    SPVAL_ALLOC(this%domc_vr_col(begc:endc, lbj:ubj))

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_carbonstate_recv_type)  :: this
  real(r8), intent(in) :: value_column

  if(index(bgc_type,'type1_bgc')/=0)then
    this%decomp_cpools_vr(:,:,:)=value_column
  else
    this%cwdc_vr_col(:,:) = value_column
    this%totlitc_vr_col(:,:) = value_column
    this%totsomc_vr_col(:,:) = value_column

    this%som1c_vr_col(:,:) = value_column
    this%som2c_vr_col(:,:) = value_column
    this%som3c_vr_col(:,:) = value_column
    this%domc_vr_col(:,:) = value_column
  endif
  end subroutine reset

  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz, zs)

  implicit none
  class(betr_carbonstate_recv_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)
  real(r8), intent(in) :: zs(bounds%begc:bounds%endc,lbj:ubj)

  integer :: c, j

  if(index(bgc_type,'type1_bgc')/=0)return
  this%cwdc_col(:) = 0._r8
  this%totlitc_col(:) = 0._r8
  this%totsomc_col(:) = 0._r8
  this%totlitc_1m_col(:) = 0._r8
  this%totsomc_1m_col(:) = 0._r8

  this%som1c_col(:) = 0.0_r8
  this%som2c_col(:) = 0.0_r8
  this%som3c_col(:) = 0.0_r8
  this%domc_col(:)  = 0._r8
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      this%som1c_col(c) = this%som1c_col(c) + dz(c,j)*this%som1c_vr_col(c,j)
      this%som2c_col(c) = this%som2c_col(c) + dz(c,j)*this%som2c_vr_col(c,j)
      this%som3c_col(c) = this%som3c_col(c) + dz(c,j)*this%som3c_vr_col(c,j)
      this%cwdc_col(c) = this%cwdc_col(c) + dz(c,j) * this%cwdc_vr_col(c,j)
      this%domc_col(c) = this%domc_col(c) + dz(c,j) * this%domc_vr_col(c,j)

      this%totlitc_col(c) = this%totlitc_col(c) + dz(c,j)*this%totlitc_vr_col(c,j)
      this%totsomc_vr_col(c,j) = this%som1c_vr_col(c,j) + this%som2c_vr_col(c,j) + &
         this%som3c_vr_col(c,j) + this%domc_vr_col(c,j)
      this%totsomc_col(c) = this%totsomc_col(c) + dz(c,j)*this%totsomc_vr_col(c,j)
    enddo
  enddo
  do c = bounds%begc, bounds%endc
     this%totomc(c) = this%totlitc_col(c) + this%totsomc_col(c) + &
        this%cwdc_col(c)
  enddo
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      if(zs(c,j)<1._r8)then
        if(zs(c,j+1)>1._r8)then
          this%totlitc_1m_col(c) = this%totlitc_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totlitc_vr_col(c,j)
          this%totsomc_1m_col(c) = this%totsomc_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totsomc_vr_col(c,j)
        else
          this%totlitc_1m_col(c) = this%totlitc_1m_col(c) + dz(c,j)*this%totlitc_vr_col(c,j)
          this%totsomc_1m_col(c) = this%totsomc_1m_col(c) + dz(c,j)*this%totsomc_vr_col(c,j)
        endif
      endif
    enddo
  enddo

  end subroutine summary
end module BeTR_carbonstateRecvType
