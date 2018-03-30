module LinearAlgebraMod

!DESCRIPTION
!code to do linear algebra

#include "bshr_assert.h"
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  implicit none
  private
  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: sparse_gemv
contains

  subroutine sparse_gemv(transp, nx, ny, a, nb, b, nz, dxdt)

  implicit none
  character(len=1), intent(in) :: transp
  integer , intent(in) :: nx, ny, nb, nz
  real(r8), intent(in) :: a(1:nx, 1:ny)
  real(r8), intent(in) :: b(1:ny)
  real(r8), intent(out) :: dxdt(1:nz)
  real(r8), parameter :: tiny_val= 1.e-12_r8

  integer :: ii, jj
  if(transp == 'N')then
    !compute dxdt = a * b
    SHR_ASSERT_ALL_EXT(((/nx,ny/) == (/nz,nb/)),   errMsg(mod_filename,__LINE__))
    dxdt(:) = 0._r8
    do jj = 1, ny
      do ii = 1, nx
        if(abs(a(ii,jj))>tiny_val .and. abs(b(jj))>tiny_val)then
          dxdt(ii) = dxdt(ii) + a(ii,jj) * b(jj)
        endif
      enddo
    enddo
  else
    !compute dxdt = a' * b
    SHR_ASSERT_ALL_EXT(((/ny,nx/) == (/nz,nb/)),   errMsg(mod_filename,__LINE__))
    dxdt(:) = 0._r8
    do jj = 1, ny
      do ii = 1, nx
        if(abs(a(ii,jj))>tiny_val .and. abs(b(ii))>tiny_val)then
          dxdt(jj) = dxdt(jj) + a(ii,jj) * b(ii)
        endif
      enddo
    enddo
  endif

  end subroutine sparse_gemv
end module LinearAlgebraMod
