module LinearAlgebraMod

!DESCRIPTION
!code to do linear algebra

#include "bshr_assert.h"
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  implicit none
  private
  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: sparse_gemv
  public :: taxpy
  interface taxpy
    module procedure taxpy_v,taxpy_m
  end interface taxpy

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
#if (defined SBETR)
    SHR_ASSERT_ALL_EXT(((/nx,ny/) == (/nz,nb/)),   errMsg(mod_filename,__LINE__))
#endif
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
#if (defined SBETR)
    SHR_ASSERT_ALL_EXT(((/ny,nx/) == (/nz,nb/)),   errMsg(mod_filename,__LINE__))
#endif
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

!-------------------------------------------------------------------------------
 subroutine taxpy_v(N,DA,DX,INCX,DY,INCY)
 implicit none
 integer, intent(in) :: N
 real(r8), intent(in) :: da
 real(r8), intent(in) :: dx(:)
 integer, intent(in) :: incx
 integer, intent(in) :: incy
 real(r8), intent(inout):: dy(:)

 integer :: i, m, mp1, ix, iy
  if (n <= 0) return
  if (da == 0._r8) return
  if (incx == 1 .AND. incy == 1) then
    m = mod(n,5)
    if (m /= 0 )then
      do i = 1,m
        dy(i) = dy(i) + da*dx(i)
      enddo
    endif
    if (n < 5) return
    mp1 = m + 1
    do i = mp1,n,5
      dy(i) = dy(i) + da*dx(i)
      dy(i+1) = dy(i+1) + da*dx(i+1)
      dy(i+2) = dy(i+2) + da*dx(i+2)
      dy(i+3) = dy(i+3) + da*dx(i+3)
      dy(i+4) = dy(i+4) + da*dx(i+4)
    enddo
  else
    ix = 1
    iy = 1
    if (incx < 0) ix = (-n+1)*incx + 1
    if (incy < 0) iy = (-n+1)*incy + 1
    do i = 1,n
      dy(iy) = dy(iy) + da*dx(ix)
      ix = ix + incx
      iy = iy + incy
    enddo
  end if

 end subroutine taxpy_v

!-------------------------------------------------------------------------------
 subroutine taxpy_m(N,DA,DX,INCX,DY,INCY)
 implicit none
 integer, intent(in) :: N
 real(r8), intent(in) :: da
 real(r8), intent(in) :: dx(:,:)
 integer, intent(in) :: incx
 integer, intent(in) :: incy
 real(r8), intent(inout):: dy(:,:)

 integer :: szy, jj, nn

 szy = size(dx,2); nn = size(dx,1)
 do jj = 1, szy
   call taxpy_v(nn,da,dx(:,jj),incx,dy(:,jj),incy)
 enddo
 end subroutine taxpy_m


end module LinearAlgebraMod
