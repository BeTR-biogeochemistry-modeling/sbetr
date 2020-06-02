module InterpolationMod
#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  ! subroutines to do polynomial interpolation
  ! author: Jinyun Tang, Sep, 2014

  ! !USES:
  use bshr_assert_mod, only : shr_assert
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use bshr_assert_mod, only : shr_assert_any
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg

  implicit none

  private

  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: Lagrange_interp
  public :: mono_Linear_interp_trj
  public :: pchip_polycc
  public :: pchip_interp
  public :: cmass_interp
  public :: loc_x, loc_xj
  public :: bmass_interp, mass_interp
  public :: layer_adjust
contains

  !-------------------------------------------------------------------------------
  subroutine mono_Linear_interp_trj(nx, x, y, nxi, xi, yi, bstatus)
    !
    ! !DESCRIPTION:
    ! do order pn lagrangian interpolation
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    integer,    intent(in)  :: nx    !order of interpolation
    real(r8),   intent(in)  :: x(1:nx)   !location of data
    real(r8),   intent(in)  :: y(1:nx)   !value of data
    integer,    intent(in)  :: nxi
    real(r8),   intent(in)  :: xi(1:nxi)  !target points to be evaluated
    real(r8),   intent(out) :: yi(1:nxi) !target values
    type(betr_status_type), intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer :: k
    integer :: pos, disp, disp1


    call bstatus%reset()
    do k = 1, nxi
       ! find the position of z in array x
       !pos = find_idx(x, xi(k), pos)
       if(k==1)then
         pos=loc_x(nx, x, xi(1), bstatus)
       else
         pos=loc_xj(nx, x, xi(k), pos)
       endif

       if(pos==0)then
         yi(k)=y(1)
       elseif(pos==nx)then
         yi(k)=y(nx)
       else
         ! call linear interpolation
         if(pos<0)then
           yi(k) = x(1)
         else
           yi(k)=twopoint_linterp(x(pos),x(pos+1),Y(pos),Y(pos+1), xi(k))
         endif
       endif
    enddo

  end subroutine mono_Linear_interp_trj

  !-------------------------------------------------------------------------------
  subroutine Lagrange_interp(pn, x, y, xi, yi, bstatus)
    !
    ! !DESCRIPTION:
    ! do order pn lagrangian interpolation
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    integer,                intent(in)  :: pn    !order of interpolation
    real(r8), dimension(:), intent(in)  :: x   !location of data
    real(r8), dimension(:), intent(in)  :: y   !value of data
    real(r8), dimension(:), intent(in)  :: xi  !target points to be evaluated
    real(r8), dimension(:), intent(out) :: yi !target values
    type(betr_status_type), intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer :: k, ni, nx
    integer :: pos, disp, disp1

    call bstatus%reset()
    SHR_ASSERT_ALL((ubound(x) == ubound(y)),   errMsg(mod_filename,__LINE__), bstatus)

    SHR_ASSERT_ALL((ubound(xi) == ubound(yi)), errMsg(mod_filename,__LINE__), bstatus)

    SHR_ASSERT_ALL((ubound(x) >= pn+1),        errMsg(mod_filename,__LINE__), bstatus)

    ni  = size(xi)
    nx = size(x)
    disp=int((pn+1)*0.5_r8+1.e-8_r8)
    !get the half size of the local window
    if(mod(pn,2)==0)then
       disp1=disp
    else
       disp1=disp-1
    endif
    pos=0
    do k = 1, ni
       ! find the position of z in array x
       pos = find_idx(x, xi(k), pos)
       if(pos == -100) then
          !left boundary
          yi(k) = y(1)
       elseif(pos == -200) then
          !right boundary
          yi(k) = y(nx)
       else
          ! call function Lagrange
          if (pos <= disp1) then
             yi(k) = Lagrange_poly(pn, x(1:pn+1), y(1:pn+1), xi(k), bstatus)
             if(bstatus%check_status())return
          else if (pos >= nx-disp) then
             yi(k) = Lagrange_poly(pn, x(nx-pn:nx), y(nx-pn:nx), xi(k), bstatus)
             if(bstatus%check_status())return
          else
             yi(k) = Lagrange_poly(pn, x(pos-disp1:pos+disp), y(pos-disp1:pos+disp), xi(k), bstatus)
             if(bstatus%check_status())return
          end if
       endif
    enddo

  end subroutine Lagrange_interp

  !-------------------------------------------------------------------------------
  function Lagrange_poly(pn, xvect, yvect, z, bstatus)result(Pz)
    !
    ! !DESCRIPTION:
    ! do lagrangian interpolation at order pn
    !
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    integer,                intent(in) :: pn   ! Order of Interpolation Polynomial
    real(r8), dimension(:), intent(in) :: xvect, yvect  ! vectors of known data: x,y-values
    real(r8),               intent(in) :: z   ! the target point "z"
    type(betr_status_type), intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer   :: i, j, n
    real(r8)  :: L(pn+1)    ! Lagrange cardinal function
    real(r8)  :: Pz    ! target value

    call bstatus%reset()
    SHR_ASSERT_ALL((size(xvect) == size(yvect)), errMsg(mod_filename,__LINE__), bstatus)

    SHR_ASSERT_ALL((size(xvect) == pn+1),        errMsg(mod_filename,__LINE__), bstatus)

    ! n = number of data points:length of each data vector
    n = size(xvect)
    ! Initializations of Pz and L
    Pz = 0._r8 ! initializing the polynomia value at z
    L(:)  = 1._r8  ! initalizing the vector of cardinal functions to 1
    ! Performing the interpolation
    do i = 1, n
       do j = 1, n
          if (i /= j) then
             ! part of L(i)
             L(i) = ( (z - xvect(j)) / (xvect(i) - xvect(j)) )* L(i)
          end if
       end do
       Pz = Pz + L(i)*yvect(i) ! update Pz ~ f(z)
    end do
  end function Lagrange_poly
  !------------------------------------------------------------
  function find_idx(xvect, x, k0)result(k)
    !
    ! !DESCRIPTION:
    ! locate the position of x in xvect
    !
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: xvect       ! vector of x-values
    real(r8),               intent(in) :: x
    integer ,               intent(in) :: k0
    integer :: i, k, n

    ! array dimension
    n = size(xvect)

    if(x<xvect(1))then
       k = -100  !beyond left boundary
    elseif(x>xvect(n))then
       k = -200  !beyond right boundary
    elseif(x==xvect(1))then
       k=1
    elseif(x==xvect(n))then
       k=n-1
    else
       ! find index k so that x[k] <= x < x[k+1]
       do i = max(k0,1), n-1
          if ((xvect(i) <= x) .and. (x < xvect(i+1))) then
             k = i
             exit
          end if
       end do
    endif


  end function find_idx
  !------------------------------------------------------------
  subroutine pchip_polycc(x, fx, di, bstatus, region)
    !
    ! DESCRIPTION
    ! Given the data, generate the coefficients of the monotonic cubic
    ! polynomials
    ! Ref, Fritsch and Carlson, 1980
    !
    ! !USES:
    use MathfuncMod      , only : diff
    use BetrStatusType   , only : betr_status_type
    use betr_constants   , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: x
    real(r8), dimension(:), intent(in) :: fx
    real(r8), dimension(:), intent(out):: di
    integer,  optional,     intent(in) :: region
    type(betr_status_type), intent(out):: bstatus

    ! !LOCAL VARIABLES:
    real(r8), allocatable :: h(:)
    real(r8), allocatable :: df(:)
    real(r8), allocatable :: slp(:)
    real(r8) :: alpha, beta, tao, rr
    integer :: region_loc
    integer :: n, j
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    SHR_ASSERT_ALL((size(x) == size(fx)), errMsg(mod_filename,__LINE__), bstatus)

    SHR_ASSERT_ALL((size(x) == size(di)), errMsg(mod_filename,__LINE__), bstatus)

    region_loc=2
    if(present(region))region_loc=region

    n = size(x)
    allocate(h(n-1))
    allocate(df(n-1))
    allocate(slp(n-1))
    !get interval length
    call diff(x, h, bstatus)
    if(bstatus%check_status())return
    !get function step
    call diff(fx, df, bstatus)
    if(bstatus%check_status())return
    !get slope
    do j = 1, n-1
       slp(j)=df(j)/h(j)
    enddo

    !get di
    di(:) = 0._r8

    j = 1
    di(j)=(fx(j+1)+fx(j+2)-2*fx(1))/(2*h(j)+h(j+1))
    do j = 2, n-1
       di(j)=(fx(j+1)-fx(j-1))/(h(j)+h(j-1))
    enddo
    j = n
    di(j)=(2._r8*fx(j)-(fx(j-1)+fx(j-2)))/(2._r8*h(j-1)+h(j-2))

    !enforce the sign condition
    if(slp(1)*di(1)<=0._r8)then
       di(1)=0._r8
    endif

    if(slp(n-1)*di(n)<=0)then
       di(n)=0._r8
    endif

    !enforce the range 2 constraint

    do j = 1, n-1
       if(abs(slp(j))<=1.e-16_r8)then
          di(j)=0._r8
          di(j+1)=0._r8
       else
          alpha=di(j)/slp(j)
          beta =di(j+1)/slp(j)
          select case (region_loc)
          case (1)
             rr=beta/alpha
             if(rr>1._r8)then
                if(beta>3._r8)then
                   beta=3._r8
                   alpha=beta/rr
                   di(j)=slp(j)*alpha
                   di(j+1)=slp(j)*beta
                endif
             else
                if(alpha>3._r8)then
                   alpha=3._r8
                   beta=rr*alpha
                   di(j)=slp(j)*alpha
                   di(j+1)=slp(j)*beta
                endif
             endif
          case (2)
             tao=3._r8/sqrt(alpha*alpha+beta*beta)
             if(tao<1._r8)then
                di(j)=tao*di(j)
                di(j+1)=tao*di(j+1)
             endif
          case (3)
             if(alpha+beta>3._r8)then
                if(alpha>0._r8)then
                   rr=beta/alpha
                   alpha=3._r8/(1._r8+rr)
                   beta=alpha*rr
                   di(j)=slp(j)*alpha
                   di(j+1)=slp(j)*beta
                else
                   beta=3._r8
                   di(j+1)=slp(j)*beta
                endif
             endif
          case (4)
             if(alpha>0._r8)then
                rr=beta/alpha
                if(rr>=1._r8)then
                   if(2._r8*alpha+beta>3._r8)then
                      alpha=3._r8/(2._r8+rr)
                      beta=alpha*rr
                      di(j)=slp(j)*alpha
                      di(j+1)=slp(j)*beta;
                   endif
                else
                   if(alpha+2._r8*beta>3._r8)then
                      alpha=3._r8/(1._r8+2._r8*rr)
                      beta=alpha*rr
                      di(j)=slp(j)*alpha
                      di(j+1)=slp(j)*beta
                   endif
                endif
             else
                if(beta>3._r8)then
                   beta=3._r8
                   di(j+1)=slp(j)*beta
                endif
             endif
          case default
             msg='an constraint region must be specified for pchip_polycc '//errMsg(mod_filename, __LINE__)
             call bstatus%set_msg(msg=msg, err=-1)
             return
          end select
       endif
    enddo

    deallocate(h)
    deallocate(df)
    deallocate(slp)
  end subroutine pchip_polycc
  !------------------------------------------------------------

  subroutine pchip_interp(x, fx, di, xi, yi, bstatus)

    ! !DESCRIPTION:
    ! do monotonic cubic spline interpolation
    use BetrStatusType   , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: x
    real(r8), dimension(:), intent(in) :: fx
    real(r8), dimension(:), intent(in) :: di
    real(r8), dimension(:), intent(in) :: xi
    real(r8), dimension(:), intent(out) :: yi
    type(betr_status_type), intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: h, t1, t2
    real(r8) :: h1x, h2x, h3x, h4x
    integer  :: n, j
    integer  :: id

    call bstatus%reset()
    SHR_ASSERT_ALL((size(x) == size(fx)),  errMsg(mod_filename,__LINE__), bstatus)

    SHR_ASSERT_ALL((size(x) == size(di)),  errMsg(mod_filename,__LINE__), bstatus)

    SHR_ASSERT_ALL((size(xi) == size(yi)), errMsg(mod_filename,__LINE__), bstatus)

    n=size(xi)  !total number of points to be interpolated

    yi(:)=0._r8
    id=0
    do j = 1, n

       id=find_idx(x,xi(j),id)
       h=x(id+1)-x(id)
       t1=(x(id+1)-xi(j))/h
       t2=(xi(j)-x(id))/h

       h1x=phi(t1)
       h2x=phi(t2)
       h3x=-h*psi(t1)
       h4x=h*psi(t2)
       yi(j)=fx(id)*h1x+fx(id+1)*h2x+di(id)*h3x+di(id+1)*h4x
    enddo

  contains

    function phi(t) result(fval)
      implicit none
      real(r8), intent(in) :: t

      real(r8) :: fval

      fval=(3._r8-2._r8*t)*t*t

    end function phi

    function psi(t) result(fval)
      implicit none
      real(r8), intent(in) :: t

      real(r8) :: fval
      fval=t*t*(t-1._r8)
    end function psi
  end subroutine pchip_interp

  !------------------------------------------------------------
  function loc_x(m, x, xi, bstatus)result(i_loc)

  ! !DESCRIPTION:
  ! locate xi within the input vector x
  use BetrStatusType   , only : betr_status_type
  use betr_constants   , only : betr_errmsg_len
  implicit none
  integer , intent(in) :: m
  real(r8), intent(in) :: x(1:m)
  real(r8), intent(in) :: xi
  type(betr_status_type), intent(out) :: bstatus

  !local variables
  integer :: jj
  integer :: i_loc
  character(len=betr_errmsg_len) :: msg

  call bstatus%reset()

  i_loc = -1
  if (xi < x(1) .or. xi > x(m)) then
    print*,'xi,xl,xr',xi,x(1),x(m)
    msg='the target value is outside the range of given x '//errMsg(mod_filename, __LINE__)
    call bstatus%set_msg(msg=msg, err=-1)
    return
  else
    if(xi==x(1))then
      i_loc=1
      return
    elseif(xi==x(m))then
      i_loc=m
      return
    endif

    do jj = 1, m-1
      if ( xi >= x(jj) .and. xi<x(jj+1)) then
        i_loc=jj
        exit
      endif
    enddo

  endif
  return
  end function loc_x

  !------------------------------------------------------------
  function loc_xj(m, x, xi, j0)result(i_loc)

  ! !DESCRIPTION:
  ! locate xi within the input vector x
  use BetrStatusType   , only : betr_status_type
  use betr_constants   , only : betr_errmsg_len
  implicit none
  integer , intent(in) :: m
  real(r8), intent(in) :: x(1:m)
  real(r8), intent(in) :: xi
  integer , intent(in) :: j0

  !local variables
  integer :: jj
  integer :: i_loc

  i_loc=-1
  do jj = max(j0-1,1), m-1
    if ( xi >= x(jj) .and. xi<x(jj+1)) then
      i_loc=jj
      exit
    endif
  enddo
  if(i_loc==-1)then
    if(xi==x(m))i_loc=m-1
  endif
  return
  end function loc_xj


  !------------------------------------------------------------
  function twopoint_linterp(x1,x2,f1,f2, xi)result(fy)

  ! !DESCRIPTION:
  ! two point linear interpolation
  implicit none
  real(r8), intent(in) :: x1, x2
  real(r8), intent(in) :: f1, f2
  real(r8), intent(in) :: xi

  !temporary variables
  real(r8) :: fy
  real(r8) :: d
  real(r8), parameter :: tiny_val=1.e-20_r8

  d=x2-x1
  if(abs(d)<tiny_val)then
    fy=(f1+f2)*0.5_r8
  else
    fy=f2*(xi-x1)/d+f1*(x2-xi)/d
  endif

  end function twopoint_linterp

  !------------------------------------------------------------

  subroutine cmass_interp(nx, x, ny, Y, nxi, xi, YI, bstatus)

  ! DESCRIPTION
  ! linear monotonic interpolation based on cumulative mass
  ! locate xi within the input vector x
  use BetrStatusType   , only : betr_status_type
  implicit none
  integer , intent(in) :: nx
  real(r8), intent(in) :: x(1:nx)
  integer , intent(in) :: ny
  real(r8), intent(in) :: Y(1:nx,1:ny)
  integer , intent(in) :: nxi
  real(r8), intent(in) :: xi(1:nxi)
  real(r8), intent(out):: Yi(1:nxi,1:ny)
  type(betr_status_type), intent(out) :: bstatus

  !local variables
  integer :: j0, j1, j2

  !check for end point
  j0=loc_x(nx, x, xi(nxi), bstatus)
  if(bstatus%check_status())return

  !check for start point
  j0=loc_x(nx, x, xi(1), bstatus)
  if(bstatus%check_status())return

  do j1 = 1, nxi
    if (j1/=1) then
      j0=loc_xj(nx, x, xi(j1), j0)
    endif
    if(j0==nx)then
      do j2 = 1, ny
        yi(j1,j2)=y(nx,j2)
      enddo
    else
      do j2 = 1, ny
        yi(j1,j2)=twopoint_linterp(x(j0),x(j0+1),Y(j0,j2),Y(j0+1,j2), xi(j1))
      enddo
    endif
  enddo
  end subroutine cmass_interp

  !------------------------------------------------------------
  subroutine mass_interp(zh, mass_curve,zl, zr, mass_new, bstatus)

  !DESCRIPTION
  !compute the tracer mass encompassed by left and right boundaries
  !this interpolation is for inner node
  use BetrStatusType   , only : betr_status_type
  use MathfuncMod      , only : cumsum
  implicit none
  real(r8), dimension(:), intent(in) :: zh
  real(r8), dimension(:), intent(in) :: mass_curve
  real(r8),               intent(in) :: zl
  real(r8),               intent(in) :: zr
  real(r8),               intent(out):: mass_new
  type(betr_status_type), intent(out)   :: bstatus

  real(r8) :: massl
  real(r8) :: massr
  real(r8) :: cmass_curve(size(zh))
  integer  :: nn

  call bstatus%reset()
  SHR_ASSERT_ALL((size(zh)  == size(mass_curve)),   errMsg(mod_filename,__LINE__),bstatus)

  nn = size(zh)

  call cumsum(bstatus, mass_curve, cmass_curve)
  if(bstatus%check_status())return

  massl=twopoint_linterp(zh(1),zh(2),cmass_curve(1),cmass_curve(2), zl)
  massr=twopoint_linterp(zh(nn),zh(nn-1),cmass_curve(nn),cmass_curve(nn-1), zr)

  mass_new=max(massr-massl,0._r8)

  end subroutine mass_interp

  !------------------------------------------------------------
  subroutine bmass_interp(zh, mass_curve,zl, zr, mass_new, bstatus)
  !DESCRIPTION
  !compute the tracer mass encompassed by left and right boundaries
  !this interpolation is for boundary node
  use BetrStatusType   , only : betr_status_type
  use MathfuncMod      , only : cumsum
  implicit none
  real(r8), dimension(:), intent(in) :: zh
  real(r8), dimension(:), intent(in) :: mass_curve
  real(r8),               intent(in) :: zl
  real(r8),               intent(in) :: zr
  real(r8),               intent(out):: mass_new
  type(betr_status_type), intent(out)   :: bstatus

  real(r8) :: massl
  real(r8) :: massr
  real(r8) :: cmass_curve(size(zh))
  integer  :: nn
  integer  :: j0

  call bstatus%reset()
  SHR_ASSERT_ALL((size(zh)  == size(mass_curve)),   errMsg(mod_filename,__LINE__),bstatus)

  nn = size(zh)

  call cumsum(bstatus, mass_curve, cmass_curve)
  if(bstatus%check_status())return

  j0=loc_x(nn, zh, zl, bstatus)
  if(bstatus%check_status())return

  massl=twopoint_linterp(zh(j0),zh(j0+1),cmass_curve(j0),cmass_curve(j0+1), zl)

  j0=loc_xj(nn, zh, zr, j0)
  massr=twopoint_linterp(zh(j0),zh(j0+1),cmass_curve(j0),cmass_curve(j0+1), zr)

  mass_new=max(massr-massl,0._r8)

  end subroutine bmass_interp

  !------------------------------------------------------------
  function loc_layer(z1,zt,nlen,k0)result(k1)
  implicit none
  integer , intent(in) :: nlen
  real(r8), intent(in) :: z1(nlen)
  real(r8), intent(in) :: zt
  integer , intent(in) :: k0

  integer :: k1
  integer :: j

  k1 = 0
  do j = 1, nlen
    if (zt <= z1(j))then
      k1 = j
    endif
  enddo
  end function loc_layer
  !------------------------------------------------------------

  subroutine layer_adjust(z1, zt, len1, len2, rmat)

  implicit none
  integer , intent(in) :: len1, len2
  real(r8), intent(in) :: z1(len1)
  real(r8), intent(in) :: zt(len2)
  real(r8), intent(out) :: rmat(len1,len2)

  integer :: k0, k1, j, kk, k

  k0 = 1
  do j = 1, len2
    k1 = loc_layer(z1,zt(j),len1,k0)
    if (j == 1) then
      if (k1 == k0) then
        rmat(j,k1) = zt(j)/z1(k0)
      else
        !middle
        do k=k0, k1-1
          rmat(j,k) = 1._r8
        enddo
        rmat(j,k1)=(zt(j)-z1(k1-1))/(z1(k1)-z1(k1-1))
      endif
    else
      if (k0 == k1)then
        if (k0==1)then
          rmat(j,k1) = (zt(j)-zt(j-1))/z1(k0)
        else
          rmat(j,k1) = (zt(j)-zt(j-1))/(z1(k0)-z1(k0-1))
        endif
      else
        if (k0 == 1)then
          rmat(j,k0) = (z1(k0)-zt(j-1))/z1(k0)
        else
          rmat(j,k0) = (z1(k0)-zt(j-1))/(z1(k0)-z1(k0-1))
        endif
        do kk = k0+1,k1-1
          rmat(j,kk) = 0._r8
        enddo
        rmat(j,k1) = (zt(j)-z1(k1-1))/(z1(k1)-z1(k-1))
      endif
    endif
  enddo

  end subroutine layer_adjust
end module InterpolationMod
