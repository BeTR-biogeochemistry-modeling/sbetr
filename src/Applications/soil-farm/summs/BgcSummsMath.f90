module BgcSummsMath
  !DESCRIPTION
!module defines math functions needed
use bshr_kind_mod       , only : r8 => shr_kind_r8
use BgcSummsDebType     , only : debs

implicit none

contains
  subroutine interp1( xData, yData, xVal, yVal )
  ! Inputs: xData = a vector of the x-values of the data to be interpolated
  !         yData = a vector of the y-values of the data to be interpolated
  !         xVal  = x-value where interpolation should be performed
  ! Output: yVal  = interpolated value

    implicit none
    real(r8), intent(in) :: yData(:), xVal
    integer, intent(in) :: xData(:)
    real(r8), intent(out) :: yVal
    integer :: dataIndex, minXdata, maxXdata, xRange
    real(r8) :: weight

    ! Possible checks on inputs could go here
    ! Things you may want to check:
    !   monotonically increasing xData
    !   size(xData) == size(yData)
    !   size(xVal) == size(yVal)

    minXdata = xData(1)
    maxXdata = xData(size(xData))
    xRange = maxXdata - minXdata

        dataIndex = floor(xVal) - minXdata

        weight = (xVal - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
        yVal = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1)

  end subroutine interp1

  subroutine brent(x, x1,x2,f1, f2, macheps, tol, func_data, func) ! func_data are any other arguments to the function

    !!DESCRIPTION:
    !Use Brent's method to find the root to a single variable function func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol.

    !!REVISION HISTORY:
    !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
    !
    !!USES:
    use BgcSummsDebType    , only : debs
    !
    !!ARGUMENTS:
    implicit none
    real(r8), intent(in) :: x1, x2, f1, f2    !minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in) :: macheps           !machine precision
    real(r8), intent(in) :: tol               !the error tolerance
    type(debs), intent(in) :: func_data ! data passed to subroutine func
    real(r8), intent(out):: x

    interface
       subroutine func(x, func_data, f)
         use bshr_kind_mod        , only : r8 => shr_kind_r8
         use BgcSummsDebType      , only : debs
         implicit none
         real(r8), intent(in)  :: x
         type(debs), intent(in) :: func_data ! data passed to subroutine func
         real(r8), intent(out) :: f
       end subroutine func
    end interface

    ! !CALLED FROM:
    ! whenever it is needed
    integer, parameter :: ITMAX = 100            !maximum number of iterations
    integer :: iter
    real(r8)  :: a,b,c,d,e,fa,fb,fc,p,q,r,s,xm,tol1

    a=x1
    b=x2
    fa=f1
    fb=f2
    if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8))then
       !write(msg,*) 'root must be bracketed for brent', new_line('A')//'a=',a,' b=',b,' fa=',fa,' fb=',fb
       print *, 'root must be bracketed for brent'
    endif
    c=b
    fc=fb
    iter = 0
    do
       if(iter==ITMAX)exit
       iter=iter+1
       if((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8))then
          c=a   !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if( abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2._r8*macheps*abs(b)+0.5_r8*tol  !Convergence check.
       xm=0.5_r8*(c-b)
       if(abs(xm) <= tol1 .or. fb == 0._r8)then
          x=b
          return
       endif
       if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa !Attempt inverse quadratic interpolation.
          if(a == c) then
             p=2._r8*xm*s
             q=1._r8-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2._r8*xm*q*(q-r)-(b-a)*(r-1._r8))
             q=(q-1._r8)*(r-1._r8)*(s-1._r8)
          endif
          if(p > 0._r8) q=-q !Check whether in bounds.
          p=abs(p)
          if(2._r8*p < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e=d !Accept interpolation.
             d=p/q
          else
             d=xm  !Interpolation failed, use bisection.
             e=d
          endif
       else !Bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b !Move last best guess to a.
       fa=fb
       if(abs(d) > tol1) then !Evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       call func(b, func_data, fb)
       if(fb==0._r8)exit
    enddo
    if(iter==ITMAX)then
      print *, 'brent exceeding maximum iterations', b, fb
    endif
    x=b

  end subroutine brent


end module BgcSummsMath
