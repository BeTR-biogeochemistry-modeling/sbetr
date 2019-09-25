module eChemMod

#include "bshr_assert.h"
  use FindRootMod     , only : imaxloc, cubic_newtonraphson
  use FindRootMod, only : quadrootbnd
  use BetrStatusType  , only : betr_status_type
  use eChemConstMod   , only : r_echem_AKD, r_echem_ABKDE
  use eChemConstMod   , only : r_echem_ABKD,r_echem_DsKAB
  use eChemConstMod   , only : r_echem_DsKA2B,r_echem_DsKAB2
  use eChemConstMod   , only : r_echem_DsK_generic,r_echem_K_generic
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  use FindRootMod     , only : quadrootbnd
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

contains

  function echem_DKA(ca,cd,Kr)result(dA)
  ! A <-> D
  ! (cd+dA)/(ca-dA)=Kr
  implicit none
  real(r8), intent(in) :: ca, cd, Kr
  real(r8) :: dA

  dA=(ca*Kr-cd)/(1._r8+Kr)

  end function echem_DKA

  !-------------------------------------------------------------------------------

  function echem_ABKDE(ca,cb,cd,ce, Kr)result(dA)
  ! A + B <-> D + E
  ! (cd+dA)*(ce+dA)/((ca-dA)*(cb-dA))=Kr
  use FindRootMod, only : quadrootbnd
  implicit none
  real(r8), intent(in) :: ca, cb, cd, ce, Kr
  real(r8) :: dA
  real(r8) :: a, b, c, xl, xr
  type(betr_status_type) :: bstatus

  a = -1._r8+Kr
  b = -(Kr*(ca+cb)+(cd+ce))
  c = Kr*ca*cb-cd*ce
  xr =  min(ca, cb)
  xl = -min(cd, ce)
  call quadrootbnd(a,b,c, xl, xr, dA, bstatus)

  end function echem_ABKDE
  !-------------------------------------------------------------------------------

  function echem_ABKD(ca,cb,cd, Kr)result(dA)
  ! A + B <-> D
  ! (cd+dA)/((ca-dA)*(cb-dA))=Kr
  implicit none
  real(r8), intent(in) :: ca, cb, cd, Kr
  real(r8) :: dA
  real(r8) :: a, b, c, xl, xr
  type(betr_status_type) :: bstatus

  a = Kr
  b = -(Kr*(ca+cb)+1._r8)
  c = Kr*ca*cb-cd
  xr = min(ca, cb)
  xl = -cd
  call quadrootbnd(a,b,c, xl, xr, dA, bstatus)

  end function echem_ABKD
  !-------------------------------------------------------------------------------

  function echem_DsKAB(cd,ca,cb,Kr)result(dA)
  ! D(s) <-> A+B
  ! (ca+dA)*(cb+dA)=Kr
  implicit none
  real(r8), intent(in) :: cd, ca, cb, Kr
  real(r8) :: dA
  real(r8) :: a, b, c, xl, xr
  type(betr_status_type) :: bstatus

  if ((ca+cd)*(cb+cd)<=Kr)then
    !all dissolve
    dA=cd
  else
    a = 1._r8
    b = (ca+cb)
    c = ca*cb-Kr
    xr = cd
    xr = -min(ca,cb)
    call quadrootbnd(a,b,c, xl, xr, dA, bstatus)
  endif
  end function echem_DsKAB
  !-------------------------------------------------------------------------------

  function echem_DsKA2B(cd, ca, cb, Kr)result(dA)
  !D(s) <-> 2A + B
  !(ca+2dA)^2*(cb+dA) = Kr
  implicit none
  real(r8), intent(in) :: cd, ca, cb, Kr
  real(r8) :: q, r, s
  type(betr_status_type) :: bstatus

  real(r8) :: dA
  if ((ca+2._r8*cd)**2._r8*(cb+cd)<=Kr)then
    !all dissolved
    dA=cd
  else
    q = (ca+cb)
    r = (ca*ca+4._r8*ca*cb)*0.25_r8
    s = (ca*ca*cb-kr)*0.25_r8
    dA=cubic_newtonraphson(q,r,s,0._r8)

    dA= max(min(cd, dA),-min(0.5_r8*ca,cb))
  endif
  end function echem_DsKA2B

  !-------------------------------------------------------------------------------
  subroutine eChemSolver(Nr, Ns, rlenmax, ystates1, log10K0, react_type, spm_stoich,tspm_stoich,ystates2)
  use SparseMatMod, only : sparseMat_type

  implicit none
  integer, intent(in) :: Nr
  integer, intent(in) :: Ns
  integer, intent(in) :: rlenmax
  real(r8),intent(in) :: ystates1(1:Ns)             !initial state variables
  real(r8),intent(in) :: log10K0(1:Nr)                !log10 of the equilibrium constants
  integer, intent(in) :: react_type(1:Nr)           !reaction type
  class(sparseMat_type), intent(in) :: spm_stoich   !the stoichiometry matrix, nsxnr
  class(sparseMat_type), intent(in) :: tspm_stoich  !transpose of the stoichiometry matrix, nrxns
  real(r8), intent(out):: ystates2(1:Nr)

  real(r8) :: Fr(1:Nr)
  logical :: update_reaction(1:Nr)
  integer :: j, j1, id
  integer :: imax
  integer :: rlen
  real(r8):: ys(1:rlenmax)  !state variable of the reaction imax
  real(r8):: vs(1:rlenmax)  !stoichiometry coefficient of the reaction imax
  integer :: ias(1:rlenmax)  !location of the state variable of the reaction imax
  type(betr_status_type) :: bstatus
  real(r8), parameter :: err_tol = 1.e-3_r8
  !compute all equilibrium deviations
  do j = 1, nr
    !do j-th reaction
    Fr(j)=-log10K0(j)
    do j1 =1, tspm_stoich%ncol(j)
      id = tspm_stoich%pB(j)+j1 -1
      Fr(j)=Fr(j)+tspm_stoich%val(id)*log10(ystates1(tspm_stoich%icol(id)))
    enddo
    Fr(j)=abs(Fr(j))
  enddo
  ystates2(:)=ystates1(:)
  do
     imax=imaxloc(Fr)
     if(Fr(imax)<err_tol)exit
     !solve the equilibrium for reaction imax
     !assemble the reaction
     rlen=tspm_stoich%ncol(imax)
     do j1 =1, tspm_stoich%ncol(imax)
       id = tspm_stoich%pB(imax)+j1 -1
       ys(j1)=ystates2(tspm_stoich%icol(id))
       vs(j1)=tspm_stoich%val(id)
       ias(j1)=tspm_stoich%icol(id)
     enddo

     !solve the reaction imax
     call eChemReaction(rlen, ys, vs, react_type(imax), log10K0(imax), bstatus)
     if(bstatus%check_status())return
     !update state variable
     do j1 =1, tspm_stoich%ncol(imax)
       id = tspm_stoich%pB(imax)+j1 -1
       ystates2(tspm_stoich%icol(id))=ys(j1)
     enddo

     update_reaction(:)=.false.
     !update relevant Fj
     do j = 1, rlen
       !reaction associated with state variable ias(j)
       if(update_reaction(ias(j)))cycle
       Fr(ias(j))=-log10K0(ias(j))
       do j1 = 1, spm_stoich%ncol(ias(j))
         id = tspm_stoich%pB(ias(j))+j1 -1
         Fr(ias(j))=Fr(ias(j))+tspm_stoich%val(id)*log10(ystates2(tspm_stoich%icol(id)))
       enddo
       Fr(ias(j))=abs(Fr(ias(j)))
       update_reaction(ias(j))=.true.
     enddo
  enddo

  end subroutine eChemSolver
  !-------------------------------------------------------------------------------
  subroutine eChemReaction(rlen, ys, vs, react_type, log10K0, betr_status)

  use MathfuncMod, only : asc_sorti_vec
  implicit none
  integer, intent(in)    :: rlen
  real(r8), intent(inout):: ys(1:rlen)
  real(r8), intent(in)   :: vs(1:rlen)
  integer, intent(in)    :: react_type
  real(r8), intent(in)   :: log10K0
  type(betr_status_type)           , intent(out)   :: betr_status

  real(r8) :: ca, cb, cd, ce, Kr
  real(r8) :: dA
  real(r8) :: ysm(1:rlen)
  real(r8) :: vsm(1:rlen)
  integer :: vsm_index(1:rlen)
  integer :: j

  call betr_status%reset()

  Kr = 10**(log10K0)
  ysm(:) = ys(:)
  vsm(:) = vsm(:)
  call asc_sorti_vec(vsm, vsm_index)
  do j = 1, rlen
    ysm(j)=ys(vsm_index(j))
  enddo
  select case (react_type)
    case (r_echem_AKD)
      SHR_ASSERT_ALL((rlen ==2), errMsg(mod_filename,__LINE__),betr_status);if(betr_status%check_status())return
      ca = ysm(1); cd=ysm(2)
      dA=echem_DKA(ca,cd,Kr)

    case (r_echem_ABKDE)
      SHR_ASSERT_ALL((rlen ==4), errMsg(mod_filename,__LINE__),betr_status);if(betr_status%check_status())return
      ca = ysm(1); cb = ysm(2)
      cd = ysm(3); ce = ysm(4)
      dA=echem_ABKDE(ca,cb,cd,ce, Kr)

    case (r_echem_ABKD)
      SHR_ASSERT_ALL((rlen ==3), errMsg(mod_filename,__LINE__),betr_status);if(betr_status%check_status())return
      ca = ysm(1); cb = ysm(2); cd = ysm(3)
      dA=echem_ABKD(ca,cb,cd, Kr)

    case (r_echem_K_generic)
      dA=echem_K_generic(rlen, ys,vs, Kr)

      !those below are dissolution reactions
    case (r_echem_DsKAB)
      SHR_ASSERT_ALL((rlen ==3), errMsg(mod_filename,__LINE__),betr_status);if(betr_status%check_status())return
      cd = ysm(1); ca=ysm(2); cb=ysm(3)
      dA=echem_DsKAB(cd,ca,cb,Kr)

    case (r_echem_DsKA2B, r_echem_DsKAB2)
      SHR_ASSERT_ALL((rlen ==3), errMsg(mod_filename,__LINE__),betr_status);if(betr_status%check_status())return
      cd = ysm(1); ca=ysm(3); cb=ysm(2)
      dA=echem_DsKA2B(cd, ca, cb, Kr)

    case (r_echem_DsK_generic)
      cd = ysm(1)
      dA=echem_DsK_generic(cd, vsm(1), rlen-1, ysm(2:rlen), vsm(2:rlen), Kr)

  end select
  !update variables
  do j = 1, rlen
    ys(j) = ys(j) + vs(j)*dA
  enddo
  end subroutine eChemReaction
  !-------------------------------------------------------------------------------
  function echem_DsK_generic(cd, vsd, rlen, ysm, vsm, Kr)result(dA)

  ! D(s)<-> vsm{j} * ysm{j}
  ! (ysm+vsm*dA)^vsm = Kr
  implicit none
  real(r8), intent(in) :: cd
  real(r8), intent(in) :: vsd
  integer,  intent(in) :: rlen
  real(r8), intent(in) :: ysm(rlen)
  real(r8), intent(in) :: vsm(rlen)
  real(r8), intent(in) :: Kr

  real(r8) :: dA, dA1, dA_min, dA_max
  real(r8) :: dt
  real(r8) :: fx,ff
  real(r8) :: Kt, logK
  integer  :: j
  logical  :: first
  real(r8) :: err_tol = 1.e-3_r8

  !test existence of solid phase
  Kt=1._r8
  do j = 1, rlen
    Kt = Kt * (ysm(j)+vsm(j)*cd)**vsm(j)
  enddo
  if (Kt <= Kr)then
    dA=cd
  else
    logK=log(Kr)
    dA = 0._r8
    dt = 0.5_r8
    dA_min=0._r8
    first = .true.
    do j = 1, rlen
      if (vsm(j) > 0._r8)then
        if (first)then
          dA_min=-ysm(j)/vsm(j)
          first = .false.
        else
          dA_min = max(dA_min,-ysm(j)/vsm(j))
        endif
      endif
    enddo
    dA_max=cd/vsd
    !find the root using the nudging method
    do
      fx=0._r8
      ff=0._r8
      do j = 1, rlen
        fx=fx+vsm(j)*vsm(j)/(ysm(j)+vsm(j)*dA)
        ff=ff+vsm(j)*log(ysm(j)+vsm(j)*dA)
      enddo
      dA1=dA-dt*(ff-logK)
      dA1=min(max(dA1,dA_min),dA_max)
      if(abs((dA1-dA)/(dA1+dA))<err_tol)then
        dA=dA1
        exit
      endif
      dA=dA1
    enddo
  endif
  end function echem_DsK_generic
  !-------------------------------------------------------------------------------
  function echem_K_generic(rlen, ys,vs, Kr)result(dA)
  ! vsm{j}*ysm{j}
  ! (ysm+vsm*dA)^vsm=Kr
  implicit none
  integer,  intent(in) :: rlen
  real(r8), intent(in) :: ys(rlen)
  real(r8), intent(in) :: vs(rlen)
  real(r8), intent(in) :: Kr

  real(r8) :: dA, dA1, dA_min
  real(r8) :: dt
  real(r8) :: fx,ff
  real(r8) :: logK
  integer  :: j
  logical  :: first
  real(r8) :: err_tol = 1.e-3_r8
  logK=log(Kr)
  dA = 0._r8
  dt = 0.5_r8
  dA_min=0._r8
  first = .true.
  do j = 1, rlen
    if (vs(j) > 0._r8)then
      if (first)then
        dA_min=-ys(j)/vs(j)
        first = .false.
      else
        dA_min = max(dA_min,-ys(j)/vs(j))
      endif
    endif
  enddo
  !find the root using the nudging method
  do
    fx=0._r8
    ff=0._r8
    do j = 1, rlen
      fx=fx+vs(j)*vs(j)/(ys(j)+vs(j)*dA)
      ff=ff+vs(j)*log(ys(j)+vs(j)*dA)
    enddo
    dA1=dA-dt*(ff-logK)
    dA1=max(dA1,dA_min)
    if(abs((dA1-dA)/(dA1+dA))<err_tol)then
      dA=dA1
      exit
    endif
    dA=dA1
  enddo
  end function echem_K_generic


end module eChemMod
