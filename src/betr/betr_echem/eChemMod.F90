module eChemMod
  use FindRootMod
  use BetrStatusType  , only : betr_status_type

implicit none

contains


  function echem_DKA(ca,cd,Kr)result(dA)
  implicit none
  real(r8), intent(in) :: ca, cd, Kr
  real(r8) :: dA

  dA=(cd-ca*Kr)/(1._r8+Kr)

  end function echem_DKA


  function echem_ABKDE(ca,cb,cd,ce, Kr)result(dA)

  implicit none
  real(r8), intent(in) :: ca, cb, cd, ce, Kr
  real(r8) :: dA
  real(r8) :: a, b, c, xl, xr
  type(betr_status_type) :: bstatus

  a = 1._r8+Kr
  b = Kr*(ca+cb)+(cd+ce)
  c = Kr*ca*cb-cd*ce
  xl = -min(ca, cb)
  xr = min(cd, ce)
  call quadrootbnd(a,b,c, xl, xr, dA, bstatus)

  end function echem_ABKDE

  function echem_ABKD(ca,cb,cd, Kr)result(dA)
  implicit none
  real(r8), intent(in) :: ca, cb, cd, Kr
  real(r8) :: dA
  real(r8) :: a, b, c, xl, xr
  type(betr_status_type) :: bstatus

  a = Kr
  b = Kr*(ca+cb)+1._r8
  c = Kr*ca*cb-cd
  xl = -min(ca, cb)
  xr = cd
  call quadrootbnd(a,b,c, xl, xr, dA, bstatus)

  end function echem_ABKD

  function echem_DsKAB(cd,ca,cb,Kr)result(dA)

  implicit none
  real(r8), intent(in) :: cd, ca, cb, Kr
  real(r8) :: dA
  real(r8) :: a, b, c, xl, xr
  type(betr_status_type) :: bstatus

  a = 1._r8
  b = -(ca+cb)
  c = ca*cb-Kr
  xl = 0._r8
  xr = min(ca,cb)
  call quadrootbnd(a,b,c, xl, xr, dA, bstatus)

  end function echem_DsKAB

  function echem_DsKA2B(cd, ca, cb, Kr)result(dA)
  implicit none
  real(r8), intent(in) :: cd, ca, cb, Kr
  real(r8) :: q, r, s
  type(betr_status_type) :: bstatus

  q = -(2._r8*ca+cb)
  r = (ca*ca+2*ca*cb)
  s = -ca*ca*cb+kr
  dA=cubic_newtonraphson(q,r,s,0._r8)

  end function echem_DsKA2B


end module eChemMod
