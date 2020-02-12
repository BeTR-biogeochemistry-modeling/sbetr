module anaforc

!
! DESCRIPTION
! subroutines output environmental forcing from analytical functions
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod , only : errMsg => shr_log_errMsg

implicit none
  public
  character(len=*), private, parameter :: mod_filename=__FILE__
contains

  subroutine anal_tempf(temp, tsec, Tave, Ay, Ad)

  use shr_const_mod, only : pi=> SHR_CONST_PI
  implicit none
  real(r8), intent(out) :: temp
  real(r8), intent(in) :: tsec
  real(r8), optional, intent(in) :: tave !annual average
  real(r8), optional, intent(in) :: Ay   !seasonal amplitude
  real(r8), optional, intent(in) :: Ad   !diurnal amplitude

  real(r8) :: Tave_loc, Ay_loc, Ad_loc
  real(r8) :: tday
  Tave_loc=290._r8
  Ay_loc = 10._r8
  Ad_loc= 8._r8
  if(present(Tave))Tave_loc=Tave
  if(present(Ay))Ay_loc=Ay
  if(present(Ad))Ad_loc=Ad
  tday=tsec/86400._r8
  temp = Tave_loc-Ay_loc*cos(2._r8*pi/365._r8*tday)+Ad_loc*sin(2_r8*pi*tday)

  end subroutine anal_tempf
end module anaforc
