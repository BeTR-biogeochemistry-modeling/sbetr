module EcosysMicDynParamMod

!
! DESCRIPTION
! functions for microbial dynamics parameterization extracted
! from the ecosys model
  use bshr_kind_mod             , only : r8 => shr_kind_r8
implicit none
  private

  character(len=*), private, parameter :: mod_filename = &
     __FILE__
  public :: calc_tm_factor
  public :: Kb_smmodifier
  interface Kb_smmodifier
    module procedure Kb_smmodifier_single, Kb_smmodifier_dual
  end interface Kb_smmodifier
contains
  !-------------------------------------------------------------------------------
  function f_actv(stk, rtk)result(actv)
  use bshr_kind_mod             , only : r8 => shr_kind_r8
  implicit none
  real(r8), intent(in) :: stk  !product between entropy and temperature
  real(r8), intent(in) :: rtk  !product between universal gas constant and temperature
  !local variables
  real(r8) :: actv
  actv = 1._r8 + exp((195000._r8-stk)/rtk)+exp((stk-222500._r8)/rtk)
  end function f_actv
  !-------------------------------------------------------------------------------
  function f_tfnx(rtk, actv)result(tfnx)
  !
  ! temperature function for substrate uptake
  implicit none
  real(r8), intent(in) :: rtk
  real(r8), intent(in) :: actv
  !local variables
  real(r8) :: tfnx

  tfnx = exp(25.227_r8 - 62500._r8/rtk)/actv
  end function f_tfnx
  !-------------------------------------------------------------------------------
  function f_tfny(tcs, offset)result(tfny)
  implicit none
  real(r8), intent(in) :: tcs !temperaure in celcius
  real(r8), intent(in) :: offset
  !local variables
  real(r8) :: tfny
  real(r8) :: tcsr
  tcsr = min(45._r8, tcs) + offset
  tfny = exp(0.0742_r8*(tcsr-25._r8))
  end function f_tfny

  !-------------------------------------------------------------------------------
  function f_wfng(psism)result(wfng)
  implicit none
  real(r8), intent(in) :: psism !soil micropore matric potential, MPa, negative
  real(r8) :: wfng

  wfng = exp(0.05_r8 * psism)
  end function f_wfng
  !-------------------------------------------------------------------------------
  subroutine calc_tm_factor(tks, psism, tfng, tfnr, offset)
  implicit none
  real(r8), intent(in) :: tks   !soil temperature in kelvin
  real(r8), intent(in) :: psism
  real(r8), intent(out) :: tfng  !temperature x moisture factor for uptake
  real(r8), intent(out) :: tfnr  !temperature x moisture factor for maintenance
  real(r8), optional, intent(in) :: offset  !this is a place holder for accounting for temperature acclimation

  real(r8) :: actv
  real(r8) :: offset_local
  real(r8) :: tkso, stk, rtk, tcs
  real(r8) :: tfnx, tfny, wfng

  if(present(offset))then
    offset_local = offset
  else
    offset_local = 0._r8
  endif
  tcs = tks - 273.15_r8
  tkso = tks + offset_local
  stk = 710._r8 * tkso
  rtk = 8.3143_r8 * tkso
  actv = f_actv(stk, rtk)
  tfnx = f_tfnx(rtk, actv)
  tfny = f_tfny(tcs, offset_local)
  wfng = f_wfng(psism)
  tfng = tfnx * wfng
  tfnr = tfny * wfng
  end subroutine calc_tm_factor

  !-------------------------------------------------------------------------------
  function Kb_smmodifier_single(tauw, h2osoi_liqvol, filmt)result(ans)
  !
  ! DESCRIPTION
  ! compute the moisture modifier to half saturation constant to microbes
  ! reference: Tang and Riley, 2017, GMD
  implicit none
  real(r8), intent(in) :: tauw  !tortuosity
  real(r8), intent(in) :: filmt       !water film thickness, m
  real(r8), intent(in) :: h2osoi_liqvol

  real(r8), parameter :: c1=0.48_r8
  real(r8), parameter :: rc=1.e-6_r8
  real(r8), parameter :: rm=3.e-6_r8
  real(r8) :: ans

  ans = h2osoi_liqvol + c1 * rc/(rm+filmt) * (filmt/rm*h2osoi_liqvol+1._r8/max(tauw,1.e-10))

  end function Kb_smmodifier_single

  !-------------------------------------------------------------------------------
  function Kb_smmodifier_dual(tauw, h2osoi_liqvol, filmt, air_vol, taug, bunsencef)result(ans)
  !
  ! DESCRIPTION
  ! compute the moisture modifier to half saturation constant to microbes
  ! dual phase tracer
  implicit none
  real(r8), intent(in) :: taug, tauw  !tortuosity
  real(r8), intent(in) :: filmt       !water film thickness
  real(r8), intent(in) :: h2osoi_liqvol
  real(r8), intent(in) :: air_vol
  real(r8), intent(in) :: bunsencef
  real(r8), parameter :: c1=0.48_r8
  real(r8), parameter :: rc=1.e-6_r8
  real(r8), parameter :: rm=3.e-6_r8
  real(r8), parameter :: dif_ratio_g2w = 1.e4_r8  ! this is a very crude approximation
  real(r8) :: ans


  ans = (h2osoi_liqvol+air_vol/bunsencef)*(1._r8+c1*rc/(rm+filmt) * &
   (filmt/rm + 1._r8/(tauw*h2osoi_liqvol+taug*air_vol/bunsencef*dif_ratio_g2w)))

  end function Kb_smmodifier_dual

end module EcosysMicDynParamMod
