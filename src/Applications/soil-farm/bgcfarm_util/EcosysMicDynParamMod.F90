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


  public :: get_film_thickness
  public :: get_soil_bacteria_Keff_gas
  public :: get_soil_bacteria_Keff_solute
  public :: get_microbe_ftn
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

  !-------------------------------------------------------------------------------
  function soil_matrix_mic_resistance(rm, vm, filmthk, diffusb)result(ans)

  use bshr_const_mod, only : pi => SHR_CONST_PI
  implicit none
  real(r8), intent(in) :: rm      !microsite radius, m
  real(r8), intent(in) :: vm      !microsite volume, m
  real(r8), intent(in) :: filmthk !water film thickness, m
  real(r8), intent(in) :: diffusb !bulk diffusivity, m^2 s^-1
  real(r8) :: ans

  ans = 4._r8 * pi/(diffusb*(rm+filmthk))
  return
  end function soil_matrix_mic_resistance


  !-------------------------------------------------------------------------------
  function microsite_mic_resistance(rm, vm, filmthk, diffusw)result(ans)

  use bshr_const_mod, only : pi => SHR_CONST_PI
  implicit none
  real(r8), intent(in) :: rm      !microsite radius, m
  real(r8), intent(in) :: vm      !microsite volume, m
  real(r8), intent(in) :: filmthk !water film thickness, m
  real(r8), intent(in) :: diffusw !aqueous diffusivity, m^2 s^-1
  real(r8) :: ans

  ans = 4._r8 * pi * vm * filmthk/(diffusw*rm*(rm+filmthk))
  return
  end function microsite_mic_resistance

  !-------------------------------------------------------------------------------

  subroutine ref_Km_k1w_bacteria(k2p, Np, rp, rc, diffusw, K0, k1w)
  use bshr_const_mod, only : Na=> SHR_CONST_AVOGAD
  use bshr_const_mod, only : pi => SHR_CONST_PI
  implicit none
  real(r8), intent(in) :: k2p     !transporter specific substrate processing rate 1/s
  real(r8), intent(in) :: Np      !number of transporters per bacteria
  real(r8), intent(in) :: rp      !radius of transporter
  real(r8), intent(in) :: rc      !radisu of bacteria
  real(r8), intent(in) :: diffusw !aqueous diffusivity
  real(r8), intent(out):: K0
  real(r8), intent(out):: k1w


  k1w = 4._r8*pi*diffusw*rc*Np*rp/(Np*rp+pi*rc)

  K0 = k2p*(Np*rp+pi*rc)/(4._r8*pi*diffusw*rc*rp*Na)


  end subroutine ref_km_k1w_bacteria
  !-------------------------------------------------------------------------------
  function get_bacteria_gamma(Ncell, k1w, rm,  filmthk, diffusb, diffusw)result(ans)
  use bshr_const_mod, only : pi => SHR_CONST_PI
  implicit none
  real(r8), intent(in) :: Ncell   !number of cells per microsite
  real(r8), intent(in) :: k1w     !substrate approaching rate at cell surface
  real(r8), intent(in) :: rm      !microsite radius, m

  real(r8), intent(in) :: filmthk !water film thickness, m
  real(r8), intent(in) :: diffusw !aqueous diffusivity, m^2 s^-1
  real(r8), intent(in) :: diffusb
  real(r8) :: ans

  real(r8) :: ikappa
  real(r8) :: vm
  real(r8) :: Rs_soil, Rs_mic
  vm = pi*4._r8/3._r8* rm **3._r8
  Rs_soil = soil_matrix_mic_resistance(rm, vm, filmthk, diffusb)
  Rs_mic  = microsite_mic_resistance(rm, vm, filmthk, diffusw)
  ikappa = Rs_soil + Rs_mic
  ans = 1._r8 + k1w * Ncell * ikappa
  return
  end function get_bacteria_gamma
  !-------------------------------------------------------------------------------

  subroutine get_bacteria_K0(ftn, diffusw, K0, k1w)

  implicit none
  real(r8), intent(in) :: ftn
  real(r8), intent(in) :: diffusw
  real(r8), intent(out):: K0
  real(r8), intent(out):: k1w
  real(r8), parameter :: k2p = 1.e3_r8
  real(r8), parameter :: Np = 3000._r8
  real(r8), parameter :: rp = 1.e-9_r8
  real(r8), parameter :: rc = 1.e-6_r8

  real(r8) :: Np_eff

  Np_eff = Np * ftn
  call ref_Km_k1w_bacteria(k2p, Np_eff, rp, rc, diffusw, K0, k1w)
  end subroutine get_bacteria_K0


  !-------------------------------------------------------------------------------
  function get_film_thickness(psiMPa)result(ans)
  !
  !DESCRIPTION
  ! compute water film thickness based on soil matrix pressure
  ! Tang and Riley, 2018, jgr-biogeosci.
  implicit none
  real(r8), intent(in) :: psiMPa
  real(r8) :: ans

  ans = exp(-13.65_r8 - 0.857_r8 * log(-psiMPa))

  ans = max(1.e-8_r8, psiMPa)
  end function get_film_thickness

  !-------------------------------------------------------------------------------
  function get_soil_bacteria_Keff_solute(ft, diffusw, diffusw0, filmthk)result(ans)
  implicit none
  real(r8), intent(in) :: ft
  real(r8), intent(in) :: filmthk
  real(r8), intent(in) :: diffusw
  real(r8), intent(in) :: diffusw0
  real(r8), parameter :: rc = 1.e-6_r8
  real(r8), parameter :: Ncell=10._r8
  real(r8), parameter :: rm = rc * (80._r8*Ncell)**(1._r8/3._r8)

  real(r8), parameter :: agg_tor = 0.5_r8
  real(r8) :: k0, k1w
  real(r8) :: gamma
  real(r8) :: ans
  real(r8) :: diffus

  diffus = diffusw0 * agg_tor

  call get_bacteria_K0(ft, diffus, K0, k1w)

  gamma=get_bacteria_gamma(Ncell, k1w, rm,  filmthk, diffusw, diffusw)

  ans=gamma * K0

  return
  end function get_soil_bacteria_Keff_solute

  !-------------------------------------------------------------------------------
  function get_soil_bacteria_Keff_gas(ftn, diffusb, diffusw, diffusw0, filmthk)result(ans)
  implicit none
  real(r8), intent(in) :: ftn      !the stoichiometry and temperature factor on interception probability
  real(r8), intent(in) :: diffusb
  real(r8), intent(in) :: diffusw
  real(r8), intent(in) :: diffusw0
  real(r8), intent(in) :: filmthk

  real(r8), parameter :: rc = 1.e-6_r8
  real(r8), parameter :: Ncell=10._r8
  real(r8), parameter :: rm = rc * (80._r8*Ncell)**(1._r8/3._r8)
  real(r8), parameter :: agg_tor = 0.5_r8

  real(r8) :: k0, k1w
  real(r8) :: gamma
  real(r8) :: diffus
  real(r8) :: ans

  diffus = diffusw0 * agg_tor

  call get_bacteria_K0(ftn, diffus,  K0, k1w)

  gamma = get_bacteria_gamma(Ncell, k1w, rm,  filmthk, diffusb, diffusw)
  ans=gamma * K0

  return
  end function get_soil_bacteria_Keff_gas
  !-------------------------------------------------------------------------------

  function get_microbe_ftn(tks, offset)result(ans)
  implicit none
  real(r8), intent(in) :: tks
  real(r8), intent(in) :: offset

  real(r8) :: tkso
  real(r8) :: stk
  real(r8) :: rtk
  real(r8) :: actv
  real(r8) :: ans

  tkso = tks + offset
  stk = 710._r8 * tkso
  rtk = 8.3143_r8 * tkso
  actv = f_actv(stk, rtk)

  ans = 1._r8/actv
  return
  end function get_microbe_ftn
end module EcosysMicDynParamMod
