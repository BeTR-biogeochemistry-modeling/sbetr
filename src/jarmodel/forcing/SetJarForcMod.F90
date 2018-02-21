
module SetJarForcMod

  use JarBgcForcType, only : JarBGC_forc_type
  use shr_kind_mod  , only : r8 => shr_kind_r8

implicit none
private
  public :: SetJarForc
  public :: set_soil_hydro_property
contains

  function ppm2molv(patm_pascal, ppmv, temp_kelvin)result(molv)
  use shr_const_mod, only : Rgas_kmol => SHR_CONST_RGAS
  implicit none
  real(r8), intent(in) :: Patm_pascal
  real(r8), intent(in) :: ppmv
  real(r8), intent(in) :: temp_kelvin
  real(r8) :: molv

  molv = patm_pascal/(Rgas_kmol*1.e-3*temp_kelvin)
  end function ppm2molv


  !-----------------------------------------------------------------------
  subroutine set_soil_hydro_property(sand, clay, om_frac, zsoi, bd, watsat, bsw, sucsat, hksat_min)

  use FuncPedotransferMod, only : pedotransf, init_pedof, get_ipedof
  implicit none
  real(r8), intent(in) :: sand
  real(r8), intent(in) :: clay
  real(r8), intent(in) :: zsoi
  real(r8), intent(in) :: om_frac
  real(r8), intent(out):: bd
  real(r8), intent(out):: watsat
  real(r8), intent(out):: bsw
  real(r8), intent(out):: sucsat
  real(r8), intent(out):: hksat_min

  integer :: ipedof
  real(r8):: om_watsat
  real(r8):: xksat
  real(r8):: om_b, om_sucsat, om_hksat
  real(r8), parameter :: zsapric   = 0.5_r8       ! depth (m) that organic matter takes on characteristics of sapric peat

  call init_pedof
  ipedof = get_ipedof(0)

  call pedotransf(ipedof, sand, clay, watsat, bsw, sucsat, xksat)

  om_watsat  = max(0.93_r8 - 0.1_r8   *(zsoi/zsapric), 0.83_r8)
  om_b       = min(2.7_r8  + 9.3_r8   *(zsoi/zsapric), 12.0_r8)
  om_sucsat  = min(10.3_r8 - 0.2_r8   *(zsoi/zsapric), 10.1_r8)
  om_hksat   = max(0.28_r8 - 0.2799_r8*(zsoi/zsapric), 0.0001_r8)

  bd  = (1._r8 - watsat)*2.7e3_r8
  watsat    = (1._r8 - om_frac) * watsat + om_watsat*om_frac

  bsw       = (1._r8-om_frac) * (2.91_r8 + 0.159_r8*clay) + om_frac*om_b
  sucsat    = (1._r8-om_frac) * sucsat + om_sucsat*om_frac
  hksat_min = xksat

  end subroutine set_soil_hydro_property

  !-----------------------------------------------------------------------
  subroutine soil_hk(hksat, imped, s, bsw, hk, dhkds)
    !
    ! !DESCRIPTION:
    ! Compute hydraulic conductivity
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hksat    !saturated hydraulic conductivity [mm/s]
    real(r8), intent(in) :: imped    !ice impedance
    real(r8), intent(in) :: s        !reletive saturation, [0, 1]
    real(r8), intent(in) :: bsw      !shape parameter
    real(r8), intent(out):: hk       !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out):: dhkds    !d[hk]/ds   [mm/s]
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'soil_hk'
    !-----------------------------------------------------------------------

    !compute hydraulic conductivity
    hk=imped*hksat*s**(2._r8*bsw+3._r8)

    !compute the derivative
    if(present(dhkds))then
       dhkds=(2._r8*bsw+3._r8)*hk/s
    endif

  end subroutine soil_hk

  !-----------------------------------------------------------------------
  subroutine soil_suction(smpsat, s, bsw, smp, dsmpds)
    !
    ! !DESCRIPTION:
    ! Compute soil suction potential
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: smpsat   !minimum soil suction, positive [mm]
    real(r8), intent(in)            :: s        !reletive saturation, [0, 1]
    real(r8), intent(in)            :: bsw      !shape parameter
    real(r8), intent(out)           :: smp      !soil suction, negative, [mm]
    real(r8), optional, intent(out) :: dsmpds   !d[smp]/ds, [mm]
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'soil_suction'
    !-----------------------------------------------------------------------

    !compute soil suction potential, negative
    smp = -smpsat*s**(-bsw)

    !compute derivative
    if(present(dsmpds))then
       dsmpds=-bsw*smp/s
    endif

  end subroutine soil_suction

  !--------------------------------------------------------------------
  subroutine SetJarForc(jar_forc)

  implicit none
  type(JarBGC_forc_type), intent(inout) :: jar_forc

  real(r8) :: patm_pascal = 1.01325e5 !Pa
  !all mass fluxes are in the unit of g C/m2/s

  jar_forc%cflx_input_litr_met = 1.e-7_r8
  jar_forc%cflx_input_litr_cel = 1.e-7_r8
  jar_forc%cflx_input_litr_lig = 1.e-7_r8
  jar_forc%cflx_input_litr_cwd = 0._r8
  jar_forc%cflx_input_litr_fwd = 0._r8
  jar_forc%cflx_input_litr_lwd = 0._r8

  jar_forc%nflx_input_litr_met = jar_forc%cflx_input_litr_met/90._r8
  jar_forc%nflx_input_litr_cel = jar_forc%cflx_input_litr_cel/90._r8
  jar_forc%nflx_input_litr_lig = jar_forc%cflx_input_litr_lig/90._r8
  jar_forc%nflx_input_litr_cwd = jar_forc%cflx_input_litr_cwd/90._r8
  jar_forc%nflx_input_litr_fwd = jar_forc%cflx_input_litr_fwd/90._r8
  jar_forc%nflx_input_litr_lwd = jar_forc%cflx_input_litr_lwd/90._r8

  jar_forc%pflx_input_litr_met = jar_forc%cflx_input_litr_met/1600._r8
  jar_forc%pflx_input_litr_cel = jar_forc%cflx_input_litr_cel/2000._r8
  jar_forc%pflx_input_litr_lig = jar_forc%cflx_input_litr_lig/2500._r8
  jar_forc%pflx_input_litr_cwd = jar_forc%cflx_input_litr_cwd/4500._r8
  jar_forc%pflx_input_litr_fwd = jar_forc%cflx_input_litr_fwd/4500._r8
  jar_forc%pflx_input_litr_lwd = jar_forc%cflx_input_litr_lwd/4500._r8

  jar_forc%air_temp = 298.15_r8
  jar_forc%temp   =290._r8            !temperature
  jar_forc%depz   =0.05_r8             !depth of the soil
  jar_forc%dzsoi  =0.1_r8             !soil thickness
  jar_forc%sucsat =0._r8             ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
  jar_forc%soilpsi=0._r8             ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
  jar_forc%bsw    =0._r8
  jar_forc%bd     =0._r8             !bulk density
  jar_forc%pct_sand=0._r8
  jar_forc%pct_clay=0._r8
  jar_forc%h2osoi_vol=0._r8
  jar_forc%h2osoi_liq=0._r8
  jar_forc%air_vol=0._r8
  jar_forc%finundated=0._r8
  jar_forc%watsat=0._r8
  jar_forc%watfc=0._r8
  jar_forc%cellorg=0._r8
  jar_forc%pH=0._r8

  jar_forc%ch4_g2b=0._r8
  jar_forc%co2_g2b=0._r8
  jar_forc%o2_g2b =0._r8
  jar_forc%o2_w2b =0._r8       !conversion parameter for o2 from aqueous to bulk conc
  jar_forc%n2_g2b =0._r8
  jar_forc%ar_g2b =0._r8
  jar_forc%n2o_g2b=0._r8

  jar_forc%conc_atm_n2 =ppm2molv(patm_pascal, 78.e4_r8, jar_forc%air_temp)   !n2 concentration in atmosphere, mol n2/m3
  jar_forc%conc_atm_n2o=ppm2molv(patm_pascal, 0._r8, jar_forc%air_temp)
  jar_forc%conc_atm_o2=ppm2molv(patm_pascal, 21.e4_r8, jar_forc%air_temp)
  jar_forc%conc_atm_ar=ppm2molv(patm_pascal, 0.9e4_r8, jar_forc%air_temp)
  jar_forc%conc_atm_co2=ppm2molv(patm_pascal, 400._r8, jar_forc%air_temp)
  jar_forc%conc_atm_ch4=ppm2molv(patm_pascal, 1.7_r8, jar_forc%air_temp)

  jar_forc%plant_froot_nn(:)=0._r8
  jar_forc%plant_froot_np(:)=0._r8
  jar_forc%plant_vtype(:)=0
  jar_forc%plant_ntypes=0
  jar_forc%soilorder=1
  jar_forc%msurf_nh4=0._r8
  jar_forc%msurf_minp=0._r8

  end subroutine SetJarForc
end module SetJarForcMod
