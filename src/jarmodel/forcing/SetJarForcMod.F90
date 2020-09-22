
module SetJarForcMod
#include "bshr_assert.h"
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use JarBgcForcType, only : JarBGC_forc_type
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use SoilForcType, only : soil_forc_type
  use AtmForcType , only : atm_forc_type
implicit none
private
  character(len=*), parameter :: mod_filename = &
       __FILE__
  public :: SetJarForc
  interface SetJarForc
    module procedure SetJarForc_const,SetJarForc_transient
  end interface SetJarForc
  public :: setJarStates
contains


  !--------------------------------------------------------------------
  subroutine SetJarForc_transient(jar_forc, om_forc, nut_forc, atm_forc, &
     soil_forc, dtime, jarpars)

  use OMForcType  , only : om_forc_type
  use NutForcType , only : nut_forc_type
  use UnitConvertMod, only : ppm2molv
  !begin_appadd
  use ecacnpParaType, only : ecacnp_para_type
  use kecaParaType, only : keca_para_type
  use v1ecaParaType, only : v1eca_para_type
  !end_appadd
  use BiogeoConType    , only : BiogeoCon_type
  implicit none
  type(JarBGC_forc_type), intent(inout) :: jar_forc
  type(nut_forc_type)   , intent(in) :: nut_forc
  type(om_forc_type)    , intent(in) :: om_forc
  type(atm_forc_type)   , intent(in) :: atm_forc
  type(soil_forc_type)  , intent(in) :: soil_forc
  real(r8)              , intent(in) :: dtime
  class(BiogeoCon_type),  intent(inout) :: jarpars
  !all mass fluxes are in the unit of g C/m2/s

  jar_forc%cflx_input_litr_met = om_forc%cflx_input_litr_met
  jar_forc%cflx_input_litr_cel = om_forc%cflx_input_litr_cel
  jar_forc%cflx_input_litr_lig = om_forc%cflx_input_litr_lig
  jar_forc%cflx_input_litr_cwd = om_forc%cflx_input_litr_cwd
  jar_forc%cflx_input_litr_fwd = om_forc%cflx_input_litr_fwd
  jar_forc%cflx_input_litr_lwd = om_forc%cflx_input_litr_lwd

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

  jar_forc%sflx_minn_input_nh4 = nut_forc%sflx_minn_input_nh4
  jar_forc%sflx_minn_input_no3 = nut_forc%sflx_minn_input_no3
  jar_forc%sflx_minp_input_po4 = nut_forc%sflx_minp_input_po4


  jar_forc%air_temp   = atm_forc%air_temp
  jar_forc%temp       = soil_forc%temp            !temperature
  jar_forc%h2osoi_vol = soil_forc%h2osoi_vol
  jar_forc%h2osoi_liq = soil_forc%h2osoi_liq
  jar_forc%air_vol    = soil_forc%air_vol
  jar_forc%finundated = soil_forc%finundated
  jar_forc%soilpsi    = soil_forc%soilpsi         ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]

  jar_forc%conc_atm_n2 =ppm2molv(atm_forc%patm_pascal, atm_forc%n2_ppmv, jar_forc%air_temp)   !n2 concentration in atmosphere, mol n2/m3
  jar_forc%conc_atm_n2o=ppm2molv(atm_forc%patm_pascal, atm_forc%n2o_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_o2=ppm2molv(atm_forc%patm_pascal, atm_forc%o2_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_ar=ppm2molv(atm_forc%patm_pascal, atm_forc%ar_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_co2=ppm2molv(atm_forc%patm_pascal, atm_forc%co2_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_ch4=ppm2molv(atm_forc%patm_pascal, atm_forc%ch4_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_nh3=ppm2molv(atm_forc%patm_pascal, atm_forc%nh3_ppmv, jar_forc%air_temp)

  call set_phase_convert_coeff(atm_forc, soil_forc, jar_forc)

  select type(jarpars)
  class is (keca_para_type)
    jar_forc%tmic_opt = jar_forc%tmic_opt + &
      dtime/jarpars%tau30*(jar_forc%temp-jar_forc%tmic_opt)
  class default

  end select

  end subroutine SetJarForc_transient

  !--------------------------------------------------------------------
  subroutine SetJarForc_const(jar_forc, soil_forc)
  use tracer_varcon, only : catomw
  implicit none
  type(JarBGC_forc_type), intent(inout) :: jar_forc
  type(soil_forc_type), intent(in) :: soil_forc

  real(r8) :: patm_pascal = 1.01325e5 !Pa
  !all mass fluxes are in the unit of g C/m2/s

  jar_forc%depz    = soil_forc%depz             !depth of the soil
  jar_forc%dzsoi   = soil_forc%dzsoi             !soil thickness
  jar_forc%sucsat  = soil_forc%sucsat             ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
  jar_forc%bsw     = soil_forc%bsw
  jar_forc%bd      = soil_forc%bd             !bulk density
  jar_forc%pct_sand= soil_forc%pct_sand
  jar_forc%pct_clay= soil_forc%pct_clay
  jar_forc%watsat  = soil_forc%watsat
  jar_forc%watfc   = soil_forc%watfc
  jar_forc%cellorg = soil_forc%cellorg
  jar_forc%pH      = soil_forc%pH

  jar_forc%plant_froot_nn(:)=0._r8
  jar_forc%plant_froot_np(:)=0._r8
  jar_forc%plant_vtype(:)=0
  jar_forc%plant_ntypes=0
  jar_forc%soilorder=1
  jar_forc%msurf_nh4=0._r8
  jar_forc%msurf_minp=0._r8
  jar_forc%msurf_OM= 500._r8/catomw
  jar_forc%KM_OM_ref=10._r8/catomw
  end subroutine SetJarForc_const

  !--------------------------------------------------------------------
  subroutine setJarStates(jar_forc, ystates)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  implicit none
  type(JarBGC_forc_type), intent(inout) :: jar_forc
  real(r8), intent(in) :: ystates(:)

  SHR_ASSERT_ALL_EXT((size(jar_forc%ystates(:)) == size(ystates)),  errMsg(mod_filename,__LINE__))
  jar_forc%ystates(:) = ystates(:)

  end subroutine setJarStates
  !--------------------------------------------------------------------
  subroutine set_phase_convert_coeff(atm_forc, soil_forc, jar_forc)
  implicit none
  type(atm_forc_type), intent(in) :: atm_forc
  type(soil_forc_type), intent(in) :: soil_forc
  type(JarBGC_forc_type), intent(inout) :: jar_forc

  real(r8) :: henrycef
  real(r8) :: bunsencef_n2o, bunsencef_n2
  real(r8) :: Dw, Dg
  associate(             &
    temp   => soil_forc%temp   , &
    air_vol=> soil_forc%air_vol, &
    h2osoi_liqvol=>soil_forc%h2osoi_liqvol, &
    Diff_Darcy  => soil_forc%Diff_Darcy, &
    tauaqu => soil_forc%tauaqu, &
    taugas => soil_forc%taugas, &
    dzsoi => soil_forc%dzsoi, &
    ra => atm_forc%ra, &
    o2_w2b => jar_forc%o2_w2b, &
    o2_g2b => jar_forc%o2_g2b, &
    n2_g2b => jar_forc%n2_g2b, &
    n2o_g2b=> jar_forc%n2o_g2b, &
    bunsencef_o2=> jar_forc%bunsen_o2, &
    aren_cond_n2o=>jar_forc%aren_cond_n2o, &
    aren_cond_o2 => jar_forc%aren_cond_o2, &
    aren_cond_n2 => jar_forc%aren_cond_n2  &
  )

  jar_forc%h2osoi_liqvol=h2osoi_liqvol
  jar_forc%air_vol = air_vol
  !compute henry's constant
  henrycef=1.3e-3_r8*exp(-1500._r8*(1._r8/temp-1._r8/298.15_r8))

  !compute bunsen coefficient
  bunsencef_o2= henrycef*temp/12.2_r8

  !compute phase conversion parameters
  o2_w2b=air_vol/bunsencef_o2+h2osoi_liqvol
  o2_g2b=air_vol+h2osoi_liqvol*bunsencef_o2

  !for n2o
  !compute henry's constant
  henrycef=2.5e-2_r8*exp(-2600._r8*(1._r8/temp-1._r8/298.15_r8))

  !compute bunsen coefficient
  bunsencef_n2o= henrycef*temp/12.2_r8

  !compute phase conversion parameters
  n2o_g2b=air_vol+h2osoi_liqvol*bunsencef_n2o

  !for n2
  !compute henry's constant
  henrycef=6.1e-4_r8*exp(-1300._r8*(1._r8/temp-1._r8/298.15_r8))

  !compute bunsen coefficient
  bunsencef_n2= henrycef*temp/12.2_r8

  !compute phase conversion parameters
  n2_g2b=air_vol+h2osoi_liqvol*bunsencef_n2

  !set conductivity coefficient
  !according to Tang and Riley (2013), HESS
  !the conductance is
  !1/(ra+0.5*dzsoi/(Da*air_vol+Dw*h2osoi_liqvol*bunsen))
  !here ra is aerodynamic resistance
  !Da is gas diffusivity in soil, including tortuosity
  !Dw is aqueous diffusivity plus Darcy diffusivity

  !compute tortuosity

  !n2
  Dg=1.93e-5_r8*(temp/273.0_r8)**1.82_r8 * taugas
  Dw=2.57e-9_r8*(temp/273.0_r8) * tauaqu
  aren_cond_n2=1._r8/(ra+dzsoi*0.5_r8/(Dg*air_vol+(Dw+Diff_Darcy)*h2osoi_liqvol*bunsencef_n2))
  aren_cond_n2=aren_cond_n2/dzsoi

  !o2
  Dg=1.8e-5_r8*(temp/273.0_r8)**1.82_r8 * taugas
  Dw=2.4e-9_r8*temp/298.0_r8 * tauaqu
  aren_cond_o2=1._r8/(ra+dzsoi*0.5_r8/(Dg*air_vol+(Dw+Diff_Darcy)*h2osoi_liqvol*bunsencef_o2))
  aren_cond_o2=aren_cond_o2/dzsoi

  !n2o
  Dg=0.159e-4_r8*(temp/293.15)**1.5_r8 * taugas
  Dw=2.6e-9_r8*temp/298.15_r8 * tauaqu
  aren_cond_o2=1._r8/(ra+dzsoi*0.5_r8/(Dg*air_vol+(Dw+Diff_Darcy)*h2osoi_liqvol*bunsencef_n2o))
  aren_cond_o2=aren_cond_o2/dzsoi

  !compute diffusivity
  jar_forc%diffusw0_nh4 =1.64e-9_r8*temp/298.15_r8
  jar_forc%diffusw_nh4   = jar_forc%diffusw0_nh4 * tauaqu

  jar_forc%diffusw0_no3  = 2.6e-9_r8*temp/298.15_r8
  jar_forc%diffusw_no3   = jar_forc%diffusw0_no3 * tauaqu

  jar_forc%diffusw0_minp = 0.7e-9_r8*temp/293.15_r8
  jar_forc%diffusw_minp  = jar_forc%diffusw0_minp * tauaqu

  end associate
  end subroutine set_phase_convert_coeff

end module SetJarForcMod
