
module SetJarForcMod

  use JarBgcForcType, only : JarBGC_forc_type
  use shr_kind_mod  , only : r8 => shr_kind_r8

implicit none
private
  public :: SetJarForc
  interface SetJarForc
    module procedure SetJarForc_const,SetJarForc_transient
  end interface SetJarForc

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



  !--------------------------------------------------------------------
  subroutine SetJarForc_transient(jar_forc, om_forc, atm_forc, soil_forc)

  use SoilForcType, only : soil_forc_type
  use AtmForcType, only : atm_forc_type
  use OMForcType, only : om_forc_type
  implicit none
  type(JarBGC_forc_type), intent(inout) :: jar_forc
  type(om_forc_type), intent(in) :: om_forc
  type(atm_forc_type), intent(in) :: atm_forc
  type(soil_forc_type), intent(in) :: soil_forc

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

  jar_forc%air_temp   = atm_forc%air_temp
  jar_forc%temp       = soil_forc%temp            !temperature
  jar_forc%h2osoi_vol = soil_forc%h2osoi_vol
  jar_forc%h2osoi_liq = soil_forc%h2osoi_liq
  jar_forc%air_vol    = soil_forc%air_vol
  jar_forc%finundated = soil_forc%finundated
  jar_forc%soilpsi= soil_forc%soilpsi             ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]

  jar_forc%conc_atm_n2 =ppm2molv(atm_forc%patm_pascal, atm_forc%n2_ppmv, jar_forc%air_temp)   !n2 concentration in atmosphere, mol n2/m3
  jar_forc%conc_atm_n2o=ppm2molv(atm_forc%patm_pascal, atm_forc%n2o_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_o2=ppm2molv(atm_forc%patm_pascal, atm_forc%o2_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_ar=ppm2molv(atm_forc%patm_pascal, atm_forc%ar_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_co2=ppm2molv(atm_forc%patm_pascal, atm_forc%co2_ppmv, jar_forc%air_temp)
  jar_forc%conc_atm_ch4=ppm2molv(atm_forc%patm_pascal, atm_forc%ch4_ppmv, jar_forc%air_temp)

  end subroutine SetJarForc_transient


  !--------------------------------------------------------------------
  subroutine SetJarForc_const(jar_forc, soil_forc)
  use SoilForcType, only : soil_forc_type
  implicit none
  type(JarBGC_forc_type), intent(inout) :: jar_forc
  type(soil_forc_type), intent(in) :: soil_forc

  real(r8) :: patm_pascal = 1.01325e5 !Pa
  !all mass fluxes are in the unit of g C/m2/s

  jar_forc%depz   = soil_forc%depz             !depth of the soil
  jar_forc%dzsoi  = soil_forc%dzsoi             !soil thickness
  jar_forc%sucsat = soil_forc%sucsat             ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
  jar_forc%bsw    = soil_forc%bsw
  jar_forc%bd     = soil_forc%bd             !bulk density
  jar_forc%pct_sand= soil_forc%pct_sand
  jar_forc%pct_clay= soil_forc%pct_clay
  jar_forc%watsat= soil_forc%watsat
  jar_forc%watfc= soil_forc%watfc
  jar_forc%cellorg= soil_forc%cellorg
  jar_forc%pH= soil_forc%pH

  jar_forc%plant_froot_nn(:)=0._r8
  jar_forc%plant_froot_np(:)=0._r8
  jar_forc%plant_vtype(:)=0
  jar_forc%plant_ntypes=0
  jar_forc%soilorder=1
  jar_forc%msurf_nh4=0._r8
  jar_forc%msurf_minp=0._r8

  end subroutine SetJarForc_const

end module SetJarForcMod
