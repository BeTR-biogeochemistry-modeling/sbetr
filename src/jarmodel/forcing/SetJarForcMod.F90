
module SetJarForcMod

  use JarBgcForcType, only : JarBGC_forc_type
  use shr_kind_mod  , only : r8 => shr_kind_r8

implicit none
private
  public :: SetJarForc
contains
  subroutine SetJarForc(jar_forc)

  implicit none
  type(JarBGC_forc_type), intent(inout) :: jar_forc

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

  end subroutine SetJarForc
end module SetJarForcMod
