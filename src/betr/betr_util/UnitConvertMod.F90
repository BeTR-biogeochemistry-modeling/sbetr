module UnitConvertMod

  use bshr_kind_mod  , only : r8 => shr_kind_r8
implicit none


contains
  function ppm2molv(patm_pascal, ppmv, temp_kelvin)result(molv)
  !
  !DESCRIPTIONS
  !convert gas concentration from ppmv to mol m-3
  use bshr_const_mod, only : Rgas_kmol => SHR_CONST_RGAS
  implicit none
  real(r8), intent(in) :: Patm_pascal
  real(r8), intent(in) :: ppmv
  real(r8), intent(in) :: temp_kelvin
  real(r8) :: molv

  molv = patm_pascal/(Rgas_kmol*1.e-3_r8*temp_kelvin) * ppmv * 1.e-6_r8

  end function ppm2molv


end module UnitConvertMod
