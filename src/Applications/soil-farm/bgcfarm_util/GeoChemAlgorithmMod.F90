module GeoChemAlgorithmMod

  use bshr_kind_mod         , only : r8 => shr_kind_r8
  use betr_decompMod , only : bounds_type => betr_bounds_type
implicit none

  public :: calc_P_weathering_flux
contains
  subroutine calc_P_weathering_flux(bounds, biophysforc, bgc_con, P_weather_flx)
  !
  !DESCRIPTION
  !calculuate P-weathering flux using Hartmann et al. 2014's model.
  !
  !!USES
  use BiogeoConType, only : BiogeoCon_type
  use BeTR_biogeophysInputType, only : betr_biogeophys_input_type
  use tracer_varcon, only : patomw
  implicit none
  type(bounds_type)                       , intent(in) :: bounds
  type(betr_biogeophys_input_type)        , intent(in) :: biophysforc
  type(BiogeoCon_type)                    , intent(in) :: bgc_con
  real(r8)                                , intent(out):: P_weather_flx(bounds%begc:bounds%endc) ! mol/m2/s

  integer :: c
  real(r8) :: ft
  associate(                            &
   begc => bounds%begc                , &
   endc => bounds%endc                , &
   E_weath=> bgc_con%E_weath          , &
   T_ref_weath=> bgc_con%T_ref_weath  , &
   b_weath  => bgc_con%b_weath        , &
   f_shield => bgc_con%f_shield       , &
   P_weip   => bgc_con%P_weip         , &
   t_soi_10cm=> biophysforc%t_soi_10cm, &
   qflx_runoff_col=> biophysforc%qflx_runoff_col, &
   lithotype => biophysforc%lithotype_col &
  )
  do c = begc, endc
    ft = exp(-E_weath(lithotype(c))*(1._r8/t_soi_10cm(c)-1._r8/T_ref_weath))

    P_weather_flx(c) = b_weath(lithotype(c)) * qflx_runoff_col(c) * ft * &
      f_shield(lithotype(c)) * P_weip(lithotype(c)) / patomw
  enddo
  end associate
  end subroutine calc_P_weathering_flux


end module GeoChemAlgorithmMod
