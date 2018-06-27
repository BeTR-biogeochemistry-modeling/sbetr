module GeoChemAlgorithmMod

  use bshr_kind_mod         , only : r8 => shr_kind_r8
  use betr_decompMod , only : bounds_type => betr_bounds_type
implicit none

  public :: calc_P_weathering_flux
  public :: calc_om_sorption_para
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
  class(BiogeoCon_type)                   , intent(in) :: bgc_con
  real(r8)                                , intent(out):: P_weather_flx(bounds%begc:bounds%endc) ! gP/m2/s

  integer :: c,ltotype
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
    ltotype=lithotype(c)
    if(lithotype(c)<0)ltotype=8
    ft = exp(-E_weath(ltotype)*(1._r8/t_soi_10cm(c)-1._r8/T_ref_weath))

    P_weather_flx(c) = b_weath(ltotype) * qflx_runoff_col(c) * ft * &
      f_shield(ltotype) * P_weip(ltotype) 
    
  enddo
  end associate
  end subroutine calc_P_weathering_flux
!------------------------------------------------------------------------------
  subroutine calc_om_sorption_para(clay, toc, bd, pH, CEC, Qmax, Kaff)
  !
  ! DESCRIPTION
  !
  !Below is a simple model for computing the Qmax and Kaff parameters for
  !the Langmuir isotherm
  !
  !USES
  use tracer_varcon, only : catomw
  implicit none
  real(r8), intent(in) :: clay ! % clay
  real(r8), intent(in) :: toc  ! gC/m3
  real(r8), intent(in) :: bd   ! bulk density, kg/m3
  real(r8), intent(in) :: pH   !
  real(r8), intent(in) :: CEC  !cation exchange capacity
  real(r8), intent(out):: Qmax !maximum sorption capacity, mol C/m3
  real(r8), intent(out):: Kaff !sorption affinity parameter

  Qmax=exp(0.483*log(clay)+2.328)*bd*1.e-3_r8/catomw
  Kaff= 1._r8/exp(-0.186_r8*pH - 0.216_r8)/catomw

  end subroutine calc_om_sorption_para
end module GeoChemAlgorithmMod
