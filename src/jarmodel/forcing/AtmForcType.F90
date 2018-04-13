module AtmForcType

  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__


  type, public :: atm_forc_type

    real(r8) :: air_temp
    real(r8) :: patm_pascal
    real(r8) :: n2_ppmv
    real(r8) :: o2_ppmv
    real(r8) :: ar_ppmv
    real(r8) :: co2_ppmv
    real(r8) :: ch4_ppmv
    real(r8) :: n2o_ppmv
    real(r8) :: no_ppmv
    real(r8) :: nh3_ppmv
  end type atm_forc_type
end module AtmForcType
