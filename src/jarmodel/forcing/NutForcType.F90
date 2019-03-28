module NutForcType

  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__


  type, public :: nut_forc_type
    real(r8) :: sflx_minn_input_nh4
    real(r8) :: sflx_minn_input_no3
    real(r8) :: sflx_minp_input_po4
  end type nut_forc_type
end module NutForcType
