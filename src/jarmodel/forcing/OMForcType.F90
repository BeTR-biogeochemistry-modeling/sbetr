module OMForcType

  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__


  type, public :: om_forc_type
    real(r8) :: cflx_input_litr_met
    real(r8) :: cflx_input_litr_cel
    real(r8) :: cflx_input_litr_lig
    real(r8) :: cflx_input_litr_cwd
    real(r8) :: cflx_input_litr_fwd
    real(r8) :: cflx_input_litr_lwd

    real(r8) :: nflx_input_litr_met
    real(r8) :: nflx_input_litr_cel
    real(r8) :: nflx_input_litr_lig
    real(r8) :: nflx_input_litr_cwd
    real(r8) :: nflx_input_litr_fwd
    real(r8) :: nflx_input_litr_lwd

    real(r8) :: pflx_input_litr_met
    real(r8) :: pflx_input_litr_cel
    real(r8) :: pflx_input_litr_lig
    real(r8) :: pflx_input_litr_cwd
    real(r8) :: pflx_input_litr_fwd
    real(r8) :: pflx_input_litr_lwd


  end type om_forc_type
end module OMForcType
