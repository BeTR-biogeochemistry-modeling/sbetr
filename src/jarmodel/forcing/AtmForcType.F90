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
    real(r8) :: ra           !aerodynamic resistance
  contains
    procedure, public :: init
  end type atm_forc_type

contains
  subroutine init(this)

  implicit none
  class(atm_forc_type), intent(inout) :: this

  this%n2_ppmv = 7.8e5_r8
  this%o2_ppmv = 2.1e5_r8
  this%ar_ppmv = 0.9e5_r8
  this%co2_ppmv= 400._r8
  this%ch4_ppmv= 1.8_r8
  this%n2o_ppmv= 250.e-3_r8
  this%nh3_ppmv= 0._r8
  this%no_ppmv = 0._r8
  this%ra      = 0._r8
  end subroutine init
end module AtmForcType
