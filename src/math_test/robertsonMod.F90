module robertsonMod

! subroutines for the robertson model
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  public :: get_robertson_rates
  public :: get_robertson_matrix
contains

  !-------------------------------------------------------------------------------
  subroutine get_robertson_matrix(robmat)

  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  implicit none
  real(r8), intent(out) :: robmat(3,3)

  robmat(1,:)=(/-0.04_r8,0._r8,1.e4_r8/)
  robmat(2,:)=(/0.04_r8,-3.e7_r8,-1.e4_r8/)
  robmat(3,:)=(/0._r8,3.e7_r8,0._r8/)

  end subroutine get_robertson_matrix

  !-------------------------------------------------------------------------------
  subroutine get_robertson_rates(ystates, rrates)

  implicit none
  real(r8), intent(in) :: ystates(3)
  real(r8), intent(out):: rrates(3)

  rrates(1) = ystates(1)
  rrates(2) = ystates(2)**2._r8
  rrates(3) = ystates(2)*ystates(3)

  end subroutine get_robertson_rates


end module robertsonMod
