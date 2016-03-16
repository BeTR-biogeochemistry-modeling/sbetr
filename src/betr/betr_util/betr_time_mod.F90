module BeTR_TimeMod

  use bshr_kind_mod        , only : r8 => shr_kind_r8

  implicit none
  
  type, public:: betr_time_type
    real(r8) :: time_end
    real(r8) :: time
    real(r8) :: restart_dtime
    integer  :: tstep
 end type betr_time_type
 
end module BeTR_TimeMod
