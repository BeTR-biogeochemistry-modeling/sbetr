  module DebGrowMod
  !  History: created by Jinyun Tang, 2013

  ! !USES:
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use BgcSummsDebType         , only : debs
  use MathfuncMod             , only : safe_div

  implicit none
  
contains
  !-------------------------------------------------------------------------------
subroutine deb_grow(gbt,deb,residual)
    !use func_data_type_mod      , only : func_data_type
    use BgcSummsDebType         , only : debs
    use MathfuncMod             , only : safe_div
    implicit none

    ! !ARGUMENTS:
    real(r8), intent(in) :: gbt
    !type(func_data_type), intent(inout) :: deb
    type(debs), intent(in) :: deb
    real(r8), intent(out) :: residual

    real(r8) :: dc, pot_growth, pot_enz
    ! real(r8), pointer :: gmax_mic
    ! real(r8), pointer  :: yld_mic
    ! real(r8), pointer  :: yld_enz
    ! real(r8), pointer  :: mr_mic
    ! real(r8), pointer  :: pmax_enz
    ! real(r8), pointer  :: je
    ! real(r8), pointer  :: ec
    ! real(r8), pointer  :: gB
    ! real(r8), pointer  :: pE

      associate(                                & !
                gmax_mic =>     deb%gmax_mic  , &
                yld_mic  =>     deb%yld_mic   , &
                yld_enz  =>     deb%yld_enz   , &
                mr_mic   =>     deb%mr_mic    , &
                pmax_enz =>     deb%pmax_enz  , &
                je       =>     deb%je        , & 
                ec       =>     deb%ec        , &
                gB       =>     deb%gB        , &
                pE       =>     deb%pE          &
                )

    ! Compute the extra carbon after maintenance requirement
    dc = (je - gbt*ec) - mr_mic
    
    if (dc>0._r8) then
    ! There is carbon to support growth activity    
        ! Compute the potential growth rate for different processes
        pot_growth = safe_div(1._r8 , safe_div(1._r8, gmax_mic) + safe_div(1._r8 , dc*yld_mic) )
        pot_enz = safe_div( 1._r8, safe_div(1._r8, pmax_enz) + safe_div(1._r8, dc*yld_enz) )

        ! Compute the actual growth rate
        gB = pot_growth*min( safe_div(dc, safe_div(pot_growth,yld_mic) + pot_enz*yld_enz) , 1._r8 )
        pE = pot_enz*min( safe_div(dc, safe_div(pot_growth,yld_mic) + pot_enz*yld_enz) , 1._r8 )

        ! Calculate residual due to growth
        residual = je - gB*ec - mr_mic - ( safe_div(gB,yld_mic) + safe_div(pE,yld_enz) )

    else 
        gB=0._r8
        pE=0._r8

        residual = je - gB*ec - mr_mic

    end if
end associate
end subroutine deb_grow
  !-------------------------------------------------------------------------------
end module DebGrowMod
