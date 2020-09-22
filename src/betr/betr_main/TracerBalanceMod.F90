module TracerBalanceMod

!
! !DESCRIPTION:
! module contains subroutines to do
! tracer mass balance check

  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use BeTR_decompMod  , only : bounds_type  => betr_bounds_type
  use BeTRTracerType  , only : betrtracer_type
  use TracerFluxType  , only : TracerFlux_type
  use TracerStateType , only : TracerState_type
  use betr_ctrl       , only : iulog  => biulog

  implicit none

  private

  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: begin_betr_tracer_massbalance
  public :: betr_tracer_massbalance_check

  contains



    !--------------------------------------------------------------------------------
    subroutine begin_betr_tracer_massbalance(bounds, col, numf, filter, &
         numfp, filterp, betrtracer_vars, tracerstate_vars, tracerflux_vars, betr_status)
      !
      ! !DESCRIPTION:
      ! Preparing for tracer mass balance check
      !
      ! !USES:
      use BetrStatusType  , only : betr_status_type
      use betr_columnType , only : betr_column_type
      implicit none
      ! !ARGUMENTS:
      type(bounds_type)      , intent(in)    :: bounds
      type(betr_column_type) , intent(in)    :: col
      integer                , intent(in)    :: numf                        ! number of columns in column filter
      integer                , intent(in)    :: filter(:)                   ! column filter
      integer                , intent(in)    :: numfp
      integer                , intent(in)    :: filterp(:)
      type(BeTRtracer_type)  , intent(in)    :: betrtracer_vars
      type(TracerFlux_type)  , intent(inout) :: tracerflux_vars
      type(TracerState_type) , intent(inout) :: tracerState_vars
      type(betr_status_type) , intent(out)   :: betr_status

      ! !LOCAL VARIABLES:
      character(len=256) :: subname='begin_betr_tracer_massbalance'
      integer :: fc, c
      integer :: lbj, ubj

      call betr_status%reset()
      lbj = bounds%lbj;  ubj = bounds%ubj
      call tracerflux_vars%Reset(bounds, numf, filter, numfp, filterp)

      call betr_tracer_mass_summary(bounds, col, lbj, ubj, numf, filter, &
           betrtracer_vars, tracerstate_vars, &
           tracerstate_vars%beg_tracer_molarmass_col, betr_status)

    end subroutine begin_betr_tracer_massbalance

    !--------------------------------------------------------------------------------
    subroutine betr_tracer_massbalance_check(betr_time, bounds, col,  numf, filter, &
         betrtracer_vars, tracerstate_vars, tracerflux_vars, betr_status)!, ldebug)
      !
      ! !DESCRIPTION:
      ! do mass balance check for betr tracers
      !
      ! for solid phase tracers, the only source/sink is biogeochemical production/consumption
      ! and it is currently assumed no solid phase input from atmospheric precipitation (either dry or wet)
      ! the equilibrium fraction is always associated with the (dual)-phase mobile tracer.
      ! However the situation is different for water isotopes, because ice is also part of the
      ! mass budget, and by assuming equilibrium partitioning, the chemical source/sink for ice is not tracked explicitly.
      !
      ! !USES:

      use betr_ctrl       , only : iulog  => biulog, do_betr_output, bgc_type
      use betr_varcon     , only : namec  => bnamec
      use tracer_varcon   , only : catomw,natomw
      use BeTR_TimeMod    , only : betr_time_type
      use BetrStatusType  , only : betr_status_type
      use betr_constants  , only : betr_errmsg_len
      use betr_columnType , only : betr_column_type
      use betr_constants  , only : betr_var_name_length
      implicit none

      ! !ARGUMENTS:
      class(betr_time_type)  , intent(in)    :: betr_time
      type(bounds_type)      , intent(in)    :: bounds
      type(betr_column_type) , intent(in)    :: col
      integer                , intent(in)    :: numf             ! number of columns in column filter
      integer                , intent(in)    :: filter(:)        ! column filter
      type(BeTRtracer_type)  , intent(inout) :: betrtracer_vars
      type(TracerFlux_type)  , intent(inout) :: tracerflux_vars
      type(TracerState_type) , intent(inout) :: tracerState_vars
      type(betr_status_type) , intent(out)   :: betr_status
!      logical                , intent(in)    :: ldebug
      ! !LOCAL VARIABLES:
      integer  :: jj, fc, c, kk
      real(r8) :: dtime
      real(r8) :: atw
      real(r8) :: err_rel, bal_beg, bal_end, bal_flx
      real(r8), parameter :: err_min = 1.e-8_r8
      real(r8), parameter :: err_min_rel=1.e-3_r8
      integer    :: lbj, ubj, jl
      character(len=betr_errmsg_len) :: msg, msg1
      character(len=betr_var_name_length) :: tracername

      call betr_status%reset()
      associate(                                                                            &
           beg_tracer_molarmass      => tracerstate_vars%beg_tracer_molarmass_col         , &
           end_tracer_molarmass      => tracerstate_vars%end_tracer_molarmass_col         , &
           tracer_flx_infl           => tracerflux_vars%tracer_flx_infl_col               , &
           tracer_flx_netpro         => tracerflux_vars%tracer_flx_netpro_col             , &
           tracer_flx_netphyloss     => tracerflux_vars%tracer_flx_netphyloss_col         , &
           is_mobile                 => betrtracer_vars%is_mobile                         , &
           is_volatile               => betrtracer_vars%is_volatile                       , &
           errtracer                 => tracerstate_vars%errtracer_col                    , &
           ngwmobile_tracers         => betrtracer_vars%ngwmobile_tracers                 , &
           do_mass_balchk            => betrtracer_vars%do_mass_balchk                    , &
           ntracers                  => betrtracer_vars%ntracers                            &
           )
      lbj = bounds%lbj
      ubj = bounds%ubj

        call betr_tracer_mass_summary(bounds, col, lbj, ubj, numf, filter, betrtracer_vars, tracerstate_vars, &
             end_tracer_molarmass, betr_status)

        dtime = betr_time%get_step_size()


        do fc = 1, numf
           c = filter(fc)
           !summarize the fluxes
           call tracerflux_vars%flux_summary(lbj, ubj, col, betr_time, c, betrtracer_vars,betr_status)
           call tracerflux_vars%Temporal_average(c,dtime)
           if(betr_status%check_status())return
           do kk = 1, ntracers
              !type1_bgc, only check for volatile tracers
              if(index(bgc_type,'type1_bgc')/=0 .and. .not. is_volatile(kk) .or. (.not. do_mass_balchk(kk)))cycle
              errtracer(c,kk) = beg_tracer_molarmass(c,kk)-end_tracer_molarmass(c,kk)  &
                   + (tracer_flx_netpro(c,kk)-tracer_flx_netphyloss(c,kk))*dtime
              if(abs(errtracer(c,kk))<err_min)then
                 err_rel=1.e-4_r8
              else
                 print*,'nstep=', betr_time%get_nstep()
                 print*,'err,  tracername,  netpro,  netphysloss'
                 print*,errtracer(c,kk),betrtracer_vars%get_tracername(kk),tracer_flx_netpro(c,kk),tracer_flx_netphyloss(c,kk)
                 err_rel = errtracer(c,kk)/max(abs(beg_tracer_molarmass(c,kk)),abs(end_tracer_molarmass(c,kk)))
              endif

              if(abs(err_rel)>err_min_rel .and. do_betr_output)then
                 tracername=betrtracer_vars%get_tracername(kk)
                 write(msg,*)'error exceeds the tolerance for tracer '//trim(tracername), &
                      ' err=',errtracer(c,kk), ' col=',c, ' trcid=',kk,&
                      ' nstep=', betr_time%get_nstep(), &
                      ' netpro=',tracer_flx_netpro(c,kk),&
                      ' netphyloss=',tracer_flx_netphyloss(c,kk),&
                      ' begm=',beg_tracer_molarmass(c,kk), &
                      ' endm=',end_tracer_molarmass(c,kk), &
                      errMsg(mod_filename, __LINE__)
                 call tracerflux_vars%flux_display(c,kk,betrtracer_vars, msg1)
                 msg = trim(msg)//new_line('A')//trim(msg1)
                 call betr_status%set_msg(msg=msg, err=-1)
                 return
              endif
           enddo
        enddo

      end associate

    end subroutine betr_tracer_massbalance_check

    !--------------------------------------------------------------------------------

    subroutine betr_tracer_mass_summary(bounds, col, lbj, ubj, numf, filter, betrtracer_vars,&
       tracerstate_vars, tracer_molarmass_col, betr_status)
      !
      ! !DESCRIPTION:
      ! summarize the column tracer mass
      !
      ! !USES:
      use tracerstatetype , only : tracerstate_type
      use BetrStatusType  , only : betr_status_type
      use betr_columnType , only : betr_column_type
      use betr_ctrl       , only : bgc_type
      implicit none
      ! !ARGUMENTS:
      type(bounds_type)       , intent(in)    :: bounds
      type(betr_column_type)  , intent(in)    :: col
      integer                 , intent(in)    :: lbj, ubj
      integer                 , intent(in)    :: numf                        ! number of columns in column filter
      integer                 , intent(in)    :: filter(:)                   ! column filter
      type(betrtracer_type)   , intent(in)    :: betrtracer_vars             ! betr configuration information
      class(tracerstate_type) , intent(inout) :: tracerstate_vars            ! tracer state variables data structure
      real(r8)                , intent(inout) :: tracer_molarmass_col(bounds%begc:bounds%endc, 1:betrtracer_vars%ntracers)
      type(betr_status_type)  , intent(out)   :: betr_status

      ! !LOCAL VARIABLES:
      integer :: jj, fc, c, kk

      call betr_status%reset()
      ! remove unused dummy args compiler warnings
      if (lbj > 0) continue
      if (ubj > 0) continue

      associate(                                                                            &
           tracer_conc_mobile        => tracerstate_vars%tracer_conc_mobile_col           , &
           tracer_conc_solid_equil   => tracerstate_vars%tracer_conc_solid_equil_col      , &
           tracer_conc_frozen        => tracerstate_vars%tracer_conc_frozen_col           , &
           dz                        => col%dz                                            , &
           ntracers                  => betrtracer_vars%ntracers                          , &
           is_adsorb                 => betrtracer_vars%is_adsorb                         , &
           adsorbid                  => betrtracer_vars%adsorbid                          , &
           is_volatile               => betrtracer_vars%is_volatile                       , &
           is_frozen                 => betrtracer_vars%is_frozen                         , &
           frozenid                  => betrtracer_vars%frozenid                            &
           )
        do jj = 1,   ntracers
           if(index(bgc_type,'type1_bgc')/=0 .and. .not. is_volatile(jj))cycle
           do fc = 1, numf
              c = filter(fc)

              tracer_molarmass_col(c,jj) = &
                 tracerstate_vars%int_mass_mobile_col(1,ubj,c,jj,dz(c,1:ubj),betr_status)
              if(betr_status%check_status())return

              if(is_frozen(jj))then
                 tracer_molarmass_col(c,jj) = tracer_molarmass_col(c,jj) + &
                      tracerstate_vars%int_mass_frozen_col(1,ubj,c,&
                      frozenid(jj),dz(c,1:ubj),betr_status)
                 if(betr_status%check_status())return
              endif

           enddo
        enddo
      end associate
    end subroutine betr_tracer_mass_summary
end module TracerBalanceMod
