module H2OIsotopeBGCReactionsType

#include "shr_assert.h"

!
! !DESCRIPTION
! This is used to do O18 isotope simulations involving H2O(18) and COO(18).
! For H2O(18), the advective part is assumed to follow the equation
! \frac{\partial Rw*vsm}{\partial t} = \frac{\partial Rw*q}{\partial z} - Transp*R
! and the evaporative part is assumed to follow the betr diffusion equation
! with evapative flux as top boundary condition.
! After diffusive and aqueous transport, an equilibration is assumed to occur simultaneously between
! solid (ice), liquid and vapor phases
! The formulation adopted by BeTR is similar as that proposed in Braud et al. (2005) for the SiSPAT-isotope model.
! However, because CLM does not consider water vapor during water movement calculation, the inclusion of water vapor
! diffusion may cause some consistency problems, even though this problem is partially fixed using
! the prescribed top boundary condition.

! HISTORY:
! Created by Jinyun Tang, Jan 15nd, 2015
! !USES

  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type
  use BGCReactionsMod       , only : bgc_reaction_type
  use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
  use abortutils            , only : endrun
implicit none

  save
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! !PUBLIC TYPES:
  public :: bgc_reaction_h2oiso_type

  type, extends(bgc_reaction_type) :: &
     bgc_reaction_h2oiso_type
     private
   contains
     procedure :: Init_betrbgc                  ! initialize betr bgc
     procedure :: set_boundary_conditions       ! set top/bottom boundary conditions for various tracers
     procedure :: calc_bgc_reaction             ! doing bgc calculation
     procedure :: init_boundary_condition_type  ! initialize type of top boundary conditions
     procedure :: do_tracer_equilibration       ! do equilibrium tracer chemistry
     procedure :: initCold
     procedure :: readParams
     procedure :: lsm_betr_flux_state_receive
     procedure :: init_betr_lsm_bgc_coupler
   end type bgc_reaction_h2oiso_type

   interface bgc_reaction_h2oiso_type
     module procedure constructor
   end interface bgc_reaction_h2oiso_type

  contains
!-------------------------------------------------------------------------------
  type(bgc_reaction_h2oiso_type) function constructor()
  !
  ! ! DESCRIPTION
  ! create an object of type bgc_reaction_h2oiso_type.
  ! Right now it is purposely left empty

  end function constructor

!-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
  !
  ! DESCRIPTIONS
  ! initialize boundary condition types

  use TracerBoundaryCondType, only : tracerboundarycond_type
  use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
  use BeTRTracerType        , only : betrtracer_type
  use BeTR_decompMod           , only : betr_bounds_type

  class(bgc_reaction_h2oiso_type)   , intent(in) :: this
  type(BeTRtracer_type )            , intent(in) :: betrtracer_vars
  type(betr_bounds_type)                 , intent(in) :: bounds
  type(tracerboundarycond_type)     , intent(in) :: tracerboundarycond_vars

  tracerboundarycond_vars%topbc_type(:) = bndcond_as_conc

  !only the water vapor is set with prescribed flux based boundary condition, Riley et al. (2002, GBC)
  !had a discussion about this.
  tracerboundarycond_vars%topbc_type(betrtracer_vars%id_trc_d_h2o)   = bndcond_as_flux
  tracerboundarycond_vars%topbc_type(betrtracer_vars%id_trc_o18_h2o) = bndcond_as_flux
  tracerboundarycond_vars%topbc_type(betrtracer_vars%id_trc_blk_h2o) = bndcond_as_flux

  end subroutine init_boundary_condition_type


!-------------------------------------------------------------------------------

  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars)
  !
  ! DESCRIPTION
  ! initialize the betrbgc
  use BeTRTracerType        , only : betrtracer_type
  use MathfuncMod           , only : addone
  use BeTR_decompMod           , only : betr_bounds_type

  class(bgc_reaction_h2oiso_type) ,  intent(in) :: this
  type(betr_bounds_type)          ,  intent(in) :: bounds
  integer                         ,  intent(in) :: lbj, ubj
  type(BeTRtracer_type )          ,intent(inout) :: betrtracer_vars
  !local variables
  character(len=*)       , parameter :: subname ='Init_betrbgc'
  integer :: itemp_gwm
  integer :: itemp_g
  integer :: itemp_s
  integer :: itemp_gwm_grp
  integer :: dum
  integer :: itemp_grp, itemp_v, itemp_vgrp, itemp_adsgrp
  integer :: itemp_frz

  itemp_gwm     = 0
  itemp_g       = 0
  itemp_s       = 0
  itemp_gwm_grp = 0
  itemp_adsgrp  = 0
  betrtracer_vars%id_trc_n2  = addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);
  betrtracer_vars%id_trc_o2  = addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);
  betrtracer_vars%id_trc_ar  = addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);
  betrtracer_vars%id_trc_co2x= addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);
  betrtracer_vars%id_trc_ch4 = addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);
  betrtracer_vars%id_trc_blk_h2o = addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);
  betrtracer_vars%id_trc_o18_h2o = addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);
  betrtracer_vars%id_trc_d_h2o = addone(itemp_gwm); itemp_gwm_grp = addone(itemp_gwm_grp);


  betrtracer_vars%ngwmobile_tracers      = itemp_gwm;   betrtracer_vars%ngwmobile_tracer_groups= itemp_gwm_grp;
  betrtracer_vars%nvolatile_tracer_groups= itemp_gwm;
  betrtracer_vars%nmem_max               = 1


  call betrtracer_vars%Init()


  itemp_grp = 0    !group id
  itemp_v = 0      !volatile id
  itemp_vgrp = 0   !volatile group
  itemp_frz = 0    !frozen tracer id
  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_n2, trc_name='N2'  ,      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
       trc_group_mem= 1,  is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))

  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_o2, trc_name='O2'  ,      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))

  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_ar, trc_name='AR'  ,      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))

  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_co2x, trc_name='CO2x',    &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp)  , &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &
       trc_volatile_group_id = addone(itemp_vgrp) )

  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_ch4, trc_name='CH4',      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))

  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_blk_h2o, trc_name='BLK_H2O' ,   &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp)      ,   &
       trc_group_mem = 1, is_trc_diffusive =.false., is_trc_volatile=.true.                   ,   &
       trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp)          ,   &
       is_trc_h2o=.true., trc_vtrans_scal=1._r8, is_trc_frozen=.true.,  trc_frozenid = addone(itemp_frz))

  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_o18_h2o, trc_name='O18_H2O' ,   &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id =  addone(itemp_grp)     ,   &
       trc_group_mem = 1, is_trc_diffusive =.false., is_trc_volatile=.true.                   ,   &
       trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp)          ,   &
       is_trc_h2o=.true., trc_vtrans_scal=1._r8, is_trc_frozen=.true.,  trc_frozenid = addone(itemp_frz))

  call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_d_h2o, trc_name='D_H2O' ,   &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id =  addone(itemp_grp)     ,   &
       trc_group_mem = 1, is_trc_diffusive =.false., is_trc_volatile=.true.                   ,   &
       trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp)          ,   &
       is_trc_h2o=.true., trc_vtrans_scal=1._r8, is_trc_frozen=.true.,  trc_frozenid = addone(itemp_frz))

  end subroutine Init_betrbgc


!-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, dz_top, betrtracer_vars, &
       waterflux_vars, tracerboundarycond_vars)
  !
  ! DESCRIPTION
  ! set up boundary conditions for tracer movement
  !
  use clm_varctl            , only : iulog
  use TracerBoundaryCondType, only : tracerboundarycond_type
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use BeTR_decompMod        , only : betr_bounds_type
  use BeTRTracerType        , only : betrtracer_type
  use BeTR_WaterfluxType    , only : betr_waterflux_type
  use clm_varcon            , only : denh2o


  class(bgc_reaction_h2oiso_type), intent(in) :: this
  type(betr_bounds_type)         , intent(in) :: bounds                     !
  integer                           , intent(in) :: num_soilc                  ! number of columns in column filter_soilc
  integer                           , intent(in) :: filter_soilc(:)            ! column filter_soilc
  type(betrtracer_type)             , intent(in) :: betrtracer_vars            !
  real(r8)                          , intent(in) :: dz_top(bounds%begc: )      !
  type(betr_waterflux_type), intent(in) :: waterflux_vars             !
  type(tracerboundarycond_type)     , intent(inout) :: tracerboundarycond_vars !

  !local variables
  integer :: fc, c
  character(len=255) :: subname = 'set_boundary_conditions'

  SHR_ASSERT_ALL((ubound(dz_top)                == (/bounds%endc/)),   errMsg(mod_filename,__LINE__))


  associate(                                                           &
    qflx_gross_evap_soil => waterflux_vars%qflx_gross_evap_soil_col    &
  )
  !eventually, the following code will be implemented using polymorphism
  !for simplicity, all gases other than water vapor are set with fixed concentration based boundary conditions
  !
  do fc = 1, num_soilc
    c = filter_soilc(fc)
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2)      = 32.8_r8                         !mol m-3, contant boundary condition, as concentration
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o2)      = 8.78_r8                         !mol m-3, contant boundary condition, as concentration
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ar)      = 0.3924_r8                       !mol m-3, contant boundary condition, as concentration
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_co2x)    = 0.0168_r8                       !mol m-3, contant boundary condition, as concentration
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ch4)     = 6.939e-5_r8                     !mol m-3, contant boundary condition, as concentration
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_blk_h2o) = -qflx_gross_evap_soil(c)        !kg m-2-s, not diffusive water vapor transport
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o18_h2o) = -qflx_gross_evap_soil(c)        !kg m-2-s, not diffusive water vapor transport
    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_d_h2o)   = -qflx_gross_evap_soil(c)        !kg m-2-s, not diffusive water vapor transport

    tracerboundarycond_vars%bot_concflux_col(c,1,:)                         = 0._r8                           !zero flux boundary condition
    tracerboundarycond_vars%condc_toplay_col(c,betrtracer_vars%id_trc_n2)   = 2._r8*1.837e-5_r8/dz_top(c)     !m/s surface conductance, this will be represented with Tang-Riley scheme (HESS, 2013)
    tracerboundarycond_vars%condc_toplay_col(c,betrtracer_vars%id_trc_o2)   = 2._r8*1.713e-5_r8/dz_top(c)     !m/s surface conductance, this will be represented with Tang-Riley scheme (HESS, 2013)
    tracerboundarycond_vars%condc_toplay_col(c,betrtracer_vars%id_trc_ar)   = 2._r8*1.532e-5_r8/dz_top(c)     !m/s surface conductance, this will be represented with Tang-Riley scheme (HESS, 2013)
    tracerboundarycond_vars%condc_toplay_col(c,betrtracer_vars%id_trc_co2x) = 2._r8*1.399e-5_r8/dz_top(c)     !m/s surface conductance, this will be represented with Tang-Riley scheme (HESS, 2013)
    tracerboundarycond_vars%condc_toplay_col(c,betrtracer_vars%id_trc_ch4)  = 2._r8*1.808e-5_r8/dz_top(c)     !m/s surface conductance, this will be represented with Tang-Riley scheme (HESS, 2013)
  enddo
  end associate
  end subroutine set_boundary_conditions

!-------------------------------------------------------------------------------

  subroutine calc_bgc_reaction(this, bounds, lbj, ubj, num_soilc, filter_soilc, &
       num_soilp,filter_soilp, jtops, dtime, betrtracer_vars, tracercoeff_vars,  cnstate_vars,    &
       tracerstate_vars, tracerflux_vars, tracerboundarycond_vars, plant_soilbgc)

  !
  ! do bgc reaction
  ! eventually this will be an abstract subroutine, but now I use the select case approach for a quick and dirty implementation.
  !USES
  !
  ! !USES:
  use BeTR_decompMod           , only : betr_bounds_type
  use tracerfluxType           , only : tracerflux_type
  use tracerstatetype          , only : tracerstate_type
  use tracercoeffType          , only : tracercoeff_type
  use BetrTracerType           , only : betrtracer_type
  use BeTR_CNStateType         , only : betr_cnstate_type
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use TracerBoundaryCondType   , only : tracerboundarycond_type
  use ColumnType               , only : col
  !ARGUMENTS
  class(bgc_reaction_h2oiso_type)     , intent(in)    :: this                       !
  type(betr_bounds_type), intent(in) :: bounds ! bounds
  integer                             , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
  integer                             , intent(in)    :: filter_soilc(:)            ! column filter_soilc
  integer                             , intent(in)    :: num_soilp                  !
  integer                             , intent(in)    :: filter_soilp(:)            ! pft filter
  integer                             , intent(in)    :: jtops(bounds%begc: )       ! top index of each column
  integer                             , intent(in)    :: lbj, ubj                   ! lower and upper bounds, make sure they are > 0
  real(r8)                            , intent(in)    :: dtime                      ! model time step
  type(betrtracer_type)               , intent(in)    :: betrtracer_vars            ! betr configuration information
  type(betr_cnstate_type)             , intent(inout) :: cnstate_vars
  type(tracercoeff_type)              , intent(in)    :: tracercoeff_vars           !
  type(tracerstate_type)              , intent(inout) :: tracerstate_vars           !
  type(tracerflux_type)               , intent(inout) :: tracerflux_vars            !
  type(tracerboundarycond_type)       , intent(inout) :: tracerboundarycond_vars !
  class(plant_soilbgc_type)           , intent(inout) :: plant_soilbgc
  character(len=*)                    , parameter     :: subname ='calc_bgc_reaction'

  integer :: jj, c, fc, ll
  integer, parameter :: nh2o_trcs=3
  integer :: jjs(nh2o_trcs), kk
  real(r8):: tot0, tot1

  associate(                                                                                  &
    tracer_mobile_phase            => tracerstate_vars%tracer_conc_mobile_col               , &
    tracer_gwdif_concflux_top_col  => tracerboundarycond_vars%tracer_gwdif_concflux_top_col , &
    volatileid                     =>  betrtracer_vars%volatileid                           , & !
    tracer_flx_dif                 =>  tracerflux_vars%tracer_flx_dif_col                   , & !
    id_trc_blk_h2o                 => betrtracer_vars%id_trc_blk_h2o                        , &
    id_trc_o18_h2o                 => betrtracer_vars%id_trc_o18_h2o                        , &
    id_trc_d_h2o                   => betrtracer_vars%id_trc_d_h2o                            &
  )



  !apply the evaporation to the water tracer, the following is a hack to avoid the
  !inconsistency between water vapor transport in betr and the hydrology code
  !in the future, when the hdyrology code is corrected, the following will be gone, jyt Feb, 17, 2016
  jjs = (/id_trc_blk_h2o,id_trc_o18_h2o, id_trc_d_h2o/)

  do kk = 1, nh2o_trcs
    jj = jjs(kk)
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      tracer_mobile_phase(c,1,jj) = tracer_mobile_phase(c,1,jj) + tracer_gwdif_concflux_top_col(c,1,jj)*dtime/col%dz(c,1)
      if(tracer_mobile_phase(c,1,jj) < 0._r8)then
        do ll = 1, 2
          tot0 = tracer_mobile_phase(c,ll,jj)*col%dz(c,ll)
          tot1 = tracer_mobile_phase(c,ll+1,jj)*col%dz(c,ll+1)
          tot1 = tot1 + tot0
          tracer_mobile_phase(c,ll,jj) = 0._r8
          tracer_mobile_phase(c,ll+1,jj) = tot1/col%dz(c,ll+1)
          if(tot1>0._r8)exit
        enddo
        !the following should rarely occur, so when it occur, end with a warning
        if(tot1<0._r8)then
          call endrun('negative H2O tracer '//errMsg(mod_filename, __LINE__))
        endif
      endif
      tracer_flx_dif(c,volatileid(jj)) = tracer_flx_dif(c,volatileid(jj))- tracer_gwdif_concflux_top_col(c,1,jj) * dtime
    enddo
  enddo

  end associate
  end subroutine calc_bgc_reaction

!-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars)
 !
  ! DESCRIPTIONS
  ! requilibrate tracers that has solid and mobile phases
  ! using the theory of mass action. When the redox-ladder is on, this
  ! subroutine will update the change of pH due to tracer transport, or
  ! USES
  !
  use tracerstatetype       , only : tracerstate_type
  use tracercoeffType       , only : tracercoeff_type
  use BeTRTracerType        , only : betrtracer_type
  use BeTR_decompMod        , only : betr_bounds_type


  class(bgc_reaction_h2oiso_type),    intent(in) :: this

  type(betr_bounds_type),      intent(in) :: bounds
  integer,                intent(in) :: lbj, ubj
  integer,                intent(in) :: jtops(bounds%begc: )        ! top label of each column
  integer,                intent(in) :: num_soilc
  integer,                intent(in) :: filter_soilc(:)
  type(betrtracer_type),  intent(in) :: betrtracer_vars
  type(tracercoeff_type), intent(in) :: tracercoeff_vars
  type(tracerstate_type), intent(inout) :: tracerstate_vars

  character(len=255) :: subname = 'do_tracer_equilibration'
  integer   :: j, fc, c
  integer   :: trc_id1, trc_id2


  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__))

  associate(                                                                         &
    aqu2equilscef                   => tracercoeff_vars%aqu2equilsolidcef_col          , &
    aqu2bulkcef_mobile              => tracercoeff_vars%aqu2bulkcef_mobile_col     , &
    tracer_solid_phase_equil        => tracerstate_vars%tracer_conc_solid_equil_col, &
    tracer_mobile_phase             => tracerstate_vars%tracer_conc_mobile_col       &
  )
  !depending on the simulation type, an implementation of aqueous chemistry will be
  !employed to separate out the adsorbed phase
  !It should be noted that this formulation excludes the use of linear isotherm, which
  !can be integrated through the retardation factor
  !assuming equilibrium fractionation between ice/water/vapor, calculate the equilibrium solid phase concentrations
  !this might introduce some bias, because soil moisture profile is updated from phase change and advective transport, while
  !the equilibration adjusts continously as water flows or phase change occurs



  trc_id1 = betrtracer_vars%id_trc_o18_h2o
  trc_id2 = betrtracer_vars%id_trc_o18_h2o_ice

  !the following code is replaced with diagnose_dtracer_freeze_thaw in TracerParamsMod

!  if(trc_id1>0)then
!    call do_h2o_isotope_equilibration(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
!      aqu2bulkcef_mobile(bounds%begc:bounds%endc, lbj:ubj, trc_id1)        , &
!      aqu2equilscef(bounds%begc:bounds%endc, lbj:ubj, trc_id2)             , &
!      tracer_solid_phase_equil(bounds%begc:bounds%endc, lbj:ubj, trc_id2)  , &
!      tracer_mobile_phase(bounds%begc:bounds%endc, lbj:ubj, trc_id1))

!  endif

  end associate


  end subroutine do_tracer_equilibration


!-------------------------------------------------------------------------------

  subroutine do_h2o_isotope_equilibration(bounds, lbj, ubj, jtops, numf, filter, aqu2bulkcef, &
    aqu2equilscef, tracer_solid_phase_equil, tracer_mobile_phase)
  !
  ! Diagnose solid phase tracer
  !
  use BeTR_decompMod           , only : betr_bounds_type
  implicit none
  type(betr_bounds_type),      intent(in) :: bounds
  integer,                intent(in) :: lbj, ubj
  integer,                intent(in) :: jtops(bounds%begc: )        ! top label of each column
  integer,                intent(in) :: numf
  integer,                intent(in) :: filter(:)
  real(r8),               intent(in) :: aqu2equilscef(bounds%begc: , lbj: )
  real(r8),               intent(in) :: aqu2bulkcef(bounds%begc: , lbj: )
  real(r8),               intent(in) :: tracer_mobile_phase(bounds%begc: , lbj: )
  real(r8),            intent(inout) :: tracer_solid_phase_equil(bounds%begc: ,lbj: )

  real(r8)  :: frac
  real(r8)  :: tracer_conc
  integer   :: c, fc, j

  SHR_ASSERT_ALL((ubound(aqu2equilscef) == (/bounds%endc, ubj/)), errMsg(mod_filename,__LINE__))
  SHR_ASSERT_ALL((ubound(aqu2bulkcef)   == (/bounds%endc, ubj/)), errMsg(mod_filename,__LINE__))
  SHR_ASSERT_ALL((ubound(tracer_solid_phase_equil) == (/bounds%endc, ubj/)), errMsg(mod_filename,__LINE__))
  SHR_ASSERT_ALL((ubound(tracer_mobile_phase) == (/bounds%endc, ubj/)), errMsg(mod_filename,__LINE__))

  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(j>=jtops(c))then
        !obtains total concentration
        tracer_conc= tracer_mobile_phase(c,j) + tracer_solid_phase_equil(c,j)
        !obtain the equilibrium conversion parameter
        frac = aqu2bulkcef(c,j) / (aqu2bulkcef(c,j) + aqu2equilscef(c, j))

        !tracer_mobile_phase(c,j) = tracer_conc * frac
        tracer_solid_phase_equil(c,j) = tracer_conc - tracer_mobile_phase(c,j)
      endif
    enddo
  enddo

  end subroutine do_h2o_isotope_equilibration

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, betrtracer_vars, waterstate_vars, tracerstate_vars)
    !
    ! !USES:
    !
    use BeTR_decompMod           , only : betr_bounds_type
    use BeTRTracerType           , only : BeTRTracer_Type
    use tracerstatetype          , only : tracerstate_type
    use BeTR_WaterstateType      , only : betr_waterstate_type
    use BeTR_PatchType           , only : pft  => betr_pft
    use clm_varcon               , only : spval, ispval
    use landunit_varcon          , only : istsoil, istcrop
    use LandunitType             , only : lun
    use ColumnType               , only : col
    use clm_varcon               , only : denh2o
    use clm_varpar               , only : nlevtrc_soil
    ! !ARGUMENTS:
    class(bgc_reaction_h2oiso_type) , intent(in)    :: this
    type(betr_bounds_type)                 , intent(in)    :: bounds
    type(BeTRTracer_Type)             , intent(in)    :: betrtracer_vars
    type(betr_waterstate_type)             , intent(in)    :: waterstate_vars
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars

    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, k, j
    integer :: fc                                        ! filter index
    integer               :: begc, endc
    integer               :: begg, endg
    integer               :: trcid
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !-----------------------------------------------------------------------


    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
         if(betrtracer_vars%ngwmobile_tracers>0)then
           tracerstate_vars%tracer_conc_mobile_col(c,:,:)        = spval
           tracerstate_vars%tracer_conc_surfwater_col(c,:)       = spval
           tracerstate_vars%tracer_conc_aquifer_col(c,:)         = spval
           tracerstate_vars%tracer_conc_grndwater_col(c,:)       = spval
         endif
         if(betrtracer_vars%ntracers > betrtracer_vars%ngwmobile_tracers)then
           tracerstate_vars%tracer_conc_solid_passive_col(c,:,:) = spval
         endif
         if(betrtracer_vars%nsolid_equil_tracers>0)then
           tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = spval
         endif
       endif
       tracerstate_vars%tracer_soi_molarmass_col(c,:)            = spval

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
         !dual phase tracers

         tracerstate_vars%tracer_conc_mobile_col(c,:, :)          = 0._r8
         tracerstate_vars%tracer_conc_surfwater_col(c,:)          = 0._r8
         tracerstate_vars%tracer_conc_aquifer_col(c,:)            = 0._r8
         tracerstate_vars%tracer_conc_grndwater_col(c,:)          = 0._r8

         !solid tracers
         if(betrtracer_vars%ngwmobile_tracers < betrtracer_vars%ntracers)then
           tracerstate_vars%tracer_conc_solid_passive_col(c,:,:) = 0._r8
         endif

         if(betrtracer_vars%nsolid_equil_tracers>0)then
           tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = 0._r8
         endif
         tracerstate_vars%tracer_soi_molarmass_col(c,:)          = 0._r8

         !set for o18_h2o, assuming no fractionation, which is equivalent to assuming concentration equals 1
         trcid = betrtracer_vars%id_trc_blk_h2o
         tracerstate_vars%tracer_conc_grndwater_col(c,trcid) = denh2o
         do j = 1, nlevtrc_soil
           tracerstate_vars%tracer_conc_mobile_col(c,j,trcid) = 1._r8 * waterstate_vars%h2osoi_liq_col(c,j)/col%dz(c,j)
           tracerstate_vars%tracer_conc_frozen_col(c,j,betrtracer_vars%frozenid(trcid)) = &
                1._r8 * waterstate_vars%h2osoi_ice_col(c,j)/col%dz(c,j)
         enddo

         trcid = betrtracer_vars%id_trc_o18_h2o
         tracerstate_vars%tracer_conc_grndwater_col(c,trcid) = denh2o
         do j = 1, nlevtrc_soil
           tracerstate_vars%tracer_conc_mobile_col(c,j,trcid) = 1._r8 * waterstate_vars%h2osoi_liq_col(c,j)/col%dz(c,j)
           tracerstate_vars%tracer_conc_frozen_col(c,j,betrtracer_vars%frozenid(trcid)) = &
                1._r8 * waterstate_vars%h2osoi_ice_col(c,j)/col%dz(c,j)
         enddo

         trcid = betrtracer_vars%id_trc_d_h2o
         tracerstate_vars%tracer_conc_grndwater_col(c,trcid) = denh2o
         do j = 1, nlevtrc_soil
           tracerstate_vars%tracer_conc_mobile_col(c,j,trcid) = 1._r8 * waterstate_vars%h2osoi_liq_col(c,j)/col%dz(c,j)
           tracerstate_vars%tracer_conc_frozen_col(c,j,betrtracer_vars%frozenid(trcid)) = &
                1._r8 * waterstate_vars%h2osoi_ice_col(c,j)/col%dz(c,j)
         enddo

       endif
     enddo

  end subroutine InitCold



  !-----------------------------------------------------------------------
  subroutine readParams(this, ncid, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! read in module specific parameters
    !
    ! !USES:

    use ncdio_pio                , only : file_desc_t
    use BeTRTracerType           , only : BeTRTracer_Type
    ! !ARGUMENTS:
    class(bgc_reaction_h2oiso_type)   , intent(in)    :: this
    type(BeTRTracer_Type)             , intent(inout) :: betrtracer_vars
    type(file_desc_t)                 , intent(inout) :: ncid  ! pio netCDF file id

    !do nothing here for the moment, but contents will eventually be filled in here

  end subroutine readParams

  !-------------------------------------------------------------------------------
  subroutine lsm_betr_flux_state_receive(this, bounds, num_soilc, filter_soilc,  &
       tracerstate_vars, tracerflux_vars,  betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! do flux and state variable change between betr and lsm.
    !
    ! !USES:
    use shr_kind_mod             , only : r8 => shr_kind_r8
    use tracerfluxType           , only : tracerflux_type
    use tracerstatetype          , only : tracerstate_type
    use decompMod                , only : bounds_type
    use BeTRTracerType           , only : BeTRTracer_Type
    use BeTR_decompMod           , only : betr_bounds_type


    ! !ARGUMENTS:
    class(bgc_reaction_h2oiso_type) , intent(in)    :: this               !
    type(betr_bounds_type)                 , intent(in)    :: bounds             ! bounds
    integer                           , intent(in)    :: num_soilc
    integer                           , intent(in)    :: filter_soilc(:)
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars    ! betr configuration information
    type(tracerstate_type)            , intent(in)    :: tracerstate_vars   !
    type(tracerflux_type)             , intent(in)    :: tracerflux_vars    !




  end subroutine lsm_betr_flux_state_receive


  !-------------------------------------------------------------------------------


  subroutine init_betr_lsm_bgc_coupler(this, bounds,                     &
       plant_soilbgc, betrtracer_vars, tracerstate_vars, cnstate_vars,  ecophyscon_vars)

    ! !DESCRIPTION:
    ! initialize the bgc coupling between betr and lsm
    !
    ! !USES:
    !use clm_instMod
    use tracerstatetype          , only : tracerstate_type
    use BetrTracerType           , only : betrtracer_type
    use clm_varpar               , only : nlevtrc_soil
    use landunit_varcon          , only : istsoil, istcrop
    use PlantSoilBGCMod          , only : plant_soilbgc_type
    use EcophysConType           , only : ecophyscon_type
    use BeTR_CNStateType         , only : betr_cnstate_type
    use BeTR_decompMod           , only : betr_bounds_type


    ! !ARGUMENTS:
    class(bgc_reaction_h2oiso_type)    , intent(in)    :: this               !
    type(betr_bounds_type)                  , intent(in)    :: bounds             !
    type(tracerstate_type)             , intent(inout) :: tracerstate_vars   !
    type(betrtracer_type)              , intent(in)    :: betrtracer_vars    ! betr configuration information
    class(plant_soilbgc_type)          , intent(inout) :: plant_soilbgc !
    type(ecophyscon_type)              , intent(in)    :: ecophyscon_vars
    type(betr_cnstate_type)            , intent(in)    :: cnstate_vars


  end subroutine init_betr_lsm_bgc_coupler

end module H2OIsotopeBGCReactionsType
