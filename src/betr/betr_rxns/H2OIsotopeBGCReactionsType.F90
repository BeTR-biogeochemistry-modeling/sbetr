module H2OIsotopeBGCReactionsType

#include "bshr_assert.h"

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
  use bshr_assert_mod, only : shr_assert
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use bshr_assert_mod, only : shr_assert_any
  use bshr_log_mod             , only : errMsg => shr_log_errMsg
  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use BGCReactionsMod          , only : bgc_reaction_type
  use tracer_varcon            , only : bndcond_as_conc, bndcond_as_flux
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BetrStatusType           , only : betr_status_type
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! !PUBLIC TYPES:
  public :: bgc_reaction_h2oiso_type

  type, extends(bgc_reaction_type) :: &
     bgc_reaction_h2oiso_type
     private
    integer, public :: parcol

   contains
     procedure :: Init_betrbgc                  ! initialize betr bgc
     procedure :: set_boundary_conditions       ! set top/bottom boundary conditions for various tracers
     procedure :: calc_bgc_reaction             ! doing bgc calculation
     procedure :: init_boundary_condition_type  ! initialize type of top boundary conditions
     procedure :: do_tracer_equilibration       ! do equilibrium tracer chemistry
     procedure :: initCold
     procedure :: retrieve_biogeoflux
     procedure :: set_kinetics_par
     procedure :: retrieve_lnd2atm
     procedure, private :: readParams
     procedure :: retrieve_biostates
     procedure :: debug_info
     procedure :: set_bgc_spinup
     procedure :: UpdateParas
     procedure :: init_iP_prof
     procedure :: reset_biostates
     procedure :: SetParCols
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

    type(bgc_reaction_h2oiso_type), allocatable :: bgc
    allocate(bgc)
    constructor = bgc

  end function constructor
  !-------------------------------------------------------------------------------
  subroutine UpdateParas(this, bounds, lbj, ubj, bstatus)
  use BeTR_decompMod         , only : betr_bounds_type
  implicit none
  class(bgc_reaction_h2oiso_type)         , intent(inout)    :: this
  type(betr_bounds_type)                    , intent(in)    :: bounds
  integer                              , intent(in)    :: lbj, ubj        ! lower and upper bounds, make sure they are > 0
  type(betr_status_type)           , intent(out)   :: bstatus
  integer :: c, j
  if (this%dummy_compiler_warning) continue
  call bstatus%reset()
  !do nothing
  end subroutine UpdateParas
  !----------------------------------------------------------------------
  subroutine init_iP_prof(this, bounds, lbj, ubj, biophysforc, tracers, tracerstate_vars)
  !
  !DESCRIPTION
  ! set up initial inorganic P profile
  use tracer_varcon, only : patomw
  use tracerstatetype        , only : tracerstate_type
  use BeTRTracerType         , only : betrtracer_type
  use BeTR_decompMod         , only : betr_bounds_type
  implicit none
  ! !ARGUMENTS:
  class(bgc_reaction_h2oiso_type)         , intent(inout)    :: this
  type(betr_bounds_type)                        , intent(in) :: bounds
  integer                                  , intent(in) :: lbj, ubj
  type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc
  type(BeTRtracer_type)                    , intent(inout) :: tracers
  type(tracerstate_type)                   , intent(inout) :: tracerstate_vars


  if (this%dummy_compiler_warning) continue
  if (bounds%begc > 0) continue

  end subroutine init_iP_prof
  !----------------------------------------------------------------------
  subroutine set_kinetics_par(this, lbj, ubj,nactpft, plantNutkinetics, tracers, tracercoeff_vars)
  use PlantNutKineticsMod, only : PlantNutKinetics_type
  use tracercoeffType          , only : tracercoeff_type
  use BeTRTracerType           , only : betrtracer_type
  implicit none
  ! !ARGUMENTS:
  class(bgc_reaction_h2oiso_type)         , intent(inout)    :: this                       !
  class(PlantNutKinetics_type), intent(in) :: plantNutkinetics
  type(betrtracer_type), intent(in) :: tracers
  type(tracercoeff_type), intent(inout) :: tracercoeff_vars
  integer, intent(in) :: lbj, ubj
  integer, intent(in) :: nactpft
  if (this%dummy_compiler_warning) continue

  end subroutine set_kinetics_par

  !-------------------------------------------------------------------------------
  subroutine set_bgc_spinup(this, bounds, lbj, ubj,  biophysforc, &
  tracers, tracerstate_vars)
  use tracerstatetype        , only : tracerstate_type
  use BeTRTracerType         , only : betrtracer_type
  use BeTR_decompMod         , only : betr_bounds_type
  implicit none
    class(bgc_reaction_h2oiso_type)         , intent(inout)    :: this                       !
    type(betr_bounds_type)                       , intent(in) :: bounds
    integer                                 , intent(in) :: lbj, ubj
    type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc
    type(BeTRtracer_type)                   , intent(inout) :: tracers
    type(tracerstate_type)                  , intent(inout) :: tracerstate_vars

    if (this%dummy_compiler_warning) continue
    if (bounds%begc > 0) continue

  end subroutine set_bgc_spinup

!-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
  !
  ! DESCRIPTIONS
  ! initialize boundary condition types
  ! USES
  use TracerBoundaryCondType , only : tracerboundarycond_type
  use tracer_varcon          , only : bndcond_as_conc, bndcond_as_flux
  use BeTRTracerType         , only : betrtracer_type
  use BeTR_decompMod         , only : betr_bounds_type
  implicit none
  !arguments
  class(bgc_reaction_h2oiso_type) , intent(inout) :: this
  type(BeTRtracer_type )          , intent(in) :: betrtracer_vars
  type(betr_bounds_type)          , intent(in) :: bounds
  type(tracerboundarycond_type)   , intent(in) :: tracerboundarycond_vars

  !by default top boundary conditions are prescribed as concentrations
  tracerboundarycond_vars%topbc_type(:) = bndcond_as_conc

  !only the water vapor is set with prescribed flux based boundary condition, Riley et al. (2002, GBC)
  !had a discussion about this.
  tracerboundarycond_vars%topbc_type(betrtracer_vars%id_trc_d_h2o)   = bndcond_as_flux
  tracerboundarycond_vars%topbc_type(betrtracer_vars%id_trc_o18_h2o) = bndcond_as_flux
  tracerboundarycond_vars%topbc_type(betrtracer_vars%id_trc_blk_h2o) = bndcond_as_flux

  end subroutine init_boundary_condition_type


   !----------------------------------------------------------------------
   subroutine retrieve_lnd2atm(this, bounds, num_soilc, filter_soilc, tracerflux_vars, &
   betrtracer_vars, biogeo_flux)

   use tracerfluxType           , only : tracerflux_type
   use BeTR_decompMod           , only : betr_bounds_type
   use BeTRTracerType           , only : BeTRTracer_Type
   use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
   implicit none
   class(bgc_reaction_h2oiso_type) , intent(inout) :: this
   type(betr_bounds_type)           , intent(in)    :: bounds                      ! bounds
   integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
   integer                          , intent(in)    :: filter_soilc(:)             ! column filter
   type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
   type(tracerflux_type)            , intent(in)    :: tracerflux_vars
   type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux


   if (this%dummy_compiler_warning) continue
   if (bounds%begc > 0)             continue

   end subroutine retrieve_lnd2atm

!-------------------------------------------------------------------------------

  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars, namelist_buffer, bstatus)
  !
  ! DESCRIPTION
  ! initialize the betrbgc
  ! USES
  use BeTRTracerType , only : betrtracer_type
  use MathfuncMod    , only : addone
  use BeTR_decompMod , only : betr_bounds_type
  use BetrStatusType , only : betr_status_type
  use gbetrType      , only : gbetr_type
  implicit none
  class(bgc_reaction_h2oiso_type) , intent(inout)    :: this
  type(betr_bounds_type)          , intent(in)    :: bounds
  integer                         , intent(in)    :: lbj, ubj
  type(BeTRtracer_type )          , intent(inout) :: betrtracer_vars
  character(len=*)                , intent(in)    :: namelist_buffer
  type(betr_status_type)          , intent(out)   :: bstatus

  !local variables
  character(len=*)       , parameter :: subname ='Init_betrbgc'
  integer :: itemp, itemp_trc

  integer :: dum
  integer :: itemp_grp, itemp_v, itemp_vgrp, itemp_adsgrp
  integer :: itemp_frz

  call bstatus%reset()
  this%parcol=1
  betrtracer_vars%is_tagged_h2o=.true.

  !volatile tracers
  itemp = 0; itemp_trc=0
  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_n2, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_n2,&
      trc_grp_end=betrtracer_vars%id_trc_end_n2, is_trc_gw=.true., is_trc_volatile = .true.)

  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_o2, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_o2, &
      trc_grp_end=betrtracer_vars%id_trc_end_o2, &
      is_trc_gw=.true., is_trc_volatile = .true.)

  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_ar,&
      trc_grp_beg= betrtracer_vars%id_trc_beg_ar, &
      trc_grp_end= betrtracer_vars%id_trc_end_ar, &
      is_trc_gw=.true., is_trc_volatile = .true.)

  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_co2x, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_co2x, &
      trc_grp_end=betrtracer_vars%id_trc_end_co2x, &
      is_trc_gw=.true., is_trc_volatile = .true.)

  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp= betrtracer_vars%id_trc_ch4, &
      trc_grp_beg= betrtracer_vars%id_trc_beg_ch4, &
      trc_grp_end= betrtracer_vars%id_trc_end_ch4, &
      is_trc_gw=.true., is_trc_volatile = .true.)

  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_blk_h2o, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_blk_h2o, &
      trc_grp_end=betrtracer_vars%id_trc_end_blk_h2o, &
      is_trc_gw=.true., is_trc_volatile = .true.)

  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_o18_h2o, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_o18_h2o, &
      trc_grp_end=betrtracer_vars%id_trc_end_o18_h2o, &
      is_trc_gw=.true., is_trc_volatile = .true.)

  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_d_h2o, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_d_h2o, &
      trc_grp_end=betrtracer_vars%id_trc_end_d_h2o, &
      is_trc_gw=.true., is_trc_volatile = .true.)

  betrtracer_vars%nmem_max               = 1

  call betrtracer_vars%Init()

  itemp_v = 0      !volatile id
  itemp_vgrp = 0   !volatile group
  itemp_frz = 0    !frozen tracer id
  call betrtracer_vars%set_tracer(bstatus=bstatus, trc_id = betrtracer_vars%id_trc_n2, trc_name='N2'  ,      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = betrtracer_vars%id_trc_n2,   &
       trc_group_mem= 1,  is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))
  if(bstatus%check_status())return

  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_o2, trc_name='O2'  ,      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = betrtracer_vars%id_trc_o2,   &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))
  if(bstatus%check_status())return

  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ar, trc_name='AR'  ,      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = betrtracer_vars%id_trc_ar,   &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))
  if(bstatus%check_status())return

  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_co2x, trc_name='CO2x',    &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = betrtracer_vars%id_trc_co2x  , &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &
       trc_volatile_group_id = addone(itemp_vgrp))
  if(bstatus%check_status())return

  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ch4, trc_name='CH4',      &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = betrtracer_vars%id_trc_ch4,   &
       trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
       trc_volatile_group_id = addone(itemp_vgrp))
  if(bstatus%check_status())return

  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_blk_h2o, trc_name='BLK_H2O' ,   &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = betrtracer_vars%id_trc_blk_h2o      ,   &
       trc_group_mem = 1, is_trc_diffusive =.false., is_trc_volatile=.true.                   ,   &
       trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp)          ,   &
       is_trc_h2o=.true., trc_vtrans_scal=1._r8, is_trc_frozen=.true.,  &
       trc_frozenid = addone(itemp_frz))
  if(bstatus%check_status())return

  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_o18_h2o, trc_name='O18_H2O' ,   &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = betrtracer_vars%id_trc_o18_h2o     ,   &
       trc_group_mem = 1, is_trc_diffusive =.false., is_trc_volatile=.true.                   ,   &
       trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp)          ,   &
       is_trc_h2o=.true., trc_vtrans_scal=1._r8, is_trc_frozen=.true.,  &
       trc_frozenid = addone(itemp_frz))
  if(bstatus%check_status())return

  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_d_h2o, trc_name='D_H2O' ,   &
       is_trc_mobile=.true., is_trc_advective = .true., trc_group_id =  betrtracer_vars%id_trc_d_h2o  ,   &
       trc_group_mem = 1, is_trc_diffusive =.false., is_trc_volatile=.true.                   ,   &
       trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp)          ,   &
       is_trc_h2o=.true., trc_vtrans_scal=1._r8, is_trc_frozen=.true., &
       trc_frozenid = addone(itemp_frz))
  if(bstatus%check_status())return

  !call betrtracer_vars%disp_betr_tracer()
  end subroutine Init_betrbgc


!-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, jtops, dz_top, betr_time, &
       betrtracer_vars, biophysforc, biogeo_flux, tracercoeff_vars, tracerboundarycond_vars, betr_status)
  !
  ! DESCRIPTION
  ! set up boundary conditions for tracer movement
  !
  ! USES
  use betr_ctrl              , only : iulog => biulog
  use TracerBoundaryCondType , only : tracerboundarycond_type
  use bshr_log_mod           , only : errMsg => shr_log_errMsg
  use BeTR_decompMod         , only : betr_bounds_type
  use BeTRTracerType         , only : betrtracer_type
  use betr_varcon            , only : denh2o  => bdenh2o
  use betr_varcon            , only : rgas => brgas
  use BeTR_biogeoFluxType    , only : betr_biogeo_flux_type
  use BetrStatusType         , only : betr_status_type
  use TracerCoeffType        , only : tracercoeff_type
  use UnitConvertMod         , only : ppm2molv
  use BeTR_TimeMod           , only : betr_time_type
  implicit none
  !ARGUMENTS
  class(bgc_reaction_h2oiso_type)  , intent(inout)    :: this
  type(betr_bounds_type)           , intent(in)    :: bounds                     !
  integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
  integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
  type(betrtracer_type)            , intent(in)    :: betrtracer_vars            !
  integer                          , intent(in)    :: jtops(bounds%begc: )
  real(r8)                         , intent(in)    :: dz_top(bounds%begc: )      !
  type(betr_time_type)             , intent(in)    :: betr_time
  type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
  type(betr_biogeo_flux_type)      , intent(in)    :: biogeo_flux
  type(tracercoeff_type)           , intent(in)    :: tracercoeff_vars
  type(tracerboundarycond_type)    , intent(inout) :: tracerboundarycond_vars !
  type(betr_status_type)           , intent(out)   :: betr_status
  !local variables
  integer :: fc, c, kk
  character(len=255) :: subname = 'set_boundary_conditions'
  real(r8) :: irt   !the inverse of R*T

  call betr_status%reset()
  SHR_ASSERT_ALL((ubound(dz_top)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__), betr_status)

  ! remove compiler warnings for unused dummy args
  if (this%dummy_compiler_warning) continue

  associate(                                                          &
    forc_pbot            => biophysforc%forc_pbot_downscaled_col    , &
    groupid              => betrtracer_vars%groupid                 , &
    ngwmobile_tracers    => betrtracer_vars%ngwmobile_tracers       , &
    is_volatile          => betrtracer_vars%is_volatile             , &
    volatilegroupid      => betrtracer_vars%volatilegroupid         , &
    snowres_col          => tracercoeff_vars%snowres_col            , &
    condc_toplay_col     => tracerboundarycond_vars%condc_toplay_col, &
    n2_ppmv              => biophysforc%n2_ppmv_col                 , &
    o2_ppmv              => biophysforc%o2_ppmv_col                 , &
    ar_ppmv              => biophysforc%ar_ppmv_col                 , &
    co2_ppmv             => biophysforc%co2_ppmv_col                , &
    ch4_ppmv             => biophysforc%ch4_ppmv_col                , &
    n2o_ppmv             => biophysforc%n2o_ppmv_col                , &
    nh3_ppmv             => biophysforc%nh3_ppmv_col                , &
    pbot_pa              => biophysforc%forc_pbot_downscaled_col    , &
    tair                 => biophysforc%forc_t_downscaled_col       , &
    diffblkm_topsoi_col  => tracercoeff_vars%diffblkm_topsoi_col    , &
    qflx_gross_evap_soil => biogeo_flux%qflx_gross_evap_soil_col      &
  )
  !eventually, the following code will be implemented using polymorphism
  !for simplicity, all gases other than water vapor are set with fixed concentration based boundary conditions
  !now the following gas composition does not make into 100% at the moment, it is about 99.15%
    !irt = 1._r8/(forc_tbot(c)*rgas)
  do fc = 1, num_soilc
    c = filter_soilc(fc)

   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2)   = ppm2molv(pbot_pa(c), n2_ppmv(c), tair(c))  !mol m-3, contant boundary condition, as concentration
   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o2)   = ppm2molv(pbot_pa(c), o2_ppmv(c), tair(c)) !mol m-3, contant boundary condition, as concentration
   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ar)   = ppm2molv(pbot_pa(c), ar_ppmv(c), tair(c)) !mol m-3, contant boundary condition, as concentration
   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_co2x) = ppm2molv(pbot_pa(c), co2_ppmv(c), tair(c))  !mol m-3, contant boundary condition, as concentration
   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ch4)  = ppm2molv(pbot_pa(c), ch4_ppmv(c), tair(c))  !mol m-3, contant boundary condition, as concentration

   !this is used for diffusion, however, because diffusion is off, these fluxes are imposed through reaction
   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_blk_h2o) = -qflx_gross_evap_soil(c)     !kg m-2-s, not diffusive water vapor transport
   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o18_h2o) = -qflx_gross_evap_soil(c)     !kg m-2-s, not diffusive water vapor transport
   tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_d_h2o)   = -qflx_gross_evap_soil(c)     !kg m-2-s, not diffusive water vapor transport

    tracerboundarycond_vars%bot_concflux_col(c,1,:)                         = 0._r8                       !zero flux boundary condition
    do kk = 1, ngwmobile_tracers
      if(.not. is_volatile(kk))cycle
         condc_toplay_col(c,groupid(kk))    = 1._r8/(snowres_col(c,volatilegroupid(kk))+&
           0.5_r8*dz_top(c)/diffblkm_topsoi_col(c,volatilegroupid(kk))) !m/s surface conductance
       enddo
  enddo
  end associate
  end subroutine set_boundary_conditions

!-------------------------------------------------------------------------------

  subroutine calc_bgc_reaction(this, bounds, col, lbj, ubj, num_soilc, filter_soilc,              &
       num_soilp,filter_soilp, jtops, betr_time, betrtracer_vars, tracercoeff_vars, biophysforc, &
       tracerstate_vars, tracerflux_vars, tracerboundarycond_vars, plant_soilbgc, &
       biogeo_flux, biogeo_state, betr_status)

  !
  ! do bgc reaction
  ! eventually this will be an abstract subroutine, but now I use the select case approach for a quick and dirty implementation.
  !USES
  !
  ! !USES:
  use BeTR_decompMod         , only : betr_bounds_type
  use tracerfluxType         , only : tracerflux_type
  use tracerstatetype        , only : tracerstate_type
  use tracercoeffType        , only : tracercoeff_type
  use BetrTracerType         , only : betrtracer_type
  use PlantSoilBGCMod        , only : plant_soilbgc_type
  use TracerBoundaryCondType , only : tracerboundarycond_type
  use BetrStatusType         , only : betr_status_type
  use betr_constants         , only : betr_errmsg_len
  use betr_columnType        , only : betr_column_type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  use BeTR_biogeoStateType     , only : betr_biogeo_state_type
  use BeTR_TimeMod             , only : betr_time_type
  use betr_varcon            , only : denh2o  => bdenh2o
  !ARGUMENTS
  class(bgc_reaction_h2oiso_type)  , intent(inout) :: this                       !
  type(betr_bounds_type)           , intent(in)    :: bounds ! bounds
  type(betr_column_type)           , intent(in)    :: col
  integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
  integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
  integer                          , intent(in)    :: num_soilp                  !
  integer                          , intent(in)    :: filter_soilp(:)            ! pft filter
  integer                          , intent(in)    :: jtops(bounds%begc: )       ! top index of each column
  integer                          , intent(in)    :: lbj, ubj                   ! lower and upper bounds, make sure they are > 0
  type(betr_time_type)             , intent(in)    :: betr_time
  type(betrtracer_type)            , intent(in)    :: betrtracer_vars            ! betr configuration information
  type(betr_biogeophys_input_type) , intent(inout) :: biophysforc
  type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars           !
  type(tracerstate_type)           , intent(inout) :: tracerstate_vars           !
  type(tracerflux_type)            , intent(inout) :: tracerflux_vars            !
  type(tracerboundarycond_type)    , intent(inout) :: tracerboundarycond_vars !
  class(plant_soilbgc_type)        , intent(inout) :: plant_soilbgc
  type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux
  type(betr_biogeo_state_type)     , intent(inout) :: biogeo_state
  type(betr_status_type)           , intent(out)   :: betr_status

  !local variables
  character(len=*)  , parameter :: subname ='calc_bgc_reaction'
  integer                       :: jj, c, fc, ll
  integer,            parameter :: nh2o_trcs=3
  integer                       :: jjs(nh2o_trcs), kk
  real(r8)                      :: tot0, tot1, frac1, tevap
  real(r8), parameter :: tiny_val=1.e-9_r8
  character(len=betr_errmsg_len) :: msg

  call betr_status%reset()

    associate(                                                                                &
    tracer_mobile_phase            => tracerstate_vars%tracer_conc_mobile_col               , &
    tracer_gwdif_concflux_top_col  => tracerboundarycond_vars%tracer_gwdif_concflux_top_col , &
    tracer_conc_frozen_vr          => tracerstate_vars%tracer_conc_frozen_col               , &
    frozenid                       => betrtracer_vars%frozenid                              , & !
    volatileid                     => betrtracer_vars%volatileid                            , & !
    tracer_flx_dif                 => tracerflux_vars%tracer_flx_dif_col                    , & !
    id_trc_blk_h2o                 => betrtracer_vars%id_trc_blk_h2o                        , &
    id_trc_o18_h2o                 => betrtracer_vars%id_trc_o18_h2o                        , &
    id_trc_d_h2o                   => betrtracer_vars%id_trc_d_h2o                          , &
    h2osoi_liq_vr                  => biophysforc%h2osoi_liq_col                              &
  )

  !apply the evaporation to the water tracer, the following is a hack to avoid the
  !inconsistency between water vapor transport in betr and the hydrology code
  !in the future, when the hdyrology code is corrected, the following will be gone, jyt Feb, 17, 2016
  jjs = (/id_trc_blk_h2o,id_trc_o18_h2o, id_trc_d_h2o/)

  do kk = 1, nh2o_trcs
    jj = jjs(kk)
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      !print*,'rac',tracer_mobile_phase(c,1,jj),tracer_gwdif_concflux_top_col(c,1,jj),h2osoi_liq_vr(c,1)/biophysforc%dz(c,1)
      tevap=tracer_gwdif_concflux_top_col(c,1,jj)*betr_time%delta_time/biophysforc%dz(c,1)
      if(tracer_conc_frozen_vr(c,1,frozenid(jj))<tiny_val)then
        tracer_mobile_phase(c,1,jj) = tracer_mobile_phase(c,1,jj) + tevap
      else
        frac1=tracer_mobile_phase(c,1,jj)/(tracer_mobile_phase(c,1,jj)+tracer_conc_frozen_vr(c,1,frozenid(jj)))
        tracer_mobile_phase(c,1,jj) = tracer_mobile_phase(c,1,jj) + tevap * frac1
        tracer_conc_frozen_vr(c,1,frozenid(jj))=tracer_conc_frozen_vr(c,1,frozenid(jj)) + tevap*(1._r8-frac1)
      endif
      if(tracer_mobile_phase(c,1,jj) < 0._r8)then
        do ll = 1, 3
          tot0 = tracer_mobile_phase(c,ll,jj)*biophysforc%dz(c,ll)
          tot1 = tracer_mobile_phase(c,ll+1,jj)*biophysforc%dz(c,ll+1)
          tot1 = tot1 + tot0
          tracer_mobile_phase(c,ll,jj) = 0._r8
          tracer_mobile_phase(c,ll+1,jj) = tot1/biophysforc%dz(c,ll+1)
          if(tot1>0._r8)exit
        enddo
        !the following should rarely occur, so when it occur, end with a warning
        if(tot1<0._r8)then
          write(msg,*)tracer_mobile_phase(c,1:2,jj),tot1
          msg=trim(msg)//new_line('A')//'negative H2O tracer '//errMsg(mod_filename, __LINE__)
          call betr_status%set_msg(msg=msg, err=-1)
        endif
      endif
      tracer_flx_dif(c,jj) = tracer_flx_dif(c,jj)- &
        tracer_gwdif_concflux_top_col(c,1,jj) * betr_time%delta_time
    enddo
  enddo

  end associate
  end subroutine calc_bgc_reaction

!-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars, betr_status)
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
  use BetrStatusType        , only : betr_status_type
  implicit none
  !ARGUMENTS
  class(bgc_reaction_h2oiso_type) , intent(inout) :: this
  type(betr_bounds_type)          , intent(in) :: bounds
  integer                         , intent(in) :: lbj, ubj
  integer                         , intent(in) :: jtops(bounds%begc: )        ! top label of each column
  integer                         , intent(in) :: num_soilc
  integer                         , intent(in) :: filter_soilc(:)
  type(betrtracer_type)           , intent(in) :: betrtracer_vars
  type(tracercoeff_type)          , intent(in) :: tracercoeff_vars
  type(tracerstate_type)          , intent(inout) :: tracerstate_vars
  type(betr_status_type)          , intent(out)   :: betr_status

  !local variables
  character(len=255) :: subname = 'do_tracer_equilibration'
  integer   :: j, fc, c
  integer   :: trc_id1, trc_id2

  call betr_status%reset()
  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__), betr_status)

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)                          continue
    if (bounds%begc > 0)                                      continue
    if (lbj > 0)                                              continue
    if (ubj > 0)                                              continue
    if (size(jtops) > 0)                                      continue
    if (num_soilc > 0)                                        continue
    if (size(filter_soilc) > 0)                               continue
    if (size(tracerstate_vars%tracer_conc_surfwater_col) > 0) continue
!    if (size(tracercoeff_vars%annsum_counter_col) > 0)        continue

  associate(                                                                         &
    aqu2equilscef                   => tracercoeff_vars%aqu2equilsolidcef_col      , &
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

  if (.false.) then
  !  if(trc_id1>0)then
    call do_h2o_isotope_equilibration(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
      aqu2bulkcef_mobile(bounds%begc:bounds%endc, lbj:ubj, trc_id1)        ,            &
      aqu2equilscef(bounds%begc:bounds%endc, lbj:ubj, trc_id2)             ,            &
      tracer_solid_phase_equil(bounds%begc:bounds%endc, lbj:ubj, trc_id2)  ,            &
      tracer_mobile_phase(bounds%begc:bounds%endc, lbj:ubj, trc_id1), betr_status)

  endif
  end associate
  end subroutine do_tracer_equilibration


!-------------------------------------------------------------------------------

  subroutine do_h2o_isotope_equilibration(bounds, lbj, ubj, jtops, numf, filter, aqu2bulkcef, &
    aqu2equilscef, tracer_solid_phase_equil, tracer_mobile_phase, bstatus)
  !
  ! Diagnose solid phase tracer
  !
  use BeTR_decompMod        , only : betr_bounds_type
  use BetrStatusType        , only : betr_status_type
  implicit none
  type(betr_bounds_type) , intent(in)    :: bounds
  integer                , intent(in)    :: lbj, ubj
  integer                , intent(in)    :: jtops(bounds%begc: )        ! top label of each column
  integer                , intent(in)    :: numf
  integer                , intent(in)    :: filter(:)
  real(r8)               , intent(in)    :: aqu2equilscef(bounds%begc: , lbj: )
  real(r8)               , intent(in)    :: aqu2bulkcef(bounds%begc: , lbj: )
  real(r8)               , intent(in)    :: tracer_mobile_phase(bounds%begc: , lbj: )
  real(r8)               , intent(inout) :: tracer_solid_phase_equil(bounds%begc: ,lbj: )
  type(betr_status_type) , intent(out)   :: bstatus
  real(r8)  :: frac
  real(r8)  :: tracer_conc
  integer   :: c, fc, j

  call bstatus%reset()
  SHR_ASSERT_ALL((ubound(aqu2equilscef,1)            == bounds%endc), errMsg(mod_filename,__LINE__),bstatus)

  SHR_ASSERT_ALL((ubound(aqu2equilscef,2)            == ubj),         errMsg(mod_filename,__LINE__),bstatus)

  SHR_ASSERT_ALL((ubound(aqu2bulkcef,1)              == bounds%endc), errMsg(mod_filename,__LINE__),bstatus)

  SHR_ASSERT_ALL((ubound(aqu2bulkcef,2)              == ubj),         errMsg(mod_filename,__LINE__),bstatus)

  SHR_ASSERT_ALL((ubound(tracer_solid_phase_equil,1) == bounds%endc), errMsg(mod_filename,__LINE__),bstatus)

  SHR_ASSERT_ALL((ubound(tracer_solid_phase_equil,2) == ubj),         errMsg(mod_filename,__LINE__),bstatus)

  SHR_ASSERT_ALL((ubound(tracer_mobile_phase,1)      == bounds%endc), errMsg(mod_filename,__LINE__),bstatus)

  SHR_ASSERT_ALL((ubound(tracer_mobile_phase,2)      == ubj),         errMsg(mod_filename,__LINE__),bstatus)

  ! remove compiler warnings for unused dummy args
  if (bounds%begc > 0) continue
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
  subroutine InitCold(this, bounds, col, betrtracer_vars, biophysforc, tracerstate_vars)
    !
    ! !USES:
    !
    use BeTR_decompMod    , only : betr_bounds_type
    use BeTRTracerType    , only : BeTRTracer_Type
    use tracerstatetype   , only : tracerstate_type
    use betr_varcon       , only : spval => bspval, ispval => bispval
    use betr_varcon       , only : denh2o => bdenh2o
    use betr_columnType   , only : betr_column_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_h2oiso_type)  , intent(inout)    :: this
    type(betr_bounds_type)           , intent(in)    :: bounds
    type(betr_column_type)           , intent(in)    :: col
    type(BeTRTracer_Type)            , intent(in)    :: betrtracer_vars
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, k, j
    integer :: fc                                        ! filter index
    integer               :: begc, endc
    integer               :: begg, endg
    integer               :: trcid
    !-----------------------------------------------------------------------

    associate(                                                                   &
      frozenid                  => betrtracer_vars%frozenid                    , &
      nsolid_equil_tracers      => betrtracer_vars%nsolid_equil_tracers        , &
      tracer_conc_mobile_vr     => tracerstate_vars%tracer_conc_mobile_col     , &
      tracer_conc_solid_equil_vr=> tracerstate_vars%tracer_conc_solid_equil_col, &
      tracer_conc_frozen_vr     => tracerstate_vars%tracer_conc_frozen_col     , &
      tracer_conc_surfwater_col => tracerstate_vars%tracer_conc_surfwater_col  , &
      tracer_conc_aquifer_col   => tracerstate_vars%tracer_conc_aquifer_col    , &
      tracer_conc_grndwater_col => tracerstate_vars%tracer_conc_grndwater_col  , &
      tracer_soi_molarmass_col  => tracerstate_vars%tracer_soi_molarmass_col   , &
      h2osoi_liq_col            => biophysforc%h2osoi_liq_col                  , &
      dz                        => biophysforc%dz                              , &
      h2osoi_ice_col            => biophysforc%h2osoi_ice_col                    &
    )
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !-----------------------------------------------------------------------
    do c = bounds%begc, bounds%endc

      !dual phase tracers
      tracer_conc_mobile_vr(c,:, :)          = 0._r8
      tracer_conc_surfwater_col(c,:)          = 0._r8
      tracer_conc_aquifer_col(c,:)            = 0._r8
      tracer_conc_grndwater_col(c,:)          = 0._r8

      if(nsolid_equil_tracers>0)then
        tracer_conc_solid_equil_vr(c, :, :) = 0._r8
      endif
      tracer_soi_molarmass_col(c,:)          = 0._r8
      !set for o18_h2o, assuming no fractionation, which is equivalent to assuming concentration equals 1
      trcid = betrtracer_vars%id_trc_blk_h2o
      tracer_conc_grndwater_col(c,trcid) = denh2o

      do j = 1, bounds%ubj
        tracer_conc_mobile_vr(c,j,trcid) = 1._r8 * h2osoi_liq_col(c,j)/dz(c,j)
        tracer_conc_frozen_vr(c,j,frozenid(trcid)) = 1._r8 * h2osoi_ice_col(c,j)/dz(c,j)
      enddo

      trcid = betrtracer_vars%id_trc_o18_h2o
      tracer_conc_grndwater_col(c,trcid) = denh2o
      do j = 1, bounds%ubj
        tracer_conc_mobile_vr(c,j,trcid) = 1._r8 * h2osoi_liq_col(c,j)/dz(c,j)
        tracer_conc_frozen_vr(c,j,frozenid(trcid)) = 1._r8 * h2osoi_ice_col(c,j)/dz(c,j)
      enddo
      trcid = betrtracer_vars%id_trc_d_h2o
      tracer_conc_grndwater_col(c,trcid) = denh2o
      do j = 1, bounds%ubj
        tracer_conc_mobile_vr(c,j,trcid) = 1._r8 * h2osoi_liq_col(c,j)/dz(c,j)
        tracer_conc_frozen_vr(c,j,frozenid(trcid)) = 1._r8 * h2osoi_ice_col(c,j)/dz(c,j)
      enddo
   enddo
   end associate
  end subroutine InitCold



  !-----------------------------------------------------------------------
  subroutine readParams(this, name_list_buffer, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! read in module specific parameters
    !
    ! !USES:
    use BeTRTracerType , only : BeTRTracer_Type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_h2oiso_type) , intent(inout)    :: this
    type(BeTRTracer_Type)           , intent(inout) :: betrtracer_vars
    character(len=*)                  , intent(in)  :: name_list_buffer

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning) continue
    if (len(betrtracer_vars%betr_simname) > 0) continue
    !do nothing here for the moment, but contents will eventually be filled in here

  end subroutine readParams

  !-------------------------------------------------------------------------------
  subroutine retrieve_biogeoflux(this, num_soilc, filter_soilc, tracerflux_vars, &
  betrtracer_vars, biogeo_flux)

  use tracerfluxType           , only : tracerflux_type
  use BeTR_decompMod           , only : betr_bounds_type
  use BeTRTracerType           , only : BeTRTracer_Type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  implicit none
   class(bgc_reaction_h2oiso_type) , intent(inout) :: this     !!
  integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
  integer                          , intent(in)    :: filter_soilc(:)             ! column filter
  type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
  type(tracerflux_type)            , intent(in)    :: tracerflux_vars
  type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)                          continue
    if (num_soilc > 0)                                        continue
    if (size(filter_soilc) > 0)                               continue
    if (len(betrtracer_vars%betr_simname) > 0)                continue
    if (size(tracerflux_vars%tracer_flx_top_soil_col) > 0)    continue

  end subroutine retrieve_biogeoflux
   !-------------------------------------------------------------------------------
  subroutine debug_info(this, bounds, num_soilc, filter_soilc, dzsoi, betrtracer_vars, tracerstate_vars, header, betr_status)

   use BeTRTracerType           , only : BeTRTracer_Type
   use tracerstatetype          , only : tracerstate_type
   use BeTR_decompMod           , only : betr_bounds_type
     ! !ARGUMENTS:
    implicit none
   class(bgc_reaction_h2oiso_type) , intent(inout) :: this     !
   type(betr_bounds_type)               , intent(in) :: bounds                      ! bounds
   integer                              , intent(in) :: num_soilc                   ! number of columns in column filter
   integer                              , intent(in) :: filter_soilc(:)             ! column filter
   real(r8)                             , intent(in) :: dzsoi(bounds%begc: ,bounds%lbj: )
   type(betrtracer_type)                , intent(in) :: betrtracer_vars             ! betr configuration information
   type(tracerstate_type)               , intent(in) :: tracerstate_vars
   character(len=*)                     , intent(in) :: header
   type(betr_status_type)               , intent(out):: betr_status
   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(dzsoi)  == (/bounds%endc, bounds%ubj/)),   errMsg(mod_filename,__LINE__), betr_status)

   if (this%dummy_compiler_warning) continue
     end subroutine debug_info
   !----------------------------------------------------------------------
   subroutine retrieve_biostates(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
     betrtracer_vars, tracerstate_vars, biogeo_state, betr_status)
   !
   !retrieve state variables for lsm mass balance check
   use tracer_varcon, only : catomw, natomw, patomw, c13atomw, c14atomw
   use BeTRTracerType           , only : BeTRTracer_Type
   use tracerstatetype          , only : tracerstate_type
   use BeTR_biogeoStateType     , only : betr_biogeo_state_type
   use BeTR_decompMod           , only : betr_bounds_type
   implicit none
   class(bgc_reaction_h2oiso_type) , intent(inout) :: this               !
   type(betr_bounds_type)               , intent(in)  :: bounds                      ! bounds
   integer                              , intent(in) :: lbj, ubj
   integer                              , intent(in) :: jtops(bounds%begc: )
   integer                              , intent(in)    :: num_soilc                   ! number of columns in column filter
   integer                              , intent(in)    :: filter_soilc(:)             ! column filter
   type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
   type(tracerstate_type)               , intent(inout) :: tracerstate_vars
   type(betr_biogeo_state_type)         , intent(inout) :: biogeo_state
   type(betr_status_type)               , intent(out)   :: betr_status

   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__), betr_status)

   if (this%dummy_compiler_warning) continue
   if (bounds%begc > 0)             continue

   end subroutine retrieve_biostates



   !----------------------------------------------------------------------
   subroutine reset_biostates(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, biophysforc,  tracerstate_vars, betr_status)

       ! !USES:
       use bshr_kind_mod            , only : r8 => shr_kind_r8
       use tracerstatetype          , only : tracerstate_type
       use BeTR_decompMod           , only : betr_bounds_type
       use BeTRTracerType           , only : BeTRTracer_Type
       use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
       use BetrStatusType           , only : betr_status_type
       use betr_columnType          , only : betr_column_type

       ! !ARGUMENTS:
     implicit none
       class(bgc_reaction_h2oiso_type) , intent(inout) :: this               !
       type(betr_bounds_type)           , intent(in)    :: bounds                      ! bounds
       integer                          , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
       integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
       integer                          , intent(in)    :: filter_soilc(:)             ! column filter
       integer                          , intent(in)    :: jtops( : )                  ! top index of each column
       type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
       type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
       type(tracerstate_type)           , intent(inout) :: tracerstate_vars
       type(betr_status_type)           , intent(out)   :: betr_status

   end subroutine reset_biostates
   !----------------------------------------------------------------------
   subroutine SetParCols(this, parcol)
     implicit none
       class(bgc_reaction_h2oiso_type) , intent(inout) :: this            !
     integer, intent(in) :: parcol

     this%parcol = parcol
   end subroutine SetParCols
end module H2OIsotopeBGCReactionsType
