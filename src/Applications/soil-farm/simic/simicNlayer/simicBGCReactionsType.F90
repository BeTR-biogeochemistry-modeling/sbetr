module SimicBGCReactionsType

#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  ! This is an example on how to use polymorphism to create your own bgc modules that will be run with BeTR
  !
  ! HISTORY:
  ! Created by Jinyun Tang, Oct 2nd, 2014
  ! !USES:
  use bshr_assert_mod, only : shr_assert
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use bshr_assert_mod, only : shr_assert_any
  use bshr_log_mod             , only : errMsg => shr_log_errMsg
  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use BeTR_decompMod           , only : bounds_type  => betr_bounds_type
  use BGCReactionsMod          , only : bgc_reaction_type
  use tracer_varcon            , only : bndcond_as_conc, bndcond_as_flux
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BetrStatusType           , only : betr_status_type
  use simicParaType            , only : simic_para
  use simicBGCIndexType        , only : simic_bgc_index_type
  use simicBGCType             , only : simic_bgc_type
  use JarBgcForcType           , only : JarBGC_forc_type
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: simic_bgc_reaction_type

  type, extends(bgc_reaction_type) :: &
       simic_bgc_reaction_type
  private
    integer :: parcol

    type(simic_bgc_type), pointer :: simic_bgc(:,:)
    type(JarBGC_forc_type), pointer :: simic_forc(:,:)
    type(simic_bgc_index_type) :: simic_bgc_index
    logical :: use_c13
    logical :: use_c14
    logical :: nop_limit
    logical :: non_limit
    integer :: nactpft               ! number of active pfts
  contains
     procedure :: Init_betrbgc                          ! initialize betr bgc
     procedure :: set_boundary_conditions               ! set top/bottom boundary conditions for various tracers
     procedure :: calc_bgc_reaction                     ! doing bgc calculation
     procedure :: init_boundary_condition_type          ! initialize type of top boundary conditions
     procedure :: do_tracer_equilibration               ! do equilibrium tracer chemistry
     procedure :: InitCold                              ! do cold initialization
     procedure :: retrieve_biogeoflux           !
     procedure :: set_kinetics_par
     procedure :: retrieve_lnd2atm
     procedure :: readParams                   ! read in parameters
     procedure :: retrieve_biostates
     procedure :: debug_info
     procedure :: set_bgc_spinup
     procedure :: UpdateParas
     procedure :: init_iP_prof
     procedure :: reset_biostates
     procedure :: SetParCols
     procedure, private :: set_tracer
     procedure, private :: InitAllocate
     procedure, private :: retrieve_output
     procedure, private :: set_bgc_forc
     procedure, private :: update_sorpphase_coeff
  end type simic_bgc_reaction_type

  interface simic_bgc_reaction_type
     module procedure constructor
  end interface simic_bgc_reaction_type

contains
  !-------------------------------------------------------------------------------
  type(simic_bgc_reaction_type) function constructor()
  !
  ! !DESCRIPTION:
  ! create an object of type simic_bgc_reaction_type.
  ! Right now it is purposely empty
    type(simic_bgc_reaction_type), allocatable :: bgc
    allocate(bgc)
    constructor = bgc
  end function constructor

  !-------------------------------------------------------------------------------
  subroutine UpdateParas(this, bounds, lbj, ubj, bstatus)
  implicit none
  class(simic_bgc_reaction_type)         , intent(inout)    :: this
  type(bounds_type)                    , intent(in)    :: bounds
  integer                              , intent(in)    :: lbj, ubj        ! lower and upper bounds, make sure they are > 0
  type(betr_status_type)               , intent(out)   :: bstatus
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
  use BeTRTracerType        , only : betrtracer_type
  implicit none
  ! !ARGUMENTS:
  class(simic_bgc_reaction_type)         , intent(inout)    :: this
  type(bounds_type)                        , intent(in) :: bounds
  integer                                  , intent(in) :: lbj, ubj
  type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc
  type(BeTRtracer_type)                    , intent(inout) :: tracers
  type(tracerstate_type)                   , intent(inout) :: tracerstate_vars


  if (this%dummy_compiler_warning) continue
  if (bounds%begc > 0) continue

  end subroutine init_iP_prof
  !----------------------------------------------------------------------
  subroutine set_kinetics_par(this, lbj, ubj, nactpft, plantNutkinetics, tracers, tracercoeff_vars)
  use PlantNutKineticsMod, only : PlantNutKinetics_type
  use tracercoeffType          , only : tracercoeff_type
  use BeTRTracerType           , only : betrtracer_type
  ! !ARGUMENTS:
  class(simic_bgc_reaction_type)         , intent(inout)    :: this
  class(PlantNutKinetics_type), intent(in) :: plantNutkinetics
  type(betrtracer_type)       , intent(in) :: tracers
  type(tracercoeff_type), intent(inout) :: tracercoeff_vars
  integer, intent(in) :: lbj, ubj
  integer, intent(in) :: nactpft

  integer :: c_l, p, j
  !in the following, only one column is assumed for the bgc
  c_l = 1
  this%nactpft = nactpft
  do j = lbj, ubj
    !effective p competing decomposers
    this%simic_forc(c_l,j)%Msurf_OM  = plantNutkinetics%minsurf_dom_compet_vr_col(c_l,j)
    this%simic_forc(c_l,j)%KM_OM_ref = plantNutkinetics%km_minsurf_dom_vr_col(c_l,j)
  enddo

  end subroutine set_kinetics_par
  !-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
    !
    ! !DESCRIPTION:
    ! initialize boundary condition types
    !
    ! !USES:
    use BeTRTracerType        , only : betrtracer_type
    use TracerBoundaryCondType, only : tracerboundarycond_type
    use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
    use BeTRTracerType        , only : betrtracer_type

    ! !ARGUMENTS:
    class(simic_bgc_reaction_type), intent(inout) :: this
    type(BeTRtracer_type),             intent(in) :: betrtracer_vars
    type(bounds_type),                 intent(in) :: bounds
    type(tracerboundarycond_type),     intent(in) :: tracerboundarycond_vars

    ! !LOCAL VARIABLES:

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning) continue
    if (bounds%begc > 0) continue
    if (len(betrtracer_vars%betr_simname) > 0) continue
    tracerboundarycond_vars%topbc_type(:) = bndcond_as_conc

    !when bottom BC is not given, it is specified as constant flux
    ! FIXME(bja, 201604) Don't we need a bottom BC?
    !X!tracerboundarycond_vars%botbc_type(:) = bndcond_as_flux

  end subroutine init_boundary_condition_type
  !-------------------------------------------------------------------------------
  subroutine set_bgc_spinup(this, bounds, lbj, ubj, biophysforc, &
  tracers, tracerstate_vars)

  use tracerstatetype        , only : tracerstate_type
  use BeTRTracerType         , only : betrtracer_type
  use BeTR_decompMod         , only : betr_bounds_type

  implicit none
    class(simic_bgc_reaction_type), intent(inout)    :: this
    type(betr_bounds_type)                       , intent(in) :: bounds
    integer                                 , intent(in) :: lbj, ubj
    type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc
    type(BeTRtracer_type)                   , intent(inout) :: tracers
    type(tracerstate_type)                  , intent(inout) :: tracerstate_vars


    if (this%dummy_compiler_warning) continue
    if (bounds%begc > 0) continue

  end subroutine set_bgc_spinup
  !-------------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, lbj, ubj, bstatus)

  use betr_varcon                      , only : betr_maxpatch_pft
  implicit none
  class(simic_bgc_reaction_type), intent(inout)    :: this
  type(betr_status_type)           , intent(out)   :: bstatus
  type(bounds_type)                , intent(in)    :: bounds
  integer                          , intent(in)    :: lbj, ubj

  integer :: j, c
  logical :: batch_mode

  batch_mode =.false.
  this%nactpft = 0

  call this%simic_bgc_index%Init(simic_para%use_c13, simic_para%use_c14, &
     simic_para%non_limit, simic_para%nop_limit, betr_maxpatch_pft)

  if(bstatus%check_status())return

  !create the models
  allocate(this%simic_bgc(bounds%begc:bounds%endc,lbj:ubj))

  !create model specific forcing data structure
  allocate(this%simic_forc(bounds%begc:bounds%endc,lbj:ubj))

  !initialize
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      call this%simic_bgc(c,j)%Init(simic_para, batch_mode, bstatus)
      if(bstatus%check_status())return
        call this%simic_forc(c,j)%Init(this%simic_bgc_index%nstvars)
    enddo
  enddo
  end subroutine InitAllocate
  !-------------------------------------------------------------------------------
  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars, namelist_buffer, bstatus)
    !
    ! DESCRIPTION:
    ! initialize the betrbgc
    !
    ! !USES:
    use BeTRTracerType , only : betrtracer_type
    use BetrStatusType , only : betr_status_type
    use gbetrType      , only : gbetr_type
    implicit none
    ! !ARGUMENTS:
    class(simic_bgc_reaction_type), intent(inout)    :: this
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: lbj, ubj
    type(BeTRtracer_type )           , intent(inout) :: betrtracer_vars
    character(len=*)                 , intent(in)    :: namelist_buffer
    type(betr_status_type)           , intent(out)   :: bstatus
    character(len=*), parameter                      :: subname ='Init_betrbgc'


    call bstatus%reset()
    this%parcol=1
    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)           continue
    if (bounds%begc > 0)                       continue
    if (ubj > lbj)                             continue
    if (len(betrtracer_vars%betr_simname) > 0) continue

    call this%InitAllocate(bounds, lbj, ubj, bstatus)
    if(bstatus%check_status())return

    this%use_c13 = simic_para%use_c13
    this%use_c14 = simic_para%use_c14
    this%nop_limit=simic_para%nop_limit
    this%non_limit=simic_para%non_limit

    call this%set_tracer(betrtracer_vars, bstatus)

  end subroutine Init_betrbgc
  !-------------------------------------------------------------------------------
  subroutine set_tracer(this, betrtracer_vars, bstatus)

  use BeTRTracerType  , only : betrtracer_type
  use MathfuncMod     , only : addone
  implicit none
  class(simic_bgc_reaction_type), intent(inout)    :: this
  type(BeTRtracer_type )           , intent(inout) :: betrtracer_vars
  type(betr_status_type)           , intent(out)   :: bstatus

  integer :: itemp_gwm
  integer :: itemp_g
  integer :: itemp_s
  integer :: itemp_gwm_grp
  integer :: dum, itemp
  integer :: itemp_grp, itemp_v, itemp_vgrp, itemp_trc, itemp_ads, itemp_ads_grp
  integer :: litr_cnt, wood_cnt, Bm_cnt, trcid, itemp_mem, ngroupmems


    associate(                           &
    nelm    => this%simic_bgc_index%nelms,   &
    c_loc   => this%simic_bgc_index%c_loc,   &
    c13_loc => this%simic_bgc_index%c13_loc, &
    c14_loc => this%simic_bgc_index%c14_loc, &
    e_loc   => this%simic_bgc_index%e_loc    &
    )
    itemp_gwm     = 0;
    itemp_g       = 0 ;
    itemp_s       = 0;
    itemp_gwm_grp = 0
    itemp_ads_grp =0!counter for sorptive groups
    itemp_ads=0     !counter for sorptive tracers

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

    if(this%use_c13)then
      call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
        trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_c13_co2x, &
        trc_grp_beg=betrtracer_vars%id_trc_beg_c13_co2x, &
        trc_grp_end=betrtracer_vars%id_trc_end_c13_co2x, &
        is_trc_gw=.true., is_trc_volatile = .true.)
    endif
    if(this%use_c14)then
      call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
        trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_c14_co2x, &
        trc_grp_beg=betrtracer_vars%id_trc_beg_c14_co2x, &
        trc_grp_end=betrtracer_vars%id_trc_end_c14_co2x, &
        is_trc_gw=.true., is_trc_volatile = .true.)
    endif

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp= betrtracer_vars%id_trc_ch4, &
      trc_grp_beg= betrtracer_vars%id_trc_beg_ch4, &
      trc_grp_end= betrtracer_vars%id_trc_end_ch4, &
      is_trc_gw=.true., is_trc_volatile = .true.)

    ngroupmems=nelm+1  !dom, element + energy
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_dom, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_dom, &
      trc_grp_end=betrtracer_vars%id_trc_end_dom, &
      is_trc_gw=.true., is_trc_volatile = .false.)

    ngroupmems=nelm+1  !pom, element + energy
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_pom, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_pom, &
      trc_grp_end=betrtracer_vars%id_trc_end_pom, &
      is_trc_passive=.true.)

    !three litter groups

    ngroupmems = 3*nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_litr, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_litr, &
      trc_grp_end=betrtracer_vars%id_trc_end_litr, &
      is_trc_passive=.true.)

    ngroupmems = nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), &
      is_trc_passive=.true., mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_wood, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_wood, &
      trc_grp_end=betrtracer_vars%id_trc_end_wood)

    ngroupmems = 2*nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_Bm, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_Bm, &
      trc_grp_end=betrtracer_vars%id_trc_end_Bm, &
      is_trc_passive=.true.)

    betrtracer_vars%nmem_max               = nelm*3

    call betrtracer_vars%Init()

    itemp_grp = 0    !group id
    itemp_v = 0      !volatile id
    itemp_vgrp = 0   !volatile group

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_n2, trc_name='N2'  ,      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem= 1,  is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_o2, trc_name='O2'  ,      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ar, trc_name='AR'  ,      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_co2x, trc_name='CO2x',    &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp)  , &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ch4, trc_name='CH4',      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    if(this%use_c13)then
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_c13_co2x, &
         trc_name='13CO2x', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_c13_co2x, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp), &
         trc_family_name='CO2x')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_c14_co2x, &
         trc_name='14CO2x', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_c14_co2x, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp), &
         trc_family_name='CO2x')
      if(bstatus%check_status())return
    endif

    itemp_mem=0
    litr_cnt = 0
    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C_C13' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C_C14' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
      if(bstatus%check_status())return
    endif
    litr_cnt = litr_cnt + 1

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C'  ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C_C13' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C_C14' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
      if(bstatus%check_status())return
    endif
    litr_cnt=litr_cnt+1

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C_C13' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C_C14' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
      if(bstatus%check_status())return
    endif

    wood_cnt=0
    itemp_mem=0
    !coarse root woody components, equivalent to default cwd
    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDC' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='CWD')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDC_C13' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='CWD')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDC_C14' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='CWD')
      if(bstatus%check_status())return
    endif

    Bm_cnt=0;itemp_mem=0
    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_live', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='MB')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_live_C13',&
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='MB')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_live_C14', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm,  trc_group_mem = addone(itemp_mem), &
         trc_family_name='MB')
      if(bstatus%check_status())return
    endif
    Bm_cnt=Bm_cnt+1
    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_dead', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='MB')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_dead_C13',&
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='MB')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_dead_C14', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm,  trc_group_mem = addone(itemp_mem), &
         trc_family_name='MB')
      if(bstatus%check_status())return
    endif

    itemp_mem=0
    trcid =  betrtracer_vars%id_trc_beg_dom+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                       , &
         trc_name='DOM_C', is_trc_mobile=.true., is_trc_advective = .true. , &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &
         trc_adsorbgroupid=addone(itemp_ads_grp), trc_sorpisotherm='LANGMUIR'            , &
         is_trc_dom=.true.,trc_family_name='DOM')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_dom+e_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                       , &
         trc_name='DOM_e', is_trc_mobile=.true., is_trc_advective = .true.  , &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &
         trc_adsorbgroupid=itemp_ads_grp, trc_sorpisotherm='LANGMUIR'                    , &
         is_trc_dom=.true., trc_family_name='DOM')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_dom+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                     , &
         trc_name='DOM_C13', is_trc_mobile=.true., is_trc_advective = .true. , &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &
         trc_adsorbgroupid=itemp_ads_grp,trc_sorpisotherm='LANGMUIR'                     , &
         is_trc_dom=.true., trc_family_name='DOM')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid=betrtracer_vars%id_trc_beg_dom+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                     , &
         trc_name='DOM_C14', is_trc_mobile=.true., is_trc_advective = .true. , &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &
         trc_adsorbgroupid=itemp_ads_grp, trc_sorpisotherm='LANGMUIR'                    , &
         is_trc_dom=.true., trc_family_name='DOM')
      if(bstatus%check_status())return
    endif

    itemp_mem=0
    trcid =  betrtracer_vars%id_trc_beg_pom+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_C' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &
         trc_family_name='POM')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_pom+e_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_e' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &
         trc_family_name='POM')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_pom+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_C13' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &
         trc_family_name='POM')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid=betrtracer_vars%id_trc_beg_pom+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_C14' ,  &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &
         trc_family_name='POM')
      if(bstatus%check_status())return
    endif

  end associate
  end subroutine set_tracer
  !-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, jtops, dz_top, betr_time, &
       betrtracer_vars, biophysforc, biogeo_flux, tracercoeff_vars, tracerboundarycond_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! set up boundary conditions for tracer movement
    !
    ! !USES:
    use betr_ctrl              , only : iulog  => biulog
    use TracerBoundaryCondType , only : tracerboundarycond_type
    use bshr_log_mod           , only : errMsg => shr_log_errMsg
    use BeTRTracerType         , only : betrtracer_type
    use betr_varcon            , only : rgas => brgas
    use BeTR_biogeoFluxType    , only : betr_biogeo_flux_type
    use BetrStatusType         , only : betr_status_type
    use UnitConvertMod         , only : ppm2molv
    use TracerCoeffType        , only : tracercoeff_type
    use BeTR_TimeMod             , only : betr_time_type
    implicit none
    ! !ARGUMENTS:
    class(simic_bgc_reaction_type) , intent(inout)    :: this                       !
    type(bounds_type)                 , intent(in)    :: bounds                     !
    integer                           , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                           , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars            !
    real(r8)                          , intent(in)    :: dz_top(bounds%begc: )      !
    integer                           , intent(in)    :: jtops(bounds%begc: )
    type(betr_time_type)              , intent(in)    :: betr_time
    type(betr_biogeophys_input_type)  , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)       , intent(in)    :: biogeo_flux
    type(tracercoeff_type)             , intent(in)   :: tracercoeff_vars
    type(tracerboundarycond_type)     , intent(inout) :: tracerboundarycond_vars !
    type(betr_status_type)            , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    integer            :: fc, c, kk
    character(len=255) :: subname = 'set_boundary_conditions'
    real(r8) :: irt   !the inverse of R*T

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(dz_top)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__), betr_status)

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)      continue
    if (size(biogeo_flux%qflx_adv_col)>0) continue
    associate(                                                           &
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
       pbot_pa              => biophysforc%forc_pbot_downscaled_col    , &
       tair                 => biophysforc%forc_t_downscaled_col       , &
       id_trc_beg_dom       => betrtracer_vars%id_trc_beg_dom          , &
       id_trc_end_dom       => betrtracer_vars%id_trc_end_dom          , &
       diffblkm_topsoi_col  => tracercoeff_vars%diffblkm_topsoi_col      &
      )

      do fc = 1, num_soilc
         c = filter_soilc(fc)
         !eventually, the following code will be implemented using polymorphism
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,:) = 0._r8
         tracerboundarycond_vars%bot_concflux_col(c,1,:)  = 0._r8                        !zero flux boundary condition for diffusion
         tracerboundarycond_vars%condc_toplay_col(c,:)    = 0._r8
         !set value for specific groups
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2)    =ppm2molv(pbot_pa(c), n2_ppmv(c), tair(c))    !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o2)    =ppm2molv(pbot_pa(c), o2_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ar)    =ppm2molv(pbot_pa(c), ar_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_co2x)  =ppm2molv(pbot_pa(c), co2_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ch4)   =ppm2molv(pbot_pa(c), ch4_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,id_trc_beg_dom:id_trc_end_dom)= 0._r8            !mol m-3, contant boundary condition, as concentration

         do kk = 1, ngwmobile_tracers
           if(.not. is_volatile(kk))cycle
           condc_toplay_col(c,groupid(kk))    = 1._r8/(snowres_col(c,volatilegroupid(kk))+&
             0.5_r8*dz_top(c)/diffblkm_topsoi_col(c,volatilegroupid(kk))) !m/s surface conductance
         enddo
      enddo

    end associate
  end subroutine set_boundary_conditions

  !-------------------------------------------------------------------------------
  subroutine calc_bgc_reaction(this, bounds, col, lbj, ubj, num_soilc, filter_soilc,               &
       num_soilp,filter_soilp, jtops, betr_time, betrtracer_vars, tracercoeff_vars,  biophysforc, &
       tracerstate_vars, tracerflux_vars, tracerboundarycond_vars, plant_soilbgc, &
       biogeo_flux, biogeo_state,  betr_status)
    !
    ! !DESCRIPTION:
    ! do bgc reaction
    !
    ! !USES:
    use TracerBoundaryCondType , only : tracerboundarycond_type
    use tracerfluxType         , only : tracerflux_type
    use tracerstatetype        , only : tracerstate_type
    use tracercoeffType        , only : tracercoeff_type
    use BetrTracerType         , only : betrtracer_type
    use PlantSoilBGCMod        , only : plant_soilbgc_type
    use BetrStatusType         , only : betr_status_type
    use betr_columnType        , only : betr_column_type
    use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
    use BeTR_biogeoStateType     , only : betr_biogeo_state_type
    use BeTR_TimeMod             , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(simic_bgc_reaction_type) , intent(inout) :: this                       !
    type(bounds_type)                 , intent(in)    :: bounds                      ! bounds
    type(betr_column_type)            , intent(in)    :: col
    integer                           , intent(in)    :: num_soilc                   ! number of columns in column filter
    integer                           , intent(in)    :: filter_soilc(:)             ! column filter
    integer                           , intent(in)    :: num_soilp
    integer                           , intent(in)    :: filter_soilp(:)
    integer                           , intent(in)    :: jtops( : )                  ! top index of each column
    integer                           , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
    type(betr_time_type)              , intent(in)    :: betr_time                       ! model time step
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars             ! betr configuration information
    type(betr_biogeophys_input_type)  , intent(inout) :: biophysforc
    type(tracercoeff_type)            , intent(inout) :: tracercoeff_vars
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars
    type(tracerflux_type)             , intent(inout) :: tracerflux_vars
    type(tracerboundarycond_type)     , intent(inout) :: tracerboundarycond_vars !
    class(plant_soilbgc_type)         , intent(inout) :: plant_soilbgc
    type(betr_biogeo_flux_type)       , intent(inout) :: biogeo_flux
    type(betr_biogeo_state_type)      , intent(inout) :: biogeo_state
    type(betr_status_type)            , intent(out)   :: betr_status
    character(len=*)                 , parameter     :: subname ='calc_bgc_reaction'

    integer :: c, fc, j
    character(len=5) :: laystr
    logical :: is_surflit  !surface litter layer?
    integer :: nstates
    real(r8), allocatable :: ystates0(:)
    real(r8), allocatable :: ystatesf(:)

    associate(                                                                    &
    tracer_mobile_phase            => tracerstate_vars%tracer_conc_mobile_col  ,  &
    tracer_flx_netpro_vr           => tracerflux_vars%tracer_flx_netpro_vr_col    &
    )

    call betr_status%reset()

    nstates = this%simic_bgc_index%nstvars
    allocate(ystates0(nstates))
    allocate(ystatesf(nstates))

    call this%set_bgc_forc(bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
        biophysforc, plant_soilbgc, betrtracer_vars, tracercoeff_vars, tracerstate_vars,betr_status)

    do j = lbj, ubj
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        if(j<jtops(c))cycle
        is_surflit=(j<=0)
        !do bgc reaction
        call this%simic_bgc(c,j)%runbgc(is_surflit, betr_time%delta_time, this%simic_forc(c,j), &
           nstates, ystates0, ystatesf, betr_status)

        if(betr_status%check_status())then
          write(laystr,'(I2.2)')j
          betr_status%msg=trim(betr_status%msg)//' lay '//trim(laystr)
          return
        endif

        call this%retrieve_output(c, j, nstates, ystates0, ystatesf, betr_time%delta_time, &
           betrtracer_vars, tracerflux_vars, tracerstate_vars, plant_soilbgc, biogeo_flux)

      enddo
    enddo

    !update phase change coefficients for tracers involved in sorptive reactions
    call this%update_sorpphase_coeff(bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
        betrtracer_vars, tracerstate_vars, tracercoeff_vars)

    deallocate(ystates0)
    deallocate(ystatesf)
   end associate
  end subroutine calc_bgc_reaction

  !-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars, betr_status)
    !
    ! DESCRIPTION:
    ! requilibrate tracers that has solid and mobile phases
    ! using the theory of mass action.
    !
    ! !USES:
    !
    use tracerstatetype , only : tracerstate_type
    use tracercoeffType , only : tracercoeff_type
    use BeTRTracerType  , only : betrtracer_type
    use BetrStatusType  , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(simic_bgc_reaction_type) , intent(inout)    :: this
    type(bounds_type)                 , intent(in)    :: bounds
    integer                           , intent(in)    :: lbj, ubj
    integer                           , intent(in)    :: jtops(bounds%begc: )        ! top label of each column
    integer                           , intent(in)    :: num_soilc
    integer                           , intent(in)    :: filter_soilc(:)
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars
    type(tracercoeff_type)            , intent(in)    :: tracercoeff_vars
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars
    type(betr_status_type)            , intent(out)   :: betr_status
    !local variables
    character(len=255) :: subname = 'do_tracer_equilibration'

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__), betr_status)

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)                          continue
    if (bounds%begc > 0)                                      continue
    if (ubj > lbj)                                            continue
    if (size(jtops) > 0)                                      continue
    if (num_soilc > 0)                                        continue
    if (size(filter_soilc) > 0)                               continue
    if (len(betrtracer_vars%betr_simname) > 0)                continue
    if (size(tracerstate_vars%tracer_conc_surfwater_col) > 0) continue
  !  if (size(tracercoeff_vars%annsum_counter_col) > 0)        continue

    !continue on the simulation type, an implementation of aqueous chemistry will be
    !employed to separate out the adsorbed phase
    !It should be noted that this formulation excludes the use of linear isotherm, which
    !can be integrated through the retardation factor

  end subroutine do_tracer_equilibration

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, col, betrtracer_vars, biophysforc, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! do cold initialization
    !
    ! !USES:
    use BeTRTracerType      , only : BeTRTracer_Type
    use tracerstatetype     , only : tracerstate_type
    use betr_varcon         , only : spval => bspval, ispval => bispval
    use BeTR_landvarconType , only : landvarcon  => betr_landvarcon
    use betr_columnType     , only : betr_column_type
    implicit none
    ! !ARGUMENTS:
    class(simic_bgc_reaction_type) , intent(inout)    :: this
    type(bounds_type)                 , intent(in)    :: bounds
    type(betr_column_type)            , intent(in)    :: col
    type(BeTRTracer_Type)             , intent(in)    :: betrtracer_vars
    type(betr_biogeophys_input_type)  , intent(in)    :: biophysforc
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars

    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, k, j
    integer :: fc                                        ! filter_soilc index
    integer :: begc, endc
    integer :: begg, endg
    integer :: trcid
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !-----------------------------------------------------------------------
    trcid = betrtracer_vars%id_trc_beg_Bm

    do c = bounds%begc, bounds%endc

      !dual phase tracers

      tracerstate_vars%tracer_conc_mobile_col(c,:, :)          = 0._r8
      tracerstate_vars%tracer_conc_surfwater_col(c,:)          = 0._r8
      tracerstate_vars%tracer_conc_aquifer_col(c,:)            = 0._r8
      tracerstate_vars%tracer_conc_grndwater_col(c,:)          = 0._r8

      if(betrtracer_vars%nsolid_equil_tracers>0)then
        tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = 0._r8
      endif
      tracerstate_vars%tracer_soi_molarmass_col(c,:)          = 0._r8
      !initialize microbial biomass
      tracerstate_vars%tracer_conc_mobile_col(c,:,trcid) = 1.e-2_r8
    enddo

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
    class(simic_bgc_reaction_type) , intent(inout)    :: this
    type(BeTRTracer_Type)             , intent(inout) :: betrtracer_vars
    character(len=*)                  , intent(in)  :: name_list_buffer

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)           continue
    if (len(betrtracer_vars%betr_simname) > 0) continue

    !do nothing here
  end subroutine readParams

  !-------------------------------------------------------------------------------
  subroutine retrieve_biogeoflux(this, num_soilc, filter_soilc, tracerflux_vars, &
  betrtracer_vars, biogeo_flux)

  use tracerfluxType           , only : tracerflux_type
  use BeTR_decompMod           , only : betr_bounds_type
  use BeTRTracerType           , only : BeTRTracer_Type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  implicit none
  class(simic_bgc_reaction_type) , intent(inout) :: this               !
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


   !----------------------------------------------------------------------
   subroutine retrieve_lnd2atm(this, bounds, num_soilc, filter_soilc, tracerflux_vars, &
   betrtracer_vars, biogeo_flux)

   use tracerfluxType           , only : tracerflux_type
   use BeTR_decompMod           , only : betr_bounds_type
   use BeTRTracerType           , only : BeTRTracer_Type
   use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
   implicit none
   class(simic_bgc_reaction_type) , intent(inout) :: this               !
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
   subroutine debug_info(this, bounds, num_soilc, filter_soilc, dzsoi, betrtracer_vars, tracerstate_vars, header, betr_status)

   use BeTRTracerType           , only : BeTRTracer_Type
   use tracerstatetype          , only : tracerstate_type
   use BeTR_decompMod           , only : betr_bounds_type
     ! !ARGUMENTS:
    implicit none
   class(simic_bgc_reaction_type) , intent(inout) :: this      !
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
   subroutine retrieve_biostates(this, bounds, lbj, ubj, jtops,num_soilc, filter_soilc, &
      betrtracer_vars, tracerstate_vars, biogeo_state, betr_status)
   !
   !retrieve state variables for lsm mass balance check
   use tracer_varcon, only : catomw, natomw, patomw, c13atomw, c14atomw
   use BeTR_decompMod           , only : betr_bounds_type
   use BeTRTracerType           , only : BeTRTracer_Type
   use tracerstatetype          , only : tracerstate_type
    use BeTR_biogeoStateType     , only : betr_biogeo_state_type
   implicit none
   class(simic_bgc_reaction_type) , intent(inout) :: this               !
   type(betr_bounds_type)               , intent(in)  :: bounds                      ! bounds
   integer                              , intent(in) :: lbj, ubj
   integer                              , intent(in) :: jtops(bounds%begc: )
   integer                              , intent(in)    :: num_soilc                   ! number of columns in column filter
   integer                              , intent(in)    :: filter_soilc(:)             ! column filter
   type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
   type(tracerstate_type)               , intent(inout) :: tracerstate_vars
   type(betr_biogeo_state_type)         , intent(inout) :: biogeo_state
   type(betr_status_type)               , intent(out):: betr_status
   !local variables

   integer :: nelm
   integer :: c_loc, c13_loc, c14_loc
   integer :: c, fc, j, kk


   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__), betr_status)

   c_loc=this%simic_bgc_index%c_loc
   c13_loc=this%simic_bgc_index%c13_loc
   c14_loc=this%simic_bgc_index%c14_loc
   nelm =this%simic_bgc_index%nelms

   do j = lbj, ubj
     do fc = 1, num_soilc
        c = filter_soilc(fc)
        if(j<jtops(c))cycle

        !add litter
        do kk = betrtracer_vars%id_trc_beg_litr, betrtracer_vars%id_trc_end_litr, nelm
          biogeo_state%c12state_vars%totlitc_vr_col(c,j) = biogeo_state%c12state_vars%totlitc_vr_col(c,j) + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)
!          print*,'totlit',c,j,biogeo_state%c12state_vars%totlitc_vr_col(c,j)
          if(this%use_c13)then
            biogeo_state%c13state_vars%totlitc_vr_col(c,j) = biogeo_state%c13state_vars%totlitc_vr_col(c,j) + &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%totlitc_vr_col(c,j) = biogeo_state%c14state_vars%totlitc_vr_col(c,j) + &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif

        enddo

        !add cwd
        do kk = betrtracer_vars%id_trc_beg_wood, betrtracer_vars%id_trc_end_wood, nelm
          biogeo_state%c12state_vars%cwdc_vr_col(c,j) = biogeo_state%c12state_vars%cwdc_vr_col(c,j) + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%cwdc_vr_col(c,j) = biogeo_state%c13state_vars%cwdc_vr_col(c,j) + &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%cwdc_vr_col(c,j) = biogeo_state%c14state_vars%cwdc_vr_col(c,j) + &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        !Microbial biomass
        do kk = betrtracer_vars%id_trc_beg_Bm, betrtracer_vars%id_trc_end_Bm-nelm, nelm
          biogeo_state%c12state_vars%som1c_vr_col(c,j) =  &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%som1c_vr_col(c,j) = &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%som1c_vr_col(c,j) =  &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        do kk = betrtracer_vars%id_trc_beg_pom, betrtracer_vars%id_trc_end_pom, nelm+1
          biogeo_state%c12state_vars%som2c_vr_col(c,j) = &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%som2c_vr_col(c,j) =  &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%som2c_vr_col(c,j) =  &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        do kk = betrtracer_vars%id_trc_beg_dom, betrtracer_vars%id_trc_end_dom, nelm+1
          biogeo_state%c12state_vars%domc_vr_col(c,j) = &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%domc_vr_col(c,j) =  &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%domc_vr_col(c,j) =  &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        do kk = betrtracer_vars%id_trc_beg_Bm+nelm, betrtracer_vars%id_trc_end_Bm, nelm
          biogeo_state%c12state_vars%som3c_vr_col(c,j) =  &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%som3c_vr_col(c,j) = &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%som3c_vr_col(c,j) =  &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo
     enddo
   enddo

   end subroutine retrieve_biostates


  !------------------------------------------------------------------------------
  subroutine retrieve_output(this, c, j, nstates, ystates0, ystatesf, dtime, betrtracer_vars, tracerflux_vars,&
     tracerstate_vars, plant_soilbgc, biogeo_flux)
  !DESCRIPTION
  !retrieve flux and state variables after evolving the bgc calculation
  !
  !USES
  use BetrTracerType           , only : betrtracer_type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  use tracerfluxType           , only : tracerflux_type
  use tracerstatetype          , only : tracerstate_type
  use betr_ctrl                , only : betr_spinup_state
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use simicPlantSoilBGCType    , only : simic_plant_soilbgc_type
  use tracer_varcon            , only : catomw, natomw, patomw, fix_ip
  implicit none
  class(simic_bgc_reaction_type) , intent(inout)    :: this
  integer                              , intent(in) :: c, j
  integer                              , intent(in) :: nstates
  real(r8)                             , intent(in) :: ystates0(nstates)
  real(r8)                             , intent(inout) :: ystatesf(nstates)
  real(r8)                             , intent(in) :: dtime
  type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
  type(tracerstate_type)               , intent(inout) :: tracerstate_vars
  type(tracerflux_type)                , intent(inout) :: tracerflux_vars
  class(plant_soilbgc_type)            , intent(inout) :: plant_soilbgc
  type(betr_biogeo_flux_type)          , intent(inout) :: biogeo_flux

  integer :: k, k1, k2, jj, p
  integer :: trcid

  associate( &
     litr_beg =>  this%simic_bgc_index%litr_beg  , &
     litr_end =>  this%simic_bgc_index%litr_end  , &
     wood_beg =>  this%simic_bgc_index%wood_beg  , &
     wood_end =>  this%simic_bgc_index%wood_end  , &
     dom_beg =>  this%simic_bgc_index%dom_beg    , &
     dom_end =>  this%simic_bgc_index%dom_end    , &
     pom_beg =>  this%simic_bgc_index%pom_beg    , &
     pom_end =>  this%simic_bgc_index%pom_end    , &
     Bm_beg  =>  this%simic_bgc_index%Bm_beg     , &
     Bm_end  =>  this%simic_bgc_index%Bm_end     , &
     volatileid            => betrtracer_vars%volatileid                   , &
     tracer_flx_netpro_vr  => tracerflux_vars%tracer_flx_netpro_vr_col     , & !
     tracer_flx_parchm_vr  => tracerflux_vars%tracer_flx_parchm_vr_col     , & !
     ngwmobile_tracers     => betrtracer_vars%ngwmobile_tracers              & !
  )

    !tracer states
    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_litr:betrtracer_vars%id_trc_end_litr) = &
        ystatesf(litr_beg:litr_end)
!    print*,'c,j',c,j,ystatesf(litr_beg:litr_end)
    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_wood:betrtracer_vars%id_trc_end_wood) = &
        ystatesf(wood_beg:wood_end)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_Bm:betrtracer_vars%id_trc_end_Bm) = &
        ystatesf(Bm_beg:Bm_end)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_dom:betrtracer_vars%id_trc_end_dom) = &
        ystatesf(dom_beg:dom_end)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_pom:betrtracer_vars%id_trc_end_pom) = &
        ystatesf(pom_beg:pom_end)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2) = &
        ystatesf(this%simic_bgc_index%lid_n2)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2) = &
        ystatesf(this%simic_bgc_index%lid_o2)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar) = &
        ystatesf(this%simic_bgc_index%lid_ar)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x) = &
        ystatesf(this%simic_bgc_index%lid_co2)

    if(this%use_c13)then
      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x) = &
          ystatesf(this%simic_bgc_index%lid_c13_co2)
    endif

    if(this%use_c14)then
      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x) = &
          ystatesf(this%simic_bgc_index%lid_c14_co2)
    endif

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4) = &
        ystatesf(this%simic_bgc_index%lid_ch4)

    !fluxes
    tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_o2) ) = &
         ystatesf(this%simic_bgc_index%lid_o2_paere )  - &
         ystates0(this%simic_bgc_index%lid_o2_paere)

    if ( betr_spinup_state == 0 ) then
      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2)  ) = &
          ystatesf(this%simic_bgc_index%lid_n2_paere)  - &
          ystates0(this%simic_bgc_index%lid_n2_paere)

      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ar)  ) = &
          ystatesf(this%simic_bgc_index%lid_ar_paere)  - &
          ystates0(this%simic_bgc_index%lid_ar_paere)

      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_co2x)) = &
          ystatesf(this%simic_bgc_index%lid_co2_paere)  - &
          ystates0(this%simic_bgc_index%lid_co2_paere)

      if(this%use_c13)then
        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c13_co2x)) = &
            ystatesf(this%simic_bgc_index%lid_c13_co2_paere)  - &
            ystates0(this%simic_bgc_index%lid_c13_co2_paere)
      endif

      if(this%use_c14)then
        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c14_co2x)) = &
            ystatesf(this%simic_bgc_index%lid_c14_co2_paere)  - &
            ystates0(this%simic_bgc_index%lid_c14_co2_paere)
      endif

      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ch4) ) = &
          ystatesf(this%simic_bgc_index%lid_ch4_paere)  - &
          ystates0(this%simic_bgc_index%lid_ch4_paere)

    endif

    !get net production for om pools
    do k = 1, litr_end-litr_beg + 1
      k1 = litr_beg+k-1; k2 = betrtracer_vars%id_trc_beg_litr + k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo
    do k = 1, wood_end-wood_beg + 1
      k1 = wood_beg+k-1; k2 = betrtracer_vars%id_trc_beg_wood + k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo
    do k = 1, Bm_end-Bm_beg + 1
      k1 = Bm_beg+k-1; k2 = betrtracer_vars%id_trc_beg_Bm+ k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo

    do k = 1, dom_end-dom_beg + 1
      k1 = dom_beg+k-1; k2 = betrtracer_vars%id_trc_beg_dom+ k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo

    do k = 1, pom_end-pom_beg + 1
      k1 = pom_beg+k-1; k2 = betrtracer_vars%id_trc_beg_pom+ k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2) = &
        ystatesf(this%simic_bgc_index%lid_n2) - &
        ystates0(this%simic_bgc_index%lid_n2)

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_co2x ) = &
        ystatesf(this%simic_bgc_index%lid_co2) - &
        ystates0(this%simic_bgc_index%lid_co2)

    if(this%use_c13)then
      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c13_co2x ) = &
          ystatesf(this%simic_bgc_index%lid_c13_co2) - &
          ystates0(this%simic_bgc_index%lid_c13_co2)
    endif

    if(this%use_c14)then
      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c14_co2x ) = &
          ystatesf(this%simic_bgc_index%lid_c14_co2) - &
          ystates0(this%simic_bgc_index%lid_c14_co2)
    endif

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_o2   ) = &
        ystatesf(this%simic_bgc_index%lid_o2) - &
        ystates0(this%simic_bgc_index%lid_o2)

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ch4  ) = &
        ystatesf(this%simic_bgc_index%lid_ch4) - &
        ystates0(this%simic_bgc_index%lid_ch4)

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ar) = &
        ystatesf(this%simic_bgc_index%lid_ar) - &
        ystates0(this%simic_bgc_index%lid_ar)

    !get net production for om pools
    do k = 1, litr_end-litr_beg + 1
      k1 = litr_beg+k-1; k2 = betrtracer_vars%id_trc_beg_litr + k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo

    do k = 1, wood_end-wood_beg + 1
      k1 = wood_beg+k-1; k2 = betrtracer_vars%id_trc_beg_wood + k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo

    do k = 1, Bm_end-Bm_beg + 1
      k1 = Bm_beg+k-1; k2 = betrtracer_vars%id_trc_beg_Bm+ k-1
      tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
    enddo

    biogeo_flux%c12flux_vars%hr_vr_col(c,j) = &
        (ystatesf(this%simic_bgc_index%lid_co2_hr) - &
        ystates0(this%simic_bgc_index%lid_co2_hr))*catomw/dtime
  end associate
  end subroutine retrieve_output


!------------------------------------------------------------------------------
  subroutine set_bgc_forc(this, bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
      biophysforc, plant_soilbgc, betrtracer_vars, tracercoeff_vars, tracerstate_vars, betr_status)

  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use tracerstatetype          , only : tracerstate_type
  use betr_decompMod           , only : betr_bounds_type
  use tracercoeffType          , only : tracercoeff_type
  use betr_columnType          , only : betr_column_type
  use BetrTracerType           , only : betrtracer_type
  use simicPlantSoilBGCType    , only : simic_plant_soilbgc_type
  use MathfuncMod              , only : fpmax
  use betr_varcon              , only : grav => bgrav
  implicit none
  class(simic_bgc_reaction_type) , intent(inout) :: this                       !
  type(bounds_type)                    , intent(in) :: bounds                         ! bounds
  type(betr_column_type)               , intent(in) :: col
  integer                              , intent(in) :: jtops(bounds%begc: ) ! top index of each column
  integer                              , intent(in) :: lbj, ubj                       ! lower and upper bounds, make sure they are > 0
  integer                              , intent(in) :: num_soilc       ! number of columns in column filter
  integer                              , intent(in) :: filter_soilc(:) ! column filter
  type(betr_biogeophys_input_type)     , intent(in) :: biophysforc
  class(plant_soilbgc_type)            , intent(in) :: plant_soilbgc
  type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
  type(tracerstate_type)               , intent(in) :: tracerstate_vars
  type(tracercoeff_type)               , intent(in) :: tracercoeff_vars
  type(betr_status_type)               , intent(out)   :: betr_status

  integer :: j, fc, c
  integer :: k1, k2
  real(r8), parameter :: tiny_cval =1.e-16_r8
  associate( &
     litr_beg =>  this%simic_bgc_index%litr_beg  , &
     litr_end =>  this%simic_bgc_index%litr_end  , &
     wood_beg =>  this%simic_bgc_index%wood_beg  , &
     wood_end =>  this%simic_bgc_index%wood_end  , &
     dom_beg =>  this%simic_bgc_index%dom_beg    , &
     dom_end =>  this%simic_bgc_index%dom_end    , &
     pom_beg =>  this%simic_bgc_index%pom_beg    , &
     pom_end =>  this%simic_bgc_index%pom_end    , &
     Bm_beg  =>  this%simic_bgc_index%Bm_beg     , &
     Bm_end  =>  this%simic_bgc_index%Bm_end       &
  )
  call betr_status%reset()
  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__),betr_status)

  do j = lbj, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      if(j<jtops(c))cycle
      this%simic_forc(c,j)%plant_ntypes = this%nactpft
      this%simic_forc(c,j)%ystates(:) = 0._r8

      !litter
      this%simic_forc(c,j)%ystates(litr_beg:litr_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_litr:betrtracer_vars%id_trc_end_litr)

      !wood
      this%simic_forc(c,j)%ystates(wood_beg:wood_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_wood:betrtracer_vars%id_trc_end_wood)

      !dom
      this%simic_forc(c,j)%ystates(dom_beg:dom_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_dom:betrtracer_vars%id_trc_end_dom)
      if(this%simic_forc(c,j)%ystates(dom_beg)<=tiny_cval)this%simic_forc(c,j)%ystates(dom_beg:dom_end)=0._r8

      !pom
      this%simic_forc(c,j)%ystates(pom_beg:pom_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_pom:betrtracer_vars%id_trc_end_pom)
      if(this%simic_forc(c,j)%ystates(pom_beg)<=tiny_cval)this%simic_forc(c,j)%ystates(pom_beg:pom_end)=0._r8


      !microbial biomass
      this%simic_forc(c,j)%ystates(Bm_beg:Bm_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_Bm:betrtracer_vars%id_trc_end_Bm)
      if(this%simic_forc(c,j)%ystates(Bm_beg)<=tiny_cval)this%simic_forc(c,j)%ystates(Bm_beg:Bm_end)=0._r8
      !non-soluble phase of mineral p

      this%simic_forc(c,j)%ystates(this%simic_bgc_index%lid_n2) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2))

      this%simic_forc(c,j)%ystates(this%simic_bgc_index%lid_o2) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2))

      this%simic_forc(c,j)%ystates(this%simic_bgc_index%lid_ar) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar))

      this%simic_forc(c,j)%ystates(this%simic_bgc_index%lid_co2)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x))

      if(this%use_c13)then
        this%simic_forc(c,j)%ystates(this%simic_bgc_index%lid_c13_co2)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%simic_forc(c,j)%ystates(this%simic_bgc_index%lid_c14_co2)= &
          fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x))
      endif

      this%simic_forc(c,j)%ystates(this%simic_bgc_index%lid_ch4)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4))


      !input
      this%simic_forc(c,j)%cflx_input_litr_met = biophysforc%c12flx%cflx_input_litr_met_vr_col(c,j)
      this%simic_forc(c,j)%cflx_input_litr_cel = biophysforc%c12flx%cflx_input_litr_cel_vr_col(c,j)
      this%simic_forc(c,j)%cflx_input_litr_lig = biophysforc%c12flx%cflx_input_litr_lig_vr_col(c,j)
      this%simic_forc(c,j)%cflx_input_litr_cwd = biophysforc%c12flx%cflx_input_litr_cwd_vr_col(c,j)
      this%simic_forc(c,j)%cflx_input_litr_lwd = biophysforc%c12flx%cflx_input_litr_lwd_vr_col(c,j)
      this%simic_forc(c,j)%cflx_input_litr_fwd = biophysforc%c12flx%cflx_input_litr_fwd_vr_col(c,j)

      !environmental variables
      this%simic_forc(c,j)%temp   = biophysforc%t_soisno_col(c,j)            !temperature
      this%simic_forc(c,j)%depz   = biophysforc%z(c,j)            !depth of the soil
      this%simic_forc(c,j)%dzsoi  = biophysforc%dz(c,j)            !soil thickness
      this%simic_forc(c,j)%sucsat  = biophysforc%sucsat_col(c,j)            ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
      this%simic_forc(c,j)%soilpsi = max(biophysforc%smp_l_col(c,j)*grav*1.e-6_r8,-15._r8)    ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
      this%simic_forc(c,j)%bsw = biophysforc%bsw_col(c,j)
      this%simic_forc(c,j)%bd   = biophysforc%bd_col(c,j)              !bulk density
      this%simic_forc(c,j)%pct_sand = biophysforc%cellsand_col(c,j)
      this%simic_forc(c,j)%pct_clay = biophysforc%cellclay_col(c,j)
      this%simic_forc(c,j)%h2osoi_vol = biophysforc%h2osoi_vol_col(c,j)
      this%simic_forc(c,j)%h2osoi_liq = biophysforc%h2osoi_liq_col(c,j)
      this%simic_forc(c,j)%air_vol = biophysforc%air_vol_col(c,j)
      this%simic_forc(c,j)%finundated = biophysforc%finundated_col(c)
      this%simic_forc(c,j)%watsat = biophysforc%watsat_col(c,j)
      this%simic_forc(c,j)%watfc = biophysforc%watfc_col(c,j)
      this%simic_forc(c,j)%cellorg = biophysforc%cellorg_col(c,j)
      this%simic_forc(c,j)%pH = biophysforc%soil_pH(c,j)

      !conductivity for plant-aided gas transport
      this%simic_forc(c,j)%aren_cond_n2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2))
      this%simic_forc(c,j)%aren_cond_o2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2))
      this%simic_forc(c,j)%aren_cond_co2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%simic_forc(c,j)%aren_cond_co2_c13 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c13_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%simic_forc(c,j)%aren_cond_co2_c14 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c14_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%simic_forc(c,j)%aren_cond_ar = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar))
      this%simic_forc(c,j)%aren_cond_ch4 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4))
      !phase conversion parameter
      this%simic_forc(c,j)%ch4_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4))
      this%simic_forc(c,j)%co2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x))
      this%simic_forc(c,j)%o2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2))
      this%simic_forc(c,j)%n2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2))
      this%simic_forc(c,j)%ar_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar))
      this%simic_forc(c,j)%o2_w2b = &
          tracercoeff_vars%aqu2bulkcef_mobile_col(c,j,betrtracer_vars%groupid(betrtracer_vars%id_trc_o2))

      !atmospheric pressure (mol/m3) for gas ventilation.
      this%simic_forc(c,j)%conc_atm_n2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_n2))
      this%simic_forc(c,j)%conc_atm_o2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_o2))
      this%simic_forc(c,j)%conc_atm_ar = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_ar))
      this%simic_forc(c,j)%conc_atm_co2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%simic_forc(c,j)%conc_atm_co2_c13 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%simic_forc(c,j)%conc_atm_co2_c14 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%simic_forc(c,j)%conc_atm_ch4 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_ch4))

      this%simic_forc(c,j)%soilorder = biophysforc%isoilorder(c)

    enddo
  enddo
  end associate
  end subroutine set_bgc_forc

!------------------------------------------------------------------------------
  subroutine update_sorpphase_coeff(this, bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
      betrtracer_vars, tracerstate_vars, tracercoeff_vars)
  !
  !DESCRIPTION
  !update sorption related phase conversion parameter
  !in this formulation, the amount sorption surface is assumed to scale linearly with moisture
  !content. Microbes are always assumed in sufficiently moist status, though the activity could be
  !smaller
  use tracerstatetype          , only : tracerstate_type
  use betr_decompMod           , only : betr_bounds_type
  use tracercoeffType          , only : tracercoeff_type
  use betr_columnType          , only : betr_column_type
  use BetrTracerType           , only : betrtracer_type
  implicit none
  class(simic_bgc_reaction_type) , intent(inout) :: this                       !
  type(bounds_type)                    , intent(in) :: bounds                         ! bounds
  type(betr_column_type)               , intent(in) :: col
  integer                              , intent(in) :: jtops(bounds%begc: ) ! top index of each column
  integer                              , intent(in) :: lbj, ubj                       ! lower and upper bounds, make sure they are > 0
  integer                              , intent(in) :: num_soilc       ! number of columns in column filter
  integer                              , intent(in) :: filter_soilc(:) ! column filter
  type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
  type(tracerstate_type)               , intent(in) :: tracerstate_vars
  type(tracercoeff_type)            , intent(inout) :: tracercoeff_vars

  real(r8) :: KM_CM, Msurf, KM_EM, BMT
  real(r8) :: denorm1, denorm0, denorm2,beta
  integer :: c_l, j
  associate(                                                           &
    aqu2bulkcef_mobile   => tracercoeff_vars%aqu2bulkcef_mobile_col  , & !Output:[real(r8)(:,:)], phase conversion coeff
    id_trc_dom           => betrtracer_vars%id_trc_dom               , &
    trcid_Bm             => betrtracer_vars%id_trc_beg_Bm            , &
    trcid_dom            => betrtracer_vars%id_trc_beg_dom           , &
    trcid_pom            => betrtracer_vars%id_trc_beg_pom           , &
    id_trc_end_dom       => betrtracer_vars%id_trc_end_dom           , &
    tracer_conc_mobile   => tracerstate_vars%tracer_conc_mobile_col  , &
    Kaff_CM              => simic_para%Kaff_CM                       , &
    Kaff_EM              => simic_para%Kaff_EM                       , &
    Kaff_BC              => simic_para%Kaff_BC                       , &
    alpha_B2E            => simic_para%alpha_B2E                     , &
    alpha_B2T            => simic_para%alpha_B2T                     , &
    nelms                =>  this%simic_bgc_index%nelms                &
  )
  c_l=1
  do j = 1, ubj
    KM_CM=aqu2bulkcef_mobile(c_l,j,id_trc_dom)*this%simic_forc(c_l,j)%KM_OM_ref*Kaff_CM
    BMT=(tracer_conc_mobile(c_l,j,trcid_Bm)+tracer_conc_mobile(c_l,j,trcid_Bm+nelms))*alpha_B2T
    Msurf=this%simic_forc(c_l,j)%Msurf_OM-tracer_conc_mobile(c_l,j,trcid_pom)
    denorm0=1._r8+Msurf/KM_CM+BMT/Kaff_BC
    denorm1=denorm0+tracer_conc_mobile(c_l,j,trcid_dom)/KM_CM
    denorm2=denorm0+tracer_conc_mobile(c_l,j,trcid_dom)/Kaff_BC
    beta=1._r8/(1._r8-Msurf/KM_CM/denorm1-BMT/Kaff_BC/denorm2)
    aqu2bulkcef_mobile(c_l,j,id_trc_dom) = aqu2bulkcef_mobile(c_l,j,id_trc_dom)*beta
    !print*,'c,j',1._r8/(1._r8-Msurf/KM_CM/denorm1-BMT/Kaff_BC/denorm2)
  enddo
  end associate
  end subroutine update_sorpphase_coeff


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
       class(simic_bgc_reaction_type) , intent(inout) :: this
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
       class(simic_bgc_reaction_type) , intent(inout) :: this
     integer, intent(in) :: parcol

     this%parcol = parcol
   end subroutine SetParCols
end module SimicBGCReactionsType
