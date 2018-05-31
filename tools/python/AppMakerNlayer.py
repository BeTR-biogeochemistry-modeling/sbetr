#!/usr/bin/env python2
#in app_nameNlayer write three files
#CMakeLists.txt
#app_nameBGCReactionsType.F90
#app_namePlantSoilBGCType.F90

def MakeNlayer(sfarm_dir, app_name):
    print "create file "+sfarm_dir+'/'+app_name+'Nlayer'+"/CMakeLists.txt"

    print "set("+app_name.upper()+"NLAYER_SOURCES"
    print "  "+app_name+"BGCReactionsType.F90"
    print "  "+app_name+"PlantSoilBGCType.F90"
    print ")"
    print ""
    print "include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/"+app_name+"/"+app_name+"Para)"
    print "include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/"+app_name+"/"+app_name+"1layer)"
    print "include(add_betr_library)"
    print "add_betr_library("+app_name+"Nlayer ${"+app_name.upper()+"NLAYER_SOURCES})"

    print "set(BETR_LIBRARIES "+app_name+"Nlayer;${BETR_LIBRARIES} PARENT_SCOPE)"
    print "set(BETR_LIBRARIES "+app_name+"Nlayer;${BETR_LIBRARIES})"

#X#add_subdirectory(tests)

    print 'if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")'
    print "  install(TARGETS "+app_name+"Nlayer DESTINATION lib)"
    print "  file(GLOB HEADERS *.h)"
    print "  install(FILES ${HEADERS} DESTINATION include/soil-farm/"+app_name+"/"+app_name+"Nlayer)"
    print "endif()"

#app_nameBGCReactionsType.F90
    print "create file "+sfarm_dir+'/'+app_name+'Nlayer/'+app_name+"BGCReactionsType.F90"
    print "module "+app_name+"BGCReactionsType"

    print '#include "bshr_assert.h"'
    print "!"
    print "! !DESCRIPTION:"
    print "! This is an example generated from AppMaker for use in betr"
    print "!"
    print "! HISTORY:"
    print "! Created by Jinyun Tang, Oct 2nd, 2014"
    print "! !USES:"
    print "  use bshr_log_mod             , only : errMsg => shr_log_errMsg"
    print "  use bshr_kind_mod            , only : r8 => shr_kind_r8"
    print "  use BeTR_decompMod           , only : bounds_type  => betr_bounds_type"
    print "  use BGCReactionsMod          , only : bgc_reaction_type"
    print "  use tracer_varcon            , only : bndcond_as_conc, bndcond_as_flux"
    print "  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type"
    print "  use BetrStatusType           , only : betr_status_type"
    print "  use "+app_name+"ParaType            , only : "+app_name+"_para"
    print "  use "+app_name+"BGCIndexType        , only : "+app_name+"_bgc_index_type"
    print "  use "+app_name+"BGCType             , only : "+app_name+"_bgc_type"
    print "  use JarBgcForcType           , only : JarBGC_forc_type"
    print "implicit none"
    print ""
    print "  private"

    print "  character(len=*), parameter :: mod_filename = &"
    print "    __FILE__"

    print "  public :: "+app_name+"_bgc_reaction_type"

    print "  type, extends(bgc_reaction_type) :: &"
    print "    "+app_name+"_bgc_reaction_type"
    print "  private"
    print "    type("+app_name+"_bgc_type), pointer :: "+app_name+"_bgc(:,:)"
    print "    type(JarBGC_forc_type), pointer :: "+app_name+"_forc(:,:)"
    print "    type("+app_name+"_bgc_index_type) :: "+app_name+"_bgc_index"
    print "    logical :: use_c13"
    print "    logical :: use_c14"
    print "    logical :: nop_limit"
    print "    logical :: non_limit"
    print "    integer :: nactpft               ! number of active pfts"
    print "  contains"
    print "    procedure :: Init_betrbgc                          ! initialize betr bgc"
    print "    procedure :: set_boundary_conditions               ! set top/bottom boundary conditions for various tracers"
    print "    procedure :: calc_bgc_reaction                     ! doing bgc calculation"
    print "    procedure :: init_boundary_condition_type          ! initialize type of top boundary conditions"
    print "    procedure :: do_tracer_equilibration               ! do equilibrium tracer chemistry"
    print "    procedure :: InitCold                              ! do cold initialization"
    print "    procedure :: retrieve_biogeoflux           !"
    print "    procedure :: set_kinetics_par"
    print "    procedure :: retrieve_lnd2atm"
    print "    procedure :: readParams                   ! read in parameters"
    print "    procedure :: retrieve_biostates"
    print "    procedure :: debug_info"
    print "    procedure :: set_bgc_spinup"
    print "    procedure :: UpdateParas"
    print "    procedure :: init_iP_prof"
    print "    procedure, private :: set_tracer"
    print "    procedure, private :: InitAllocate"
    print "    procedure, private :: retrieve_output"
    print "    procedure, private :: set_bgc_forc"
    print "    procedure, private :: update_sorpphase_coeff"
    print "  end type "+app_name+"_bgc_reaction_type"
    print ""
    print "  interface "+app_name+"_bgc_reaction_type"
    print "    module procedure constructor"
    print "  end interface "+app_name+"_bgc_reaction_type"
    print ""
    print "contains"
    print "!-------------------------------------------------------------------------------"
    print "  type("+app_name+"_bgc_reaction_type) function constructor()"
    print "  !"
    print "  ! !DESCRIPTION:"
    print "  ! create an object of type "+app_name+"_bgc_reaction_type."
    print ""
    print "  type("+app_name+"_bgc_reaction_type), allocatable :: bgc"
    print "  allocate(bgc)"
    print "  constructor = bgc"
    print "  end function constructor"
    print "!-------------------------------------------------------------------------------"
    print "  subroutine UpdateParas(this, bounds, lbj, ubj, bstatus)"
    print "  implicit none"
    print "  class("+app_name+"_bgc_reaction_type)         , intent(inout)    :: this"
    print "  type(bounds_type)                    , intent(in)    :: bounds"
    print "  integer                              , intent(in)    :: lbj, ubj        ! lower and upper bounds, make sure they are > 0"
    print "  type(betr_status_type)               , intent(out)   :: bstatus"
    print "  integer :: c, j"
    print ""
    print "  if (this%dummy_compiler_warning) continue"
    print "  call bstatus%reset()"
    print '  !do nothing'
    print "  end subroutine UpdateParas"
    print '!----------------------------------------------------------------------"
    print "  subroutine init_iP_prof(this, bounds, lbj, ubj, biophysforc, tracers, tracerstate_vars)"
    print "  !"
    print "  !DESCRIPTION"
    print "  ! set up initial inorganic P profile"
    print "  use tracer_varcon, only : patomw"
    print "  use tracerstatetype        , only : tracerstate_type"
    print "  use BeTRTracerType         , only : betrtracer_type"
    print "  implicit none
    print "  ! !ARGUMENTS:"
    print "  class("+app_name+"_bgc_reaction_type)         , intent(inout)    :: this"
    print "  type(bounds_type)                        , intent(in) :: bounds"
    print "  integer                                  , intent(in) :: lbj, ubj"
    print "  type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc"
    print "  type(BeTRtracer_type)                    , intent(inout) :: tracers"
    print "  type(tracerstate_type)                   , intent(inout) :: tracerstate_vars"
    print ""
    print "  if (this%dummy_compiler_warning) continue"
    print "  if (bounds%begc > 0) continue"
    print ""
    print "  end subroutine init_iP_prof"
    print "!----------------------------------------------------------------------"
    print "  subroutine set_kinetics_par(this, lbj, ubj, nactpft, plantNutkinetics)"
    print "  !"
    print "  !DESCRIPTION"
    print "  !set up kinetic parameters needed by the model"
    print "  use PlantNutKineticsMod, only : PlantNutKinetics_type"
    print "  ! !ARGUMENTS:"
    print "  class("+app_name+"_bgc_reaction_type)         , intent(inout)    :: this"
    print "  class(PlantNutKinetics_type), intent(in) :: plantNutkinetics"
    print "  integer, intent(in) :: lbj, ubj"
    print "  integer, intent(in) :: nactpft"
    print "  integer :: c_l, p, j"
    print "  !in the following, only one column is assumed for the bgc"
    print "  c_l = 1"
    print "  this%nactpft = nactpft"
    print "  do j = lbj, ubj"
    print "    !effective p competing decomposers"
    print "    this%"+app_name+"_forc(c_l,j)%Msurf_OM  = plantNutkinetics%minsurf_dom_compet_vr_col(c_l,j)"
    print "    this%"+app_name+"_forc(c_l,j)%KM_OM_ref = plantNutkinetics%km_minsurf_dom_vr_col(c_l,j)"
    print "  enddo"
    print "  end subroutine set_kinetics_par"
    print "!-------------------------------------------------------------------------------"
    print "  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )"
    print "  !"
    print "  ! !DESCRIPTION:"
    print "  ! initialize boundary condition types"
    print "  !"
    print "  ! !USES:"
    print "  use BeTRTracerType        , only : betrtracer_type"
    print "  use TracerBoundaryCondType, only : tracerboundarycond_type"
    print "  use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux"
    print "  use BeTRTracerType        , only : betrtracer_type"
    print "  ! !ARGUMENTS:"
    print "  class("+app_name+"_bgc_reaction_type), intent(inout) :: this"
    print "  type(BeTRtracer_type),             intent(in) :: betrtracer_vars"
    print "  type(bounds_type),                 intent(in) :: bounds"
    print "  type(tracerboundarycond_type),     intent(in) :: tracerboundarycond_vars"
    print "  ! !LOCAL VARIABLES:"
    print "  ! remove compiler warnings for unused dummy args"
    print "  if (this%dummy_compiler_warning) continue"
    print "  if (bounds%begc > 0) continue"
    print "  if (len(betrtracer_vars%betr_simname) > 0) continue"
    print "  tracerboundarycond_vars%topbc_type(:) = bndcond_as_conc"
    print "  !when bottom BC is not given, it is specified as constant flux"
    print "  !X!tracerboundarycond_vars%botbc_type(:) = bndcond_as_flux"
    print "  end subroutine init_boundary_condition_type"
    print "!-------------------------------------------------------------------------------"
    print "subroutine set_bgc_spinup(this, bounds, lbj, ubj, biophysforc, &"
    print "  tracers, tracerstate_vars)"
    print "  use tracerstatetype        , only : tracerstate_type"
    print "  use BeTRTracerType         , only : betrtracer_type"
    print "  use BeTR_decompMod         , only : betr_bounds_type"
    print ""
    print "  implicit none"
    print "  class("+app_name+"_bgc_reaction_type), intent(inout)    :: this"
    print "  type(betr_bounds_type)                  , intent(in) :: bounds"
    print "  integer                                 , intent(in) :: lbj, ubj"
    print "  type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc"
    print "  type(BeTRtracer_type)                   , intent(inout) :: tracers"
    print "  type(tracerstate_type)                  , intent(inout) :: tracerstate_vars"
    print ""
    print "  if (this%dummy_compiler_warning) continue"
    print "  if (bounds%begc > 0) continue"
    print "  end subroutine set_bgc_spinup"
    print "!-------------------------------------------------------------------------------"
    print "  subroutine InitAllocate(this, bounds, lbj, ubj, bstatus)"
    print "  !"
    print "  !!DESCRIPTION"
    print "  !allocate memory"
    print "  use betr_varcon                      , only : betr_maxpatch_pft"
    print "  implicit none"
    print "  class("+app_name+"_bgc_reaction_type), intent(inout)    :: this"
    print "  type(betr_status_type)           , intent(out)   :: bstatus"
    print "  type(bounds_type)                , intent(in)    :: bounds"
    print "  integer                          , intent(in)    :: lbj, ubj"
    print "  integer :: j, c"
    print "  logical :: batch_mode"
    print "  batch_mode =.false."
    print "  this%nactpft = 0"
    print ""
    print "  call this%"+app_name+"_bgc_index%Init("+app_name+"_para%use_c13, "+app_name+"_para%use_c14, &"
    print "  "+app_name+"_para%non_limit, "+app_name+"_para%nop_limit, betr_maxpatch_pft)"
    print "  if(bstatus%check_status())return"
    print "  !create the models"
    print "  allocate(this%"+app_name+"_bgc(bounds%begc:bounds%endc,lbj:ubj))"
    print "  !create model specific forcing data structure"
    print "  allocate(this%"+app_name+"_forc(bounds%begc:bounds%endc,lbj:ubj))"
    print "  !initialize"
    print "  do j = lbj, ubj"
    print "    do c = bounds%begc, bounds%endc"
    print "      call this%"+app_name+"_bgc(c,j)%Init("+app_name+"_para, batch_mode, bstatus)"
    print "      if(bstatus%check_status())return"
    print "      call this%"+app_name+"_forc(c,j)%Init(this%"+app_name+"_bgc_index%nstvars)"
    print "    enddo"
    print "  enddo"
    print "  end subroutine InitAllocate"
    print "!-------------------------------------------------------------------------------"
    print "  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars, namelist_buffer, bstatus)"
    print "!"
    print "! DESCRIPTION:"
    print "! initialize the betrbgc"
    print "!"
    print "! !USES:"
    print "  use BeTRTracerType , only : betrtracer_type"
    print "  use BetrStatusType , only : betr_status_type"
    print "  use gbetrType      , only : gbetr_type"
    print "  implicit none"
    print "  ! !ARGUMENTS:"
    print "  class("+app_name+"_bgc_reaction_type), intent(inout)    :: this"
    print "  type(bounds_type)                , intent(in)    :: bounds"
    print "  integer                          , intent(in)    :: lbj, ubj"
    print "  type(BeTRtracer_type )           , intent(inout) :: betrtracer_vars"
    print "  character(len=*)                 , intent(in)    :: namelist_buffer"
    print "  type(betr_status_type)           , intent(out)   :: bstatus"
    print "  !LOCAL VARIABLES"
    print "  character(len=*), parameter                      :: subname ='Init_betrbgc'"
    print "  call bstatus%reset()"
    print "  ! remove compiler warnings for unused dummy args"
    print "  if (this%dummy_compiler_warning)           continue"
    print "  if (bounds%begc > 0)                       continue"
    print "  if (ubj > lbj)                             continue"
    print "  if (len(betrtracer_vars%betr_simname) > 0) continue"
    print "  call this%InitAllocate(bounds, lbj, ubj, bstatus)"
    print "  if(bstatus%check_status())return"
    print "  this%use_c13 = "+app_name+"_para%use_c13"
    print "  this%use_c14 = "+app_name+"_para%use_c14"
    print "  this%nop_limit="+app_name+"_para%nop_limit"
    print "  this%non_limit="+app_name+"_para%non_limit"
    print "  call this%set_tracer(betrtracer_vars, bstatus)"
    print "  end subroutine Init_betrbgc"
    print "!-------------------------------------------------------------------------------"
    print "  subroutine set_tracer(this, betrtracer_vars, bstatus)"
    print "  use BeTRTracerType  , only : betrtracer_type"
    print "  use MathfuncMod     , only : addone"
    print "  implicit none"
    print "  class("+app_name+"_bgc_reaction_type), intent(inout)    :: this"
    print "  type(BeTRtracer_type )           , intent(inout) :: betrtracer_vars"
    print "  type(betr_status_type)           , intent(out)   :: bstatus"
    print "  !LOCAL VARIABLES"
    print "  integer :: itemp_gwm"
    print "  integer :: itemp_g"
    print "  integer :: itemp_s"
    print "  integer :: itemp_gwm_grp"
    print "  integer :: dum, itemp"
    print "  integer :: itemp_grp, itemp_v, itemp_vgrp, itemp_trc, itemp_ads, itemp_ads_grp"
    print "  integer :: litr_cnt, wood_cnt, Bm_cnt, trcid, itemp_mem, ngroupmems"
    print ""
    print "  associate(                           &"
    print "    nelm    => this%"+app_name+"_bgc_index%nelms,   &"
    print "    c_loc   => this%"+app_name+"_bgc_index%c_loc,   &"
    print "    c13_loc => this%"+app_name+"_bgc_index%c13_loc, &"
    print "    c14_loc => this%"+app_name+"_bgc_index%c14_loc, &"
    print "    e_loc   => this%"+app_name+"_bgc_index%e_loc    &"
    print "  )"
    print "  itemp_gwm     = 0"
    print "  itemp_g       = 0"
    print "  itemp_s       = 0"
    print "  itemp_gwm_grp = 0"
    print "  itemp_ads_grp =0!counter for sorptive groups"
    print "  itemp_ads=0     !counter for sorptive tracers"
    print "  !volatile tracers"
    print "  itemp = 0; itemp_trc=0"
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_n2, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_n2,&"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_n2, is_trc_gw=.true., is_trc_volatile = .true.)"
    print ""
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_o2, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_o2, &"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_o2, &"
    print "    is_trc_gw=.true., is_trc_volatile = .true.)"
    print ""
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_ar,&"
    print "    trc_grp_beg= betrtracer_vars%id_trc_beg_ar, &"
    print "    trc_grp_end= betrtracer_vars%id_trc_end_ar, &"
    print "    is_trc_gw=.true., is_trc_volatile = .true.)"
    print ""
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_co2x, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_co2x, &"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_co2x, &"
    print "    is_trc_gw=.true., is_trc_volatile = .true.)"
    print ""
    print "  if(this%use_c13)then"
    print "    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &"
    print "      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_c13_co2x, &"
    print "      trc_grp_beg=betrtracer_vars%id_trc_beg_c13_co2x, &"
    print "      trc_grp_end=betrtracer_vars%id_trc_end_c13_co2x, &"
    print "      is_trc_gw=.true., is_trc_volatile = .true.)"
    print "  endif"
    print "  if(this%use_c14)then"
    print "    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &"
    print "      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_c14_co2x, &"
    print "      trc_grp_beg=betrtracer_vars%id_trc_beg_c14_co2x, &"
    print "      trc_grp_end=betrtracer_vars%id_trc_end_c14_co2x, &"
    print "      is_trc_gw=.true., is_trc_volatile = .true.)"
    print "  endif"
    print ""
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &"
    print "    trc_cnt=itemp_trc, trc_grp= betrtracer_vars%id_trc_ch4, &"
    print "    trc_grp_beg= betrtracer_vars%id_trc_beg_ch4, &"
    print "    trc_grp_end= betrtracer_vars%id_trc_end_ch4, &"
    print "    is_trc_gw=.true., is_trc_volatile = .true.)"
    print ""
    print "  ngroupmems=nelm+1  !dom, element + energy"
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = ngroupmems, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_dom, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_dom, &"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_dom, &"
    print "    is_trc_gw=.true., is_trc_volatile = .false.)"
    print ""
    print "  ngroupmems=nelm+1  !pom, element + energy"
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = ngroupmems, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_pom, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_pom, &"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_pom, &"
    print "    is_trc_passive=.true.)"
    print ''
    print "  !three litter groups"
    print "  ngroupmems = 3*nelm"
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = ngroupmems, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_litr, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_litr, &"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_litr, &"
    print "    is_trc_passive=.true.)"
    print ""
    print "  ngroupmems = nelm"
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), &"
    print "    is_trc_passive=.true., mem = ngroupmems, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_wood, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_wood, &"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_wood)"
    print ""
    print "  ngroupmems = 2*nelm"
    print "  call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),mem = ngroupmems, &"
    print "    trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_Bm, &"
    print "    trc_grp_beg=betrtracer_vars%id_trc_beg_Bm, &"
    print "    trc_grp_end=betrtracer_vars%id_trc_end_Bm, &"
    print "    is_trc_passive=.true.)"
    print ""
    print "  betrtracer_vars%nmem_max               = nelm*3"
    print ""
    print "  call betrtracer_vars%Init()"
    print ""
    print "  itemp_grp = 0    !group id"
    print "  itemp_v = 0      !volatile id"
    print "  itemp_vgrp = 0   !volatile group"
    print ""
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_n2, trc_name='N2'  ,      &"
    print "     is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &"
    print "     trc_group_mem= 1,  is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &"
    print "     trc_volatile_group_id = addone(itemp_vgrp))"
    print "  if(bstatus%check_status())return"
    print ""
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_o2, trc_name='O2'  ,      &"
    print "     is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &"
    print "     trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &"
    print "     trc_volatile_group_id = addone(itemp_vgrp))"
    print "  if(bstatus%check_status())return"
    print ""
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ar, trc_name='AR'  ,      &"
    print "     is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &"
    print "     trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &"
    print "     trc_volatile_group_id = addone(itemp_vgrp))"
    print "  if(bstatus%check_status())return"
    print ""
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_co2x, trc_name='CO2x',    &"
    print "     is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp)  , &"
    print "     trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &"
    print "     trc_volatile_group_id = addone(itemp_vgrp))"
    print "  if(bstatus%check_status())return"
    print ""
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ch4, trc_name='CH4',      &"
    print "     is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &"
    print "     trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &"
    print "     trc_volatile_group_id = addone(itemp_vgrp))"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_c13_co2x, &"
    print "      trc_name='13CO2x', is_trc_mobile=.true., is_trc_advective = .true., &"
    print "      trc_group_id = betrtracer_vars%id_trc_c13_co2x, trc_group_mem = 1, is_trc_volatile=.true., &"
    print "      trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp), &"
    print "      trc_family_name='CO2x')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_c14_co2x, &"
    print "      trc_name='14CO2x', is_trc_mobile=.true., is_trc_advective = .true., &"
    print "      trc_group_id = betrtracer_vars%id_trc_c14_co2x, trc_group_mem = 1, is_trc_volatile=.true., &"
    print "      trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp), &"
    print "      trc_family_name='CO2x')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  itemp_mem=0"
    print "  litr_cnt = 0"
    print "  trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C' ,    &"
    print "     is_trc_mobile=.true., is_trc_advective = .false., &"
    print "     trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "     trc_family_name='LIT1')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C_C13' ,    &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='LIT1')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C_C14' ,    &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='LIT1')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print "  litr_cnt = litr_cnt + 1"
    print ""
    print "  trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C'  ,    &"
    print "     is_trc_mobile=.true., is_trc_advective = .false., &"
    print "     trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "     trc_family_name='LIT2')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C_C13' ,    &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='LIT2')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C_C14' ,    &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='LIT2')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print "  litr_cnt=litr_cnt+1"
    print ""
    print "  trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C' ,    &"
    print "     is_trc_mobile=.true., is_trc_advective = .false., &"
    print "     trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "     trc_family_name='LIT3')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C_C13' ,    &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='LIT3')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print "  if(this%use_c14)then"
    print "    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C_C14' ,    &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='LIT3')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  wood_cnt=0"
    print "  itemp_mem=0"
    print "  !coarse root woody components, equivalent to default cwd"
    print "  trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDC' ,    &"
    print "    is_trc_mobile=.false., is_trc_advective = .false., &"
    print "    trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &"
    print "    trc_family_name='CWD')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDC_C13' ,    &"
    print "      is_trc_mobile=.false., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='CWD')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDC_C14' ,    &"
    print "      is_trc_mobile=.false., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='CWD')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  Bm_cnt=0;itemp_mem=0"
    print "  trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_live', &"
    print "     is_trc_mobile=.true., is_trc_advective = .false., &"
    print "     trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &"
    print "     trc_family_name='MB')"
    print "  f(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_live_C13',&"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &"
    print "      trc_family_name='MB')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_live_C14', &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_Bm,  trc_group_mem = addone(itemp_mem), &"
    print "      trc_family_name='MB')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print "  Bm_cnt=Bm_cnt+1"
    print ""
    print "  trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_dead', &"
    print "    is_trc_mobile=.true., is_trc_advective = .false., &"
    print "    trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &"
    print "    trc_family_name='MB')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_dead_C13',&"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &"
    print "      trc_family_name='MB')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='MB_dead_C14', &"
    print "      is_trc_mobile=.true., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_Bm,  trc_group_mem = addone(itemp_mem), &"
    print "      trc_family_name='MB')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  itemp_mem=0"
    print "  trcid =  betrtracer_vars%id_trc_beg_dom+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                       , &"
    print "     trc_name='DOM_C', is_trc_mobile=.true., is_trc_advective = .true. , &"
    print "     trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &"
    print "     is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &"
    print "     trc_adsorbgroupid=addone(itemp_ads_grp), trc_sorpisotherm='LANGMUIR'            , &"
    print "     is_trc_dom=.true.,trc_family_name='DOM')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  trcid = betrtracer_vars%id_trc_beg_dom+e_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                       , &"
    print "     trc_name='DOM_e', is_trc_mobile=.true., is_trc_advective = .true.  , &"
    print "     trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &"
    print "     is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &"
    print "     trc_adsorbgroupid=itemp_ads_grp, trc_sorpisotherm='LANGMUIR'                    , &"
    print "     is_trc_dom=.true., trc_family_name='DOM')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_dom+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                     , &"
    print "      trc_name='DOM_C13', is_trc_mobile=.true., is_trc_advective = .true. , &"
    print "      trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &"
    print "      is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &"
    print "      trc_adsorbgroupid=itemp_ads_grp,trc_sorpisotherm='LANGMUIR'                     , &"
    print "      is_trc_dom=.true., trc_family_name='DOM')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    trcid=betrtracer_vars%id_trc_beg_dom+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid                     , &"
    print "      trc_name='DOM_C14', is_trc_mobile=.true., is_trc_advective = .true. , &"
    print "      trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &"
    print "      is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads) , &"
    print "      trc_adsorbgroupid=itemp_ads_grp, trc_sorpisotherm='LANGMUIR'                    , &"
    print "      is_trc_dom=.true., trc_family_name='DOM')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  itemp_mem=0"
    print "  trcid =  betrtracer_vars%id_trc_beg_pom+c_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_C' ,    &"
    print "     is_trc_mobile=.false., is_trc_advective = .false., &"
    print "     trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &"
    print "     trc_family_name='POM')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  trcid = betrtracer_vars%id_trc_beg_pom+e_loc-1"
    print "  call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_e' ,    &"
    print "     is_trc_mobile=.false., is_trc_advective = .false., &"
    print "     trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &"
    print "     trc_family_name='POM')"
    print "  if(bstatus%check_status())return"
    print ""
    print "  if(this%use_c13)then"
    print "    trcid = betrtracer_vars%id_trc_beg_pom+c13_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_C13' ,    &"
    print "      is_trc_mobile=.false., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='POM')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  if(this%use_c14)then"
    print "    trcid=betrtracer_vars%id_trc_beg_pom+c14_loc-1"
    print "    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='POM_C14' ,  &"
    print "      is_trc_mobile=.false., is_trc_advective = .false., &"
    print "      trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem= addone(itemp_mem), &"
    print "      trc_family_name='POM')"
    print "    if(bstatus%check_status())return"
    print "  endif"
    print ""
    print "  end associate"
    print "  end subroutine set_tracer"
    print "!-------------------------------------------------------------------------------"
    print "  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, dz_top, betrtracer_vars, &"
    print "   biophysforc, biogeo_flux, tracerboundarycond_vars, betr_status)"
    print "  !"
    print "  ! !DESCRIPTION:"
    print "  ! set up boundary conditions for tracer movement"
    print "  !"
    print "  ! !USES:"
    print "  use betr_ctrl              , only : iulog  => biulog"
    print "  use TracerBoundaryCondType , only : tracerboundarycond_type"
    print "  use bshr_log_mod           , only : errMsg => shr_log_errMsg"
    print "  use BeTRTracerType         , only : betrtracer_type"
    print "  use betr_varcon            , only : rgas => brgas"
    print "  use BeTR_biogeoFluxType    , only : betr_biogeo_flux_type"
    print "  use BetrStatusType         , only : betr_status_type"
    print "  use UnitConverMod          , only : ppm2molv"
    print "  implicit none"
    print "  ! !ARGUMENTS:"
    print "  class("+app_name+"_bgc_reaction_type) , intent(inout)    :: this        !"
    print "  type(bounds_type)                 , intent(in)    :: bounds             !"
    print "  integer                           , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc"
    print "  integer                           , intent(in)    :: filter_soilc(:)            ! column filter_soilc"
    print "  type(betrtracer_type)             , intent(in)    :: betrtracer_vars            !"
    print "  real(r8)                          , intent(in)    :: dz_top(bounds%begc: )      !"
    print "  type(betr_biogeophys_input_type)  , intent(in)    :: biophysforc"
    print "  type(betr_biogeo_flux_type)       , intent(in)    :: biogeo_flux"
    print "  type(tracerboundarycond_type)     , intent(inout) :: tracerboundarycond_vars !"
    print "  type(betr_status_type)            , intent(out)   :: betr_status"
    print ""
    print "  ! !LOCAL VARIABLES:"
    print "  integer            :: fc, c"
    print "  character(len=255) :: subname = 'set_boundary_conditions'"
    print "  real(r8) :: irt   !the inverse of R*T"
    print ""
    print "  call betr_status%reset()"
    print "  SHR_ASSERT_ALL((ubound(dz_top)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__), betr_status)"
    print ""
    print "  ! remove compiler warnings for unused dummy args"
    print "  if (this%dummy_compiler_warning)      continue"
    print "  if (size(biogeo_flux%qflx_adv_col)>0) continue"
    print "  associate(                                                    &"
    print "         groupid    => betrtracer_vars%groupid              , &"
    print "         n2_ppmv    => biophysforc%n2_ppmv_col              , &"
    print "         o2_ppmv    => biophysforc%o2_ppmv_col              , &"
    print "         ar_ppmv    => biophysforc%ar_ppmv_col              , &"
    print "         co2_ppmv   => biophysforc%co2_ppmv_col             , &"
    print "         ch4_ppmv   => biophysforc%ch4_ppmv_col             , &"
    print "         pbot_pa    => biophysforc%forc_pbot_downscaled_col , &"
    print "         tair       => biophysforc%forc_t_downscaled_col    , &"
    print "         id_trc_beg_dom=>betrtracer_vars%id_trc_beg_dom     , &"
    print "         id_trc_end_dom=>betrtracer_vars%id_trc_end_dom       &"
    print "     )"
    print ""
    print "  do fc = 1, num_soilc"
    print "    c = filter_soilc(fc)"
    print "    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,:) = 0._r8"
    print "    tracerboundarycond_vars%bot_concflux_col(c,1,:)  = 0._r8                        !zero flux boundary condition for diffusion"
    print "    tracerboundarycond_vars%condc_toplay_col(c,:)    = 0._r8"
    print "    !set value for specific groups"
    print "    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2)    =ppm2molv(pbot_pa(c), n2_ppmv(c), tair(c))    !mol m-3, contant boundary condition"
    print "    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o2)    =ppm2molv(pbot_pa(c), o2_ppmv(c), tair(c))!mol m-3, contant boundary condition"
    print "    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ar)    =ppm2molv(pbot_pa(c), ar_ppmv(c), tair(c))!mol m-3, contant boundary condition"
    print "    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_co2x)  =ppm2molv(pbot_pa(c), co2_ppmv(c), tair(c))!mol m-3, contant boundary condition"
    print "    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ch4)   =ppm2molv(pbot_pa(c), ch4_ppmv(c), tair(c))!mol m-3, contant boundary condition"
    print "    tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,id_trc_beg_dom:id_trc_end_dom)= 0._r8            !mol m-3, contant boundary condition, as concentration"
    print ""
    print "    tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_n2))           = 2._r8*1.837e-5_r8/dz_top(c)  !m/s surface conductance"
    print "    tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_o2))           = 2._r8*1.713e-5_r8/dz_top(c)  !m/s surface conductance"
    print "    tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ar))           = 2._r8*1.532e-5_r8/dz_top(c)  !m/s surface conductance"
    print "    tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_co2x))         = 2._r8*1.399e-5_r8/dz_top(c)  !m/s surface conductance"
    print "    tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ch4))          = 2._r8*1.808e-5_r8/dz_top(c)  !m/s surface conductance"
    print "  enddo"
    print ""
    print "  end associate"
    print "  end subroutine set_boundary_conditions"
    print "!-------------------------------------------------------------------------------"
    print "  subroutine calc_bgc_reaction(this, bounds, col, lbj, ubj, num_soilc, filter_soilc,               &"
    print "    num_soilp,filter_soilp, jtops, dtime, betrtracer_vars, tracercoeff_vars,  biophysforc, &"
    print "    tracerstate_vars, tracerflux_vars, tracerboundarycond_vars, plant_soilbgc, &"
    print "    biogeo_flux,  betr_status)"
    print "  !"
    print "  ! !DESCRIPTION:"
    print "  ! do bgc reaction"
    print "  !"
    print "  ! !USES:"
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
    implicit none
    !ARGUMENTS
    class("+app_name+"_bgc_reaction_type) , intent(inout) :: this                       !
    type(bounds_type)                 , intent(in)    :: bounds                      ! bounds
    type(betr_column_type)            , intent(in)    :: col
    integer                           , intent(in)    :: num_soilc                   ! number of columns in column filter
    integer                           , intent(in)    :: filter_soilc(:)             ! column filter
    integer                           , intent(in)    :: num_soilp
    integer                           , intent(in)    :: filter_soilp(:)
    integer                           , intent(in)    :: jtops( : )                  ! top index of each column
    integer                           , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
    real(r8)                          , intent(in)    :: dtime                       ! model time step
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars             ! betr configuration information
    type(betr_biogeophys_input_type)  , intent(inout) :: biophysforc
    type(tracercoeff_type)            , intent(inout) :: tracercoeff_vars
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars
    type(tracerflux_type)             , intent(inout) :: tracerflux_vars
    type(tracerboundarycond_type)     , intent(inout) :: tracerboundarycond_vars !
    class(plant_soilbgc_type)         , intent(inout) :: plant_soilbgc
    type(betr_biogeo_flux_type)       , intent(inout) :: biogeo_flux
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

    nstates = this%"+app_name+"_bgc_index%nstvars
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
        call this%"+app_name+"_bgc(c,j)%runbgc(is_surflit, dtime, this%"+app_name+"_forc(c,j), nstates, ystates0, ystatesf, betr_status)

        if(betr_status%check_status())then
          write(laystr,'(I2.2)')j
          betr_status%msg=trim(betr_status%msg)//' lay '//trim(laystr)
          return
        endif

        call this%retrieve_output(c, j, nstates, ystates0, ystatesf, dtime, betrtracer_vars, tracerflux_vars,&
           tracerstate_vars, plant_soilbgc, biogeo_flux)

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
    class("+app_name+"_bgc_reaction_type) , intent(inout)    :: this
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
    if (size(tracercoeff_vars%annsum_counter_col) > 0)        continue

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
    class("+app_name+"_bgc_reaction_type) , intent(inout)    :: this
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
    class("+app_name+"_bgc_reaction_type) , intent(inout)    :: this
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
  class("+app_name+"_bgc_reaction_type) , intent(inout) :: this               !
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
   class("+app_name+"_bgc_reaction_type) , intent(inout) :: this               !
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
   class("+app_name+"_bgc_reaction_type) , intent(inout) :: this      !
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
   class("+app_name+"_bgc_reaction_type) , intent(inout) :: this               !
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

   c_loc=this%"+app_name+"_bgc_index%c_loc
   c13_loc=this%"+app_name+"_bgc_index%c13_loc
   c14_loc=this%"+app_name+"_bgc_index%c14_loc
   nelm =this%"+app_name+"_bgc_index%nelms

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
  use "+app_name+"PlantSoilBGCType    , only : "+app_name+"_plant_soilbgc_type
  use tracer_varcon            , only : catomw, natomw, patomw, fix_ip
  implicit none
  class("+app_name+"_bgc_reaction_type) , intent(inout)    :: this
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
     litr_beg =>  this%"+app_name+"_bgc_index%litr_beg  , &
     litr_end =>  this%"+app_name+"_bgc_index%litr_end  , &
     wood_beg =>  this%"+app_name+"_bgc_index%wood_beg  , &
     wood_end =>  this%"+app_name+"_bgc_index%wood_end  , &
     dom_beg =>  this%"+app_name+"_bgc_index%dom_beg    , &
     dom_end =>  this%"+app_name+"_bgc_index%dom_end    , &
     pom_beg =>  this%"+app_name+"_bgc_index%pom_beg    , &
     pom_end =>  this%"+app_name+"_bgc_index%pom_end    , &
     Bm_beg  =>  this%"+app_name+"_bgc_index%Bm_beg     , &
     Bm_end  =>  this%"+app_name+"_bgc_index%Bm_end     , &
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
        ystatesf(this%"+app_name+"_bgc_index%lid_n2)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_o2)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_ar)

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_co2)

    if(this%use_c13)then
      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_c13_co2)
    endif

    if(this%use_c14)then
      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_c14_co2)
    endif

    tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_ch4)

    !fluxes
    tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_o2) ) = &
         ystatesf(this%"+app_name+"_bgc_index%lid_o2_paere )  - &
         ystates0(this%"+app_name+"_bgc_index%lid_o2_paere)

    if ( betr_spinup_state == 0 ) then
      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2)  ) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_n2_paere)  - &
          ystates0(this%"+app_name+"_bgc_index%lid_n2_paere)

      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ar)  ) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_ar_paere)  - &
          ystates0(this%"+app_name+"_bgc_index%lid_ar_paere)

      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_co2x)) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_co2_paere)  - &
          ystates0(this%"+app_name+"_bgc_index%lid_co2_paere)

      if(this%use_c13)then
        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c13_co2x)) = &
            ystatesf(this%"+app_name+"_bgc_index%lid_c13_co2_paere)  - &
            ystates0(this%"+app_name+"_bgc_index%lid_c13_co2_paere)
      endif

      if(this%use_c14)then
        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c14_co2x)) = &
            ystatesf(this%"+app_name+"_bgc_index%lid_c14_co2_paere)  - &
            ystates0(this%"+app_name+"_bgc_index%lid_c14_co2_paere)
      endif

      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ch4) ) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_ch4_paere)  - &
          ystates0(this%"+app_name+"_bgc_index%lid_ch4_paere)

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
        ystatesf(this%"+app_name+"_bgc_index%lid_n2) - &
        ystates0(this%"+app_name+"_bgc_index%lid_n2)

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_co2x ) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_co2) - &
        ystates0(this%"+app_name+"_bgc_index%lid_co2)

    if(this%use_c13)then
      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c13_co2x ) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_c13_co2) - &
          ystates0(this%"+app_name+"_bgc_index%lid_c13_co2)
    endif

    if(this%use_c14)then
      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c14_co2x ) = &
          ystatesf(this%"+app_name+"_bgc_index%lid_c14_co2) - &
          ystates0(this%"+app_name+"_bgc_index%lid_c14_co2)
    endif

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_o2   ) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_o2) - &
        ystates0(this%"+app_name+"_bgc_index%lid_o2)

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ch4  ) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_ch4) - &
        ystates0(this%"+app_name+"_bgc_index%lid_ch4)

    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ar) = &
        ystatesf(this%"+app_name+"_bgc_index%lid_ar) - &
        ystates0(this%"+app_name+"_bgc_index%lid_ar)

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
        (ystatesf(this%"+app_name+"_bgc_index%lid_co2_hr) - &
        ystates0(this%"+app_name+"_bgc_index%lid_co2_hr))*catomw/dtime
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
  use "+app_name+"PlantSoilBGCType    , only : "+app_name+"_plant_soilbgc_type
  use MathfuncMod              , only : fpmax
  use betr_varcon              , only : grav => bgrav
  implicit none
  class("+app_name+"_bgc_reaction_type) , intent(inout) :: this                       !
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
     litr_beg =>  this%"+app_name+"_bgc_index%litr_beg  , &
     litr_end =>  this%"+app_name+"_bgc_index%litr_end  , &
     wood_beg =>  this%"+app_name+"_bgc_index%wood_beg  , &
     wood_end =>  this%"+app_name+"_bgc_index%wood_end  , &
     dom_beg =>  this%"+app_name+"_bgc_index%dom_beg    , &
     dom_end =>  this%"+app_name+"_bgc_index%dom_end    , &
     pom_beg =>  this%"+app_name+"_bgc_index%pom_beg    , &
     pom_end =>  this%"+app_name+"_bgc_index%pom_end    , &
     Bm_beg  =>  this%"+app_name+"_bgc_index%Bm_beg     , &
     Bm_end  =>  this%"+app_name+"_bgc_index%Bm_end       &
  )
  call betr_status%reset()
  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__),betr_status)

  do j = lbj, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      if(j<jtops(c))cycle
      this%"+app_name+"_forc(c,j)%plant_ntypes = this%nactpft
      this%"+app_name+"_forc(c,j)%ystates(:) = 0._r8

      !litter
      this%"+app_name+"_forc(c,j)%ystates(litr_beg:litr_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_litr:betrtracer_vars%id_trc_end_litr)

      !wood
      this%"+app_name+"_forc(c,j)%ystates(wood_beg:wood_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_wood:betrtracer_vars%id_trc_end_wood)

      !dom
      this%"+app_name+"_forc(c,j)%ystates(dom_beg:dom_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_dom:betrtracer_vars%id_trc_end_dom)
      if(this%"+app_name+"_forc(c,j)%ystates(dom_beg)<=tiny_cval)this%"+app_name+"_forc(c,j)%ystates(dom_beg:dom_end)=0._r8

      !pom
      this%"+app_name+"_forc(c,j)%ystates(pom_beg:pom_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_pom:betrtracer_vars%id_trc_end_pom)
      if(this%"+app_name+"_forc(c,j)%ystates(pom_beg)<=tiny_cval)this%"+app_name+"_forc(c,j)%ystates(pom_beg:pom_end)=0._r8


      !microbial biomass
      this%"+app_name+"_forc(c,j)%ystates(Bm_beg:Bm_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_Bm:betrtracer_vars%id_trc_end_Bm)
      if(this%"+app_name+"_forc(c,j)%ystates(Bm_beg)<=tiny_cval)this%"+app_name+"_forc(c,j)%ystates(Bm_beg:Bm_end)=0._r8
      !non-soluble phase of mineral p

      this%"+app_name+"_forc(c,j)%ystates(this%"+app_name+"_bgc_index%lid_n2) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2))

      this%"+app_name+"_forc(c,j)%ystates(this%"+app_name+"_bgc_index%lid_o2) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2))

      this%"+app_name+"_forc(c,j)%ystates(this%"+app_name+"_bgc_index%lid_ar) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar))

      this%"+app_name+"_forc(c,j)%ystates(this%"+app_name+"_bgc_index%lid_co2)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x))

      if(this%use_c13)then
        this%"+app_name+"_forc(c,j)%ystates(this%"+app_name+"_bgc_index%lid_c13_co2)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%"+app_name+"_forc(c,j)%ystates(this%"+app_name+"_bgc_index%lid_c14_co2)= &
          fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x))
      endif

      this%"+app_name+"_forc(c,j)%ystates(this%"+app_name+"_bgc_index%lid_ch4)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4))


      !input
      this%"+app_name+"_forc(c,j)%cflx_input_litr_met = biophysforc%c12flx%cflx_input_litr_met_vr_col(c,j)
      this%"+app_name+"_forc(c,j)%cflx_input_litr_cel = biophysforc%c12flx%cflx_input_litr_cel_vr_col(c,j)
      this%"+app_name+"_forc(c,j)%cflx_input_litr_lig = biophysforc%c12flx%cflx_input_litr_lig_vr_col(c,j)
      this%"+app_name+"_forc(c,j)%cflx_input_litr_cwd = biophysforc%c12flx%cflx_input_litr_cwd_vr_col(c,j)
      this%"+app_name+"_forc(c,j)%cflx_input_litr_lwd = biophysforc%c12flx%cflx_input_litr_lwd_vr_col(c,j)
      this%"+app_name+"_forc(c,j)%cflx_input_litr_fwd = biophysforc%c12flx%cflx_input_litr_fwd_vr_col(c,j)

      !environmental variables
      this%"+app_name+"_forc(c,j)%temp   = biophysforc%t_soisno_col(c,j)            !temperature
      this%"+app_name+"_forc(c,j)%depz   = col%z(c,j)            !depth of the soil
      this%"+app_name+"_forc(c,j)%dzsoi  = col%dz(c,j)            !soil thickness
      this%"+app_name+"_forc(c,j)%sucsat  = biophysforc%sucsat_col(c,j)            ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
      this%"+app_name+"_forc(c,j)%soilpsi = max(biophysforc%smp_l_col(c,j)*grav*1.e-6_r8,-15._r8)    ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
      this%"+app_name+"_forc(c,j)%bsw = biophysforc%bsw_col(c,j)
      this%"+app_name+"_forc(c,j)%bd   = biophysforc%bd_col(c,j)              !bulk density
      this%"+app_name+"_forc(c,j)%pct_sand = biophysforc%cellsand_col(c,j)
      this%"+app_name+"_forc(c,j)%pct_clay = biophysforc%cellclay_col(c,j)
      this%"+app_name+"_forc(c,j)%h2osoi_vol = biophysforc%h2osoi_vol_col(c,j)
      this%"+app_name+"_forc(c,j)%h2osoi_liq = biophysforc%h2osoi_liq_col(c,j)
      this%"+app_name+"_forc(c,j)%air_vol = biophysforc%air_vol_col(c,j)
      this%"+app_name+"_forc(c,j)%finundated = biophysforc%finundated_col(c)
      this%"+app_name+"_forc(c,j)%watsat = biophysforc%watsat_col(c,j)
      this%"+app_name+"_forc(c,j)%watfc = biophysforc%watfc_col(c,j)
      this%"+app_name+"_forc(c,j)%cellorg = biophysforc%cellorg_col(c,j)
      this%"+app_name+"_forc(c,j)%pH = biophysforc%soil_pH(c,j)

      !conductivity for plant-aided gas transport
      this%"+app_name+"_forc(c,j)%aren_cond_n2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2))
      this%"+app_name+"_forc(c,j)%aren_cond_o2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2))
      this%"+app_name+"_forc(c,j)%aren_cond_co2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%"+app_name+"_forc(c,j)%aren_cond_co2_c13 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c13_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%"+app_name+"_forc(c,j)%aren_cond_co2_c14 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c14_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%"+app_name+"_forc(c,j)%aren_cond_ar = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar))
      this%"+app_name+"_forc(c,j)%aren_cond_ch4 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4))
      !phase conversion parameter
      this%"+app_name+"_forc(c,j)%ch4_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4))
      this%"+app_name+"_forc(c,j)%co2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x))
      this%"+app_name+"_forc(c,j)%o2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2))
      this%"+app_name+"_forc(c,j)%n2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2))
      this%"+app_name+"_forc(c,j)%ar_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar))
      this%"+app_name+"_forc(c,j)%o2_w2b = &
          tracercoeff_vars%aqu2bulkcef_mobile_col(c,j,betrtracer_vars%groupid(betrtracer_vars%id_trc_o2))

      !atmospheric pressure (mol/m3) for gas ventilation.
      this%"+app_name+"_forc(c,j)%conc_atm_n2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_n2))
      this%"+app_name+"_forc(c,j)%conc_atm_o2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_o2))
      this%"+app_name+"_forc(c,j)%conc_atm_ar = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_ar))
      this%"+app_name+"_forc(c,j)%conc_atm_co2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%"+app_name+"_forc(c,j)%conc_atm_co2_c13 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%"+app_name+"_forc(c,j)%conc_atm_co2_c14 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%"+app_name+"_forc(c,j)%conc_atm_ch4 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_ch4))

      this%"+app_name+"_forc(c,j)%soilorder = biophysforc%isoilorder(c)

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
  class("+app_name+"_bgc_reaction_type) , intent(inout) :: this                       !
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
    Kaff_CM              => "+app_name+"_para%Kaff_CM                       , &
    Kaff_EM              => "+app_name+"_para%Kaff_EM                       , &
    Kaff_BC              => "+app_name+"_para%Kaff_BC                       , &
    alpha_B2E            => "+app_name+"_para%alpha_B2E                     , &
    alpha_B2T            => "+app_name+"_para%alpha_B2T                     , &
    nelms                =>  this%"+app_name+"_bgc_index%nelms                &
  )
  c_l=1
  do j = 1, ubj
    KM_CM=aqu2bulkcef_mobile(c_l,j,id_trc_dom)*this%"+app_name+"_forc(c_l,j)%KM_OM_ref*Kaff_CM
    BMT=(tracer_conc_mobile(c_l,j,trcid_Bm)+tracer_conc_mobile(c_l,j,trcid_Bm+nelms))*alpha_B2T
    Msurf=this%"+app_name+"_forc(c_l,j)%Msurf_OM-tracer_conc_mobile(c_l,j,trcid_pom)
    denorm0=1._r8+Msurf/KM_CM+BMT/Kaff_BC
    denorm1=denorm0+tracer_conc_mobile(c_l,j,trcid_dom)/KM_CM
    denorm2=denorm0+tracer_conc_mobile(c_l,j,trcid_dom)/Kaff_BC
    beta=1._r8/(1._r8-Msurf/KM_CM/denorm1-BMT/Kaff_BC/denorm2)
    aqu2bulkcef_mobile(c_l,j,id_trc_dom) = aqu2bulkcef_mobile(c_l,j,id_trc_dom)*beta
    !print*,'c,j',1._r8/(1._r8-Msurf/KM_CM/denorm1-BMT/Kaff_BC/denorm2)
  enddo
  end associate
  end subroutine update_sorpphase_coeff

end module "+app_name+"BGCReactionsType

    print "create file "+sfarm_dir+'/'+app_name+'Nlayer/'+app_name+"PlantSoilBGCType.F90"
