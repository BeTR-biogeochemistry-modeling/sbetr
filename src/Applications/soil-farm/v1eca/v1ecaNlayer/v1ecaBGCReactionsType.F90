module v1ecaBGCReactionsType

#include "bshr_assert.h"
#include "bshr_assign.h"
  !
  ! !DESCRIPTION:
  ! Do ECA based nitrogen/phosphorus competition within betr.
  ! This code uses the operator automated down-regulation scheme
  ! Ideally, all belowground biogeochemistry processes should be solved
  ! simultaneously using the ODE solver below. Yet because of the potential
  ! conflict of interest (in coding) with others, the belowground BGC
  ! considers the nutrient interaction between decomposers, nitrifiers, denitrifiers
  ! and plants (uptake). The P cycle does not include P demand from processes
  ! other than aerobic decomposition and plant growth, therefore nitrifiers and denitrifiers
  ! are never P limited.
  !
  ! Also, because I'm solving enzymatic P extraction simultaneously with decomposition
  ! and plant P uptake, each OM pool (execpt CWD) is assigned a targeting CP ratio to
  ! impose the P limitation of decomposition. This treatment is equivalent to assume
  ! the decomposers are having fixed stoichiometry. In contrast, other implementations
  ! in ACME LAND treats enzymatic P extraction and decomposition as two separate processes,
  ! which causes another ordering ambiguity, in that if one switches the decomposition
  ! and P extraction, the model will potentially make very different predictions.
  !
  ! Further, because I am solving the inorganic P dynamics using the ECA formulation
  ! the labile P pool is implicitly represented. Also, it is assumed the secondary pool
  ! are competing for adsorption space with the labile P, so there is an adsoprtion
  ! saturation effect, which is apparently missing form other ACME implementations.
  !
  ! HISTORY:
  ! Created by Jinyun Tang, Nov 20th, 2015
  ! Note: ECA parameters are note tuned.
  !
  ! !USES:
  !
  use bshr_assert_mod, only : shr_assert
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use bshr_assert_mod, only : shr_assert_any
  use bshr_log_mod          , only : errMsg => shr_log_errMsg
  use bshr_kind_mod         , only : r8 => shr_kind_r8
  use bshr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use BeTR_decompMod        , only : bounds_type  => betr_bounds_type
  use BGCReactionsMod       , only : bgc_reaction_type
  use betr_varcon           , only : spval => bspval, ispval => bispval
  use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
  use tracer_varcon         , only : natomw, patomw, catomw,c13atomw,c14atomw
  use v1ecaBGCType         , only : v1eca_bgc_type
  use JarBgcForcType        , only : JarBGC_forc_type
  use BetrStatusType        , only : betr_status_type
  use v1ecaBGCIndexType    , only : v1eca_bgc_index_type
  use v1ecaParaType        , only : v1eca_para
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BeTR_decompMod    , only : betr_bounds_type
  use BeTRTracerType              , only : betrtracer_type
  use tracerstatetype        , only : tracerstate_type
  implicit none

  save
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! !PUBLIC TYPES:

  logical :: ldebug
  !integer, private :: lpr
  type, public, extends(bgc_reaction_type) :: &
    v1eca_bgc_reaction_type
     private
    integer :: parcol
    type(v1eca_bgc_type), pointer :: v1eca(:,:)
    type(JarBGC_forc_type), pointer :: v1eca_forc(:,:)
    type(v1eca_bgc_index_type) :: v1eca_bgc_index
    logical :: use_c13
    logical :: use_c14
    logical :: nop_limit
    logical :: non_limit
    integer :: nactpft               ! number of active pfts
  contains
    procedure :: Init_betrbgc                 ! initialize betr bgc
    procedure :: set_boundary_conditions      ! set top/bottom boundary conditions for various tracers
    procedure :: calc_bgc_reaction            ! doing bgc calculation
    procedure :: init_boundary_condition_type ! initialize type of top boundary conditions
    procedure :: do_tracer_equilibration      ! do equilibrium tracer chemistry
    procedure :: initCold
    procedure :: readParams
    procedure :: retrieve_biogeoflux
    procedure :: set_kinetics_par
    procedure :: retrieve_lnd2atm
    procedure :: retrieve_biostates
    procedure :: debug_info
    procedure :: set_bgc_spinup
    procedure :: UpdateParas
    procedure :: init_iP_prof
    procedure :: reset_biostates
    procedure :: SetParCols
    procedure, private :: set_bgc_forc
    procedure, private :: retrieve_output
    procedure, private :: precision_filter
    procedure, private :: update_atm_conc
  end type v1eca_bgc_reaction_type

  interface v1eca_bgc_reaction_type
     module procedure constructor
  end interface v1eca_bgc_reaction_type

contains

  !-------------------------------------------------------------------------------
  type(v1eca_bgc_reaction_type) function constructor()
    !
    ! ! DESCRIPTION:
    !
    ! create an object of type v1eca_bgc_reaction_type.
    ! Right now it is purposely empty
   type(v1eca_bgc_reaction_type), allocatable :: bgc

   allocate(bgc)
   constructor = bgc
  end function constructor

  !-------------------------------------------------------------------------------
  subroutine UpdateParas(this, bounds, lbj, ubj, bstatus)
  implicit none
  class(v1eca_bgc_reaction_type), intent(inout) :: this
  type(bounds_type)                    , intent(in)    :: bounds
  integer                              , intent(in)    :: lbj, ubj        ! lower and upper bounds, make sure they are > 0
  type(betr_status_type)               , intent(out)   :: bstatus
  integer :: c, j
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      call this%v1eca(c,j)%UpdateParas(v1eca_para, bstatus)
      if(bstatus%check_status())return
    enddo
  enddo

  end subroutine UpdateParas
  !-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
    !
    ! DESCRIPTION:
    ! initialize boundary condition types
    ! !USES:
    use TracerBoundaryCondType      , only : tracerboundarycond_type


    ! !ARGUMENTS:
    class(v1eca_bgc_reaction_type), intent(inout) :: this
    type(bounds_type)               , intent(in) :: bounds
    type(BeTRtracer_type )          ,  intent(in) :: betrtracer_vars
    type(tracerboundarycond_type)   ,  intent(in) :: tracerboundarycond_vars

    ! !LOCAL VARIABLES:
    integer :: c

    associate(                               &
         groupid  => betrtracer_vars%groupid &
         )
    if (this%dummy_compiler_warning) continue
    tracerboundarycond_vars%topbc_type(1:betrtracer_vars%ngwmobile_tracer_groups) = bndcond_as_conc
    tracerboundarycond_vars%topbc_type(groupid(betrtracer_vars%id_trc_no3x)) = bndcond_as_flux
    tracerboundarycond_vars%topbc_type(groupid(betrtracer_vars%id_trc_p_sol)) = bndcond_as_flux

    tracerboundarycond_vars%topbc_type(betrtracer_vars%ngwmobile_tracer_groups+1:betrtracer_vars%ntracer_groups) = bndcond_as_flux

    end associate
  end subroutine init_boundary_condition_type

  !-------------------------------------------------------------------------------
  subroutine set_bgc_spinup(this, bounds, lbj, ubj, biophysforc, &
    tracers, tracerstate_vars)

  !
  !DESCRIPTION
  !set initial conditions for regular or spinup runs. It makes two assumptions
  !1. the initial conditions are defined
  !2. spinup scalar was defiend with sufficient temporal average.

  use BeTRTracerType         , only : betrtracer_type
  use betr_ctrl              , only : exit_spinup, enter_spinup, betr_spinup_state
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use tracer_varcon          , only : patomw
  implicit none
  class(v1eca_bgc_reaction_type) , intent(inout)    :: this      !
  type(bounds_type)                       , intent(in) :: bounds
  integer                                 , intent(in) :: lbj, ubj
  type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc
  type(BeTRtracer_type)                   , intent(inout) :: tracers
  type(tracerstate_type)                  , intent(inout) :: tracerstate_vars
  integer :: kk, c, j, c_l
  real(r8):: latacc

  return
  end subroutine set_bgc_spinup

  !----------------------------------------------------------------------
  subroutine init_iP_prof(this, bounds, lbj, ubj, biophysforc, tracers, tracerstate_vars)
  !
  !DESCRIPTION
  ! set up initial inorganic P profile
  use tracer_varcon, only : patomw
  use BeTR_biogeophysInputType         , only : betr_biogeophys_input_type
  use BeTRTracerType                   , only : betrtracer_type
  use tracerstatetype                  , only : tracerstate_type
  implicit none
  ! !ARGUMENTS:
  class(v1eca_bgc_reaction_type)  , intent(inout)    :: this
  type(bounds_type)                        , intent(in) :: bounds
  integer                                  , intent(in) :: lbj, ubj
  type(betr_biogeophys_input_type)        , intent(inout) :: biophysforc
  type(BeTRtracer_type)                    , intent(inout) :: tracers
  type(tracerstate_type)                   , intent(inout) :: tracerstate_vars

  integer :: c, j

  return
  end subroutine init_iP_prof
  !----------------------------------------------------------------------
  subroutine set_kinetics_par(this, lbj, ubj, nactpft, plantNutkinetics, tracers, tracercoeff_vars)
  use PlantNutKineticsMod, only : PlantNutKinetics_type
  use tracercoeffType          , only : tracercoeff_type
  use BeTRTracerType           , only : betrtracer_type
  use tracer_varcon            , only : lbcalib
  ! !ARGUMENTS:
  class(v1eca_bgc_reaction_type)         , intent(inout)    :: this                       !
  class(PlantNutKinetics_type), intent(in) :: plantNutkinetics
  type(betrtracer_type)       , intent(in) :: tracers
  type(tracercoeff_type), intent(inout) :: tracercoeff_vars
  integer, intent(in) :: lbj, ubj
  integer, intent(in) :: nactpft  !number of active pfts

  integer :: c_l, p, j, trcid, gid

  associate(                                                      &
     k_sorbsurf    => tracercoeff_vars%k_sorbsurf_col           , &
     Q_sorbsurf    => tracercoeff_vars%Q_sorbsurf_col           , &
     tracer_group_memid => tracers%tracer_group_memid           , &
     id_trc_nh3x   => tracers%id_trc_nh3x                       , &
     id_trc_p_sol  => tracers%id_trc_p_sol                      , &
     adsorbgroupid => tracers%adsorbgroupid                       &
  )

  !in the following, only one column is assumed for the bgc
  c_l = 1
  this%nactpft = nactpft
  do j = lbj, ubj
    do p = 1, nactpft
      this%v1eca(c_l,j)%competECA%vmax_minn_nh4_plant(p) = plantNutkinetics%plant_nh4_vmax_vr_patch(p,j)
      this%v1eca(c_l,j)%competECA%vmax_minn_no3_plant(p) = plantNutkinetics%plant_no3_vmax_vr_patch(p,j)
      this%v1eca(c_l,j)%competECA%vmax_minp_plant(p) = plantNutkinetics%plant_p_vmax_vr_patch(p,j)
      this%v1eca(c_l,j)%competECA%kaff_minn_no3_plant(p)= plantNutkinetics%plant_no3_km_vr_patch(p,j)
      this%v1eca(c_l,j)%competECA%kaff_minn_nh4_plant(p)= plantNutkinetics%plant_nh4_km_vr_patch(p,j)
      this%v1eca(c_l,j)%competECA%kaff_minp_plant(p)   = plantNutkinetics%plant_p_km_vr_patch(p,j)
      this%v1eca(c_l,j)%competECA%plant_froot_nn(p) = plantNutkinetics%plant_eff_ncompet_b_vr_patch(p,j)
      this%v1eca(c_l,j)%competECA%plant_froot_np(p) = plantNutkinetics%plant_eff_pcompet_b_vr_patch(p,j)
    enddo
    if(lbcalib)then
      this%v1eca(c_l,j)%competECA%kaff_minn_nh4_mic =plantNutkinetics%km_decomp_nh4_vr_col(c_l,j)
      this%v1eca(c_l,j)%competECA%kaff_minn_no3_mic =plantNutkinetics%km_decomp_no3_vr_col(c_l,j)
      this%v1eca(c_l,j)%competECA%kaff_minp_mic = plantNutkinetics%km_decomp_p_vr_col(c_l,j)
      this%v1eca(c_l,j)%competECA%kaff_minn_nh4_nit = plantNutkinetics%km_nit_nh4_vr_col(c_l,j)
      this%v1eca(c_l,j)%competECA%kaff_minn_no3_den = plantNutkinetics%km_den_no3_vr_col(c_l,j)
    endif
    !mineral surfaces
    this%v1eca(c_l,j)%competECA%kaff_minn_nh4_msurf= plantNutkinetics%km_minsurf_nh4_vr_col(c_l,j)   !this is ignored at this moment
    this%v1eca(c_l,j)%competECA%kaff_minp_msurf= plantNutkinetics%km_minsurf_p_vr_col(c_l,j)

    !effective n competing decomposers
    this%v1eca(c_l,j)%competECA%compet_bn_mic=plantNutkinetics%decomp_eff_ncompet_b_vr_col(c_l,j)
    !effective p competing decomposers
    this%v1eca(c_l,j)%competECA%compet_bp_mic=plantNutkinetics%decomp_eff_pcompet_b_vr_col(c_l,j)

    this%v1eca(c_l,j)%competECA%compet_bn_nit=plantNutkinetics%nit_eff_ncompet_b_vr_col(c_l,j)
    this%v1eca(c_l,j)%competECA%compet_bn_den=plantNutkinetics%den_eff_ncompet_b_vr_col(c_l,j)
    this%v1eca(c_l,j)%competECA%vmax_minsurf_p = plantNutkinetics%vmax_minsurf_p_vr_col(c_l,j)
    this%v1eca_forc(c_l,j)%msurf_nh4 = plantNutkinetics%minsurf_nh4_compet_vr_col(c_l,j)   !this  number needs update
    this%v1eca_forc(c_l,j)%msurf_minp= plantNutkinetics%minsurf_p_compet_vr_col(c_l,j)    !this  number needs update
    this%v1eca(c_l,j)%competECA%dsolutionp_dt = plantNutkinetics%dsolutionp_dt_vr_col(c_l,j)
    this%v1eca(c_l,j)%competECA%dlabp_dt = plantNutkinetics%dlabp_dt_vr_col(c_l,j)
    !be carefule about the following lines
    trcid=tracer_group_memid(id_trc_p_sol,1)
    !gid = adsorbgroupid(trcid)
    !k_sorbsurf(c_l,j,gid) = plantNutkinetics%km_minsurf_p_vr_col(c_l,j)
    !Q_sorbsurf(c_l,j,gid) = plantNutkinetics%minsurf_p_compet_vr_col(c_l,j)

!    trcid=tracer_group_memid(id_trc_nh3x,1); gid = adsorbgroupid(trcid)
!    k_sorbsurf(c_l,j,gid) = plantNutkinetics%km_minsurf_nh4_vr_col(c_l,j)
!    Q_sorbsurf(c_l,j,gid) = plantNutkinetics%minsurf_nh4_compet_vr_col(c_l,j)

  enddo
  end associate
  end subroutine set_kinetics_par
  !-------------------------------------------------------------------------------

  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars, namelist_buffer, bstatus)
    !
    ! DESCRIPTION:
    ! initialize the betrbgc
    !
    ! !USES:
    use BeTRTracerType                   , only : betrtracer_type
    use MathfuncMod                      , only : addone
    use betr_varcon                      , only : betr_maxpatch_pft
    use betr_constants                   , only : betr_namelist_buffer_size_ext
    use tracer_varcon                    , only : fix_ip
    implicit none
    ! !ARGUMENTS:
    class(v1eca_bgc_reaction_type) , intent(inout)    :: this
    type(bounds_type)                    , intent(in)    :: bounds
    integer                              , intent(in)    :: lbj, ubj        ! lower and upper bounds, make sure they are > 0
    type(BeTRtracer_type )               , intent(inout) :: betrtracer_vars !
    character(len=*)                     , intent(in) :: namelist_buffer
    type(betr_status_type)               , intent(out)   :: bstatus

    ! !LOCAL VARIABLES:
    character(len=32), parameter                         :: subname ='Init_betrbgc'
    integer   :: jj
    integer   :: nelm, itemp_mem
    integer   :: itemp, itemp_vgrp, itemp_v,itemp_trc
    integer   :: c_loc, n_loc, p_loc, trcid, c13_loc, c14_loc
    integer   :: c, j, litr_cnt, wood_cnt, Bm_cnt, pom_cnt, som_cnt, itemp_ads,itemp_ads_grp
    integer   :: ngroupmems
    logical   :: batch_mode
    logical   :: som_move

    call bstatus%reset()
    this%parcol=1
    batch_mode = .false.
    som_move=.false.
    this%nactpft = 0

    call this%v1eca_bgc_index%Init(v1eca_para%use_c13, v1eca_para%use_c14, &
       v1eca_para%non_limit, v1eca_para%nop_limit, betr_maxpatch_pft)

    if(bstatus%check_status())return

    !create the models
    allocate(this%v1eca(bounds%begc:bounds%endc,lbj:ubj))

    !create model specific forcing data structure
    allocate(this%v1eca_forc(bounds%begc:bounds%endc,lbj:ubj))

    !initialize
    do j = lbj, ubj
      do c = bounds%begc, bounds%endc
        call this%v1eca(c,j)%InitVL(v1eca_para, batch_mode, this%v1eca_bgc_index, bstatus)
        if(bstatus%check_status())return

        call this%v1eca_forc(c,j)%Init(this%v1eca_bgc_index%nstvars, this%v1eca_bgc_index%nom_pools)
      enddo
    enddo
    this%use_c13 = v1eca_para%use_c13
    this%use_c14 = v1eca_para%use_c14
    this%nop_limit=v1eca_para%nop_limit
    this%non_limit=v1eca_para%non_limit

    !set up betr
    nelm =this%v1eca_bgc_index%nelms
    c_loc=this%v1eca_bgc_index%c_loc
    n_loc=this%v1eca_bgc_index%n_loc
    p_loc=this%v1eca_bgc_index%p_loc
    c13_loc=this%v1eca_bgc_index%c13_loc
    c14_loc=this%v1eca_bgc_index%c14_loc
    !volatile tracers
    itemp = 0; itemp_trc=0

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_n2, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_n2, &
      trc_grp_end=betrtracer_vars%id_trc_end_n2, &
      is_trc_gw=.true., is_trc_volatile = .true.)

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
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_n2o, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_n2o, &
      trc_grp_end=betrtracer_vars%id_trc_end_n2o, &
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

    !dissolved nh3x, no volatilization is allowed at this moment.
!    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
!      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_nh3x, &
!      trc_grp_beg=betrtracer_vars%id_trc_beg_nh3x, &
!      trc_grp_end=betrtracer_vars%id_trc_end_nh3x, &
!      is_trc_gw=.true., is_trc_volatile = .true., is_trc_adsorb=.true.)

    !non-volatile tracers
    !nitrate
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),  mem = 1, &
      trc_cnt = itemp_trc, trc_grp=betrtracer_vars%id_trc_no3x, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_no3x, &
      trc_grp_end=betrtracer_vars%id_trc_end_no3x, &
      is_trc_gw=.true., is_trc_volatile = .false.)

    !soluble phosphate
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_p_sol, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_p_sol, &
      trc_grp_end=betrtracer_vars%id_trc_end_p_sol, &
      is_trc_gw=.true., is_trc_volatile = .false.,is_trc_adsorb=.true.)

    !three litter groups
    ngroupmems = 3*nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_litr, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_litr, &
      trc_grp_end=betrtracer_vars%id_trc_end_litr, &
      is_trc_passive=.true.)

    !three woody groups
    ngroupmems = nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), &
      is_trc_passive=.true., mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_wood, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_wood, &
      trc_grp_end=betrtracer_vars%id_trc_end_wood)

    !group of microbial biomass
    ngroupmems = nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_Bm, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_Bm, &
      trc_grp_end=betrtracer_vars%id_trc_end_Bm, &
      is_trc_passive=.true.)

    !group of particulate organic matter
    ngroupmems = nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_pom, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_pom, &
      trc_grp_end=betrtracer_vars%id_trc_end_pom, &
      is_trc_passive=.true.)

    !group of som, which is not dom or microbial biomass
    ngroupmems = nelm
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),mem = ngroupmems, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_som, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_som, &
      trc_grp_end=betrtracer_vars%id_trc_end_som, &
      is_trc_passive=.true.)

    betrtracer_vars%nmem_max                     = nelm*3                       ! maximum number of group elements

    call betrtracer_vars%Init()
    betrtracer_vars%is_mobile(:) = .true.

    !-------------------------------------------------------------------------------
    !set up the tracers
    itemp_vgrp = 0  !counter for volatile groups
    itemp_v    = 0  !counter for volatile tracers
    itemp_ads_grp =0
    itemp_ads=0     !counter for sorptive tracers

    call betrtracer_vars%set_tracer(bstatus=bstatus, trc_id = betrtracer_vars%id_trc_n2, &
         trc_name='N2', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_n2, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp))
    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_o2, &
         trc_name='O2', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_o2, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp))
    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ar, &
         trc_name='AR', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_ar, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp))
    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_co2x, &
         trc_name='CO2x', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_co2x, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp))
    if(bstatus%check_status())return


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

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ch4, &
         trc_name='CH4', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_ch4, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v),  trc_volatile_group_id = addone(itemp_vgrp))
    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_n2o, &
         trc_name='N2O' , is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_n2o, trc_group_mem = 1, is_trc_volatile=.true., &
         trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp))
    if(bstatus%check_status())return

!    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_nh3x, &
!         trc_name='NH3x', is_trc_mobile=.true., is_trc_advective = .true., &
!         trc_group_id = betrtracer_vars%id_trc_nh3x, trc_group_mem = 1, is_trc_volatile=.true., &
!         trc_volatile_id = addone(itemp_v), trc_volatile_group_id = addone(itemp_vgrp), &
!         is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads), &
!         trc_adsorbgroupid=addone(itemp_ads_grp), trc_sorpisotherm='LANGMUIR', &
!         trc_vtrans_scal=1._r8)
!    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_no3x, &
         trc_name='NO3x', is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_no3x, trc_group_mem = 1, is_trc_volatile=.false., &
         trc_vtrans_scal=0._r8)
    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_p_sol, &
         trc_name='P_SOL', is_trc_mobile=.false. .and. (.not. fix_ip), &
         is_trc_advective = .false. .and. (.not. fix_ip), &
         trc_group_id = betrtracer_vars%id_trc_p_sol, trc_group_mem = 1, is_trc_volatile=.false., &
         is_trc_adsorb = .false., trc_vtrans_scal=0._r8)
    if(bstatus%check_status())return

    !------------------------------------------------------------------------------------
    !only one group passive solid litter tracers
    !define litter group
    itemp_mem=0
    litr_cnt = 0
    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1N' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1P' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C_C13' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C_C14' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
      if(bstatus%check_status())return
    endif
    litr_cnt = litr_cnt + 1
    !------------------------------------------------------------------------------------

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C'  ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2N' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2P' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false.,  &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C_C13' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C_C14' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
      if(bstatus%check_status())return
    endif
    litr_cnt=litr_cnt+1
    !------------------------------------------------------------------------------------

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3N' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3P' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
    if(bstatus%check_status())return
    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C_C13' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C_C14' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
      if(bstatus%check_status())return
    endif
    !------------------------------------------------------------------------------------
    !define woody group
    wood_cnt=0
    itemp_mem=0
    !coarse root woody components, equivalent to default cwd
    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDC' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='CWD')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDN' ,    &
         is_trc_mobile=.false., is_trc_advective = .false.,&
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='CWD')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='CWDP' ,    &
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

    !------------------------------------------------------------------------------------
    !define som group
    Bm_cnt=0;itemp_mem=0
    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1C', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1_MB')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1N', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1_MB')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1P', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1_MB')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1C_C13',&
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1_MB')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1C_C14', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_Bm,  trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1_MB')
      if(bstatus%check_status())return
    endif

    !------------------------------------------------------------------------------------
    !define som group
    pom_cnt=0;itemp_mem=0
    trcid = betrtracer_vars%id_trc_beg_pom+pom_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM2C', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM2_POM')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_pom+pom_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM2N', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM2_POM')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_pom+pom_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM2P', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM2_POM')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_pom+pom_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM2C_C13',&
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM2_POM')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_pom+pom_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM2C_C14', &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_pom,  trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM2_POM')
      if(bstatus%check_status())return
    endif
    !------------------------------------------------------------------------------------
    !new group
    som_cnt = 0; itemp_mem=0
    trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3C' ,     &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3N'  ,    &
         is_trc_mobile=som_move, is_trc_advective = .false.,  &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3P' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false.,  &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
    if(bstatus%check_status())return
    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3C_C13' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3C_C14' ,    &
         is_trc_mobile=som_move, is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
      if(bstatus%check_status())return
    endif

    !------------------------------------------------------------------------------------

  end subroutine Init_betrbgc

  !-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, jtops, dz_top, betr_time, &
       betrtracer_vars, biophysforc, biogeo_flux, tracercoeff_vars, tracerboundarycond_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! set up boundary conditions for tracer movement
    !
    ! !USES:
    use TracerBoundaryCondType, only : tracerboundarycond_type
    use BeTRTracerType        , only : betrtracer_type
    use BetrStatusType        , only : betr_status_type
    use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
    use UnitConvertMod         , only : ppm2molv
    use TracerCoeffType        , only : tracercoeff_type
    use BeTR_TimeMod             , only : betr_time_type
    implicit none
    ! !ARGUMENTS:
    class(v1eca_bgc_reaction_type) , intent(inout)    :: this
    type(bounds_type)                    , intent(in)    :: bounds
    integer                              , intent(in)    :: num_soilc               ! number of columns in column filter
    integer                              , intent(in)    :: filter_soilc(:)         ! column filter
    type(betrtracer_type)                , intent(in)    :: betrtracer_vars
    real(r8)                             , intent(in)    :: dz_top(bounds%begc: )
    integer                              , intent(in)    :: jtops(bounds%begc: )
    type(betr_time_type)                 , intent(in)    :: betr_time
    type(betr_biogeophys_input_type)     , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)          , intent(in)    :: biogeo_flux
    type(tracercoeff_type)               , intent(in)   :: tracercoeff_vars
    type(tracerboundarycond_type)        , intent(inout) :: tracerboundarycond_vars !
    type(betr_status_type)               , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'set_boundary_conditions'
    integer :: fc, c, kk

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(dz_top)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__),betr_status)

    if (this%dummy_compiler_warning) continue
    associate(                                       &
         groupid  => betrtracer_vars%groupid    ,    &
         ngwmobile_tracers => betrtracer_vars%ngwmobile_tracers, &
         volatilegroupid =>  betrtracer_vars%volatilegroupid, &
         ntracers => betrtracer_vars%ntracers        , &
         is_volatile=> betrtracer_vars%is_volatile   , &
         n2_ppmv    => biophysforc%n2_ppmv_col       , &
         o2_ppmv    => biophysforc%o2_ppmv_col       , &
         ar_ppmv    => biophysforc%ar_ppmv_col       , &
         co2_ppmv   => biophysforc%co2_ppmv_col      , &
         ch4_ppmv   => biophysforc%ch4_ppmv_col      , &
         n2o_ppmv   => biophysforc%n2o_ppmv_col      , &
         nh3_ppmv   => biophysforc%nh3_ppmv_col      , &
         pbot_pa    => biophysforc%forc_pbot_downscaled_col, &
         tair       => biophysforc%forc_t_downscaled_col , &
         tracer_gwdif_concflux_top_col=> tracerboundarycond_vars%tracer_gwdif_concflux_top_col, &
         condc_toplay_col => tracerboundarycond_vars%condc_toplay_col  ,&
         bot_concflux_col => tracerboundarycond_vars%bot_concflux_col ,&
         id_trc_n2 => betrtracer_vars%id_trc_n2, &
         id_trc_o2 => betrtracer_vars%id_trc_o2, &
         id_trc_ar=> betrtracer_vars%id_trc_ar, &
         id_trc_co2x => betrtracer_vars%id_trc_co2x, &
         id_trc_ch4 => betrtracer_vars%id_trc_ch4, &
         id_trc_n2o => betrtracer_vars%id_trc_n2o, &
         id_trc_nh3x => betrtracer_vars%id_trc_nh3x, &
         id_trc_c13_co2x => betrtracer_vars%id_trc_c13_co2x, &
         id_trc_c14_co2x => betrtracer_vars%id_trc_c14_co2x,  &
         diffblkm_topsoi_col=> tracercoeff_vars%diffblkm_topsoi_col, &
         snowres_col => tracercoeff_vars%snowres_col &
         )
      do fc = 1, num_soilc
         c = filter_soilc(fc)

         !values below will be updated with datastream
         !eventually, the following code will be implemented using polymorphism
         tracer_gwdif_concflux_top_col(c,1:2,1:ntracers)                   =0._r8                        !zero incoming flux
         tracer_gwdif_concflux_top_col(c,1:2,id_trc_n2)    =ppm2molv(pbot_pa(c), n2_ppmv(c), tair(c))    !mol m-3, contant boundary condition
         tracer_gwdif_concflux_top_col(c,1:2,id_trc_o2)    =ppm2molv(pbot_pa(c), o2_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracer_gwdif_concflux_top_col(c,1:2,id_trc_ar)    =ppm2molv(pbot_pa(c), ar_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracer_gwdif_concflux_top_col(c,1:2,id_trc_co2x)  =ppm2molv(pbot_pa(c), co2_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracer_gwdif_concflux_top_col(c,1:2,id_trc_ch4)   =ppm2molv(pbot_pa(c), ch4_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracer_gwdif_concflux_top_col(c,1:2,id_trc_n2o)   =ppm2molv(pbot_pa(c), n2o_ppmv(c), tair(c))!mol m-3, contant boundary condition
!         tracer_gwdif_concflux_top_col(c,1:2,id_trc_nh3x)  =ppm2molv(pbot_pa(c), nh3_ppmv(c), tair(c))!mol m-3, contant boundary condition
         tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_no3x)  = 0._r8
         tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_p_sol) = 0._r8

         if(this%use_c13)then
           tracer_gwdif_concflux_top_col(c,1:2,id_trc_c13_co2x)  =0.0168_r8
         endif
         if(this%use_c14)then
           tracer_gwdif_concflux_top_col(c,1:2,id_trc_c14_co2x)  =0.0168_r8
         endif

         bot_concflux_col(c,1,:)       = 0._r8      !zero flux boundary condition
         condc_toplay_col(c,:) = 0._r8              !those will be updated with snow resistance and hydraulic wicking resistance
         do kk = 1, ngwmobile_tracers
           if(.not. is_volatile(kk))cycle
           condc_toplay_col(c,groupid(kk))    = 1._r8/(snowres_col(c,volatilegroupid(kk))+&
             0.5_r8*dz_top(c)/diffblkm_topsoi_col(c,volatilegroupid(kk))) !m/s surface conductance
         enddo
      enddo
    end associate
  end subroutine set_boundary_conditions
  !-------------------------------------------------------------------------------

  subroutine calc_bgc_reaction(this, bounds, col, lbj, ubj, num_soilc, filter_soilc, &
       num_soilp,filter_soilp, jtops, betr_time, betrtracer_vars, tracercoeff_vars, biophysforc,    &
       tracerstate_vars, tracerflux_vars, tracerboundarycond_vars, plant_soilbgc, &
       biogeo_flux, biogeo_state, betr_status)

    !
    ! !DESCRIPTION:
    ! do bgc reaction
    ! this returns net carbon fluxes from decay and translocation
    ! and also update the related carbon/nitrogen/phosphorus(potentially) pools of OM
    ! note it is assumed the stoichiometry of the om pools are not changed during decomposition
    !
    ! !USES:
    !
    use tracerfluxType           , only : tracerflux_type
    use tracerstatetype          , only : tracerstate_type
    use tracercoeffType          , only : tracercoeff_type
    use BetrTracerType           , only : betrtracer_type
    use TracerBoundaryCondType   , only : tracerboundarycond_type
    use PlantSoilBGCMod          , only : plant_soilbgc_type
    use BetrStatusType           , only : betr_status_type
    use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
    use betr_columnType          , only : betr_column_type
    use BeTR_biogeoStateType     , only : betr_biogeo_state_type
    use v1ecaPlantSoilBGCType   , only : v1eca_plant_soilbgc_type
    use betr_ctrl                , only : betr_spinup_state
    use BeTR_TimeMod             , only : betr_time_type
    implicit none
    ! !ARGUMENTS
    class(v1eca_bgc_reaction_type) , intent(inout) :: this
    type(bounds_type)                    , intent(in) :: bounds                        ! bounds
    type(betr_column_type)               , intent(in) :: col
    integer                              , intent(in) :: num_soilc                     ! number of columns in column filter
    integer                              , intent(in) :: filter_soilc(:)               ! column filter
    integer                              , intent(in) :: num_soilp
    integer                              , intent(in) :: filter_soilp(:)               ! pft filter
    integer                              , intent(in) :: jtops(bounds%begc: )          ! top index of each column
    integer                              , intent(in) :: lbj, ubj                      ! lower and upper bounds, make sure they are > 0
    type(betr_time_type)                 , intent(in) :: betr_time                         ! model time step
    type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
    type(tracercoeff_type)               , intent(inout) :: tracercoeff_vars
    type(betr_biogeophys_input_type)     , intent(inout) :: biophysforc
    type(tracerboundarycond_type)        , intent(inout) :: tracerboundarycond_vars !
    type(tracerstate_type)               , intent(inout) :: tracerstate_vars
    type(tracerflux_type)                , intent(inout) :: tracerflux_vars
    class(plant_soilbgc_type)            , intent(inout) :: plant_soilbgc
    type(betr_biogeo_flux_type)          , intent(inout) :: biogeo_flux
    type(betr_biogeo_state_type)         , intent(inout) :: biogeo_state
    type(betr_status_type)               , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    character(len=32), parameter :: subname ='calc_bgc_reaction'
    integer                      :: fc, c, j, k
    logical :: is_surflit  !surface litter layer?
    integer :: nstates
    real(r8), allocatable :: ystates0(:)
    real(r8), allocatable :: ystatesf(:)
    real(r8) :: tnmass(num_soilc),n_mass, c_mass1, n_mass1, p_mass1
    character(len=5) :: laystr
    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__),betr_status)


    nstates = this%v1eca_bgc_index%nstvars
    allocate(ystates0(nstates))
    allocate(ystatesf(nstates))

    !pass in fluxes and state varaibles into the 1D soil bgc model
    call this%set_bgc_forc(bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
        biophysforc, biogeo_flux, plant_soilbgc, betrtracer_vars, tracercoeff_vars, &
        tracerstate_vars,betr_status)

    select type(plant_soilbgc)
    type is(v1eca_plant_soilbgc_type)
      plant_soilbgc%plant_minn_active_yield_flx_col(:) = 0._r8
      plant_soilbgc%plant_minp_active_yield_flx_col(:) = 0._r8
    end select
    !run simulation layer by layer
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      do j = lbj, ubj
        if(j<jtops(c))cycle
        is_surflit=(j<=0)
        !temperature adaptation
        call this%v1eca(c,j)%runbgc(is_surflit, betr_time%delta_time, this%v1eca_forc(c,j),nstates, &
            ystates0, ystatesf, betr_status)

	      if(betr_status%check_status())then
	        write(laystr,'(I2.2)')j
          betr_status%msg=trim(betr_status%msg)//' lay '//trim(laystr)
          return
        endif

        call this%precision_filter(nstates, ystatesf)

        call this%retrieve_output(c, j, nstates, ystates0, ystatesf, betr_time%delta_time, &
           betrtracer_vars, tracerflux_vars, tracerstate_vars, plant_soilbgc, biogeo_flux, biogeo_state)

        select type(plant_soilbgc)
        type is(v1eca_plant_soilbgc_type)
          plant_soilbgc%plant_minn_active_yield_flx_col(c)=plant_soilbgc%plant_minn_active_yield_flx_col(c) + &
             (plant_soilbgc%plant_minn_no3_active_yield_flx_vr_col(c,j) + &
              plant_soilbgc%plant_minn_nh4_active_yield_flx_vr_col(c,j))*biophysforc%dz(c,j)

          plant_soilbgc%plant_minp_active_yield_flx_col(c)=  plant_soilbgc%plant_minp_active_yield_flx_col(c) + &
            plant_soilbgc%plant_minp_active_yield_flx_vr_col(c,j) * biophysforc%dz(c,j)

        end select
      enddo
    enddo
    deallocate(ystates0)
    deallocate(ystatesf)

  end subroutine calc_bgc_reaction

  !-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! equilibrate tracers that has solid and mobile phases
    ! using the theory of mass action.
    !
    ! !USES:
    !
    use tracercoeffType       , only : tracercoeff_type
    use BeTRTracerType        , only : betrtracer_type
    use BetrStatusType        , only : betr_status_type
    ! !ARGUMENTS:
    class(v1eca_bgc_reaction_type), intent(inout)    :: this
    type(bounds_type),                    intent(in)    :: bounds
    integer,                              intent(in)    :: lbj, ubj
    integer,                              intent(in)    :: jtops(bounds%begc: )        ! top label of each column
    integer,                              intent(in)    :: num_soilc
    integer,                              intent(in)    :: filter_soilc(:)
    type(betrtracer_type),                intent(in)    :: betrtracer_vars
    type(tracercoeff_type),               intent(in)    :: tracercoeff_vars
    type(tracerstate_type),               intent(inout) :: tracerstate_vars
    type(betr_status_type)              , intent(out)   :: betr_status
    !
    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'do_tracer_equilibration'

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__),betr_status)

    if (this%dummy_compiler_warning) continue
  end subroutine do_tracer_equilibration

  !-----------------------------------------------------------------------
  subroutine readParams(this, namelist_buffer, betrtracer_vars )
    !
    ! !DESCRIPTION:
    ! read model parameters
    ! !USES:
    use BeTRTracerType   , only : BeTRTracer_Type

    ! !ARGUMENTS:
    class(v1eca_bgc_reaction_type) , intent(inout)    :: this
    character(len=*)                  , intent(in)  :: namelist_buffer
    type(BeTRTracer_Type)                , intent(inout) :: betrtracer_vars

   !x
    if (this%dummy_compiler_warning) continue
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, col, betrtracer_vars, biophysforc, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! cold initialization
    ! !USES:
    !

    use BeTRTracerType    , only : BeTRTracer_Type
    use betr_varcon       , only : spval => bspval, ispval => bispval
    use betr_varcon       , only : denh2o => bdenh2o
    use betr_columnType   , only : betr_column_type
    use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
    use UnitConvertMod, only : ppm2molv
    implicit none
    ! !ARGUMENTS:
    class(v1eca_bgc_reaction_type) , intent(inout)    :: this
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
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !-----------------------------------------------------------------------

    associate(                                         &
         volatileid => betrtracer_vars%volatileid    , &
         nfrozen_tracers => betrtracer_vars%nfrozen_tracers, &
         n2_ppmv    => biophysforc%n2_ppmv_col       , &
         o2_ppmv    => biophysforc%o2_ppmv_col       , &
         ar_ppmv    => biophysforc%ar_ppmv_col       , &
         co2_ppmv   => biophysforc%co2_ppmv_col      , &
         ch4_ppmv   => biophysforc%ch4_ppmv_col      , &
         n2o_ppmv   => biophysforc%n2o_ppmv_col      , &
         nh3_ppmv   => biophysforc%nh3_ppmv_col      , &
         pbot_pa    => biophysforc%forc_pbot_downscaled_col, &
         tair       => biophysforc%forc_t_downscaled_col &
         )
    do c = bounds%begc, bounds%endc

      !dual phase tracers
      tracerstate_vars%tracer_conc_mobile_col    (c,:, :                                   )  = 0._r8
      if(nfrozen_tracers>0)tracerstate_vars%tracer_conc_frozen_col    (c,:, :              )  = 0._r8
      tracerstate_vars%tracer_conc_surfwater_col (c,:                                      )  = 0._r8
      tracerstate_vars%tracer_conc_aquifer_col   (c,:                                      )  = 0._r8
      tracerstate_vars%tracer_conc_grndwater_col (c,:                                      )  = 0._r8
      tracerstate_vars%tracer_conc_mobile_col    (c,:, betrtracer_vars%id_trc_o2) = ppm2molv(pbot_pa(c), o2_ppmv(c), tair(c))
      if(betrtracer_vars%nsolid_equil_tracers>0)then
        tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = 0._r8
      endif
      tracerstate_vars%tracer_soi_molarmass_col(c,:)          = 0._r8
    enddo
    call this%update_atm_conc(bounds, betrtracer_vars, biophysforc, tracerstate_vars)
    end associate
  end subroutine InitCold

  !------------------------------------------------------------------------------
  subroutine update_atm_conc(this, bounds, betrtracer_vars, biophysforc, tracerstate_vars)
    use UnitConvertMod, only : ppm2molv
  implicit none
    ! !ARGUMENTS:
    class(v1eca_bgc_reaction_type) , intent(inout)    :: this
    type(betr_bounds_type)           , intent(in)    :: bounds
    type(BeTRTracer_Type)            , intent(in)    :: betrtracer_vars
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars

  integer :: c
    associate(                                         &
         volatileid => betrtracer_vars%volatileid    , &
         n2_ppmv    => biophysforc%n2_ppmv_col       , &
         o2_ppmv    => biophysforc%o2_ppmv_col       , &
         ar_ppmv    => biophysforc%ar_ppmv_col       , &
         co2_ppmv   => biophysforc%co2_ppmv_col      , &
         ch4_ppmv   => biophysforc%ch4_ppmv_col      , &
         n2o_ppmv   => biophysforc%n2o_ppmv_col      , &
         nh3_ppmv   => biophysforc%nh3_ppmv_col      , &
         pbot_pa    => biophysforc%forc_pbot_downscaled_col, &
         tair       => biophysforc%forc_t_downscaled_col &
         )

  do c = bounds%begc, bounds%endc
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_n2   )) = ppm2molv(pbot_pa(c), n2_ppmv(c), tair(c))
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_o2   )) = ppm2molv(pbot_pa(c), o2_ppmv(c), tair(c))
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_ar   )) = ppm2molv(pbot_pa(c), ar_ppmv(c), tair(c))
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_co2x )) = ppm2molv(pbot_pa(c), co2_ppmv(c), tair(c))
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_ch4  )) = ppm2molv(pbot_pa(c), ch4_ppmv(c), tair(c))
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_n2o  )) = ppm2molv(pbot_pa(c), n2o_ppmv(c), tair(c))
  enddo

  end associate
  end subroutine update_atm_conc
  !------------------------------------------------------------------------------
  subroutine retrieve_biogeoflux(this, num_soilc, filter_soilc, tracerflux_vars, &
  betrtracer_vars, biogeo_flux)

  use tracerfluxType           , only : tracerflux_type
  use BeTRTracerType           , only : BeTRTracer_Type
  implicit none
  class(v1eca_bgc_reaction_type) , intent(inout)    :: this !!
  integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
  integer                          , intent(in)    :: filter_soilc(:)             ! column filter
  type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
  type(tracerflux_type)            , intent(in)    :: tracerflux_vars
  type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux

   integer :: fc, c, trcid, c_loc, n_loc, p_loc

   associate(   &
      tracer_flx_leaching_col => tracerflux_vars%tracer_flx_leaching_col, &
      tracer_flx_surfrun_col  => tracerflux_vars%tracer_flx_surfrun_col, &
      tracer_flx_drain_col    => tracerflux_vars%tracer_flx_drain_col, &
      tracer_flx_surfemi_col  => tracerflux_vars%tracer_flx_surfemi_col, &
      volatileid              => betrtracer_vars%volatileid   , &
      id_trc_nh3x             => betrtracer_vars%id_trc_nh3x  , &
      id_trc_no3x             => betrtracer_vars%id_trc_no3x,  &
      id_trc_co2x             => betrtracer_vars%id_trc_co2x,  &
      id_trc_p_sol            => betrtracer_vars%id_trc_p_sol  &
   )

    c_loc=this%v1eca_bgc_index%c_loc
    n_loc=this%v1eca_bgc_index%n_loc
    p_loc=this%v1eca_bgc_index%p_loc

   !retrieve tracer losses through surface and subsurface runoffs
   !no3 leach, no3 runoff
   do fc = 1, num_soilc
     c = filter_soilc(fc)
     biogeo_flux%n14flux_vars%smin_no3_leached_col(c) = tracer_flx_leaching_col(c,id_trc_no3x) * natomw  ![gN/m2/s]
     biogeo_flux%n14flux_vars%smin_no3_runoff_col(c) = tracer_flx_surfrun_col(c,id_trc_no3x) * natomw
     biogeo_flux%n14flux_vars%smin_no3_qdrain_col(c) = tracer_flx_drain_col(c,id_trc_no3x) * natomw

     biogeo_flux%c12flux_vars%co2_soi_flx_col(c) = tracer_flx_surfemi_col(c,volatileid(id_trc_co2x))*catomw

     !return mineral p
     biogeo_flux%p31flux_vars%sminp_leached_col(c) = tracer_flx_leaching_col(c,id_trc_p_sol) * patomw
     biogeo_flux%p31flux_vars%sminp_runoff_col(c) = tracer_flx_surfrun_col(c,id_trc_p_sol) * patomw
     biogeo_flux%p31flux_vars%sminp_qdrain_col(c) = tracer_flx_drain_col(c,id_trc_p_sol) * patomw

     biogeo_flux%qflx_rofliq_qsub_dic_col(c) = tracer_flx_surfrun_col(c,id_trc_co2x) * catomw
     biogeo_flux%qflx_rofliq_qsur_dic_col(c) = tracer_flx_drain_col(c,id_trc_co2x) * catomw

   enddo

   end associate

  end subroutine retrieve_biogeoflux

  !------------------------------------------------------------------------------
  subroutine set_bgc_forc(this, bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
      biophysforc, biogeo_flux, plant_soilbgc, betrtracer_vars, tracercoeff_vars,  &
      tracerstate_vars, betr_status)
  !DESCRIPTION
  !set up forcing for running bgc
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use tracerstatetype          , only : tracerstate_type
  use betr_decompMod           , only : betr_bounds_type
  use tracercoeffType          , only : tracercoeff_type
  use betr_columnType          , only : betr_column_type
  use BetrTracerType           , only : betrtracer_type
  use v1ecaPlantSoilBGCType   , only : v1eca_plant_soilbgc_type
  use betr_ctrl                , only : betr_spinup_state
  use betr_varcon              , only : grav => bgrav
  implicit none
  class(v1eca_bgc_reaction_type) , intent(inout)    :: this
  type(bounds_type)                    , intent(in) :: bounds                         ! bounds
  type(betr_column_type)               , intent(in) :: col
  integer                              , intent(in) :: jtops(bounds%begc: ) ! top index of each column
  integer                              , intent(in) :: lbj, ubj                       ! lower and upper bounds, make sure they are > 0
  integer                              , intent(in) :: num_soilc       ! number of columns in column filter
  integer                              , intent(in) :: filter_soilc(:) ! column filter
  type(betr_biogeophys_input_type)     , intent(in) :: biophysforc
  type(betr_biogeo_flux_type)          , intent(in) :: biogeo_flux
  class(plant_soilbgc_type)            , intent(in) :: plant_soilbgc
  type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
  type(tracerstate_type)               , intent(in) :: tracerstate_vars
  type(tracercoeff_type)               , intent(in) :: tracercoeff_vars
  type(betr_status_type)               , intent(out)   :: betr_status

  integer :: j, fc, c
  integer :: k1, k2
  real(r8):: latacc(bounds%begc:bounds%endc)
  real(r8), parameter :: tiny_cval =1.e-16_r8
  associate( &
     litr_beg =>  this%v1eca_bgc_index%litr_beg  , &
     litr_end =>  this%v1eca_bgc_index%litr_end  , &
     wood_beg =>  this%v1eca_bgc_index%wood_beg  , &
     wood_end =>  this%v1eca_bgc_index%wood_end  , &
     som_beg =>  this%v1eca_bgc_index%som_beg    , &
     som_end =>  this%v1eca_bgc_index%som_end    , &
     pom_beg =>  this%v1eca_bgc_index%pom_beg    , &
     pom_end =>  this%v1eca_bgc_index%pom_end    , &
     Bm_beg  =>  this%v1eca_bgc_index%Bm_beg     , &
     Bm_end  =>  this%v1eca_bgc_index%Bm_end     , &
     ncentpools => this%v1eca_bgc_index%nom_pools, &
     groupid =>  betrtracer_vars%groupid         ,  &
     nelms   => this%v1eca_bgc_index%nelms, &
     c_loc   => this%v1eca_bgc_index%c_loc, &
     n_loc   => this%v1eca_bgc_index%n_loc, &
     p_loc   => this%v1eca_bgc_index%p_loc, &
     c13_loc => this%v1eca_bgc_index%c13_loc, &
     c14_loc => this%v1eca_bgc_index%c14_loc, &
     volatilegroupid => betrtracer_vars%volatilegroupid, &
     volatileid =>  betrtracer_vars%volatileid &
  )
  call betr_status%reset()
  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__),betr_status)

  do j = lbj, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      if(j<jtops(c))cycle
      this%v1eca_forc(c,j)%latacc = latacc(c)
      this%v1eca_forc(c,j)%plant_ntypes = this%nactpft
      this%v1eca_forc(c,j)%ystates(:) = 0._r8
      this%v1eca_forc(c,j)%decomp_k(1:ncentpools)=biogeo_flux%c12flux_vars%decomp_k(c,j,1:ncentpools)
!      print*,'sefok,j',j,this%v1eca_forc(c,j)%decomp_k(1:ncentpools)
      this%v1eca_forc(c,j)%t_scalar = biophysforc%c12flx%in_t_scalar(c,j)
      this%v1eca_forc(c,j)%w_scalar = biophysforc%c12flx%in_w_scalar(c,j)
      this%v1eca_forc(c,j)%o_scalar = biophysforc%c12flx%in_o_scalar(c,j)
      !litter C
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+c_loc-1),biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,1)/catomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+nelms+c_loc-1),biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,2)/catomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+2*nelms+c_loc-1),biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,3)/catomw)

      !wood C
      FPMAX(this%v1eca_forc(c,j)%ystates(wood_beg+c_loc-1),biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,4)/catomw)

      !som C
      FPMAX(this%v1eca_forc(c,j)%ystates(Bm_beg+c_loc-1),biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,5)/catomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(pom_beg+c_loc-1),biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,6)/catomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(som_beg+c_loc-1),biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,7)/catomw)

      if(this%use_c13)then
        !litter C
        FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+c13_loc-1),biophysforc%c13flx%in_decomp_cpools_vr_col(c,j,1))/c13atomw
        FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+nelms+c13_loc-1),biophysforc%c13flx%in_decomp_cpools_vr_col(c,j,2)/c13atomw)
        FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+2*nelms+c13_loc-1),biophysforc%c13flx%in_decomp_cpools_vr_col(c,j,3)/c13atomw)

        !wood C
        FPMAX(this%v1eca_forc(c,j)%ystates(wood_beg+c13_loc-1),biophysforc%c13flx%in_decomp_cpools_vr_col(c,j,4)/c13atomw)

        !som C
        FPMAX(this%v1eca_forc(c,j)%ystates(Bm_beg+c13_loc-1),biophysforc%c13flx%in_decomp_cpools_vr_col(c,j,5)/c13atomw)
        FPMAX(this%v1eca_forc(c,j)%ystates(pom_beg+c13_loc-1),biophysforc%c13flx%in_decomp_cpools_vr_col(c,j,6)/c13atomw)
        FPMAX(this%v1eca_forc(c,j)%ystates(som_beg+c13_loc-1),biophysforc%c13flx%in_decomp_cpools_vr_col(c,j,7)/c13atomw)
      endif


      if(this%use_c14)then
        !litter C
        FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+c14_loc-1), biophysforc%c14flx%in_decomp_cpools_vr_col(c,j,1)/c14atomw)
        FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+nelms+c14_loc-1), biophysforc%c14flx%in_decomp_cpools_vr_col(c,j,2)/c14atomw)
        FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+2*nelms+c14_loc-1), biophysforc%c14flx%in_decomp_cpools_vr_col(c,j,3)/c14atomw)

        !wood C
        FPMAX(this%v1eca_forc(c,j)%ystates(wood_beg+c_loc-1), biophysforc%c14flx%in_decomp_cpools_vr_col(c,j,4)/c14atomw)

        !som C
        FPMAX(this%v1eca_forc(c,j)%ystates(Bm_beg+c_loc-1), biophysforc%c14flx%in_decomp_cpools_vr_col(c,j,5)/c14atomw)
        FPMAX(this%v1eca_forc(c,j)%ystates(pom_beg+c_loc-1), biophysforc%c14flx%in_decomp_cpools_vr_col(c,j,6)/c14atomw)
        FPMAX(this%v1eca_forc(c,j)%ystates(som_beg+c_loc-1), biophysforc%c14flx%in_decomp_cpools_vr_col(c,j,7)/c14atomw)
      endif
      !litter N
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+n_loc-1), biophysforc%n14flx%in_decomp_npools_vr_col(c,j,1)/natomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+nelms+n_loc-1), biophysforc%n14flx%in_decomp_npools_vr_col(c,j,2)/natomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+2*nelms+n_loc-1), biophysforc%n14flx%in_decomp_npools_vr_col(c,j,3)/natomw)

      !wood N
      FPMAX(this%v1eca_forc(c,j)%ystates(wood_beg+n_loc-1), biophysforc%n14flx%in_decomp_npools_vr_col(c,j,4)/natomw)

      !som N
      FPMAX(this%v1eca_forc(c,j)%ystates(Bm_beg+n_loc-1), biophysforc%n14flx%in_decomp_npools_vr_col(c,j,5)/natomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(pom_beg+n_loc-1), biophysforc%n14flx%in_decomp_npools_vr_col(c,j,6)/natomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(som_beg+n_loc-1), biophysforc%n14flx%in_decomp_npools_vr_col(c,j,7)/natomw)

      !litter P
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+p_loc-1), biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,1)/patomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+nelms+p_loc-1), biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,2)/patomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(litr_beg+2*nelms+p_loc-1), biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,3)/patomw)

      !wood P
      FPMAX(this%v1eca_forc(c,j)%ystates(wood_beg+p_loc-1), biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,4)/patomw)

      !som P
      FPMAX(this%v1eca_forc(c,j)%ystates(Bm_beg+p_loc-1), biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,5)/patomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(pom_beg+p_loc-1), biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,6)/patomw)
      FPMAX(this%v1eca_forc(c,j)%ystates(som_beg+p_loc-1), biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,7)/patomw)

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_n2), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2))

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_o2), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2))

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_ar), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar))

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_co2), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x))

      if(this%use_c13)then
        FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_c13_co2), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_c14_co2), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x))
      endif

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_ch4), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4))

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_nh4), biophysforc%n14flx%in_sminn_nh4_vr_col(c,j)/natomw)

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_no3), biophysforc%n14flx%in_sminn_no3_vr_col(c,j)/natomw)

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_n2o), tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2o))

      FPMAX(this%v1eca_forc(c,j)%ystates(this%v1eca_bgc_index%lid_minp_soluble), biophysforc%p31flx%in_sminp_vr_col(c,j)/patomw)

      !burning fraction
      !environmental variables
      this%v1eca_forc(c,j)%temp   = biophysforc%t_soisno_col(c,j)            !temperature
      this%v1eca_forc(c,j)%depz   = biophysforc%z(c,j)            !depth of the soil
      this%v1eca_forc(c,j)%dzsoi  = biophysforc%dz(c,j)            !soil thickness
      this%v1eca_forc(c,j)%sucsat  = biophysforc%sucsat_col(c,j)            ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
      this%v1eca_forc(c,j)%soilpsi = max(biophysforc%smp_l_col(c,j)*grav*1.e-6_r8,-15._r8)    ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
      this%v1eca_forc(c,j)%bsw = biophysforc%bsw_col(c,j)
      this%v1eca_forc(c,j)%bd   = biophysforc%bd_col(c,j)              !bulk density
      this%v1eca_forc(c,j)%pct_sand = biophysforc%cellsand_col(c,j)
      this%v1eca_forc(c,j)%pct_clay = biophysforc%cellclay_col(c,j)
      this%v1eca_forc(c,j)%h2osoi_vol = biophysforc%h2osoi_vol_col(c,j)
      this%v1eca_forc(c,j)%h2osoi_liq = biophysforc%h2osoi_liq_col(c,j)
      this%v1eca_forc(c,j)%air_vol = biophysforc%air_vol_col(c,j)
      this%v1eca_forc(c,j)%finundated = biophysforc%finundated_col(c)
      this%v1eca_forc(c,j)%watsat = biophysforc%watsat_col(c,j)
      this%v1eca_forc(c,j)%watfc = biophysforc%watfc_col(c,j)
      this%v1eca_forc(c,j)%cellorg = biophysforc%cellorg_col(c,j)
      this%v1eca_forc(c,j)%pH = 6.5_r8

      !conductivity for plant-aided gas transport
      this%v1eca_forc(c,j)%aren_cond_n2 = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_n2))

      this%v1eca_forc(c,j)%aren_cond_o2 = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_o2))

      this%v1eca_forc(c,j)%aren_cond_n2o = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_n2o))

      this%v1eca_forc(c,j)%aren_cond_co2 = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_co2x))

      if(this%use_c13)then
        this%v1eca_forc(c,j)%aren_cond_co2_c13 = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%v1eca_forc(c,j)%aren_cond_co2_c14 = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%v1eca_forc(c,j)%aren_cond_ar = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_ar))
      this%v1eca_forc(c,j)%aren_cond_ch4 = &
          tracercoeff_vars%aere_cond_col(c,volatilegroupid(betrtracer_vars%id_trc_ch4))

      !phase conversion parameter
      this%v1eca_forc(c,j)%ch4_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,volatilegroupid(betrtracer_vars%id_trc_ch4))
      this%v1eca_forc(c,j)%co2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,volatilegroupid(betrtracer_vars%id_trc_co2x))
      this%v1eca_forc(c,j)%o2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,volatilegroupid(betrtracer_vars%id_trc_o2))
      this%v1eca_forc(c,j)%n2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,volatilegroupid(betrtracer_vars%id_trc_n2))
      this%v1eca_forc(c,j)%ar_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,volatilegroupid(betrtracer_vars%id_trc_ar))
      this%v1eca_forc(c,j)%n2o_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,volatilegroupid(betrtracer_vars%id_trc_n2o))
      this%v1eca_forc(c,j)%o2_w2b = &
          tracercoeff_vars%aqu2bulkcef_mobile_col(c,j,groupid(betrtracer_vars%id_trc_o2))

      !atmospheric pressure (mol/m3) for gas ventilation.
      this%v1eca_forc(c,j)%conc_atm_n2 = &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_n2))
      this%v1eca_forc(c,j)%conc_atm_n2o= &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_n2o))
      this%v1eca_forc(c,j)%conc_atm_o2 = &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_o2))
      this%v1eca_forc(c,j)%conc_atm_ar = &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_ar))
      this%v1eca_forc(c,j)%conc_atm_co2 = &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%v1eca_forc(c,j)%conc_atm_co2_c13 = &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%v1eca_forc(c,j)%conc_atm_co2_c14 = &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%v1eca_forc(c,j)%conc_atm_ch4 = &
          tracerstate_vars%tracer_conc_atm_col(c,volatileid(betrtracer_vars%id_trc_ch4))

      this%v1eca_forc(c,j)%soilorder = biophysforc%isoilorder(c)

    enddo
  enddo

  select type(plant_soilbgc)
  type is(v1eca_plant_soilbgc_type)
  do j = lbj, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      this%v1eca_forc(c,j)%rt_ar  = biophysforc%c12flx%rt_vr_col(c,j)            !root autotrophic respiration, mol CO2/m3/s
    enddo
  enddo
  end select
  end associate
  end subroutine set_bgc_forc

  !------------------------------------------------------------------------------
  subroutine precision_filter(this, nstates, ystatesf)

  !
  ! DESCRIPTION
  ! reset tiny som values to zero.
  ! when the carbon amount is below some minimum value, the relative magnitudes
  ! of N and P for SOM pool may be close to random error. Occaisonly,
  ! the P amount may be larger than N amount, causing the code to crash.
  ! This fix set C, N and P to zero when C is below a threshold.
  implicit none
  class(v1eca_bgc_reaction_type) , intent(inout)    :: this
  integer                              , intent(in) :: nstates
  real(r8)                             , intent(inout) :: ystatesf(nstates)

  real(r8), parameter :: tiny_val=1.e-13_r8
  integer :: jj
  integer :: kc, kn, kp, kc13, kc14
  associate(                              &
    nelms   => this%v1eca_bgc_index%nelms, &
    c_loc   => this%v1eca_bgc_index%c_loc, &
    n_loc   => this%v1eca_bgc_index%n_loc, &
    p_loc   => this%v1eca_bgc_index%p_loc, &
    c13_loc => this%v1eca_bgc_index%c13_loc, &
    c14_loc => this%v1eca_bgc_index%c14_loc, &
    lit2    => this%v1eca_bgc_index%lit2 , &
    lit3    => this%v1eca_bgc_index%lit3 , &
    ncentpools => this%v1eca_bgc_index%nom_pools, &
    is_cenpool_som => this%v1eca_bgc_index%is_cenpool_som &
  )
  do jj = 1, ncentpools
    kc = (jj-1) * nelms + c_loc
    kn = (jj-1) * nelms + n_loc
    kp = (jj-1) * nelms + p_loc
    if( ystatesf(kc) <= tiny_val  .and. is_cenpool_som(jj))then
      ystatesf(kc)=0._r8
      ystatesf(kn)=0._r8
      ystatesf(kp)=0._r8
      if(this%use_c14)then
        kc14 = (jj-1) * nelms + c14_loc
        ystatesf(kc14)=0._r8
      endif
      if(this%use_c13)then
        kc13 = (jj-1) * nelms + c13_loc
        ystatesf(kc13)=0._r8
      endif
    endif
  enddo


  end associate

  end subroutine precision_filter
  !------------------------------------------------------------------------------
  subroutine retrieve_output(this, c, j, nstates, ystates0, ystatesf, dtime, betrtracer_vars, tracerflux_vars,&
     tracerstate_vars, plant_soilbgc, biogeo_flux, biogeo_state)
  !DESCRIPTION
  !retrieve flux and state variables after evolving the bgc calculation
  !
  !USES
  use BetrTracerType           , only : betrtracer_type
  use tracerfluxType           , only : tracerflux_type
  use tracerstatetype          , only : tracerstate_type
  use betr_ctrl                , only : betr_spinup_state
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use v1ecaPlantSoilBGCType      , only : v1eca_plant_soilbgc_type
  use BeTR_biogeoStateType     , only : betr_biogeo_state_type
  use tracer_varcon            , only : catomw, natomw, patomw, fix_ip
  implicit none
  class(v1eca_bgc_reaction_type) , intent(inout)    :: this
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
  type(betr_biogeo_state_type)         , intent(inout) :: biogeo_state

  integer :: k, k1, k2, jj, p
  integer :: trcid

  associate(                                                                &
    nom_pools             => this%v1eca_bgc_index%nom_pools               , & !
    nelms                 => this%v1eca_bgc_index%nelms                   , & !
    litr_beg              => this%v1eca_bgc_index%litr_beg                , & !
    litr_end              => this%v1eca_bgc_index%litr_end                , & !
    wood_beg              => this%v1eca_bgc_index%wood_beg                , & !
    wood_end              => this%v1eca_bgc_index%wood_end                , & !
    som_beg               => this%v1eca_bgc_index%som_beg                 , & !
    som_end               => this%v1eca_bgc_index%som_end                 , & !
    pom_beg               => this%v1eca_bgc_index%pom_beg                 , & !
    pom_end               => this%v1eca_bgc_index%pom_end                 , & !
    Bm_beg                => this%v1eca_bgc_index%Bm_beg                  , & !
    Bm_end                => this%v1eca_bgc_index%Bm_end                  , & !
    c_loc                 => this%v1eca_bgc_index%c_loc                   , &
    n_loc                 => this%v1eca_bgc_index%n_loc                   , &
    p_loc                 => this%v1eca_bgc_index%p_loc                   , &
    c13_loc               => this%v1eca_bgc_index%c13_loc                 , &
    c14_loc               => this%v1eca_bgc_index%c14_loc                 , &
    volatileid            => betrtracer_vars%volatileid                   , &
    tracer_flx_netpro_vr  => tracerflux_vars%tracer_flx_netpro_vr_col     , & !
    tracer_flx_parchm_vr  => tracerflux_vars%tracer_flx_parchm_vr_col     , & !
    ngwmobile_tracers     => betrtracer_vars%ngwmobile_tracers              & !
  )

      !tracer states
      biogeo_state%c12state_vars%decomp_cpools_vr(c,j,1) = ystatesf(litr_beg+c_loc-1)*catomw
      biogeo_state%c12state_vars%decomp_cpools_vr(c,j,2) = ystatesf(litr_beg+nelms+c_loc-1)*catomw
      biogeo_state%c12state_vars%decomp_cpools_vr(c,j,3) = ystatesf(litr_beg+2*nelms+c_loc-1)*catomw
      biogeo_state%c12state_vars%decomp_cpools_vr(c,j,4) = ystatesf(wood_beg+c_loc-1)*catomw
      biogeo_state%c12state_vars%decomp_cpools_vr(c,j,5) = ystatesf(Bm_beg+c_loc-1)*catomw
      biogeo_state%c12state_vars%decomp_cpools_vr(c,j,6) = ystatesf(pom_beg+c_loc-1)*catomw
      biogeo_state%c12state_vars%decomp_cpools_vr(c,j,7) = ystatesf(som_beg+c_loc-1)*catomw

      biogeo_state%n14state_vars%decomp_npools_vr(c,j,1) = ystatesf(litr_beg+n_loc-1)*natomw
      biogeo_state%n14state_vars%decomp_npools_vr(c,j,2) = ystatesf(litr_beg+nelms+n_loc-1)*natomw
      biogeo_state%n14state_vars%decomp_npools_vr(c,j,3) = ystatesf(litr_beg+2*nelms+n_loc-1)*natomw
      biogeo_state%n14state_vars%decomp_npools_vr(c,j,4) = ystatesf(wood_beg+n_loc-1)*natomw
      biogeo_state%n14state_vars%decomp_npools_vr(c,j,5) = ystatesf(Bm_beg+n_loc-1)*natomw
      biogeo_state%n14state_vars%decomp_npools_vr(c,j,6) = ystatesf(pom_beg+n_loc-1)*natomw
      biogeo_state%n14state_vars%decomp_npools_vr(c,j,7) = ystatesf(som_beg+n_loc-1)*natomw

      biogeo_state%p31state_vars%decomp_ppools_vr(c,j,1) = ystatesf(litr_beg+p_loc-1)*patomw
      biogeo_state%p31state_vars%decomp_ppools_vr(c,j,2) = ystatesf(litr_beg+nelms+p_loc-1)*patomw
      biogeo_state%p31state_vars%decomp_ppools_vr(c,j,3) = ystatesf(litr_beg+2*nelms+p_loc-1)*patomw
      biogeo_state%p31state_vars%decomp_ppools_vr(c,j,4) = ystatesf(wood_beg+p_loc-1)*patomw
      biogeo_state%p31state_vars%decomp_ppools_vr(c,j,5) = ystatesf(Bm_beg  + p_loc-1)*patomw
      biogeo_state%p31state_vars%decomp_ppools_vr(c,j,6) = ystatesf(pom_beg + p_loc-1)*patomw
      biogeo_state%p31state_vars%decomp_ppools_vr(c,j,7) = ystatesf(som_beg + p_loc-1)*patomw

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2) = &
        ystatesf(this%v1eca_bgc_index%lid_n2)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2) = &
        ystatesf(this%v1eca_bgc_index%lid_o2)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar) = &
        ystatesf(this%v1eca_bgc_index%lid_ar)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x) = &
        ystatesf(this%v1eca_bgc_index%lid_co2)

      if(this%use_c13)then
        tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x) = &
          ystatesf(this%v1eca_bgc_index%lid_c13_co2)
      endif

      if(this%use_c14)then
        tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x) = &
          ystatesf(this%v1eca_bgc_index%lid_c14_co2)
      endif

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4) = &
        ystatesf(this%v1eca_bgc_index%lid_ch4)

      if(this%v1eca_bgc_index%lid_supp_minn>0)then
        biogeo_flux%n14flux_vars%supplement_to_sminn_vr_col(c,j) = &
          ystatesf(this%v1eca_bgc_index%lid_supp_minn)*natomw/dtime
      else
        biogeo_flux%n14flux_vars%supplement_to_sminn_vr_col(c,j) = 0._r8
      endif
      biogeo_state%n14state_vars%sminn_nh4_vr_col(c, j) = ystatesf(this%v1eca_bgc_index%lid_nh4)*natomw
      biogeo_state%n14state_vars%sminn_no3_vr_col(c, j)= ystatesf(this%v1eca_bgc_index%lid_no3)*natomw

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2o) = &
        ystatesf(this%v1eca_bgc_index%lid_n2o)

      !not fixing inorganic Phosphorus
      if(.not. fix_ip)then
        k2 = this%v1eca_bgc_index%lid_minp_sorb
         biogeo_flux%p31flux_vars%adsorb_to_labilep_vr_col(c,j)= ystatesf(k2)*patomw/dtime

        if(this%v1eca_bgc_index%lid_supp_minp>0)then
            !no P-limitation in this time step
            biogeo_flux%p31flux_vars%supplement_to_sminp_vr_col(c,j) = &
              ystatesf(this%v1eca_bgc_index%lid_supp_minp)*patomw/dtime
        else
            biogeo_flux%p31flux_vars%supplement_to_sminp_vr_col(c,j) = 0._r8
        endif

        biogeo_state%p31state_vars%solutionp_vr_col(c,j) = ystatesf(this%v1eca_bgc_index%lid_minp_soluble)*patomw

        !fluxes
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_p_sol) =      &
          tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_p_sol) + &
          ystatesf(this%v1eca_bgc_index%lid_minp_soluble) - &
          ystates0(this%v1eca_bgc_index%lid_minp_soluble)
        biogeo_flux%p31flux_vars%net_mineralization_p_vr_col(c,j)  = &
           -ystatesf(this%v1eca_bgc_index%lid_minp_immob)*patomw/dtime

      endif
      !tracer fluxes
!      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_nh3x) =  &
!        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_nh3x) + &
!        ystatesf(this%v1eca_bgc_index%lid_nh4) - &
!        ystates0(this%v1eca_bgc_index%lid_nh4)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_no3x)  =  &
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_no3x) + &
        ystatesf(this%v1eca_bgc_index%lid_no3) - &
        ystates0(this%v1eca_bgc_index%lid_no3)

      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_o2) ) = &
         ystatesf(this%v1eca_bgc_index%lid_o2_paere )  - &
         ystates0(this%v1eca_bgc_index%lid_o2_paere)

      if ( betr_spinup_state == 0 ) then
        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2)  ) = &
          ystatesf(this%v1eca_bgc_index%lid_n2_paere)  - &
          ystates0(this%v1eca_bgc_index%lid_n2_paere)

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ar)  ) = &
          ystatesf(this%v1eca_bgc_index%lid_ar_paere)  - &
          ystates0(this%v1eca_bgc_index%lid_ar_paere)

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_co2x)) = &
          ystatesf(this%v1eca_bgc_index%lid_co2_paere)  - &
          ystates0(this%v1eca_bgc_index%lid_co2_paere)

        if(this%use_c13)then
          tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c13_co2x)) = &
            ystatesf(this%v1eca_bgc_index%lid_c13_co2_paere)  - &
            ystates0(this%v1eca_bgc_index%lid_c13_co2_paere)
        endif

        if(this%use_c14)then
          tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c14_co2x)) = &
            ystatesf(this%v1eca_bgc_index%lid_c14_co2_paere)  - &
            ystates0(this%v1eca_bgc_index%lid_c14_co2_paere)
        endif

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ch4) ) = &
          ystatesf(this%v1eca_bgc_index%lid_ch4_paere)  - &
          ystates0(this%v1eca_bgc_index%lid_ch4_paere)

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2o) ) = &
          ystatesf(this%v1eca_bgc_index%lid_n2o_paere)  - &
          ystates0(this%v1eca_bgc_index%lid_n2o_paere)
      endif

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2) = &
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2) + &
        ystatesf(this%v1eca_bgc_index%lid_n2) - &
        ystates0(this%v1eca_bgc_index%lid_n2)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_co2x ) = &
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_co2x ) + &
        ystatesf(this%v1eca_bgc_index%lid_co2) - &
        ystates0(this%v1eca_bgc_index%lid_co2)

      if(this%use_c13)then
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c13_co2x ) = &
          tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c13_co2x ) + &
          ystatesf(this%v1eca_bgc_index%lid_c13_co2) - &
          ystates0(this%v1eca_bgc_index%lid_c13_co2)
      endif

      if(this%use_c14)then
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c14_co2x ) = &
          tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c14_co2x ) + &
          ystatesf(this%v1eca_bgc_index%lid_c14_co2) - &
          ystates0(this%v1eca_bgc_index%lid_c14_co2)
      endif

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2o  ) = &
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2o  ) + &
        ystatesf(this%v1eca_bgc_index%lid_n2o) - &
        ystates0(this%v1eca_bgc_index%lid_n2o)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_o2   ) = &
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_o2   ) + &
        ystatesf(this%v1eca_bgc_index%lid_o2) - &
        ystates0(this%v1eca_bgc_index%lid_o2)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ch4  ) = &
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ch4  ) + &
        ystatesf(this%v1eca_bgc_index%lid_ch4) - &
        ystates0(this%v1eca_bgc_index%lid_ch4)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ar) = &
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ar) + &
        ystatesf(this%v1eca_bgc_index%lid_ar) - &
        ystates0(this%v1eca_bgc_index%lid_ar)


      !get net production for om pools
      do k = 1, litr_end-litr_beg + 1
        k1 = litr_beg+k-1; k2 = betrtracer_vars%id_trc_beg_litr + k-1
        tracer_flx_netpro_vr(c,j,k2) = tracer_flx_netpro_vr(c,j,k2) + &
           ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, wood_end-wood_beg + 1
        k1 = wood_beg+k-1; k2 = betrtracer_vars%id_trc_beg_wood + k-1
        tracer_flx_netpro_vr(c,j,k2) =  tracer_flx_netpro_vr(c,j,k2) + &
          ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, Bm_end-Bm_beg + 1
        k1 = Bm_beg+k-1; k2 = betrtracer_vars%id_trc_beg_Bm+ k-1
        tracer_flx_netpro_vr(c,j,k2) = tracer_flx_netpro_vr(c,j,k2) + &
          ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, som_end-som_beg + 1
        k1 = som_beg+k-1; k2 = betrtracer_vars%id_trc_beg_som+ k-1
        tracer_flx_netpro_vr(c,j,k2) = tracer_flx_netpro_vr(c,j,k2) + &
          ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, pom_end-pom_beg + 1
        k1 = pom_beg+k-1; k2 = betrtracer_vars%id_trc_beg_pom+ k-1
        tracer_flx_netpro_vr(c,j,k2) = tracer_flx_netpro_vr(c,j,k2) + &
          ystatesf(k1) - ystates0(k1)
      enddo


      !plant soil bgc

      !biogeo_flux
      biogeo_flux%c12flux_vars%hr_vr_col(c,j) = &
        (ystatesf(this%v1eca_bgc_index%lid_co2_hr) - &
        ystates0(this%v1eca_bgc_index%lid_co2_hr))*catomw/dtime

      biogeo_flux%c12flux_vars%phr_vr_col(c,j) = &
        (ystatesf(this%v1eca_bgc_index%lid_pot_co2_hr) - &
        ystates0(this%v1eca_bgc_index%lid_pot_co2_hr))*catomw/dtime

      biogeo_flux%c12flux_vars%somhr_vr_col(c,j) = &
        (ystatesf(this%v1eca_bgc_index%lid_co2_somhr) - &
        ystates0(this%v1eca_bgc_index%lid_co2_somhr))*catomw/dtime

      biogeo_flux%c12flux_vars%lithr_vr_col(c,j) = &
        (ystatesf(this%v1eca_bgc_index%lid_co2_lithr) - &
        ystates0(this%v1eca_bgc_index%lid_co2_lithr))*catomw/dtime

      biogeo_flux%c12flux_vars%cwdhr_vr_col(c,j) = &
        (ystatesf(this%v1eca_bgc_index%lid_co2_cwdhr) - &
        ystates0(this%v1eca_bgc_index%lid_co2_cwdhr))*catomw/dtime

      biogeo_flux%c12flux_vars%o_scalar_col(c,j) = ystatesf(this%v1eca_bgc_index%lid_o_scalar)

      biogeo_flux%c12flux_vars%decomp_k(c,j,1:this%v1eca_bgc_index%nom_pools) = biogeo_flux%c12flux_vars%decomp_k(c,j,1:this%v1eca_bgc_index%nom_pools) * &
          biogeo_flux%c12flux_vars%o_scalar_col(c,j)

      biogeo_flux%n14flux_vars%f_denit_vr_col(c,j)= &
        (ystatesf(this%v1eca_bgc_index%lid_no3_den) - &
         ystates0(this%v1eca_bgc_index%lid_no3_den))*natomw/dtime

      biogeo_flux%n14flux_vars%f_n2o_denit_vr_col(c,j)= &
        (ystatesf(this%v1eca_bgc_index%lid_n2o_den) - &
         ystates0(this%v1eca_bgc_index%lid_n2o_den))*natomw/dtime

      biogeo_flux%n14flux_vars%f_nit_vr_col(c,j) = &
        (ystatesf(this%v1eca_bgc_index%lid_nh4_nit) - &
         ystates0(this%v1eca_bgc_index%lid_nh4_nit))*natomw/dtime

      biogeo_flux%n14flux_vars%f_n2o_nit_vr_col(c,j) = &
        (ystatesf(this%v1eca_bgc_index%lid_n2o_nit) - &
         ystates0(this%v1eca_bgc_index%lid_n2o_nit))*natomw/dtime

      biogeo_flux%p31flux_vars%pflx_minp_weathering_po4_vr_col(c,j) =this%v1eca_forc(c,j)%sflx_minp_weathering_po4

      biogeo_flux%n14flux_vars%smin_nh4_to_plant_vr_col(c,j) = &
          (ystatesf(this%v1eca_bgc_index%lid_plant_minn_nh4) - &
          ystates0(this%v1eca_bgc_index%lid_plant_minn_nh4))*natomw/dtime

      biogeo_flux%n14flux_vars%smin_no3_to_plant_vr_col(c,j) = &
          (ystatesf(this%v1eca_bgc_index%lid_plant_minn_no3) - &
          ystates0(this%v1eca_bgc_index%lid_plant_minn_no3))*natomw/dtime

      biogeo_flux%n14flux_vars%smin_nh4_immob_vr_col(c,j) = &
          (ystatesf(this%v1eca_bgc_index%lid_minn_nh4_immob) - &
           ystates0(this%v1eca_bgc_index%lid_minn_nh4_immob))*natomw/dtime

      biogeo_flux%n14flux_vars%smin_no3_immob_vr_col(c,j) = &
          (ystatesf(this%v1eca_bgc_index%lid_minn_no3_immob) - &
           ystates0(this%v1eca_bgc_index%lid_minn_no3_immob))*natomw/dtime

      biogeo_flux%p31flux_vars%sminp_to_plant_vr_col(c,j)  = &
          (ystatesf(this%v1eca_bgc_index%lid_plant_minp) - &
           ystates0(this%v1eca_bgc_index%lid_plant_minp))*patomw/dtime

      select type(plant_soilbgc)
      type is(v1eca_plant_soilbgc_type)
        do p = 1, this%nactpft
          plant_soilbgc%plant_minn_no3_active_yield_flx_vr_patch(p,j) = &
            (ystatesf(this%v1eca_bgc_index%lid_plant_minn_no3_pft(p)) - &
            ystates0(this%v1eca_bgc_index%lid_plant_minn_no3_pft(p)))*natomw/dtime

          plant_soilbgc%plant_minn_nh4_active_yield_flx_vr_patch(p,j) = &
            (ystatesf(this%v1eca_bgc_index%lid_plant_minn_nh4_pft(p)) - &
            ystates0(this%v1eca_bgc_index%lid_plant_minn_nh4_pft(p)))*natomw/dtime

          plant_soilbgc%plant_minp_active_yield_flx_vr_patch(p,j) = &
            (ystatesf(this%v1eca_bgc_index%lid_plant_minp_pft(p)) - &
             ystates0(this%v1eca_bgc_index%lid_plant_minp_pft(p)))*patomw/dtime
        enddo
        if(this%nactpft>0)then
          plant_soilbgc%plant_minn_no3_active_yield_flx_vr_col(c,j) = &
            (ystatesf(this%v1eca_bgc_index%lid_plant_minn_no3) - &
            ystates0(this%v1eca_bgc_index%lid_plant_minn_no3))*natomw/dtime

          plant_soilbgc%plant_minn_nh4_active_yield_flx_vr_col(c,j) = &
            (ystatesf(this%v1eca_bgc_index%lid_plant_minn_nh4) - &
            ystates0(this%v1eca_bgc_index%lid_plant_minn_nh4))*natomw/dtime

          plant_soilbgc%plant_minp_active_yield_flx_vr_col(c,j) = &
            (ystatesf(this%v1eca_bgc_index%lid_plant_minp) - &
             ystates0(this%v1eca_bgc_index%lid_plant_minp))*patomw/dtime

          biogeo_flux%p31flux_vars%col_plant_pdemand_vr(c,j) = &
            plant_soilbgc%plant_minp_active_yield_flx_vr_col(c,j)
        else
          biogeo_flux%p31flux_vars%col_plant_pdemand_vr(c,j) = 0._r8
        endif

      end select


  end associate
  end subroutine retrieve_output

   !----------------------------------------------------------------------
   subroutine retrieve_lnd2atm(this, bounds, num_soilc, filter_soilc, tracerflux_vars, &
   betrtracer_vars, biogeo_flux)

   use tracerfluxType           , only : tracerflux_type
   use BeTR_decompMod           , only : betr_bounds_type
   use BeTRTracerType           , only : BeTRTracer_Type
   use tracer_varcon            , only : natomw, patomw
   implicit none
   class(v1eca_bgc_reaction_type) , intent(inout)    :: this
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
   subroutine debug_info(this, bounds, num_soilc, filter_soilc, dzsoi, betrtracer_vars, tracerstate_vars,  header, betr_status)

   use BeTR_decompMod           , only : betr_bounds_type
   use BeTRTracerType           , only : BeTRTracer_Type
   use tracerstatetype          , only : tracerstate_type
   use tracer_varcon            , only : catomw, natomw, patomw

   ! !ARGUMENTS:
   implicit none
   class(v1eca_bgc_reaction_type) , intent(inout)    :: this !
   type(betr_bounds_type)               , intent(in) :: bounds                      ! bounds
   integer                              , intent(in) :: num_soilc                   ! number of columns in column filter
   integer                              , intent(in) :: filter_soilc(:)             ! column filter
   real(r8)                             , intent(in) :: dzsoi(bounds%begc: ,bounds%lbj: )
   type(betrtracer_type)                , intent(in) :: betrtracer_vars             ! betr configuration information
   type(tracerstate_type)               , intent(in) :: tracerstate_vars
   character(len=*)                     , intent(in) :: header
   type(betr_status_type)               , intent(out):: betr_status
   integer :: fc, c
   integer :: c_loc, n_loc, p_loc, nelm, j, kk
   real(r8):: c_mass, n_mass, p_mass, minp, min_nh4,min_no3,p_massocl

   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(dzsoi)  == (/bounds%endc, bounds%ubj/)),   errMsg(mod_filename,__LINE__),betr_status)

   return
   write(*,*)trim(header)//': debug info c n p mass'

   c_loc=this%v1eca_bgc_index%c_loc
   n_loc=this%v1eca_bgc_index%n_loc
   p_loc=this%v1eca_bgc_index%p_loc
   nelm =this%v1eca_bgc_index%nelms
   c_mass = 0._r8; n_mass = 0._r8; p_mass = 0._r8; min_nh4=0._r8; min_no3=0._r8;minp=0._r8;p_massocl=0._r8
   do j = 1, bounds%ubj
     do fc = 1, num_soilc
        c = filter_soilc(fc)

        !add litter
        do kk = betrtracer_vars%id_trc_beg_litr, betrtracer_vars%id_trc_end_litr, nelm
          c_mass = c_mass  + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc) * dzsoi(c,j)

          n_mass = n_mass  + &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc) * dzsoi(c,j)

          p_mass = p_mass  + &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc) * dzsoi(c,j)

        enddo

        !add cwd
        do kk = betrtracer_vars%id_trc_beg_wood, betrtracer_vars%id_trc_end_wood, nelm
          c_mass = c_mass + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc) * dzsoi(c,j)

          n_mass = n_mass + &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc) * dzsoi(c,j)

          p_mass = p_mass + &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc) * dzsoi(c,j)

        enddo

        !add som

        do kk = betrtracer_vars%id_trc_beg_pom, betrtracer_vars%id_trc_end_pom, nelm
          c_mass = c_mass + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc) * dzsoi(c,j)

          n_mass = n_mass + &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc) * dzsoi(c,j)

          p_mass = p_mass + &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc) * dzsoi(c,j)

        enddo

        !Microbial biomass
        do kk = betrtracer_vars%id_trc_beg_Bm, betrtracer_vars%id_trc_end_Bm, nelm
          c_mass = c_mass + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc) * dzsoi(c,j)

          n_mass = n_mass + &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc) * dzsoi(c,j)

          p_mass = p_mass + &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc) * dzsoi(c,j)

        enddo

        do kk = betrtracer_vars%id_trc_beg_som, betrtracer_vars%id_trc_end_som, nelm
          c_mass = c_mass + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc) * dzsoi(c,j)

          n_mass = n_mass + &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc) * dzsoi(c,j)

          p_mass = p_mass + &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc) * dzsoi(c,j)

        enddo

        !non occluded phosphorus, soluble and adsorbed
        p_mass = p_mass + patomw * &
           (tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_beg_minp) + &
            tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_p_sol)) * dzsoi(c,j)

        !occluded
!        p_mass = p_mass + patomw * &
!           tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_end_minp) * dzsoi(c,j)

        p_massocl = p_massocl + patomw * &
           tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_end_minp) * dzsoi(c,j)

        !mineral nitrogen
!        n_mass = n_mass + natomw * &
!           (tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_nh3x) + &
!            tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_no3x)) * dzsoi(c,j)

        minp = minp + patomw * &
           (tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_beg_minp) + &
            tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_p_sol)) * dzsoi(c,j)
!        min_nh4=min_nh4+ natomw * &
!           tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_nh3x) * dzsoi(c,j)
        min_no3=min_no3+ natomw * &
           tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_no3x) * dzsoi(c,j)
     enddo
   enddo

   write(*,*)'c_mass    =', c_mass
   write(*,*)'n_mass    =', n_mass
   write(*,*)'p_mass    =', p_mass
   write(*,*)'p_massocl =', p_massocl
   write(*,*)'min_nh4   =', min_nh4
   write(*,*)'min_no3   =', min_no3
   write(*,*)'minp      =', minp
   write(*,*)'----------------------------------------'
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
   class(v1eca_bgc_reaction_type) , intent(inout)    :: this
   type(betr_bounds_type)               , intent(in)  :: bounds                      ! bounds
   integer                              , intent(in) :: lbj, ubj
   integer                              , intent(in) :: jtops(bounds%begc: )
   integer                              , intent(in)    :: num_soilc                   ! number of columns in column filter
   integer                              , intent(in)    :: filter_soilc(:)             ! column filter
   type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
   type(tracerstate_type)               , intent(inout) :: tracerstate_vars
   type(betr_biogeo_state_type)         , intent(inout) :: biogeo_state
   type(betr_status_type)               , intent(out)   :: betr_status

   integer :: nelm
   integer :: c_loc, c13_loc, c14_loc
   integer :: n_loc, p_loc
   integer :: c, fc, j, kk, kk1

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__),betr_status)

    c_loc=this%v1eca_bgc_index%c_loc
    n_loc=this%v1eca_bgc_index%n_loc
    p_loc=this%v1eca_bgc_index%p_loc
    c13_loc=this%v1eca_bgc_index%c13_loc
    c14_loc=this%v1eca_bgc_index%c14_loc
    nelm =this%v1eca_bgc_index%nelms

   do j = lbj, ubj
     do fc = 1, num_soilc
        c = filter_soilc(fc)
        if(j<jtops(c))cycle

        !add litter
        kk1=1
        do kk = betrtracer_vars%id_trc_beg_litr, betrtracer_vars%id_trc_end_litr, nelm
          biogeo_state%c12state_vars%decomp_cpools_vr(c,j,kk1) = &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          biogeo_state%n14state_vars%decomp_npools_vr(c,j,kk1)= &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc)

          biogeo_state%p31state_vars%decomp_ppools_vr(c,j,kk1) = &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%decomp_cpools_vr(c,j,kk1) =  &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%decomp_cpools_vr(c,j,kk1) = &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
          kk1 = kk1 + 1
        enddo

        !add cwd
        kk1=4
        do kk = betrtracer_vars%id_trc_beg_wood, betrtracer_vars%id_trc_end_wood, nelm
          biogeo_state%c12state_vars%decomp_cpools_vr(c,j,kk1) =  &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          biogeo_state%n14state_vars%decomp_npools_vr(c,j,kk1) =  &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc)

          biogeo_state%p31state_vars%decomp_ppools_vr(c,j,kk1) = &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%decomp_cpools_vr(c,j,kk1) =  &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%decomp_cpools_vr(c,j,kk1) = &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        !add som

        !Microbial biomass
        do kk = betrtracer_vars%id_trc_beg_Bm, betrtracer_vars%id_trc_end_Bm, nelm
          biogeo_state%c12state_vars%decomp_cpools_vr(c,j,5) =  &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          biogeo_state%n14state_vars%decomp_npools_vr(c,j,5) =  &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc)

          biogeo_state%p31state_vars%decomp_ppools_vr(c,j,5) =  &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%decomp_cpools_vr(c,j,5) = &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%decomp_cpools_vr(c,j,5) =  &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        do kk = betrtracer_vars%id_trc_beg_som, betrtracer_vars%id_trc_end_som, nelm
          biogeo_state%c12state_vars%decomp_cpools_vr(c,j,7) =  &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          biogeo_state%n14state_vars%decomp_npools_vr(c,j,7) =  &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc)

          biogeo_state%p31state_vars%decomp_ppools_vr(c,j,7) =  &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%decomp_cpools_vr(c,j,7) =  &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%decomp_cpools_vr(c,j,7) =  &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        !POM
        do kk = betrtracer_vars%id_trc_beg_pom, betrtracer_vars%id_trc_end_pom, nelm
          biogeo_state%c12state_vars%decomp_cpools_vr(c,j,6) = &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

          biogeo_state%n14state_vars%decomp_npools_vr(c,j,6) =  &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc)

          biogeo_state%p31state_vars%decomp_ppools_vr(c,j,6) =  &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc)

          if(this%use_c13)then
            biogeo_state%c13state_vars%decomp_cpools_vr(c,j,6) =  &
              c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
          endif

          if(this%use_c14)then
            biogeo_state%c14state_vars%decomp_cpools_vr(c,j,6) =  &
              c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
          endif
        enddo

        !soluble P
        biogeo_state%p31state_vars%sminp_vr_col(c,j) = patomw * &
            tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_p_sol)

        !mineral nitrogen

        biogeo_state%n14state_vars%sminn_no3_vr_col(c,j) = natomw * &
            tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_no3x)
     enddo
   enddo
   contains
     subroutine sum_totsom(c, j, ibeg, iend, nelm)
     implicit none
     integer , intent(in) :: c, j, ibeg, iend, nelm

     integer :: kk
     do kk = ibeg, iend, nelm
        biogeo_state%c12state_vars%totsomc_vr_col(c,j) = biogeo_state%c12state_vars%totsomc_vr_col(c,j) + &
            catomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c_loc)

        biogeo_state%n14state_vars%totsomn_vr_col(c,j) = biogeo_state%n14state_vars%totsomn_vr_col(c,j) + &
            natomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+n_loc)

        biogeo_state%p31state_vars%totsomp_vr_col(c,j) = biogeo_state%p31state_vars%totsomp_vr_col(c,j) + &
            patomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+p_loc)

        if(this%use_c13)then
          biogeo_state%c13state_vars%totsomc_vr_col(c,j) = biogeo_state%c13state_vars%totsomc_vr_col(c,j) + &
            c13atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c13_loc)
        endif

        if(this%use_c14)then
          biogeo_state%c14state_vars%totsomc_vr_col(c,j) = biogeo_state%c14state_vars%totsomc_vr_col(c,j) + &
            c14atomw * tracerstate_vars%tracer_conc_mobile_col(c, j, kk-1+c14_loc)
        endif
     enddo
     end subroutine sum_totsom
   end subroutine retrieve_biostates
   !----------------------------------------------------------------------
   function calc_latacc(lat)result(ans)
   implicit none
   real(r8), intent(in) :: lat
   real(r8) :: ans

   ans = 1._r8 + 50._r8/(1._r8+exp(-0.1_r8*(abs(lat)-60._r8)))
   end function calc_latacc



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
      class(v1eca_bgc_reaction_type) , intent(inout)    :: this
       type(betr_bounds_type)           , intent(in)    :: bounds                      ! bounds
       integer                          , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
       integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
       integer                          , intent(in)    :: filter_soilc(:)             ! column filter
       integer                          , intent(in)    :: jtops( : )                  ! top index of each column
       type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
       type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
       type(tracerstate_type)           , intent(inout) :: tracerstate_vars
       type(betr_status_type)           , intent(out)   :: betr_status

   integer :: fc, c, j, kk, kc, kn, kp

   associate(                                &
    c13_loc=>  this%v1eca_bgc_index%c13_loc, &
    c14_loc=>  this%v1eca_bgc_index%c14_loc, &
    c_loc=>  this%v1eca_bgc_index%c_loc    , &
    n_loc=>  this%v1eca_bgc_index%n_loc    , &
    p_loc=>  this%v1eca_bgc_index%p_loc    , &
    nelm => this%v1eca_bgc_index%nelms , &
    id_trc_no3x => betrtracer_vars%id_trc_no3x,  &
    id_trc_p_sol => betrtracer_vars%id_trc_p_sol, &
    id_trc_beg_litr => betrtracer_vars%id_trc_beg_litr, &
    id_trc_beg_wood  => betrtracer_vars%id_trc_beg_wood, &
    id_trc_beg_Bm => betrtracer_vars%id_trc_beg_Bm, &
    id_trc_beg_pom => betrtracer_vars%id_trc_beg_pom, &
    id_trc_beg_som => betrtracer_vars%id_trc_beg_som &
   )

    do j = lbj, ubj
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        if(j<jtops(c))cycle
          tracerstate_vars%tracer_conc_mobile_col(c, j, id_trc_no3x) =biophysforc%n14flx%in_sminn_no3_vr_col(c,j)/natomw
          tracerstate_vars%tracer_conc_mobile_col(c,j,id_trc_p_sol) =biophysforc%p31flx%in_sminp_vr_col(c,j)/patomw
      enddo
    enddo


    do j = lbj, ubj
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        kc = id_trc_beg_litr+c_loc-1;kn = id_trc_beg_litr+n_loc-1; kp = id_trc_beg_litr+p_loc-1; kk=1
        tracerstate_vars%tracer_conc_mobile_col(c, j, kc) = biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,kk)/catomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kn) = biophysforc%n14flx%in_decomp_npools_vr_col(c,j,kk)/natomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kp) = biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,kk)/patomw

        kc = id_trc_beg_litr+nelm*kk+c_loc-1;kn = id_trc_beg_litr+nelm*kk+n_loc-1; kp = id_trc_beg_litr+nelm*kk+p_loc-1; kk=2
        tracerstate_vars%tracer_conc_mobile_col(c, j, kc) = biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,kk)/catomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kn) = biophysforc%n14flx%in_decomp_npools_vr_col(c,j,kk)/natomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kp) = biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,kk)/patomw

        kc = id_trc_beg_litr+nelm*kk+c_loc-1;kn = id_trc_beg_litr+nelm*kk+n_loc-1; kp = id_trc_beg_litr+nelm*kk+p_loc-1; kk=3
        tracerstate_vars%tracer_conc_mobile_col(c, j, kc) = biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,kk)/catomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kn) = biophysforc%n14flx%in_decomp_npools_vr_col(c,j,kk)/natomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kp) = biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,kk)/patomw

        kc = id_trc_beg_wood+c_loc-1;kn = id_trc_beg_wood+n_loc-1; kp = id_trc_beg_wood+p_loc-1; kk=4
        tracerstate_vars%tracer_conc_mobile_col(c, j, kc) = biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,kk)/catomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kn) = biophysforc%n14flx%in_decomp_npools_vr_col(c,j,kk)/natomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kp) = biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,kk)/patomw

        kc = id_trc_beg_Bm+c_loc-1;kn = id_trc_beg_Bm+n_loc-1; kp = id_trc_beg_Bm+p_loc-1; kk=5
        tracerstate_vars%tracer_conc_mobile_col(c, j, kc) = biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,kk)/catomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kn) = biophysforc%n14flx%in_decomp_npools_vr_col(c,j,kk)/natomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kp) = biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,kk)/patomw

        kc = id_trc_beg_pom+c_loc-1;kn = id_trc_beg_pom+n_loc-1; kp = id_trc_beg_pom+p_loc-1; kk=6
        tracerstate_vars%tracer_conc_mobile_col(c, j, kc) = biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,kk)/catomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kn) = biophysforc%n14flx%in_decomp_npools_vr_col(c,j,kk)/natomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kp) = biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,kk)/patomw

        kc = id_trc_beg_som+c_loc-1;kn = id_trc_beg_som+n_loc-1; kp = id_trc_beg_som+p_loc-1; kk=6
        tracerstate_vars%tracer_conc_mobile_col(c, j, kc) = biophysforc%c12flx%in_decomp_cpools_vr_col(c,j,kk)/catomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kn) = biophysforc%n14flx%in_decomp_npools_vr_col(c,j,kk)/natomw
        tracerstate_vars%tracer_conc_mobile_col(c, j, kp) = biophysforc%p31flx%in_decomp_ppools_vr_col(c,j,kk)/patomw

      enddo
    enddo

   end associate
   end subroutine reset_biostates

   !----------------------------------------------------------------------
   subroutine SetParCols(this, parcol)
     implicit none
      class(v1eca_bgc_reaction_type) , intent(inout)    :: this
     integer, intent(in) :: parcol

     this%parcol = parcol
   end subroutine SetParCols
end module v1ecaBGCReactionsType
