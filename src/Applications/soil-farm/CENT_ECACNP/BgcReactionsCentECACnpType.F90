module BgcReactionsCentECACnpType

#include "bshr_assert.h"

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
  use bshr_log_mod           , only : errMsg => shr_log_errMsg
  use bshr_kind_mod          , only : r8 => shr_kind_r8
  use bshr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)
  use BeTR_decompMod         , only : bounds_type  => betr_bounds_type
  use BGCReactionsMod       , only : bgc_reaction_type
  use betr_varcon           , only : spval => bspval, ispval => bispval
  use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
  use BgcCentCnpType        , only : create_centuryeca_type, centurybgceca_type
  use BgcCentCnpForcType    , only : centuryeca_forc_type, create_century_forc_type
  use BetrStatusType        , only : betr_status_type
  use BiogeoConType         , only : BiogeoCon_type
  use BgcCentCnpIndexType   , only : centurybgc_index_type
  use BiogeoConType         , only : bgc_con_eca
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
    bgc_reaction_CENTURY_ECACNP_type
    type(centurybgceca_type), private, pointer :: centuryeca(:,:)
    type(centuryeca_forc_type), private, pointer :: centuryforc(:,:)

    type(centurybgc_index_type),private :: centurybgc_index
    logical, private :: use_c13
    logical, private :: use_c14
  contains
    procedure :: Init_betrbgc                 ! initialize betr bgc
    procedure :: set_boundary_conditions      ! set top/bottom boundary conditions for various tracers
    procedure :: calc_bgc_reaction            ! doing bgc calculation
    procedure :: init_boundary_condition_type ! initialize type of top boundary conditions
    procedure :: do_tracer_equilibration      ! do equilibrium tracer chemistry
    procedure :: initCold
    procedure :: readParams
    procedure :: lsm_betr_flux_state_receive
    procedure, private:: set_century_forc
    procedure, private:: retrieve_output
  end type bgc_reaction_CENTURY_ECACNP_type

  interface bgc_reaction_CENTURY_ECACNP_type
     module procedure constructor
  end interface bgc_reaction_CENTURY_ECACNP_type

contains


  !-------------------------------------------------------------------------------
  type(bgc_reaction_CENTURY_ECACNP_type) function constructor()
    !
    ! ! DESCRIPTION:
    !
    ! create an object of type bgc_reaction_CENTURY_ECACNP_type.
    ! Right now it is purposely empty
   type(bgc_reaction_CENTURY_ECACNP_type), allocatable :: bgc

   allocate(bgc)
   constructor = bgc
  end function constructor


  !-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
    !
    ! DESCRIPTION:
    ! initialize boundary condition types
    ! !USES:
    use TracerBoundaryCondType      , only : tracerboundarycond_type
    use BeTRTracerType              , only : betrtracer_type

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECACNP_type), intent(inout) :: this
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

  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars, namelist_buffer, bstatus)
    !
    ! DESCRIPTION:
    ! initialize the betrbgc
    !
    ! !USES:
    use BeTRTracerType                   , only : betrtracer_type
    use BeTRTracerType                   , only : betrtracer_type
    use MathfuncMod                      , only : addone
    use betr_varcon         , only : betr_maxpatch_pft
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout)    :: this
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
    integer   :: c, j, litr_cnt, wood_cnt, Bm_cnt, som_cnt, itemp_ads, itemp_ads_grp
    logical   :: carbon_only = .false.

    call bstatus%reset()

    if (this%dummy_compiler_warning) continue
    !initialize model parameters
    call bgc_con_eca%Init(namelist_buffer, bstatus)

    call this%centurybgc_index%Init(bgc_con_eca%use_c13, bgc_con_eca%use_c14, betr_maxpatch_pft)

    if(bstatus%check_status())return

    !create the models
    allocate(this%centuryeca(bounds%begc:bounds%endc,lbj:ubj), source=create_centuryeca_type())

    !create model specific forcing data structure
    allocate(this%centuryforc(bounds%begc:bounds%endc,lbj:ubj), source=create_century_forc_type())

    !initialize
    do j = lbj, ubj
      do c = bounds%begc, bounds%endc
        call this%centuryeca(c,j)%Init(bgc_con_eca, bstatus)
        if(bstatus%check_status())return
        call this%centuryforc(c,j)%Init(this%centurybgc_index)

      enddo
    enddo
    this%use_c13 = bgc_con_eca%use_c13
    this%use_c14 = bgc_con_eca%use_c14
    print*,'l;dgladagausexx 13 14',this%use_c13, this%use_c14
    !set up betr
    nelm =this%centurybgc_index%nelms
    c_loc=this%centurybgc_index%c_loc
    n_loc=this%centurybgc_index%n_loc
    p_loc=this%centurybgc_index%p_loc
    c13_loc=this%centurybgc_index%c13_loc
    c14_loc=this%centurybgc_index%c14_loc
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

    !non-volatile tracers
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),  mem = 1, &
      trc_cnt = itemp_trc, trc_grp=betrtracer_vars%id_trc_no3x, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_no3x, &
      trc_grp_end=betrtracer_vars%id_trc_end_no3x, &
      is_trc_gw=.true., is_trc_volatile = .false.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_nh3x, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_nh3x, &
      trc_grp_end=betrtracer_vars%id_trc_end_nh3x, &
      is_trc_gw=.true., is_trc_volatile = .false.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_p_sol, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_p_sol, &
      trc_grp_end=betrtracer_vars%id_trc_end_p_sol, &
      is_trc_gw=.true., is_trc_volatile = .false.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = nelm, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_dom, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_dom, &
      trc_grp_end=betrtracer_vars%id_trc_end_dom, &
      is_trc_gw=.true., is_trc_volatile = .false.)

    !three litter groups
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = nelm*3, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_litr, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_litr, &
      trc_grp_end=betrtracer_vars%id_trc_end_litr, &
      is_trc_passive=.true.)

    !three woody groups
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), &
      is_trc_passive=.true., mem = nelm*3, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_wood, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_wood, &
      trc_grp_end=betrtracer_vars%id_trc_end_wood)

    !group of microbial biomass
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),mem = nelm, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_som, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_Bm, &
      trc_grp_end=betrtracer_vars%id_trc_end_Bm, &
      is_trc_passive=.true.)

    !group of som, which is not dom or microbial biomass
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp),mem = nelm, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_som, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_som, &
      trc_grp_end=betrtracer_vars%id_trc_end_som, &
      is_trc_passive=.true.)
    !group of solid phase mineral phosphorus
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 2, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_minp, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_minp, &
      trc_grp_end=betrtracer_vars%id_trc_end_minp, &
      is_trc_passive=.false.)

    betrtracer_vars%nmem_max                     = nelm*3                       ! maximum number of group elements

    call betrtracer_vars%Init()
    betrtracer_vars%is_mobile(:) = .true.

    itemp_vgrp = 0  !counter for volatile groups
    itemp_v    = 0  !counter for volatile tracers
    itemp_ads_grp =0!counter for sorptive groups
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

    print*,'use_c13xxx',this%use_c13
    print*,'use_c14xxx',this%use_c14

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

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_nh3x, &
         trc_name='NH3x', is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_nh3x, trc_group_mem = 1, is_trc_volatile=.false.)
    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_no3x, &
         trc_name='NO3x', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_no3x, trc_group_mem = 1, is_trc_volatile=.false., &
         trc_vtrans_scal=1._r8)
    if(bstatus%check_status())return

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_p_sol, &
         trc_name='P_SOL', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_p_sol, trc_group_mem = 1, is_trc_volatile=.false., &
         trc_vtrans_scal=1._r8)
    if(bstatus%check_status())return

    itemp_mem=0
    trcid =  betrtracer_vars%id_trc_beg_dom+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, &
         trc_name='SOM2C', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem),&
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads), &
         trc_adsorbgroupid=addone(itemp_ads_grp), trc_family_name='DOM')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_dom+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, &
         trc_name='SOM2N', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads), &
         trc_adsorbgroupid=itemp_ads_grp, trc_family_name='DOM')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_dom+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, &
         trc_name='SOM2P', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads), &
         trc_adsorbgroupid=itemp_ads_grp, trc_family_name='DOM')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_dom+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, &
         trc_name='SOM2C_C13', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads), &
         trc_adsorbgroupid=itemp_ads_grp, trc_family_name='DOM')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid=betrtracer_vars%id_trc_beg_dom+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, &
         trc_name='SOM2C_C14', is_trc_mobile=.true., is_trc_advective = .true., &
         trc_group_id = betrtracer_vars%id_trc_dom, trc_group_mem = addone(itemp_mem), &
         is_trc_volatile=.false., is_trc_adsorb = .true., trc_adsorbid=addone(itemp_ads), &
         trc_adsorbgroupid=itemp_ads_grp, trc_family_name='DOM')
      if(bstatus%check_status())return
    endif
    !------------------------------------------------------------------------------------
    !only one group passive solid litter tracers
    !define litter group
    itemp_mem=0
    litr_cnt = 0
    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1C' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1N' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT1')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT1P' ,    &
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
    !------------------------------------------------------------------------------------

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2C'  ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2N' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT2')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT2P' ,    &
         is_trc_mobile=.true., is_trc_advective = .false.,  &
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
    !------------------------------------------------------------------------------------

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3C' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3N' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_litr, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LIT3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_litr+litr_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LIT3P' ,    &
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
    !------------------------------------------------------------------------------------
    !define woody group
    wood_cnt=0
    itemp_mem=0
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

    wood_cnt=wood_cnt+1
    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LWDC' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LWD')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LWDN' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LWD')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LWDP' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LWD')
    if(bstatus%check_status())return
    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LWDC_C13' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LWD')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='LWDC_C14' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='LWD')
      if(bstatus%check_status())return
    endif

    wood_cnt=wood_cnt+1
    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='FWDC' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='FWD')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='FWDN' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='FWD')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='FWDP' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='FWD')
    if(bstatus%check_status())return
    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='FWDC_C13' ,    &
         is_trc_mobile=.false., is_trc_advective = .false.,&
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='FWD')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_wood+wood_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='FWDC_C14' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_wood, trc_group_mem= addone(itemp_mem), &
         trc_family_name='FWD')
      if(bstatus%check_status())return
    endif

    !------------------------------------------------------------------------------------
    !define som group
    Bm_cnt=0;itemp_mem=0
    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1C', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1')
    if(bstatus%check_status())return


    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1N', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1P', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1')
    if(bstatus%check_status())return

    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1C_C13',&
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1')
      if(bstatus%check_status())return
    endif

    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_Bm+Bm_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM1C_C14', &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som,  trc_group_mem = addone(itemp_mem), &
         trc_family_name='SOM1')
      if(bstatus%check_status())return
    endif
    !------------------------------------------------------------------------------------
    !new group
    som_cnt = 0; itemp_mem=0
    trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+c_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3C' ,     &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+n_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3N'  ,    &
         is_trc_mobile=.true., is_trc_advective = .false.,  &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+p_loc-1
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3P' ,    &
         is_trc_mobile=.true., is_trc_advective = .false.,  &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
    if(bstatus%check_status())return
    if(this%use_c13)then
      trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+c13_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3C_C13' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
      if(bstatus%check_status())return
    endif
    if(this%use_c14)then
      trcid = betrtracer_vars%id_trc_beg_som+som_cnt*nelm+c14_loc-1
      call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='SOM3C_C14' ,    &
         is_trc_mobile=.true., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_som, trc_group_mem= addone(itemp_mem), &
         trc_family_name='SOM3')
      if(bstatus%check_status())return
    endif

    !------------------------------------------------------------------------------------
    !new group
    itemp_mem=0

    trcid = betrtracer_vars%id_trc_beg_minp
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='P_2ND' ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_minp,  trc_group_mem= addone(itemp_mem), &
         trc_family_name='MINP')
    if(bstatus%check_status())return

    trcid = betrtracer_vars%id_trc_end_minp;
    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = trcid, trc_name='P_OCL'  ,    &
         is_trc_mobile=.false., is_trc_advective = .false., &
         trc_group_id = betrtracer_vars%id_trc_minp,  trc_group_mem= addone(itemp_mem), &
         trc_family_name='MINP')

  end subroutine Init_betrbgc

  !-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, dz_top, betrtracer_vars, &
       biophysforc, biogeo_flux, tracerboundarycond_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! set up boundary conditions for tracer movement
    !
    ! !USES:
    use TracerBoundaryCondType, only : tracerboundarycond_type
    use BeTRTracerType        , only : betrtracer_type
    use BeTR_biogeoFluxType   , only : betr_biogeo_flux_type
    use BetrStatusType        , only : betr_status_type
    use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout)    :: this
    type(bounds_type)                    , intent(in)    :: bounds
    integer                              , intent(in)    :: num_soilc               ! number of columns in column filter
    integer                              , intent(in)    :: filter_soilc(:)         ! column filter
    type(betrtracer_type)                , intent(in)    :: betrtracer_vars
    real(r8)                             , intent(in)    :: dz_top(bounds%begc: )
    type(betr_biogeophys_input_type)     , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)          , intent(in)    :: biogeo_flux
    type(tracerboundarycond_type)        , intent(inout) :: tracerboundarycond_vars !
    type(betr_status_type)               , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'set_boundary_conditions'
    integer :: fc, c

    call betr_status%reset()
    SHR_ASSERT_ALL(ubound(dz_top)  == (/bounds%endc/),   errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    if (this%dummy_compiler_warning) continue
    associate(                                       &
         groupid  => betrtracer_vars%groupid    ,    &
         ngwmobile_tracers => betrtracer_vars%ngwmobile_tracers, &
         ntracers => betrtracer_vars%ntracers  &
         )

      do fc = 1, num_soilc
         c = filter_soilc(fc)

         !values below will be updated with datastream
         !eventually, the following code will be implemented using polymorphism
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,1:ntracers)                   =0._r8                        !zero incoming flux
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2)    =32.8_r8                      !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o2)    =8.78_r8                      !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ar)    =0.3924_r8                    !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_co2x)  =0.0168_r8                    !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ch4)   =6.939e-5_r8                  !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2o)   =1.195e-5_r8                  !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_no3x)  = 0._r8
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_p_sol) = 0._r8

         tracerboundarycond_vars%bot_concflux_col(c,1,:)                                          = 0._r8                       !zero flux boundary condition
         tracerboundarycond_vars%condc_toplay_col(c,:) = 0._r8                                                                  !those will be updated with snow resistance and hydraulic wicking resistance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_n2))    = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_o2))    = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ar))    = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_co2x))  = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ch4))   = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_n2o))   = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance

         if(this%use_c13)then
           tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_c13_co2x)  =0.0168_r8
           tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_c13_co2x))  = 2._r8*1.267e-5_r8/dz_top(c)
         endif
         if(this%use_c14)then
           tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_c14_co2x)  =0.0168_r8
           tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_c14_co2x))  = 2._r8*1.267e-5_r8/dz_top(c)
         endif

      enddo
    end associate
  end subroutine set_boundary_conditions
  !-------------------------------------------------------------------------------

  subroutine calc_bgc_reaction(this, bounds, col, lbj, ubj, num_soilc, filter_soilc, &
       num_soilp,filter_soilp, jtops, dtime, betrtracer_vars, tracercoeff_vars, biophysforc,    &
       tracerstate_vars, tracerflux_vars, tracerboundarycond_vars, plant_soilbgc, biogeo_flux, betr_status)

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
    use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
    implicit none
    ! !ARGUMENTS
    class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout) :: this
    type(bounds_type)                    , intent(in) :: bounds                        ! bounds
    type(betr_column_type)               , intent(in) :: col
    integer                              , intent(in) :: num_soilc                     ! number of columns in column filter
    integer                              , intent(in) :: filter_soilc(:)               ! column filter
    integer                              , intent(in) :: num_soilp
    integer                              , intent(in) :: filter_soilp(:)               ! pft filter
    integer                              , intent(in) :: jtops(bounds%begc: )          ! top index of each column
    integer                              , intent(in) :: lbj, ubj                      ! lower and upper bounds, make sure they are > 0
    real(r8)                             , intent(in) :: dtime                         ! model time step
    type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
    type(tracercoeff_type)               , intent(in) :: tracercoeff_vars
    type(betr_biogeophys_input_type)     , intent(in)    :: biophysforc
    type(tracerboundarycond_type)        , intent(inout) :: tracerboundarycond_vars !
    type(tracerstate_type)               , intent(inout) :: tracerstate_vars
    type(tracerflux_type)                , intent(inout) :: tracerflux_vars
    class(plant_soilbgc_type)            , intent(inout) :: plant_soilbgc
    type(betr_biogeo_flux_type)          , intent(inout) :: biogeo_flux
    type(betr_status_type)               , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    character(len=32), parameter :: subname ='calc_bgc_reaction'
    integer                      :: fc, c, j, k
    logical :: is_surf  !surface litter layer?
    integer :: nstates
    real(r8), allocatable :: ystates0(:)
    real(r8), allocatable :: ystatesf(:)

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    if (this%dummy_compiler_warning) continue

    nstates = this%centurybgc_index%nstvars
    allocate(ystates0(nstates))
    allocate(ystatesf(nstates))

    !pass in fluxes and state varaibles into the 1D soil bgc model
    call this%set_century_forc(bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
        biophysforc, plant_soilbgc, betrtracer_vars, tracercoeff_vars, tracerstate_vars)

    !run simulation layer by layer
    do j = lbj, ubj
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        if(j<jtops(c))cycle
        is_surf=(j<=0)
        !do century eca bgc simulation
!        print*,'-----------------------------------'
!        print*,'lay',j
        call this%centuryeca(c,j)%runbgc(is_surf, dtime, this%centuryforc(c,j),nstates, ystates0, ystatesf, betr_status)
        call this%retrieve_output(c, j, nstates, ystates0, ystatesf, dtime, betrtracer_vars, tracerflux_vars,&
           tracerstate_vars, plant_soilbgc, biogeo_flux)
      enddo
    enddo
!    print*,''

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
    use tracerstatetype       , only : tracerstate_type
    use tracercoeffType       , only : tracercoeff_type
    use BeTRTracerType        , only : betrtracer_type
    use BetrStatusType        , only : betr_status_type
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECACNP_type), intent(inout)    :: this
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
    if(betr_status%check_status())return
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
    class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout)    :: this
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
    use BeTR_decompMod    , only : betr_bounds_type
    use BeTRTracerType    , only : BeTRTracer_Type
    use tracerstatetype   , only : tracerstate_type
    use betr_varcon       , only : spval => bspval, ispval => bispval
    use betr_varcon       , only : denh2o => bdenh2o
    use tracer_varcon     , only : nlevtrc_soil  => betr_nlevtrc_soil
    use betr_columnType   , only : betr_column_type
    use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout)    :: this
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

    associate(                                    &
         volatileid => betrtracer_vars%volatileid &
         )
    do c = bounds%begc, bounds%endc

      !dual phase tracers

      tracerstate_vars%tracer_conc_mobile_col    (c,:, :                                   )  = 0._r8
      tracerstate_vars%tracer_conc_frozen_col    (c,:, :                                   )  = 0._r8
      tracerstate_vars%tracer_conc_surfwater_col (c,:                                      )  = 0._r8
      tracerstate_vars%tracer_conc_aquifer_col   (c,:                                      )  = 0._r8
      tracerstate_vars%tracer_conc_grndwater_col (c,:                                      )  = 0._r8
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_n2   )) = 32.8_r8
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_o2   )) = 8.78_r8
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_ar   )) = 0.3924_r8
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_co2x )) = 0.0168_r8
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_ch4  )) = 6.939e-5_r8
      tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_n2o  )) = 1.195e-5_r8

      if(betrtracer_vars%nsolid_equil_tracers>0)then
        tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = 0._r8
      endif
      tracerstate_vars%tracer_soi_molarmass_col(c,:)          = 0._r8

    enddo
    end associate
  end subroutine InitCold

  !------------------------------------------------------------------------------
  subroutine lsm_betr_flux_state_receive(this, bounds, num_soilc, filter_soilc, &
        tracerstate_vars, tracerflux_vars,  betrtracer_vars)

    !
    ! !DESCRIPTION:
    ! do state and flux variable exchange between betr and lsm
    ! !USES:
    use bshr_kind_mod             , only : r8 => shr_kind_r8
    use tracerfluxType           , only : tracerflux_type
    use tracerstatetype          , only : tracerstate_type
    use betr_decompMod           , only : betr_bounds_type
    use BeTRTracerType           , only : BeTRTracer_Type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout)    :: this
    type(betr_bounds_type)                    , intent(in)    :: bounds          ! bounds
    integer                              , intent(in)    :: num_soilc       ! number of columns in column filter
    integer                              , intent(in)    :: filter_soilc(:) ! column filter
    type(betrtracer_type)                , intent(in)    :: betrtracer_vars ! betr configuration information
    type(tracerstate_type)               , intent(in)    :: tracerstate_vars
    type(tracerflux_type)                , intent(in)    :: tracerflux_vars

  !x  call assign_nitrogen_hydroloss(bounds, num_soilc, filter_soilc, &
  !x       tracerflux_vars, nitrogenflux_vars, phosphorusflux_vars, &
  !x       betrtracer_vars)

  !x  call assign_OM_CNpools(bounds, num_soilc, filter_soilc, &
  !x       carbonstate_vars, nitrogenstate_vars, phosphorusstate_vars, &
  !x       tracerstate_vars, betrtracer_vars, centurybgc_index)

  end subroutine lsm_betr_flux_state_receive

  !------------------------------------------------------------------------------
  subroutine set_century_forc(this, bounds, col, lbj, ubj, jtops, num_soilc, filter_soilc, &
      biophysforc, plant_soilbgc, betrtracer_vars, tracercoeff_vars, tracerstate_vars)
  !DESCRIPTION
  !set up forcing for running bgc
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use tracerstatetype          , only : tracerstate_type
  use betr_decompMod           , only : betr_bounds_type
  use tracercoeffType          , only : tracercoeff_type
  use betr_columnType          , only : betr_column_type
  use BetrTracerType           , only : betrtracer_type
  use PlantSoilBgcCnpType      , only : plant_soilbgc_cnp_type
  use MathfuncMod              , only : fpmax
  use betr_varcon              , only : grav => bgrav
  implicit none
  class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout)    :: this
  type(bounds_type)                    , intent(in) :: bounds                         ! bounds
  type(betr_column_type)               , intent(in) :: col
  integer                              , intent(in) :: jtops(bounds%begc:bounds%endc) ! top index of each column
  integer                              , intent(in) :: lbj, ubj                       ! lower and upper bounds, make sure they are > 0
  integer                              , intent(in) :: num_soilc       ! number of columns in column filter
  integer                              , intent(in) :: filter_soilc(:) ! column filter
  type(betr_biogeophys_input_type)     , intent(in) :: biophysforc
  class(plant_soilbgc_type)            , intent(in) :: plant_soilbgc
  type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
  type(tracerstate_type)               , intent(in) :: tracerstate_vars
  type(tracercoeff_type)               , intent(in) :: tracercoeff_vars

  integer :: j, fc, c
  integer :: k1, k2
  associate( &
     litr_beg =>  this%centurybgc_index%litr_beg  , &
     litr_end =>  this%centurybgc_index%litr_end  , &
     wood_beg =>  this%centurybgc_index%wood_beg  , &
     wood_end =>  this%centurybgc_index%wood_end  , &
     som_beg =>  this%centurybgc_index%som_beg    , &
     som_end =>  this%centurybgc_index%som_end    , &
     dom_beg =>  this%centurybgc_index%dom_beg    , &
     dom_end =>  this%centurybgc_index%dom_end    , &
     Bm_beg  =>  this%centurybgc_index%Bm_beg     , &
     Bm_end  =>  this%centurybgc_index%Bm_end       &
  )
  do j = lbj, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      if(j<jtops(c))cycle
      this%centuryforc(c,j)%ystates(:) = 0._r8
      !litter
      this%centuryforc(c,j)%ystates(litr_beg:litr_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_litr:betrtracer_vars%id_trc_end_litr)
!      if(j==1)then
!        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!        print*,'inital state litter'
!        print*,this%centuryforc(c,j)%ystates(litr_beg:litr_end)
!      endif
      !wood
      this%centuryforc(c,j)%ystates(wood_beg:wood_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_wood:betrtracer_vars%id_trc_end_wood)
      !microbial biomass
      this%centuryforc(c,j)%ystates(Bm_beg:Bm_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_Bm:betrtracer_vars%id_trc_end_Bm)
      !som
      this%centuryforc(c,j)%ystates(som_beg:som_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_som:betrtracer_vars%id_trc_end_som)
      !dom
      this%centuryforc(c,j)%ystates(dom_beg:dom_end)= &
          tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_dom:betrtracer_vars%id_trc_end_dom)

      !non-soluble phase of mineral p
      k1= betrtracer_vars%id_trc_beg_minp; k2 = this%centurybgc_index%lid_minp_secondary
      this%centuryforc(c,j)%ystates(k2) = fpmax(tracerstate_vars%tracer_conc_mobile_col(c,j,k1))
      k1 = betrtracer_vars%id_trc_end_minp;   k2 = this%centurybgc_index%lid_minp_occlude
      this%centuryforc(c,j)%ystates(k2) = fpmax(tracerstate_vars%tracer_conc_mobile_col(c,j,k1))

      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_n2) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2))
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_o2) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2))
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_ar) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar))
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_co2)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_c13_co2)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_c14_co2)= &
          fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x))
      endif
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_ch4)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4))
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_nh4)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_nh3x))
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_no3)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_no3x))
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_n2o)= &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2o))
      this%centuryforc(c,j)%ystates(this%centurybgc_index%lid_minp_soluble) = &
           fpmax(tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_p_sol))

      !input
      this%centuryforc(c,j)%cflx_input_litr_met = biophysforc%c12flx%cflx_input_litr_met_vr_col(c,j)
      this%centuryforc(c,j)%cflx_input_litr_cel = biophysforc%c12flx%cflx_input_litr_cel_vr_col(c,j)
      this%centuryforc(c,j)%cflx_input_litr_lig = biophysforc%c12flx%cflx_input_litr_lig_vr_col(c,j)
      this%centuryforc(c,j)%cflx_input_litr_cwd = biophysforc%c12flx%cflx_input_litr_cwd_vr_col(c,j)

      this%centuryforc(c,j)%nflx_input_litr_met = biophysforc%n14flx%nflx_input_litr_met_vr_col(c,j)
      this%centuryforc(c,j)%nflx_input_litr_cel = biophysforc%n14flx%nflx_input_litr_cel_vr_col(c,j)
      this%centuryforc(c,j)%nflx_input_litr_lig = biophysforc%n14flx%nflx_input_litr_lig_vr_col(c,j)
      this%centuryforc(c,j)%nflx_input_litr_cwd = biophysforc%n14flx%nflx_input_litr_cwd_vr_col(c,j)

      this%centuryforc(c,j)%pflx_input_litr_met = biophysforc%p31flx%pflx_input_litr_met_vr_col(c,j)
      this%centuryforc(c,j)%pflx_input_litr_cel = biophysforc%p31flx%pflx_input_litr_cel_vr_col(c,j)
      this%centuryforc(c,j)%pflx_input_litr_lig = biophysforc%p31flx%pflx_input_litr_lig_vr_col(c,j)
      this%centuryforc(c,j)%pflx_input_litr_cwd = biophysforc%p31flx%pflx_input_litr_cwd_vr_col(c,j)
      !loss
      this%centuryforc(c,j)%cflx_output_litr_met= biophysforc%c12flx%cflx_output_litr_met_vr_col(c,j)
      this%centuryforc(c,j)%cflx_output_litr_cel= biophysforc%c12flx%cflx_output_litr_cel_vr_col(c,j)
      this%centuryforc(c,j)%cflx_output_litr_lig= biophysforc%c12flx%cflx_output_litr_lig_vr_col(c,j)
      this%centuryforc(c,j)%cflx_output_litr_cwd= biophysforc%c12flx%cflx_output_litr_cwd_vr_col(c,j)
      this%centuryforc(c,j)%cflx_output_litr_fwd= biophysforc%c12flx%cflx_output_litr_fwd_vr_col(c,j)
      this%centuryforc(c,j)%cflx_output_litr_lwd= biophysforc%c12flx%cflx_output_litr_lwd_vr_col(c,j)

      this%centuryforc(c,j)%nflx_output_litr_met= biophysforc%n14flx%nflx_output_litr_met_vr_col(c,j)
      this%centuryforc(c,j)%nflx_output_litr_cel= biophysforc%n14flx%nflx_output_litr_cel_vr_col(c,j)
      this%centuryforc(c,j)%nflx_output_litr_lig= biophysforc%n14flx%nflx_output_litr_lig_vr_col(c,j)
      this%centuryforc(c,j)%nflx_output_litr_cwd= biophysforc%n14flx%nflx_output_litr_cwd_vr_col(c,j)
      this%centuryforc(c,j)%nflx_output_litr_fwd= biophysforc%n14flx%nflx_output_litr_fwd_vr_col(c,j)
      this%centuryforc(c,j)%nflx_output_litr_lwd= biophysforc%n14flx%nflx_output_litr_lwd_vr_col(c,j)

      this%centuryforc(c,j)%pflx_output_litr_met= biophysforc%p31flx%pflx_output_litr_met_vr_col(c,j)
      this%centuryforc(c,j)%pflx_output_litr_cel= biophysforc%p31flx%pflx_output_litr_cel_vr_col(c,j)
      this%centuryforc(c,j)%pflx_output_litr_lig= biophysforc%p31flx%pflx_output_litr_lig_vr_col(c,j)
      this%centuryforc(c,j)%pflx_output_litr_cwd= biophysforc%p31flx%pflx_output_litr_cwd_vr_col(c,j)
      this%centuryforc(c,j)%pflx_output_litr_fwd= biophysforc%p31flx%pflx_output_litr_fwd_vr_col(c,j)
      this%centuryforc(c,j)%pflx_output_litr_lwd= biophysforc%p31flx%pflx_output_litr_lwd_vr_col(c,j)


      this%centuryforc(c,j)%sflx_minn_input_nh4 = biophysforc%n14flx%nflx_minn_input_nh4_vr_col(c,j)     !nh4 from deposition and fertilization
      this%centuryforc(c,j)%sflx_minn_nh4_fix_nomic = biophysforc%n14flx%nflx_minn_nh4_fix_nomic_vr_col(c,j)       !nh4 from fixation
      this%centuryforc(c,j)%sflx_minp_input_po4 = biophysforc%p31flx%pflx_minp_input_po4_vr_col(c,j)     !inorganic P from deposition and fertilization
      this%centuryforc(c,j)%sflx_minp_weathering_po4 = biophysforc%p31flx%pflx_minp_weathering_po4_vr_col(c,j)

      this%centuryforc(c,j)%temp   = biophysforc%t_soisno_col(c,j)            !temperature
      this%centuryforc(c,j)%depz   = col%z(c,j)            !depth of the soil
      this%centuryforc(c,j)%dzsoi  = col%dz(c,j)            !soil thickness
      this%centuryforc(c,j)%sucsat  = biophysforc%sucsat_col(c,j)            ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
      this%centuryforc(c,j)%soilpsi = max(biophysforc%smp_l_col(c,j)*grav*1.e-6_r8,-15._r8)    ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
      this%centuryforc(c,j)%bsw = biophysforc%bsw_col(c,j)
      this%centuryforc(c,j)%bd   = biophysforc%bd_col(c,j)              !bulk density
      this%centuryforc(c,j)%pct_sand = biophysforc%cellsand_col(c,j)
      this%centuryforc(c,j)%h2osoi_vol = biophysforc%h2osoi_vol_col(c,j)
      this%centuryforc(c,j)%h2osoi_liq = biophysforc%h2osoi_liq_col(c,j)
      this%centuryforc(c,j)%air_vol = biophysforc%air_vol_col(c,j)
      this%centuryforc(c,j)%finundated = biophysforc%finundated_col(c)
      this%centuryforc(c,j)%watsat = biophysforc%watsat_col(c,j)
      this%centuryforc(c,j)%watfc = biophysforc%watfc_col(c,j)
      this%centuryforc(c,j)%cellorg = biophysforc%cellorg_col(c,j)
      this%centuryforc(c,j)%pH = biophysforc%soil_pH(c,j)

      this%centuryforc(c,j)%aren_cond_n2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2))
      this%centuryforc(c,j)%aren_cond_o2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2))
      this%centuryforc(c,j)%aren_cond_n2o = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2o)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2o))
      this%centuryforc(c,j)%aren_cond_co2 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%centuryforc(c,j)%aren_cond_co2_c13 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c13_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%centuryforc(c,j)%aren_cond_co2_c14 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c14_co2x)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%centuryforc(c,j)%aren_cond_ar = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar))
      this%centuryforc(c,j)%aren_cond_ch4 = &
          tracercoeff_vars%aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4)) * &
          tracercoeff_vars%scal_aere_cond_col(c,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4))

      this%centuryforc(c,j)%ch4_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ch4))
      this%centuryforc(c,j)%co2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_co2x))
      this%centuryforc(c,j)%o2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_o2))
      this%centuryforc(c,j)%n2_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2))
      this%centuryforc(c,j)%ar_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_ar))
      this%centuryforc(c,j)%n2o_g2b = &
          tracercoeff_vars%gas2bulkcef_mobile_col(c,j,betrtracer_vars%volatilegroupid(betrtracer_vars%id_trc_n2o))
      this%centuryforc(c,j)%o2_w2b = &
          tracercoeff_vars%aqu2bulkcef_mobile_col(c,j,betrtracer_vars%groupid(betrtracer_vars%id_trc_o2))            !conversion parameter for o2 from aqueous to bulk conc

      this%centuryforc(c,j)%conc_atm_n2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_n2))  !n2 concentration in atmosphere, mol n2/m3
      this%centuryforc(c,j)%conc_atm_n2o= &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_n2o))
      this%centuryforc(c,j)%conc_atm_o2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_o2))
      this%centuryforc(c,j)%conc_atm_ar = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_ar))
      this%centuryforc(c,j)%conc_atm_co2 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_co2x))
      if(this%use_c13)then
        this%centuryforc(c,j)%conc_atm_co2_c13 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_c13_co2x))
      endif
      if(this%use_c14)then
        this%centuryforc(c,j)%conc_atm_co2_c14 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_c14_co2x))
      endif
      this%centuryforc(c,j)%conc_atm_ch4 = &
          tracerstate_vars%tracer_conc_atm_col(c,betrtracer_vars%volatileid(betrtracer_vars%id_trc_ch4))

      this%centuryforc(c,j)%soilorder = biophysforc%isoilorder(c)
      this%centuryforc(c,j)%msurf_nh4 = 0._r8    !this  number needs update
      this%centuryforc(c,j)%msurf_minp = 0._r8   !this  number needs update
    enddo
  enddo

  select type(plant_soilbgc)
  type is(plant_soilbgc_cnp_type)
  do j = lbj, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      this%centuryforc(c,j)%rt_ar  = plant_soilbgc%rt_vr_col(c,j)            !root autotrophic respiration, mol CO2/m3/s
      this%centuryforc(c,j)%plant_froot_nn(:) = plant_soilbgc%plant_froot_nn(:,j)
      this%centuryforc(c,j)%plant_froot_np(:) = plant_soilbgc%plant_froot_np(:,j)
      this%centuryforc(c,j)%plant_vtype(:) = plant_soilbgc%plant_vtype(c,:)
      this%centuryforc(c,j)%plant_ntypes  = plant_soilbgc%plant_ntypes(c)
    enddo
  enddo
  end select
  end associate
  end subroutine set_century_forc
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
  use betr_varcon              , only : spinup_state => bspinup_state
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use PlantSoilBgcCnpType      , only : plant_soilbgc_cnp_type
  use tracer_varcon, only : catomw, natomw, patomw
  implicit none
  class(bgc_reaction_CENTURY_ECACNP_type) , intent(inout)    :: this
  integer                              , intent(in) :: c, j
  integer                              , intent(in) :: nstates
  real(r8)                             , intent(in) :: ystates0(nstates)
  real(r8)                             , intent(in) :: ystatesf(nstates)
  real(r8)                             , intent(in) :: dtime
  type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
  type(tracerstate_type)               , intent(inout) :: tracerstate_vars
  type(tracerflux_type)                , intent(inout) :: tracerflux_vars
  class(plant_soilbgc_type)            , intent(inout) :: plant_soilbgc
  type(betr_biogeo_flux_type)          , intent(inout) :: biogeo_flux
  integer :: k, k1, k2, jj, p
  integer :: trcid

  associate(                                                                &
    nom_pools             => this%centurybgc_index%nom_pools              , & !
    nelms                 => this%centurybgc_index%nelms                  , & !
    litr_beg              => this%centurybgc_index%litr_beg               , & !
    litr_end              => this%centurybgc_index%litr_end               , & !
    wood_beg              => this%centurybgc_index%wood_beg               , & !
    wood_end              => this%centurybgc_index%wood_end               , & !
    som_beg               => this%centurybgc_index%som_beg                , & !
    som_end               => this%centurybgc_index%som_end                , & !
    dom_beg               => this%centurybgc_index%dom_beg                , & !
    dom_end               => this%centurybgc_index%dom_end                , & !
    Bm_beg                => this%centurybgc_index%Bm_beg                 , & !
    Bm_end                => this%centurybgc_index%Bm_end                 , & !
    volatileid            => betrtracer_vars%volatileid                   , &
    tracer_flx_netpro_vr  => tracerflux_vars%tracer_flx_netpro_vr_col     , & !
    tracer_flx_parchm_vr  => tracerflux_vars%tracer_flx_parchm_vr_col     , & !
    ngwmobile_tracers     => betrtracer_vars%ngwmobile_tracers              & !
  )

      !tracer states
      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_litr:betrtracer_vars%id_trc_end_litr) = &
        ystatesf(litr_beg:litr_end)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_wood:betrtracer_vars%id_trc_end_wood) = &
        ystatesf(wood_beg:wood_end)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_Bm:betrtracer_vars%id_trc_end_Bm) = &
        ystatesf(Bm_beg:Bm_end)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_som:betrtracer_vars%id_trc_end_som) = &
        ystatesf(som_beg:som_end)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_beg_dom:betrtracer_vars%id_trc_end_dom) = &
        ystatesf(dom_beg:dom_end)

      k1= betrtracer_vars%id_trc_beg_minp; k2 = this%centurybgc_index%lid_minp_secondary
      tracerstate_vars%tracer_conc_mobile_col(c,j,k1) = ystatesf(k2)

      k1 = betrtracer_vars%id_trc_end_minp;   k2 = this%centurybgc_index%lid_minp_occlude
      tracerstate_vars%tracer_conc_mobile_col(c,j,k1) = ystatesf(k2)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2) = &
        ystatesf(this%centurybgc_index%lid_n2)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2) = &
        ystatesf(this%centurybgc_index%lid_o2)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar) = &
        ystatesf(this%centurybgc_index%lid_ar)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x) = &
        ystatesf(this%centurybgc_index%lid_co2)

      if(this%use_c13)then
        tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c13_co2x) = &
          ystatesf(this%centurybgc_index%lid_c13_co2)
      endif

      if(this%use_c14)then
        tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_c14_co2x) = &
          ystatesf(this%centurybgc_index%lid_c14_co2)
      endif

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4) = &
        ystatesf(this%centurybgc_index%lid_ch4)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_nh3x) = &
        ystatesf(this%centurybgc_index%lid_nh4)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_no3x) = &
        ystatesf(this%centurybgc_index%lid_no3)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2o) = &
        ystatesf(this%centurybgc_index%lid_n2o)

      tracerstate_vars%tracer_conc_mobile_col(c,j,betrtracer_vars%id_trc_p_sol) = &
        ystatesf(this%centurybgc_index%lid_minp_soluble)

      !tracer fluxes
      tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_o2) ) = &
         ystatesf(this%centurybgc_index%lid_o2_paere )  - &
         ystates0(this%centurybgc_index%lid_o2_paere)

      if ( spinup_state /= 1 ) then
        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2)  ) = &
          ystatesf(this%centurybgc_index%lid_n2_paere)  - &
          ystates0(this%centurybgc_index%lid_n2_paere)

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ar)  ) = &
          ystatesf(this%centurybgc_index%lid_ar_paere)  - &
          ystates0(this%centurybgc_index%lid_ar_paere)

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_co2x)) = &
          ystatesf(this%centurybgc_index%lid_co2_paere)  - &
          ystates0(this%centurybgc_index%lid_co2_paere)

        if(this%use_c13)then
          tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c13_co2x)) = &
            ystatesf(this%centurybgc_index%lid_c13_co2_paere)  - &
            ystates0(this%centurybgc_index%lid_c13_co2_paere)
        endif

        if(this%use_c14)then
          tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_c14_co2x)) = &
            ystatesf(this%centurybgc_index%lid_c14_co2_paere)  - &
            ystates0(this%centurybgc_index%lid_c14_co2_paere)
        endif

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ch4) ) = &
          ystatesf(this%centurybgc_index%lid_ch4_paere)  - &
          ystates0(this%centurybgc_index%lid_ch4_paere)

        tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2o) ) = &
          ystatesf(this%centurybgc_index%lid_n2o_paere)  - &
          ystates0(this%centurybgc_index%lid_n2o_paere)
      endif

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_nh3x) =  &
        ystatesf(this%centurybgc_index%lid_nh4) - &
        ystates0(this%centurybgc_index%lid_nh4)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_no3x)  =  &
        ystatesf(this%centurybgc_index%lid_no3) - &
        ystates0(this%centurybgc_index%lid_no3)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2) = &
        ystatesf(this%centurybgc_index%lid_n2) - &
        ystates0(this%centurybgc_index%lid_n2)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_co2x ) = &
        ystatesf(this%centurybgc_index%lid_co2) - &
        ystates0(this%centurybgc_index%lid_co2)

      if(this%use_c13)then
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c13_co2x ) = &
          ystatesf(this%centurybgc_index%lid_c13_co2) - &
          ystates0(this%centurybgc_index%lid_c13_co2)
      endif

      if(this%use_c14)then
        tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_c14_co2x ) = &
          ystatesf(this%centurybgc_index%lid_c14_co2) - &
          ystates0(this%centurybgc_index%lid_c14_co2)
      endif

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2o  ) = &
        ystatesf(this%centurybgc_index%lid_n2o) - &
        ystates0(this%centurybgc_index%lid_n2o)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_o2   ) = &
        ystatesf(this%centurybgc_index%lid_o2) - &
        ystates0(this%centurybgc_index%lid_o2)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ch4  ) = &
        ystatesf(this%centurybgc_index%lid_ch4) - &
        ystates0(this%centurybgc_index%lid_ch4)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ar) = &
        ystatesf(this%centurybgc_index%lid_ar) - &
        ystates0(this%centurybgc_index%lid_ar)

      tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_p_sol) =      &
        ystatesf(this%centurybgc_index%lid_minp_soluble) - &
        ystates0(this%centurybgc_index%lid_minp_soluble)

      !get net production for om pools
      do k = 1, litr_end-litr_beg + 1
        k1 = litr_beg+k-1; k2 = betrtracer_vars%id_trc_beg_litr + k-1
        tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, wood_end-wood_beg + 1
        k1 = wood_beg+k-1; k2 = betrtracer_vars%id_trc_beg_wood + k-1
        tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, som_end-som_beg + 1
        k1 = som_beg+k-1; k2 = betrtracer_vars%id_trc_beg_som+ k-1
        tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, dom_end-dom_beg + 1
        k1 = dom_beg+k-1; k2 = betrtracer_vars%id_trc_beg_dom+ k-1
        tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
      enddo

      do k = 1, Bm_end-Bm_beg + 1
        k1 = Bm_beg+k-1; k2 = betrtracer_vars%id_trc_beg_Bm+ k-1
        tracer_flx_netpro_vr(c,j,k2) =  ystatesf(k1) - ystates0(k1)
      enddo

      trcid = betrtracer_vars%id_trc_beg_minp
      tracer_flx_netpro_vr(c,j, trcid) = &
        ystatesf(this%centurybgc_index%lid_minp_secondary) - &
        ystates0(this%centurybgc_index%lid_minp_secondary)

      trcid = betrtracer_vars%id_trc_end_minp
      tracer_flx_netpro_vr(c,j, trcid) =  &
        ystatesf(this%centurybgc_index%lid_minp_occlude) - &
        ystates0(this%centurybgc_index%lid_minp_occlude)
      !plant soil bgc

      !biogeo_flux
      biogeo_flux%c12flux_vars%hr_vr_col(c,j) = &
        (ystatesf(this%centurybgc_index%lid_co2_hr) - &
        ystates0(this%centurybgc_index%lid_co2_hr))*catomw/dtime

      biogeo_flux%n14flux_vars%f_denit_vr_col(c,j)= &
        (ystatesf(this%centurybgc_index%lid_nh4_nit) - &
         ystates0(this%centurybgc_index%lid_nh4_nit))*natomw/dtime

      biogeo_flux%n14flux_vars%f_nit_vr_col(c,j) = &
        (ystatesf(this%centurybgc_index%lid_n2o_nit) - &
         ystates0(this%centurybgc_index%lid_n2o_nit))*natomw/dtime

      biogeo_flux%n14flux_vars%f_n2o_nit_vr(c,j) = &
        (ystatesf(this%centurybgc_index%lid_no3_den) - &
         ystates0(this%centurybgc_index%lid_no3_den))*natomw/dtime


  select type(plant_soilbgc)
  type is(plant_soilbgc_cnp_type)

    do jj = 1, plant_soilbgc%plant_ntypes(c)
      p = plant_soilbgc%plant_vtype(c,jj)
      plant_soilbgc%plant_minn_no3_active_yield_flx_vr_patch(p,j) = &
          (ystatesf(this%centurybgc_index%lid_plant_minn_no3_pft(jj)) - &
          ystates0(this%centurybgc_index%lid_plant_minn_no3_pft(jj)))*natomw/dtime

      plant_soilbgc%plant_minn_nh4_active_yield_flx_vr_patch(p,j) = &
          (ystatesf(this%centurybgc_index%lid_plant_minn_nh4_pft(jj)) - &
          ystates0(this%centurybgc_index%lid_plant_minn_nh4_pft(jj)))*natomw/dtime

      plant_soilbgc%plant_minp_active_yield_flx_vr_patch(p,j) = &
          (ystatesf(this%centurybgc_index%lid_plant_minp_pft(jj)) - &
           ystates0(this%centurybgc_index%lid_plant_minp_pft(jj)))*patomw/dtime
    enddo
  end select
  end associate
  end subroutine retrieve_output
end module BgcReactionsCentECACnpType
