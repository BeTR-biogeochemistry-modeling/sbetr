module BgcSummsSOMType
!
!DESCRIPTION
!module defines the summs decomposition
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !5 microbial pools, 3 litter pools, 3 woody pools
  !met, cell, lig, cwd, mic, res, mono, poly, enz
  !during decomposition, the cnp ratios of microbial pools are all fixed
  !the litter pools have their cnp ratios varying with nutrient status
  !We consider mono and enz as DOM; poly as som
  !
  integer :: nsummspools
  type, public :: SummsSom_type
    real(r8), pointer :: icn_ratios(:)  => null()
    real(r8), pointer :: icp_ratios(:) => null()
    real(r8), pointer :: icc14_ratios(:)=> null()
    real(r8), pointer :: icc13_ratios(:)=> null()
    real(r8), pointer :: def_cn(:)=> null()
    real(r8), pointer :: def_cp(:)=> null()
    real(r8), pointer :: def_cc13(:)=> null()
    real(r8), pointer :: def_cc14(:)=> null()

    !private parameters
    real(r8) :: gmax_mic !constants
    real(r8) :: yld_mic
    real(r8) :: yld_enz
    real(r8) :: pmax_enz
    real(r8) :: minsite
    real(r8) :: mic_transp
    real(r8) :: decay_mic0
    real(r8) :: decay_mic1
    real(r8) :: decay_enz
    real(r8) :: yld_res 
    real(r8) :: fenz2poly
    real(r8) :: part_mic2poly !partitioning
    real(r8) :: part_mic2enz
    real(r8) :: part_res2mono
    real(r8) :: rate_co2
    real(r8) :: cue
    real(r8) :: rf_l1s1_bgc(2)    !co2 production when metabolic carbon is decomposed
    real(r8) :: rf_l2s1_bgc(2)    !co2 production when cellulose is decomposed
    real(r8) :: rf_l3s2_bgc    !co2 production when lignin is decomposed
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8) :: lwd_fcel
    real(r8) :: lwd_flig
    real(r8) :: fwd_fcel
    real(r8) :: fwd_flig
    real(r8) :: lit_flig
    real(r8) :: k_decay_lit1(2)
    real(r8) :: k_decay_lit2(2)
    real(r8) :: k_decay_lit3(2)
    real(r8) :: k_decay_cwd  !coarse root
    real(r8) :: k_decay_lwd  !large wood
    real(r8) :: k_decay_fwd  !fine branch wood
    real(r8) :: actgB
    real(r8) :: actpE
    real(r8) :: actmr
    real(r8) :: decay_mic

    logical  :: use_c13
    logical  :: use_c14
    logical  :: use_warm
  contains
    procedure, public :: Init
    procedure, public :: run_decomp
    procedure, public :: calc_pot_min_np_flx
    procedure, private :: calc_cascade_matrix
    procedure, private :: calc_som_decay_k
    procedure, private :: calc_som_scale_k
    procedure, private :: calc_som_decay_r
    procedure, private :: calc_cnp_ratios
    procedure, private :: InitAllocate
    procedure, public  :: UpdateParas
    procedure, private :: calc_potential_aerobic_hr
    procedure, private :: apply_spinupf
   !private :: deb_grow
  end type SummsSom_type
contains

  subroutine Init(this, summsbgc_index, biogeo_con, bstatus)

  use BgcSummsIndexType   , only : summsbgc_index_type
  use SummsParaType     , only : SummsPara_type
  use BetrStatusType      , only : betr_status_type
  implicit none
  class(SummsSom_type)          , intent(inout) :: this
  type(summsbgc_index_type)     , intent(in)    :: summsbgc_index
  type(SummsPara_type)        , intent(in)    :: biogeo_con
  type(betr_status_type)        , intent(out)   :: bstatus

  call bstatus%reset()
  nsummspools = summsbgc_index%nom_pools

  call this%InitAllocate()

  end subroutine Init
!------------------------------------------
  subroutine InitAllocate (this)
  !
  !
  implicit none
  class(SummsSom_type), intent(inout) :: this

  allocate(this%icn_ratios(nsummspools));
  allocate(this%icp_ratios(nsummspools));
  allocate(this%icc14_ratios(nsummspools)); this%icc14_ratios(:) = 0._r8
  allocate(this%icc13_ratios(nsummspools)); this%icc13_ratios(:) = 0._r8
  allocate(this%def_cn(nsummspools));
  allocate(this%def_cp(nsummspools));
  allocate(this%def_cc13(nsummspools));this%def_cc13(:) = 0._r8
  allocate(this%def_cc14(nsummspools));this%def_cc14(:) = 0._r8
  end subroutine InitAllocate
!------------------------------------------
  subroutine UpdateParas(this,summsbgc_index,  biogeo_con)
  !
  ! intialize model parameters
  use SummsParaType , only : SummsPara_type
  use BgcSummsIndexType , only : summsbgc_index_type
  use tracer_varcon    , only : catomw, natomw, patomw
  implicit none
  class(SummsSom_type)  , intent(inout) :: this
  type(summsbgc_index_type) , intent(in) :: summsbgc_index
  type(SummsPara_type) , intent(in)    :: biogeo_con

  this%yld_res      =biogeo_con%yld_res
  this%fenz2poly    =biogeo_con%fenz2poly
  this%yld_enz      =biogeo_con%yld_enz
  this%gmax_mic     =biogeo_con%gmax_mic
  this%yld_mic      =biogeo_con%yld_mic
  this%minsite      =biogeo_con%minsite
  this%mic_transp   =biogeo_con%mic_transp
  this%decay_mic0   =biogeo_con%decay_mic0
  this%decay_mic1   =biogeo_con%decay_mic1
  this%decay_enz    =biogeo_con%decay_enz
  this%pmax_enz     =biogeo_con%pmax_enz

  this%rf_l1s1_bgc    = biogeo_con%rf_l1s1_bgc
  this%rf_l2s1_bgc    = biogeo_con%rf_l2s1_bgc
  this%rf_l3s2_bgc    = biogeo_con%rf_l3s2_bgc
  this%cwd_fcel   = biogeo_con%cwd_fcel_bgc
  this%cwd_flig   = biogeo_con%cwd_flig_bgc
  this%lwd_fcel   = biogeo_con%lwd_fcel_bgc
  this%lwd_flig   = biogeo_con%lwd_flig_bgc
  this%fwd_fcel   = biogeo_con%fwd_fcel_bgc
  this%fwd_flig   = biogeo_con%fwd_flig_bgc

  this%k_decay_lit1   =  biogeo_con%k_decay_lit1
  this%k_decay_lit2   =  biogeo_con%k_decay_lit2
  this%k_decay_lit3   =  biogeo_con%k_decay_lit3
  this%k_decay_cwd    =  biogeo_con%k_decay_cwd
  this%k_decay_lwd    =  biogeo_con%k_decay_lwd
  this%k_decay_fwd    =  biogeo_con%k_decay_fwd

  this%def_cn(summsbgc_index%lit1) = biogeo_con%init_cn_met * natomw/catomw
  this%def_cn(summsbgc_index%lit2) = biogeo_con%init_cn_cel * natomw/catomw
  this%def_cn(summsbgc_index%lit3) = biogeo_con%init_cn_lig * natomw/catomw
  this%def_cn(summsbgc_index%cwd)  = biogeo_con%init_cn_cwd * natomw/catomw
  this%def_cn(summsbgc_index%lwd)  = biogeo_con%init_cn_lwd * natomw/catomw
  this%def_cn(summsbgc_index%fwd)  = biogeo_con%init_cn_fwd * natomw/catomw
  
  this%def_cn(summsbgc_index%poly) = biogeo_con%init_cn_poly * natomw/catomw
  this%def_cn(summsbgc_index%mono) = biogeo_con%init_cn_mono * natomw/catomw
  this%def_cn(summsbgc_index%mic)  = biogeo_con%init_cn_mic * natomw/catomw
  this%def_cn(summsbgc_index%enz)  = biogeo_con%init_cn_enz * natomw/catomw
  this%def_cn(summsbgc_index%res)  = biogeo_con%init_cn_res * natomw/catomw

  this%def_cp(summsbgc_index%lit1) = biogeo_con%init_cp_met * patomw/catomw
  this%def_cp(summsbgc_index%lit2) = biogeo_con%init_cp_cel * patomw/catomw
  this%def_cp(summsbgc_index%lit3) = biogeo_con%init_cp_lig * patomw/catomw
  this%def_cp(summsbgc_index%cwd)  = biogeo_con%init_cp_cwd * patomw/catomw
  this%def_cp(summsbgc_index%lwd)  = biogeo_con%init_cp_lwd * patomw/catomw
  this%def_cp(summsbgc_index%fwd)  = biogeo_con%init_cp_fwd * patomw/catomw
  
  this%def_cp(summsbgc_index%poly) = biogeo_con%init_cp_poly * patomw/catomw
  this%def_cp(summsbgc_index%mono) = biogeo_con%init_cp_mono * patomw/catomw
  this%def_cp(summsbgc_index%mic)  = biogeo_con%init_cp_mic * patomw/catomw
  this%def_cp(summsbgc_index%enz)  = biogeo_con%init_cp_enz * patomw/catomw
  this%def_cp(summsbgc_index%res)  = biogeo_con%init_cp_res * patomw/catomw

  this%use_c13=biogeo_con%use_c13
  this%use_c14=biogeo_con%use_c14
  this%use_warm=biogeo_con%use_warm

  if(this%use_c13)then
    this%def_cc13(summsbgc_index%lit1) = biogeo_con%init_cc13_met
    this%def_cc13(summsbgc_index%lit2) = biogeo_con%init_cc13_cel
    this%def_cc13(summsbgc_index%lit3) = biogeo_con%init_cc13_lig
    this%def_cc13(summsbgc_index%cwd)  = biogeo_con%init_cc13_cwd
    this%def_cc13(summsbgc_index%lwd)  = biogeo_con%init_cc13_lwd
    this%def_cc13(summsbgc_index%fwd)  = biogeo_con%init_cc13_fwd
    this%def_cc13(summsbgc_index%poly) = biogeo_con%init_cc13_poly
    this%def_cc13(summsbgc_index%mono) = biogeo_con%init_cc13_mono
    this%def_cc13(summsbgc_index%mic)  = biogeo_con%init_cc13_mic
    this%def_cc13(summsbgc_index%enz)  = biogeo_con%init_cc13_enz
    this%def_cc13(summsbgc_index%res)  = biogeo_con%init_cc13_res
  endif

  if(this%use_c14)then
    this%def_cc14(summsbgc_index%lit1) = biogeo_con%init_cc14_met
    this%def_cc14(summsbgc_index%lit2) = biogeo_con%init_cc14_cel
    this%def_cc14(summsbgc_index%lit3) = biogeo_con%init_cc14_lig
    this%def_cc14(summsbgc_index%cwd)  = biogeo_con%init_cc14_cwd
    this%def_cc14(summsbgc_index%lwd)  = biogeo_con%init_cc14_lwd
    this%def_cc14(summsbgc_index%fwd)  = biogeo_con%init_cc14_fwd
    this%def_cc14(summsbgc_index%poly) = biogeo_con%init_cc14_poly
    this%def_cc14(summsbgc_index%mono) = biogeo_con%init_cc14_mono
    this%def_cc14(summsbgc_index%mic) = biogeo_con%init_cc14_mic
    this%def_cc14(summsbgc_index%enz) = biogeo_con%init_cc14_enz
    this%def_cc14(summsbgc_index%res) = biogeo_con%init_cc14_res

  endif

  end subroutine UpdateParas
!------------------------------------------

  subroutine run_decomp(this, is_surflit, summsbgc_index, dtime, ystates,&
      decompkf_eca, alpha_n, alpha_p, cascade_matrix, &
      k_decay, pot_co2_hr, bstatus)
  !
  !DESCRIPTION
  !
  use BgcSummsIndexType     , only : summsbgc_index_type
  use BgcSummsDecompType    , only : DecompSumms_type
  use BetrStatusType        , only : betr_status_type
  use betr_ctrl           , only : betr_spinup_state
  implicit none
  class(SummsSom_type)        , intent(inout) :: this
  type(summsbgc_index_type)   , intent(in) :: summsbgc_index
  real(r8)                    , intent(in) :: dtime
  real(r8)                    , intent(inout) :: ystates(1:summsbgc_index%nom_tot_elms)
  type(DecompSumms_type)      , intent(in) :: decompkf_eca
  logical                     , intent(in) :: is_surflit
  !real(r8)                    , intent(in) :: pct_sand
  !real(r8)                    , intent(in) :: pct_clay
  real(r8)                    , intent(inout) :: cascade_matrix(summsbgc_index%nstvars, summsbgc_index%nreactions)
  real(r8)                    , intent(out) :: k_decay(1:nsummspools)
  real(r8)                    , intent(out) :: pot_co2_hr
  real(r8)                    , intent(out) :: alpha_n(1:nsummspools)
  real(r8)                    , intent(out) :: alpha_p(1:nsummspools)
  type(betr_status_type)      , intent(out) :: bstatus

  !local variables
  real(r8) :: pot_om_decay_rates(1:nsummspools)
  integer :: kc, jj, lay

  associate(                                      &
    nelms => summsbgc_index%nelms,              &
    nom_tot_elms=> summsbgc_index%nom_tot_elms, &
    c_loc => summsbgc_index%c_loc               &
  )
  call bstatus%reset()

  if(is_surflit)then
    lay=1
  else
    lay=2
  endif
  call this%calc_cnp_ratios(summsbgc_index, ystates, bstatus)
  if (bstatus%check_status())return
  !calculate potential decay coefficients (1/s)
  !call this%calc_som_decay_k(lay, summsbgc_index, decompkf_eca, k_decay, ystates, bstatus)
  call this%calc_som_decay_k(lay, summsbgc_index, decompkf_eca, k_decay(1:nsummspools), ystates, bstatus)

  !scale potential decay coefficients by temp (1/s)
  call this%calc_som_scale_k(lay, summsbgc_index, decompkf_eca, k_decay(1:nsummspools))

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(summsbgc_index, dtime, k_decay(1:nsummspools), &
      ystates(1:nom_tot_elms), pot_om_decay_rates)

  !calculate custom fluxes
  ystates(summsbgc_index%lid_decomp) = pot_om_decay_rates(summsbgc_index%poly)
  ystates(summsbgc_index%lid_uptake) = pot_om_decay_rates(summsbgc_index%mono)*(1._r8-this%yld_res)
  ystates(summsbgc_index%lid_cue)    = this%cue
  ystates(summsbgc_index%lid_maint)  = ystates((summsbgc_index%mic-1) * nelms + c_loc)*this%actmr
  ystates(summsbgc_index%lid_kaffmm) = decompkf_eca%kaff_mono_msurf
  ystates(summsbgc_index%lid_kaffem) = decompkf_eca%kaff_enz_msurf
  ystates(summsbgc_index%lid_micgrow)= ystates((summsbgc_index%mic-1) * nelms + c_loc)*this%actgB
  ystates(summsbgc_index%lid_enzprod)= ystates((summsbgc_index%mic-1) * nelms + c_loc)*this%actpE
  ystates(summsbgc_index%lid_turnover)= ystates((summsbgc_index%mic-1) * nelms + c_loc)*this%decay_mic

  do jj = 1, nsummspools
    kc = (jj-1) * nelms + c_loc
    !the following avoids over-estimation of potential hr which is used for nitri-denit estimation
    pot_om_decay_rates(jj) = min(pot_om_decay_rates(jj), ystates(kc)/dtime)
  enddo

  call this%calc_cascade_matrix(lay, summsbgc_index, alpha_n, alpha_p, cascade_matrix)

  !calculate potential respiration rates by summarizing all om decomposition pathways
  call this%calc_potential_aerobic_hr(summsbgc_index, pot_om_decay_rates, &
    cascade_matrix, pot_co2_hr, bstatus)

  end associate
  end subroutine run_decomp
!------------------------------------------

  subroutine calc_cascade_matrix(this, lay, summsbgc_index, alpha_n, alpha_p, cascade_matrix)
  !subroutine calc_cascade_matrix(this, is_surf, summsbgc_index, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix)

  !
  ! DESCRIPTION
  ! calculate cascade matrix for decomposition
  ! in all the reactions, the nominal carbon oxidation status is assumed as zero, which is apparently not correct.
  ! It is also assumed the recycling of nitrogen and phosphorus during decomposition is 100%, which is likely
  ! not quite right as well.
  use BgcSummsIndexType , only : summsbgc_index_type
  use MathfuncMod         , only : safe_div, fpmax
  implicit none
  class(SummsSom_type),          intent(inout) :: this
  type(summsbgc_index_type)    , intent(in)    :: summsbgc_index
  integer                      , intent(in)    :: lay
  !real(r8)                     , intent(in)    :: pct_sand
  !real(r8)                     , intent(in)    :: pct_clay
  real(r8)                     , intent(out)   :: alpha_n(nsummspools) !indicating factor for nitrogen limitation
  real(r8)                     , intent(out)   :: alpha_p(nsummspools) !indicating factor for phosphorus limitation
  real(r8)                     , intent(inout) :: cascade_matrix(summsbgc_index%nstvars, summsbgc_index%nreactions)

  integer  :: reac,jj
  real(r8) :: f1, f2, rf_s1

    !logical  :: use_c13
    !logical  :: use_c14

  associate(                                                   &
    lit1      => summsbgc_index%lit1                       , & !
    lit2      => summsbgc_index%lit2                       , & !
    lit3      => summsbgc_index%lit3                       , & !
    poly      => summsbgc_index%poly                       , & !
    mono      => summsbgc_index%mono                       , & !
    res      => summsbgc_index%res                         , & !
    mic      => summsbgc_index%mic                         , & !
    enz      => summsbgc_index%enz                         , & !
    cwd       => summsbgc_index%cwd                        , & !
    lwd       => summsbgc_index%lwd                        , & !
    fwd       => summsbgc_index%fwd                        , & !
    c_loc     => summsbgc_index%c_loc                      , & !
    n_loc     => summsbgc_index%n_loc                      , & !
    p_loc     => summsbgc_index%p_loc                      , & !
    c13_loc   => summsbgc_index%c13_loc                    , & !
    c14_loc   => summsbgc_index%c14_loc                    , & !
    nelms     => summsbgc_index%nelms                      , & !
    lid_o2    => summsbgc_index%lid_o2                     , & !
    lid_co2   => summsbgc_index%lid_co2                    , & !
    lid_nh4   => summsbgc_index%lid_nh4                    , & !
    lid_c14_co2=> summsbgc_index%lid_c14_co2               , & !
    lid_c13_co2=> summsbgc_index%lid_c13_co2               , & !
    lid_co2_hr => summsbgc_index%lid_co2_hr                , &
    lid_minn_nh4_immob=> summsbgc_index%lid_minn_nh4_immob , &
    lid_minp_immob => summsbgc_index%lid_minp_immob        , &
    lid_minp_soluble=> summsbgc_index%lid_minp_soluble     , &
    lit1_dek_reac => summsbgc_index%lit1_dek_reac          , &
    lit2_dek_reac => summsbgc_index%lit2_dek_reac          , &
    lit3_dek_reac => summsbgc_index%lit3_dek_reac          , &
    poly_dek_reac => summsbgc_index%poly_dek_reac          , &
    mono_dek_reac => summsbgc_index%mono_dek_reac          , &
    mic_dek_reac => summsbgc_index%mic_dek_reac            , &
    enz_dek_reac => summsbgc_index%enz_dek_reac            , &
    res_dek_reac => summsbgc_index%res_dek_reac            , &
    cwd_dek_reac => summsbgc_index%cwd_dek_reac            , &
    lwd_dek_reac => summsbgc_index%lwd_dek_reac            , &
    fwd_dek_reac => summsbgc_index%fwd_dek_reac            , &
    cwd_fcel     => this%cwd_fcel                            , &
    cwd_flig     => this%cwd_flig                            , &
    lwd_fcel     => this%lwd_fcel                            , &
    lwd_flig     => this%lwd_flig                            , &
    fwd_fcel     => this%fwd_fcel                            , &
    fwd_flig     => this%fwd_flig                            , &
    rf_l2s1_bgc  => this%rf_l2s1_bgc                         , &
    rf_l3s2_bgc  => this%rf_l3s2_bgc                         , &
    rf_l1s1_bgc  => this%rf_l1s1_bgc                         , &
    yld_res      => this%yld_res                             , &
    fenz2poly    => this%fenz2poly                           , &
    part_mic2poly => this%part_mic2poly                      , &
    part_mic2enz  => this%part_mic2enz                       , &
    part_res2mono => this%part_res2mono                      , &
    rate_co2     => this%rate_co2                            , &
    cue          => this%cue                                 , &
    debug        => summsbgc_index%debug                       &
  )

    alpha_n = 0._r8; alpha_p = 0._r8
    !---------------------------------------------------------------------------------
    !reaction1, lit1 -> mono
    reac=lit1_dek_reac
    !lit1 + 0.55*o2 -> 0.45 mono + 0.55co2 + (1/cn_ratios(lit1) - 0.45/cn_ratios(mono))min_n+ (1/cp_ratios(lit1)-0.45/cp_ratios(mono))min_p
    cascade_matrix((lit1-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((lit1-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(lit1)
    cascade_matrix((lit1-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(lit1)

    cascade_matrix((mono-1)*nelms+c_loc   ,reac)  = 1._r8-rf_l1s1_bgc(lay)
    cascade_matrix((mono-1)*nelms+n_loc   ,reac)  = cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icn_ratios(mono)
    cascade_matrix((mono-1)*nelms+p_loc   ,reac)  = cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icp_ratios(mono)

    cascade_matrix(lid_co2                ,reac)  = rf_l1s1_bgc(lay)

    cascade_matrix(lid_o2                 ,reac)  = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac)  = -cascade_matrix((lit1-1)*nelms+n_loc   ,reac) - &
        cascade_matrix((mono-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble         ,reac)  = -cascade_matrix((lit1-1)*nelms+p_loc   ,reac) - &
        cascade_matrix((mono-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)  = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2           ,reac)
    cascade_matrix(lid_minp_immob         ,reac)  = -cascade_matrix(lid_minp_soluble  ,reac)

    if(this%use_c14)then
      cascade_matrix((lit1-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(lit1)
      cascade_matrix(lid_c14_co2              , reac) = rf_l1s1_bgc(lay)*this%icc14_ratios(lit1)
      cascade_matrix((mono-1)*nelms+c14_loc   , reac) = cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icc14_ratios(lit1)
    endif

    if(this%use_c13)then
      cascade_matrix((lit1-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(lit1)
      cascade_matrix(lid_c13_co2              , reac) = rf_l1s1_bgc(lay)*this%icc13_ratios(lit1)
      cascade_matrix((mono-1)*nelms+c13_loc   , reac) = cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icc13_ratios(lit1)
    endif

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8
    if(debug)then
      write(*,*)'lit1 carbon',cascade_matrix((lit1-1)*nelms+c_loc   ,reac)+cascade_matrix((mono-1)*nelms+c_loc   ,reac)+&
        cascade_matrix(lid_co2             ,reac)
      write(*,*)'lit1 nitrogen',cascade_matrix((lit1-1)*nelms+n_loc   ,reac)+cascade_matrix((mono-1)*nelms+n_loc   ,reac)+&
       cascade_matrix(lid_nh4         ,reac)
      write(*,*)'lit1 phosp',cascade_matrix((lit1-1)*nelms+p_loc   ,reac)+cascade_matrix((mono-1)*nelms+p_loc   ,reac)+&
       cascade_matrix(lid_minp_soluble        ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 2, lit2 -> mono
    reac = lit2_dek_reac
    !lit2 + 0.5 o2  -> 0.5 mono + 0.5 co2 + (1/cn_ratios(lit2)-0.5/cn_ratios(mono))min_n +(1/cp_ratios(lit2)-0.5/cp_ratios(mono))min_p
    cascade_matrix((lit2-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((lit2-1)*nelms+n_loc   ,reac)   = -this%icn_ratios(lit2)
    cascade_matrix((lit2-1)*nelms+p_loc   ,reac)   = -this%icp_ratios(lit2)

    cascade_matrix((mono-1)*nelms+c_loc   ,reac)   =  1._r8-rf_l2s1_bgc(lay)
    cascade_matrix((mono-1)*nelms+n_loc   ,reac)   =  cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icn_ratios(mono)
    cascade_matrix((mono-1)*nelms+p_loc   ,reac)   =  cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icp_ratios(mono)

    cascade_matrix(lid_co2                ,reac)   =  rf_l2s1_bgc(lay)
    cascade_matrix(lid_o2                 ,reac)   = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac)   = -cascade_matrix((lit2-1)*nelms+n_loc   ,reac) - &
                                                      cascade_matrix((mono-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac)   = -cascade_matrix((lit2-1)*nelms+p_loc   ,reac) - &
                                                       cascade_matrix((mono-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)   = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac)   = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac)   = cascade_matrix(lid_co2        ,reac)


    if(cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if(cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((lit2-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(lit2)
      cascade_matrix(lid_c14_co2              , reac) = rf_l2s1_bgc(lay)*this%icc14_ratios(lit2)
      cascade_matrix((mono-1)*nelms+c14_loc   , reac) = cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icc14_ratios(lit2)
    endif
    if(this%use_c13)then
      cascade_matrix((lit2-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(lit2)
      cascade_matrix(lid_c13_co2              , reac) =  rf_l2s1_bgc(lay)*this%icc13_ratios(lit2)
      cascade_matrix((mono-1)*nelms+c13_loc   , reac) =  cascade_matrix((mono-1)*nelms+c_loc,reac)*this%icc13_ratios(lit2)
    endif
    if(debug)then
      write(*,*)'lit2 carbon',cascade_matrix((lit2-1)*nelms+c_loc   ,reac) + cascade_matrix((mono-1)*nelms+c_loc   ,reac)+&
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'lit2 nitrogen',cascade_matrix((lit2-1)*nelms+n_loc   ,reac) + cascade_matrix((mono-1)*nelms+n_loc   ,reac)+&
        cascade_matrix(lid_nh4                ,reac)
      write(*,*)'lit2 phosp',cascade_matrix((lit2-1)*nelms+p_loc   ,reac) + cascade_matrix((mono-1)*nelms+p_loc   ,reac)+&
        cascade_matrix(lid_minp_soluble                ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 3, lit3->poly
    reac = lit3_dek_reac
    !lit3 + 0.5 o2 -> 0.5 poly + 0.5 co2 + (1/cn_ratios(lit3) - 0.5/cn_ratios(poly))min_n + (1/cp_ratios(lit3)-0.5_r8/cp_ratios(poly))minp
    cascade_matrix((lit3-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((lit3-1)*nelms+n_loc   ,reac) = -this%icn_ratios(lit3)
    cascade_matrix((lit3-1)*nelms+p_loc   ,reac) = -this%icp_ratios(lit3)

    cascade_matrix((poly-1)*nelms+c_loc   ,reac) =  1._r8-rf_l3s2_bgc
    cascade_matrix((poly-1)*nelms+n_loc   ,reac) =  cascade_matrix((poly-1)*nelms+c_loc,reac)*this%icn_ratios(poly)
    cascade_matrix((poly-1)*nelms+p_loc   ,reac) =  cascade_matrix((poly-1)*nelms+c_loc,reac)*this%icp_ratios(poly)

    cascade_matrix(lid_co2                ,reac) = rf_l3s2_bgc
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((lit3-1)*nelms+n_loc   ,reac) - &
                                                    cascade_matrix((poly-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble      ,reac) = -cascade_matrix((lit3-1)*nelms+p_loc   ,reac) - &
                                                   cascade_matrix((poly-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)


    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((lit3-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(lit3)
      cascade_matrix(lid_c14_co2              , reac) =  rf_l3s2_bgc*this%icc14_ratios(lit3)
      cascade_matrix((poly-1)*nelms+c14_loc   , reac) =  cascade_matrix((poly-1)*nelms+c_loc,reac)*this%icc14_ratios(lit3)
    endif

    if(this%use_c13)then
      cascade_matrix((lit3-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(lit3)
      cascade_matrix(lid_c13_co2              , reac) =  rf_l3s2_bgc*this%icc13_ratios(lit3)
      cascade_matrix((poly-1)*nelms+c13_loc   , reac) =  cascade_matrix((poly-1)*nelms+c_loc,reac)*this%icc13_ratios(lit3)
    endif
    if(debug)then
      write(*,*)'lit3 carbon',cascade_matrix((lit3-1)*nelms+c_loc   ,reac) + cascade_matrix((poly-1)*nelms+c_loc   ,reac) + &
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'lit3 nitrogen',cascade_matrix((lit3-1)*nelms+n_loc   ,reac) + cascade_matrix((poly-1)*nelms+n_loc   ,reac) + &
         cascade_matrix(lid_nh4                ,reac)
      write(*,*)'lit3 phosp',cascade_matrix((lit3-1)*nelms+p_loc   ,reac) + cascade_matrix((poly-1)*nelms+p_loc   ,reac) + &
         cascade_matrix(lid_minp_soluble                ,reac)
    endif
    !---------------------------------------------------------------------------------
    !double check those stoichiometry parameters
    !reaction 4, poly -> mono
    !depolymerization by ECA

    reac = poly_dek_reac

    cascade_matrix((poly-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((poly-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(poly)
    cascade_matrix((poly-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(poly)

    cascade_matrix((mono-1)*nelms+c_loc   ,reac) = 1._r8
    cascade_matrix((mono-1)*nelms+n_loc   ,reac) = 1._r8*this%icn_ratios(mono)
    cascade_matrix((mono-1)*nelms+p_loc   ,reac) = 1._r8*this%icp_ratios(mono)

    !no co2 evolution or o2 consumption, since there is no mineralization of C, some mineralization of N and P based on stiochiometry differences
    cascade_matrix(lid_co2                ,reac) = 0._r8
    cascade_matrix(lid_o2                 ,reac) = 0._r8
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((poly-1)*nelms+n_loc   ,reac) - &
        cascade_matrix((mono-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((poly-1)*nelms+p_loc   ,reac) - &
        cascade_matrix((mono-1)*nelms+p_loc   ,reac)
    cascade_matrix(lid_minn_nh4_immob     ,reac) = 0._r8
    cascade_matrix(lid_minp_immob         ,reac) = 0._r8
    cascade_matrix(lid_co2_hr             ,reac) = 0._r8

    if(this%use_c14)then
      cascade_matrix((poly-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(poly)
      cascade_matrix(lid_c14_co2              , reac) =  0._r8
      cascade_matrix((mono-1)*nelms+c14_loc   , reac) =  1._r8*this%icc14_ratios(poly)
    endif

    if(this%use_c13)then
      cascade_matrix((poly-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(poly)
      cascade_matrix(lid_c13_co2              , reac) = 0._r8
      cascade_matrix((mono-1)*nelms+c13_loc   , reac) = 1._r8*this%icc13_ratios(poly)
    endif

    if(debug)then
      write(*,*)'poly carbon',cascade_matrix((poly-1)*nelms+c_loc   ,reac) +cascade_matrix((mono-1)*nelms+c_loc   ,reac)+&
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'poly nitrogen',cascade_matrix((poly-1)*nelms+n_loc   ,reac) +cascade_matrix((mono-1)*nelms+n_loc   ,reac)+&
         cascade_matrix(lid_nh4                ,reac)
      write(*,*)'poly phosp',cascade_matrix((poly-1)*nelms+p_loc   ,reac) +cascade_matrix((mono-1)*nelms+p_loc   ,reac)+&
         cascade_matrix(lid_minp_soluble       ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 5, mono -> res + co2
    !uptake and assimilation cost
    reac = mono_dek_reac

    cascade_matrix((mono-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((mono-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(mono)
    cascade_matrix((mono-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(mono)

    cascade_matrix((res-1)*nelms+c_loc   ,reac) = yld_res
    cascade_matrix((res-1)*nelms+n_loc   ,reac) = yld_res*this%icn_ratios(res)
    cascade_matrix((res-1)*nelms+p_loc   ,reac) = yld_res*this%icp_ratios(res)

    cascade_matrix(lid_co2                ,reac) = 1._r8 - yld_res
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((mono-1)*nelms+n_loc   ,reac) - &
                                                    cascade_matrix((res-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble      ,reac) = -cascade_matrix((mono-1)*nelms+p_loc   ,reac) - &
                                                   cascade_matrix((res-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8 

    if(this%use_c14)then
      cascade_matrix((mono-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(mono)
      cascade_matrix(lid_c14_co2              , reac) =  (1._r8 - yld_res)*this%icc14_ratios(mono)
      cascade_matrix((res-1)*nelms+c14_loc    , reac) =  yld_res*this%icc14_ratios(mono)
    endif

    if(this%use_c13)then
      cascade_matrix((mono-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(mono)
      cascade_matrix(lid_c13_co2              , reac) = (1._r8 - yld_res)*this%icc13_ratios(mono)
      cascade_matrix((res-1)*nelms+c13_loc    , reac) = yld_res*this%icc13_ratios(mono)
    endif

    if(debug)then
      write(*,*)'mono carbon',cascade_matrix((mono-1)*nelms+c_loc   ,reac) +cascade_matrix((res-1)*nelms+c_loc   ,reac)+&
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'mono nitrogen',cascade_matrix((mono-1)*nelms+n_loc   ,reac) +cascade_matrix((res-1)*nelms+n_loc   ,reac)+&
         cascade_matrix(lid_nh4                ,reac)
      write(*,*)'mono phosp',cascade_matrix((mono-1)*nelms+p_loc   ,reac) +cascade_matrix((res-1)*nelms+p_loc   ,reac)+&
         cascade_matrix(lid_minp_soluble       ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 6, mic -> poly + enz + co2
    !microbial production, death, and structural turnover
    reac = mic_dek_reac

    cascade_matrix((mic-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((mic-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(mic)
    cascade_matrix((mic-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(mic)

    cascade_matrix((poly-1)*nelms+c_loc   ,reac) = part_mic2poly
    cascade_matrix((poly-1)*nelms+n_loc   ,reac) = part_mic2poly*this%icn_ratios(poly)
    cascade_matrix((poly-1)*nelms+p_loc   ,reac) = part_mic2poly*this%icp_ratios(poly)

    cascade_matrix((enz-1)*nelms+c_loc   ,reac) = part_mic2enz
    cascade_matrix((enz-1)*nelms+n_loc   ,reac) = part_mic2enz*this%icn_ratios(enz)
    cascade_matrix((enz-1)*nelms+p_loc   ,reac) = part_mic2enz*this%icp_ratios(enz)

    cascade_matrix(lid_co2                ,reac) = 1._r8-part_mic2poly-part_mic2enz
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((mic-1)*nelms+n_loc   ,reac) - &
                                                    cascade_matrix((poly-1)*nelms+n_loc  ,reac) - &
                                                    cascade_matrix((enz-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((mic-1)*nelms+p_loc   ,reac) - &
                                                    cascade_matrix((poly-1)*nelms+p_loc  ,reac) - &
                                                    cascade_matrix((enz-1)*nelms+p_loc   ,reac)
    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)  

    if(this%use_c14)then
      cascade_matrix((mic-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(mic)
      cascade_matrix(lid_c14_co2              , reac) = (1._r8-part_mic2poly-part_mic2enz)*this%icc14_ratios(mic)
      cascade_matrix((poly-1)*nelms+c14_loc   , reac) = 1._r8*this%icc14_ratios(mic)
      cascade_matrix((enz-1)*nelms+c14_loc   , reac) = 1._r8*this%icc14_ratios(mic)
    endif

    if(this%use_c13)then
      cascade_matrix((mic-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(mic)
      cascade_matrix(lid_c13_co2              , reac) = (1._r8-part_mic2poly-part_mic2enz)*this%icc13_ratios(mic)
      cascade_matrix((poly-1)*nelms+c13_loc   , reac) = 1._r8*this%icc13_ratios(mic)
      cascade_matrix((enz-1)*nelms+c13_loc   , reac) = 1._r8*this%icc13_ratios(mic)
    endif

    if(debug)then
      write(*,*)'mic carbon',cascade_matrix((mic-1)*nelms+c_loc   ,reac) +cascade_matrix((poly-1)*nelms+c_loc   ,reac)+&
      cascade_matrix((enz-1)*nelms+c_loc   ,reac)+&
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'mic nitrogen',cascade_matrix((mic-1)*nelms+n_loc   ,reac) +cascade_matrix((poly-1)*nelms+n_loc   ,reac)+&
      cascade_matrix((enz-1)*nelms+n_loc   ,reac)+&
         cascade_matrix(lid_nh4                ,reac)
      write(*,*)'mic phosp',cascade_matrix((mic-1)*nelms+p_loc   ,reac) +cascade_matrix((poly-1)*nelms+p_loc   ,reac)+&
      cascade_matrix((enz-1)*nelms+p_loc   ,reac)+&
         cascade_matrix(lid_minp_soluble       ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 7, enz -> poly + mono
    !enzyme turnover
    reac = enz_dek_reac

    cascade_matrix((enz-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((enz-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(enz)
    cascade_matrix((enz-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(enz)

    cascade_matrix((poly-1)*nelms+c_loc   ,reac) = fenz2poly
    cascade_matrix((poly-1)*nelms+n_loc   ,reac) = fenz2poly*this%icn_ratios(poly)
    cascade_matrix((poly-1)*nelms+p_loc   ,reac) = fenz2poly*this%icp_ratios(poly)

    cascade_matrix((mono-1)*nelms+c_loc   ,reac)  = 1._r8-fenz2poly
    cascade_matrix((mono-1)*nelms+n_loc   ,reac)  = (1._r8-fenz2poly)*this%icn_ratios(mono)
    cascade_matrix((mono-1)*nelms+p_loc   ,reac)  = (1._r8-fenz2poly)*this%icp_ratios(mono)

    cascade_matrix(lid_co2                ,reac) = 0._r8
    cascade_matrix(lid_o2                 ,reac) = 0._r8
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((enz-1)*nelms+n_loc   ,reac) - &
        cascade_matrix((poly-1)*nelms+n_loc   ,reac) - cascade_matrix((mono-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((enz-1)*nelms+p_loc   ,reac) - &
        cascade_matrix((poly-1)*nelms+p_loc   ,reac) - cascade_matrix((mono-1)*nelms+p_loc   ,reac)
    cascade_matrix(lid_minn_nh4_immob     ,reac) = 0._r8
    cascade_matrix(lid_minp_immob         ,reac) = 0._r8
    cascade_matrix(lid_co2_hr             ,reac) = 0._r8

    if(this%use_c14)then
      cascade_matrix((enz-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(enz)
      cascade_matrix(lid_c14_co2              , reac) =  0._r8
      cascade_matrix((poly-1)*nelms+c14_loc   , reac) =  fenz2poly*this%icc14_ratios(enz)
      cascade_matrix((mono-1)*nelms+c14_loc   , reac) =  (1._r8-fenz2poly)*this%icc14_ratios(enz)
    endif

    if(this%use_c13)then
      cascade_matrix((enz-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(enz)
      cascade_matrix(lid_c13_co2              , reac) = 0._r8
      cascade_matrix((poly-1)*nelms+c13_loc   , reac) = fenz2poly*this%icc13_ratios(enz)
      cascade_matrix((mono-1)*nelms+c13_loc   , reac) = (1._r8-fenz2poly)*this%icc13_ratios(enz)
    endif

    if(debug)then
      write(*,*)'enz carbon',cascade_matrix((enz-1)*nelms+c_loc   ,reac) +cascade_matrix((poly-1)*nelms+c_loc   ,reac)+&
         cascade_matrix((mono-1)*nelms+c_loc   ,reac) +cascade_matrix(lid_co2                ,reac)
      write(*,*)'enz nitrogen',cascade_matrix((enz-1)*nelms+n_loc   ,reac) +cascade_matrix((poly-1)*nelms+n_loc   ,reac)+&
         cascade_matrix((mono-1)*nelms+n_loc   ,reac) +cascade_matrix(lid_nh4                ,reac)
      write(*,*)'enz phosp',cascade_matrix((enz-1)*nelms+p_loc   ,reac) +cascade_matrix((poly-1)*nelms+p_loc   ,reac)+&
         cascade_matrix((mono-1)*nelms+p_loc   ,reac) +cascade_matrix(lid_minp_soluble       ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 8, res -> mono + co2
    !reserve allocation and microbial growth
    reac = res_dek_reac

    cascade_matrix((res-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((res-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(res)
    cascade_matrix((res-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(res)

    cascade_matrix((mono-1)*nelms+c_loc   ,reac) = part_res2mono
    cascade_matrix((mono-1)*nelms+n_loc   ,reac) = part_res2mono*this%icn_ratios(mono)
    cascade_matrix((mono-1)*nelms+p_loc   ,reac) = part_res2mono*this%icp_ratios(mono)

    cascade_matrix(lid_co2                ,reac) = 1._r8 - part_res2mono

    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2   ,reac)

    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((res-1)*nelms+n_loc   ,reac) - &
                                                    cascade_matrix((mono-1)*nelms+n_loc  ,reac)

    cascade_matrix(lid_minp_soluble      ,reac) = -cascade_matrix((res-1)*nelms+p_loc   ,reac) - &
                                                   cascade_matrix((mono-1)*nelms+p_loc  ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8 

    if(this%use_c14)then
      cascade_matrix((res-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(res)
      cascade_matrix(lid_c14_co2             , reac) = (1._r8 - part_res2mono)*this%icc14_ratios(res)
      cascade_matrix((mono-1)*nelms+c14_loc  , reac) =  part_res2mono*this%icc14_ratios(res)
    endif

    if(this%use_c13)then
      cascade_matrix((res-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(res)
      cascade_matrix(lid_c13_co2             , reac) = (1._r8 - part_res2mono)*this%icc13_ratios(res)
      cascade_matrix((mono-1)*nelms+c13_loc   , reac) = part_res2mono*this%icc13_ratios(res)
    endif

    if(debug)then
      write(*,*)'res carbon',cascade_matrix((res-1)*nelms+c_loc ,reac) +&
         cascade_matrix((mono-1)*nelms+c_loc   ,reac)+&
         cascade_matrix(lid_co2 ,  reac)
      write(*,*)'res nitrogen',cascade_matrix((res-1)*nelms+n_loc ,reac) +&
         cascade_matrix((mono-1)*nelms+n_loc   ,reac)+&
         cascade_matrix(lid_nh4 ,  reac)
      write(*,*)'res phosp',cascade_matrix((res-1)*nelms+p_loc ,reac) +&
         cascade_matrix((mono-1)*nelms+p_loc   ,reac)+&
         cascade_matrix(lid_minp_soluble,  reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 9, the partition cwd into poly and mono
    reac = cwd_dek_reac
    !cwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*poly+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*mono+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(poly)-f2/cn_ratios(mono))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(poly)-f2/cp_ratios(mono))
    f1 = cwd_fcel*(1._r8-rf_l2s1_bgc(lay))
    f2 = (1._r8-cwd_fcel)*(1._r8-rf_l3s2_bgc)

    call wood_decomp_cascade(fwd, reac, f1, f2)

    !---------------------------------------------------------------------------------
    !reaction 8, the partition lwd into som1 and som2
    reac = lwd_dek_reac
    !lwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = lwd_fcel*(1._r8-rf_l2s1_bgc(lay))
    f2 = (1._r8-lwd_fcel)*(1._r8-rf_l3s2_bgc)

    call wood_decomp_cascade(lwd, reac, f1, f2)

    !---------------------------------------------------------------------------------
    !reaction 9, the partition fwd into som1 and som2
    reac = fwd_dek_reac
    !fwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = fwd_fcel*(1._r8-rf_l2s1_bgc(lay))
    f2 = (1._r8-fwd_fcel)*(1._r8-rf_l3s2_bgc)

    call wood_decomp_cascade(fwd, reac, f1, f2)

  end associate
  contains

    subroutine wood_decomp_cascade(iwd, reac, f1, f2)

    implicit none
    integer , intent(in) :: iwd, reac
    real(r8), intent(in):: f1, f2
    associate(                                                   &
      c_loc     => summsbgc_index%c_loc                      , & !
      n_loc     => summsbgc_index%n_loc                      , & !
      p_loc     => summsbgc_index%p_loc                      , & !
      c13_loc   => summsbgc_index%c13_loc                    , & !
      c14_loc   => summsbgc_index%c14_loc                    , & !
      nelms     => summsbgc_index%nelms                      , & !
      lid_o2    => summsbgc_index%lid_o2                     , & !
      lid_co2   => summsbgc_index%lid_co2                    , & !
      lid_nh4   => summsbgc_index%lid_nh4                    , & !
      lid_c14_co2=> summsbgc_index%lid_c14_co2               , & !
      lid_c13_co2=> summsbgc_index%lid_c13_co2               , & !
      lid_co2_hr => summsbgc_index%lid_co2_hr                , &
      lid_minn_nh4_immob=> summsbgc_index%lid_minn_nh4_immob , &
      lid_minp_immob => summsbgc_index%lid_minp_immob        , &
      lid_minp_soluble=> summsbgc_index%lid_minp_soluble     , &
      poly      => summsbgc_index%poly                       , & !
      mono      => summsbgc_index%mono                         & !
    )
    cascade_matrix((iwd-1)*nelms+c_loc    ,reac) = -1._r8
    cascade_matrix((iwd-1)*nelms+n_loc    ,reac) = -this%icn_ratios(iwd)
    cascade_matrix((iwd-1)*nelms+p_loc    ,reac) = -this%icp_ratios(iwd)

    cascade_matrix((poly-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((poly-1)*nelms+n_loc   ,reac) = f1*this%icn_ratios(poly)
    cascade_matrix((poly-1)*nelms+p_loc   ,reac) = f1*this%icp_ratios(poly)

    cascade_matrix((mono-1)*nelms+c_loc   ,reac) = f2
    cascade_matrix((mono-1)*nelms+n_loc   ,reac) = f2*this%icn_ratios(mono)
    cascade_matrix((mono-1)*nelms+p_loc   ,reac) = f2*this%icp_ratios(mono)

    cascade_matrix(lid_co2                ,reac) = 1._r8-f1-f2

    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((iwd-1)*nelms+n_loc    ,reac)  &
                                                   -cascade_matrix((poly-1)*nelms+n_loc   ,reac)  &
                                                   -cascade_matrix((mono-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((iwd-1)*nelms+p_loc    ,reac)  &
                                                     -cascade_matrix((poly-1)*nelms+p_loc   ,reac)  &
                                                     -cascade_matrix((mono-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((iwd-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(iwd)
      cascade_matrix((poly-1)*nelms+c14_loc  , reac) =  f1*this%icc14_ratios(iwd)
      cascade_matrix((mono-1)*nelms+c14_loc  , reac) =  f2*this%icc14_ratios(iwd)
    endif

    if(this%use_c14)then
      cascade_matrix((iwd-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(iwd)
      cascade_matrix((poly-1)*nelms+c13_loc  , reac) =  f1*this%icc13_ratios(iwd)
      cascade_matrix((mono-1)*nelms+c13_loc  , reac) =  f2*this%icc13_ratios(iwd)
    endif
    ! remember to define debug in this subroutine if want use
    ! if(debug)then
    !   write(*,*)'iwd carbon',   cascade_matrix((iwd-1)*nelms+c_loc    ,reac) + cascade_matrix((mono-1)*nelms+c_loc   ,reac) + &
    !     cascade_matrix((poly-1)*nelms+c_loc   ,reac) + cascade_matrix(lid_co2                ,reac)
    !   write(*,*)'iwd nitrogen', cascade_matrix((iwd-1)*nelms+n_loc    ,reac) + cascade_matrix((mono-1)*nelms+n_loc   ,reac) + &
    !     cascade_matrix((poly-1)*nelms+n_loc   ,reac) + cascade_matrix(lid_nh4         ,reac)
    !   write(*,*)'iwd phosp',    cascade_matrix((iwd-1)*nelms+p_loc    ,reac) + cascade_matrix((mono-1)*nelms+p_loc   ,reac) + &
    !     cascade_matrix((poly-1)*nelms+p_loc   ,reac) + cascade_matrix(lid_minp_soluble         ,reac)
    ! endif
    end associate
    end subroutine wood_decomp_cascade
  end subroutine calc_cascade_matrix

  !-----------------------------------------------------------------------
  subroutine calc_potential_aerobic_hr(this, summsbgc_index, pot_decay_rates, &
    cascade_matrix, pot_co2_hr, bstatus)
    !
    ! DESCRIPTION:
    ! calculate potential aerobic heteorotrophic respiration, and potential oxygen consumption based on cascade_matrix
    ! !USES:
    use MathfuncMod         , only : dot_sum
    use MathfuncMod         , only : safe_div
    use BgcSummsIndexType   , only : summsbgc_index_type
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(SummsSom_type)    , intent(inout) :: this
    type(summsbgc_index_type) , intent(in) :: summsbgc_index
    real(r8)                , intent(in) :: pot_decay_rates(nsummspools)
    real(r8)                , intent(in) :: cascade_matrix(summsbgc_index%nstvars, summsbgc_index%nreactions)
    real(r8)                , intent(out):: pot_co2_hr
    type(betr_status_type)  , intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    real(r8) :: cascade_matrix_hr(nsummspools)
    integer  :: reac

    associate(                                           & !
         nom_pools => summsbgc_index%nom_pools        , & !
         lid_co2_hr=> summsbgc_index%lid_co2_hr       , & !
         lit1      => summsbgc_index%lit1             , & !
         lit2      => summsbgc_index%lit2             , & !
         lit3      => summsbgc_index%lit3             , & !
         mic      => summsbgc_index%mic               , & !
         res      => summsbgc_index%res               , & !
         enz      => summsbgc_index%enz               , & !
         mono      => summsbgc_index%mono             , & !
         poly      => summsbgc_index%poly             , & !
         cwd       => summsbgc_index%cwd              , & !
         lwd       => summsbgc_index%lwd              , & !
         fwd       => summsbgc_index%fwd              , & !
         lit1_dek_reac=> summsbgc_index%lit1_dek_reac , & !
         lit2_dek_reac=> summsbgc_index%lit2_dek_reac , & !
         lit3_dek_reac=> summsbgc_index%lit3_dek_reac , & !
         poly_dek_reac=> summsbgc_index%poly_dek_reac , & !
         mono_dek_reac=> summsbgc_index%mono_dek_reac , & !
         mic_dek_reac=>  summsbgc_index%mic_dek_reac  , & !
         enz_dek_reac=>  summsbgc_index%enz_dek_reac  , & !
         res_dek_reac=>  summsbgc_index%res_dek_reac  , & !
         cwd_dek_reac=> summsbgc_index%cwd_dek_reac   , & !
         lwd_dek_reac=> summsbgc_index%lwd_dek_reac   , & !
         fwd_dek_reac=> summsbgc_index%fwd_dek_reac     & !
         )

    cascade_matrix_hr = 0._r8
    reac=lit1_dek_reac; cascade_matrix_hr(lit1)=cascade_matrix(lid_co2_hr,reac)
    reac=lit2_dek_reac; cascade_matrix_hr(lit2)=cascade_matrix(lid_co2_hr,reac)
    reac=lit3_dek_reac; cascade_matrix_hr(lit3)=cascade_matrix(lid_co2_hr,reac)
    reac=poly_dek_reac; cascade_matrix_hr(poly)=cascade_matrix(lid_co2_hr,reac)
    reac=mono_dek_reac; cascade_matrix_hr(mono)=cascade_matrix(lid_co2_hr,reac)
    reac=mic_dek_reac ; cascade_matrix_hr(mic) =cascade_matrix(lid_co2_hr,reac)
    reac=enz_dek_reac ; cascade_matrix_hr(enz) =cascade_matrix(lid_co2_hr,reac)
    reac=res_dek_reac ; cascade_matrix_hr(res) =cascade_matrix(lid_co2_hr,reac)
    reac=cwd_dek_reac ; cascade_matrix_hr(cwd) =cascade_matrix(lid_co2_hr,reac)
    reac=lwd_dek_reac ; cascade_matrix_hr(lwd) =cascade_matrix(lid_co2_hr,reac)
    reac=fwd_dek_reac ; cascade_matrix_hr(fwd) =cascade_matrix(lid_co2_hr,reac)

    pot_co2_hr = dot_sum(cascade_matrix_hr, pot_decay_rates, bstatus)  !mol CO2/m3/s
    end associate
  end subroutine calc_potential_aerobic_hr

  !-----------------------------------------------------------------------
  subroutine calc_cnp_ratios(this, summsbgc_index, ystates, bstatus)
  !
  ! DESCRIPTION
  ! compute the cnp ratios for the om pools
  use BetrStatusType      , only : betr_status_type
  use MathfuncMod         , only : safe_div
  use BgcSummsIndexType       , only : summsbgc_index_type
  implicit none
  class(SummsSom_type)        , intent(inout) :: this
  type(summsbgc_index_type)   , intent(in) :: summsbgc_index
  real(r8)                    , intent(inout) :: ystates(summsbgc_index%nstvars)
  type(betr_status_type)      , intent(out) :: bstatus
  integer :: jj
  integer :: kc, kn, kp, kc13, kc14, kc1, kc2
  real(r8):: rat
  real(r8) :: difn
  real(r8) :: stoibal_ncon
  character(len=255) :: msg
  real(r8), parameter :: tiny_val=1.e-14_r8
  real(r8), parameter :: tiny_ncon = 1.e-15_r8 
  associate(                         &
    nelms => summsbgc_index%nelms, &
    c_loc => summsbgc_index%c_loc, &
    n_loc => summsbgc_index%n_loc, &
    p_loc => summsbgc_index%p_loc, &
    c13_loc => summsbgc_index%c13_loc, &
    c14_loc => summsbgc_index%c14_loc, &
    lit2  => summsbgc_index%lit2, &
    lit3  => summsbgc_index%lit3, &
    is_sumpool_som => summsbgc_index%is_sumpool_som, &
    ompoolnames => summsbgc_index%ompoolnames & 
  )
  
call bstatus%reset()

  !for om pools
  do jj = 1, nsummspools
    kc = (jj-1) * nelms + c_loc
    kn = (jj-1) * nelms + n_loc
    kp = (jj-1) * nelms + p_loc
  if(ystates(kc)<tiny_val)then
      rat = 0._r8
    else
      rat=ystates(kc)/(ystates(kc)+tiny_val)
    endif
    if(ystates(kn)<tiny_val*this%def_cn(jj) .or. ystates(kc)<tiny_val)then
      this%icn_ratios(jj)= 1._r8/this%def_cn(jj)
    else
      this%icn_ratios(jj) = 1._r8/this%def_cn(jj)*(1._r8-rat)+ystates(kn)/ystates(kc)*rat
    endif
    if(ystates(kp)<tiny_val*this%def_cp(jj) .or. ystates(kc)<tiny_val)then
      this%icp_ratios(jj)=1._r8/this%def_cp(jj)
    else
      this%icp_ratios(jj) = 1._r8/this%def_cp(jj)*(1._r8-rat)+ystates(kp)/ystates(kc)*rat
    endif
    if(summsbgc_index%debug)then
       write(*,'(A,X,I2,5(X,E20.10))')'cnp',jj,ystates(kc),ystates(kn),ystates(kp),1._r8/this%icn_ratios(jj),1._r8/this%icp_ratios(jj)
    endif
    if(is_sumpool_som(jj) .and. ystates(kc)>tiny_val)then
      stoibal_ncon = ystates(kc)*this%icn_ratios(jj)
      difn=ystates(kn)-stoibal_ncon
      if(difn<-tiny_ncon)then
        ystates(kn)=stoibal_ncon
      endif
 !     write(msg,*)'phosphorus weirdo',jj,trim(ompoolnames(jj)),ystates(kc),ystates(kn),ystates(kp), rat, this%def_cn(jj),this%def_cp(jj),&
 !        1._r8/this%icn_ratios(jj),1._r8/this%icp_ratios(jj) 
 !     print*,msg
 !     call bstatus%set_msg(msg,err=-1)
 !     return
    endif
    if(this%use_c14)then
      kc14 = (jj-1) * nelms + c14_loc
      this%icc14_ratios(jj) = 1._r8/this%def_cc14(jj)*(1._r8-rat)+ystates(kc14)/ystates(kc)
      if(summsbgc_index%debug)then
        write(*,'(A,X,I4,2(X,E20.10))') 'c14rrr som jj',jj,1._r8/this%def_cc14(jj),this%icc14_ratios(jj)
      endif
    endif
    if(this%use_c13)then
      kc13 = (jj-1) * nelms + c13_loc
      kc13 = (jj-1) * nelms + c13_loc
      this%icc13_ratios(jj) = 1._r8/this%def_cc13(jj)*(1._r8-rat)+ystates(kc13)/ystates(kc)
    endif

  enddo
  kc1 = (lit2-1)*nelms+c_loc
  kc2 = (lit3-1)*nelms+c_loc
  !lignin fraction of the structural carbon
  this%lit_flig = safe_div(ystates(kc2),ystates(kc1)+ystates(kc2))


  end associate
  end subroutine calc_cnp_ratios

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  subroutine stoichiometry_fix(this, summsbgc_index,ystates) 

  !
  ! DESCRIPTION
  ! this fixes the stoichiometric drift due to limite precision of 
  ! double precision.
  use BgcSummsIndexType         , only : summsbgc_index_type
  implicit none
  class(SummsSom_type)          , intent(inout) :: this
  type(summsbgc_index_type)     , intent(in) :: summsbgc_index
  real(r8)                      , intent(inout) :: ystates(summsbgc_index%nstvars)

  associate(                         &
    nelms => summsbgc_index%nelms, &
    c_loc => summsbgc_index%c_loc, &
    n_loc => summsbgc_index%n_loc, &
    c13_loc => summsbgc_index%c13_loc, &
    c14_loc => summsbgc_index%c14_loc, &
    poly  => summsbgc_index%poly , &
    mono  => summsbgc_index%mono , &
    mic  => summsbgc_index%mic , &
    enz  => summsbgc_index%enz , &
    res  => summsbgc_index%res , &
    is_sumpool_som => summsbgc_index%is_sumpool_som, &
    ompoolnames => summsbgc_index%ompoolnames &
  )


  !for om pools
  call stoi_fix(poly)

  call stoi_fix(mono)

  call stoi_fix(mic)

  call stoi_fix(enz)

  call stoi_fix(res)

  end associate
  contains
    subroutine stoi_fix(jj)
    implicit none
    integer, intent(in) :: jj
    real(r8) :: difn
    real(r8) :: stoibal_ncon
    integer  :: kc, kn
    real(r8), parameter :: tiny_ncon = 1.e-15_r8

    associate(                       &
      nelms => summsbgc_index%nelms, &
      c_loc => summsbgc_index%c_loc, &
      n_loc => summsbgc_index%n_loc  &
    )

    kc = (jj-1) * nelms + c_loc
    kn = (jj-1) * nelms + n_loc
    stoibal_ncon = ystates(kc)*this%icn_ratios(jj)
    difn=ystates(kn)-stoibal_ncon
    if(difn<-tiny_ncon)then
      ystates(kn)=stoibal_ncon
    endif
    end associate
    end subroutine stoi_fix
  end subroutine stoichiometry_fix

  !-------------------------------------------------------------------------------
subroutine calc_som_decay_r(this, summsbgc_index, dtime, om_k_decay, om_pools, om_decay_rates)
    !
    ! !DESCRIPTION:
    ! calculate degradation for all different pools
    !
    ! !USES:
    use BgcSummsIndexType       , only : summsbgc_index_type
   implicit none
   class(SummsSom_type)     , intent(inout) :: this
   type(summsbgc_index_type) , intent(in)    :: summsbgc_index
    real(r8)  , intent(in)    :: dtime
    real(r8)  , intent(in)    :: om_k_decay(nsummspools)
    real(r8)  , intent(in)    :: om_pools(summsbgc_index%nom_tot_elms)
    real(r8)  , intent(out)   :: om_decay_rates(nsummspools)

    ! !LOCAL VARIABLES:
    integer :: jj, fc, c, j
    integer :: kc, kn
    associate(                                        &
         nelms => summsbgc_index%nelms            , &
         nom_pools => summsbgc_index%nom_pools    , &
         c_loc => summsbgc_index%c_loc              &
    )

    !for om pools
    do jj = 1, nom_pools
      kc = (jj-1) * nelms + c_loc
      om_decay_rates(jj) = om_pools(kc) * om_k_decay(jj)
    enddo
    end associate
  end subroutine calc_som_decay_r

  !-------------------------------------------------------------------------------
 
  subroutine apply_spinupf(this, summsbgc_index, decompkf_eca, k_decay, spinup_scalar, spinup_flg)
  use BgcSummsIndexType       , only : summsbgc_index_type
  use BgcSummsDecompType      , only : DecompSumms_type
  use betr_varcon               , only : kyr_spinup
  implicit none
  class(SummsSom_type)     , intent(inout) :: this
  type(DecompSumms_type), intent(in) :: decompkf_eca
  type(summsbgc_index_type)     , intent(in)    :: summsbgc_index
  real(r8)                      , intent(inout) :: k_decay(nsummspools)
  real(r8)                      , intent(inout) :: spinup_scalar
  integer                       , intent(in)    :: spinup_flg
  integer :: jj

  associate(   &
   t_scalar       => decompkf_eca%t_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil temperature scalar for decomp
   w_scalar       => decompkf_eca%w_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil water scalar for decomp
   o_scalar       => decompkf_eca%o_scalar        , & ! Intput: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
   poly           => summsbgc_index%poly        , & !
   mono           => summsbgc_index%mono        , & !
   mic            => summsbgc_index%mic         , & !
   enz            => summsbgc_index%enz         , & !
   res            => summsbgc_index%res           & !
  )

  if(spinup_flg==2)then
    !cumulated more than 2 * kyr_spinup
    do jj = 1, nsummspools
      k_decay(jj) = k_decay(jj)/spinup_scalar
    enddo
  elseif(spinup_flg==1)then
    !cumulated more than kyr_spinup but less than 2 * kyr_spinup
    spinup_scalar = spinup_scalar + t_scalar * w_scalar * o_scalar / (365._r8 * 86400._r8 * kyr_spinup)
  endif

  end associate

  end subroutine apply_spinupf
  !-------------------------------------------------------------------------------
  subroutine calc_som_scale_k(this, lay, summsbgc_index, decompkf_eca, k_decay)

  use BgcSummsIndexType       , only : summsbgc_index_type
  use BgcSummsDecompType      , only : DecompSumms_type

  implicit none
  class(SummsSom_type)        , intent(inout) :: this
  integer                     , intent(in) :: lay
  type(DecompSumms_type)      , intent(in)    :: decompkf_eca
  type(summsbgc_index_type)   , intent(in)    :: summsbgc_index
  real(r8)                    , intent(inout)    :: k_decay(nsummspools)
  integer :: jj

  associate(   &
   t_scalar       => decompkf_eca%t_scalar        , & ! Input: [real(r8) (:,:)   ]  soil temperature scalar for decomp
   w_scalar       => decompkf_eca%w_scalar        , & ! Input: [real(r8) (:,:)   ]  soil water scalar for decomp
   o_scalar       => decompkf_eca%o_scalar        , & ! Input: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
   depth_scalar   => decompkf_eca%depth_scalar    , & ! Input: [real(r8) (:,:)   ]  rate constant for decomposition (1./sec)
   lit1           => summsbgc_index%lit1               , & !
   lit2           => summsbgc_index%lit2               , & !
   lit3           => summsbgc_index%lit3               , & !
   mic            => summsbgc_index%mic                , & !
   res            => summsbgc_index%res                , & !
   enz            => summsbgc_index%enz                , & !
   mono           => summsbgc_index%mono               , & !
   poly           => summsbgc_index%poly               , & !
   cwd            => summsbgc_index%cwd                , & !
   lwd            => summsbgc_index%lwd                , & !
   fwd            => summsbgc_index%fwd                  & !
  )

  k_decay(lit1) = k_decay(lit1) * t_scalar *  w_scalar * o_scalar * depth_scalar
  k_decay(lit2) = k_decay(lit2) * t_scalar *  w_scalar * o_scalar * depth_scalar
  k_decay(lit3) = k_decay(lit3) * t_scalar *  w_scalar * o_scalar * depth_scalar

  k_decay(poly) = k_decay(poly)  * w_scalar            
  k_decay(mic)  = k_decay(mic)   * o_scalar * w_scalar
  k_decay(res)  = k_decay(res)   * o_scalar * w_scalar  
  k_decay(enz)  = k_decay(enz)   * w_scalar            
  k_decay(mono) = k_decay(mono)  * o_scalar * w_scalar
  
  k_decay(cwd)  = k_decay(cwd)   * t_scalar *  w_scalar * o_scalar * depth_scalar
  k_decay(lwd)  = k_decay(lwd)   * t_scalar *  w_scalar * o_scalar * depth_scalar
  k_decay(fwd)  = k_decay(fwd)   * t_scalar *  w_scalar * o_scalar * depth_scalar
  !impose the ligin effect
  k_decay(cwd)  = k_decay(cwd) * exp(-3._r8*this%cwd_flig)
  k_decay(lwd)  = k_decay(lwd) * exp(-3._r8*this%lwd_flig)
  k_decay(fwd)  = k_decay(fwd) * exp(-3._r8*this%fwd_flig)
  k_decay(lit2) = k_decay(lit2)* exp(-3._r8*this%lit_flig)
  k_decay(lit3) = k_decay(lit3)* exp(-3._r8*this%lit_flig)

  end associate
  end subroutine calc_som_scale_k
  !-------------------------------------------------------------------------------
   subroutine calc_som_decay_k(this, lay, summsbgc_index, decompkf_eca, k_decay, ystates, bstatus)

  use BgcSummsIndexType       , only : summsbgc_index_type
  use BgcSummsDecompType      , only : DecompSumms_type
  use MathfuncMod             , only : safe_div
  !use FindRootMod             , only : brent
  !use func_data_type_mod      , only : func_data_type
  use BgcSummsMath            , only : brent
  use BgcSummsDebType         , only : debs
  use DebGrowMod              , only : deb_grow
  use BetrStatusType      , only : betr_status_type

  implicit none
  class(SummsSom_type)        , intent(inout)   :: this !this will update the relevant values for cascade_matrix
  integer                     , intent(in) :: lay
  type(DecompSumms_type)      , intent(in)    :: decompkf_eca
  type(summsbgc_index_type)   , intent(in)    :: summsbgc_index
  real(r8)                    , intent(in)    :: ystates(1:summsbgc_index%nom_tot_elms)
  !real(r8)                    , intent(out)   :: k_decay(nsummspools)
  real(r8)                    , intent(out)   :: k_decay(1:nsummspools)
  type(betr_status_type)        , intent(out)   :: bstatus

  !type(func_data_type) :: deb
  integer    :: jj
    real(r8) :: residual
    real(r8) :: y_mic
    real(r8) :: y_enz
    real(r8) :: y_res
    real(r8) :: y_poly
    real(r8) :: y_mono
    real(r8) :: dc
    integer  :: isgrw
    real(r8) :: scal_c
    real(r8) :: m0
    real(r8) :: g0
    real(r8) :: jeg
    real(r8) :: ev
    real(r8) :: aa
    real(r8) :: bb
    real(r8) :: cc
    real(r8) :: delta
    real(r8) :: jxx    
    real(r8), parameter :: macheps = 1.e-8_r8
    real(r8), parameter :: tol = 1.e-8_r8
    real(r8), parameter :: tiny_val=1.e-35_r8

  type(debs) :: deb
  allocate(deb%gB)
  allocate(deb%pE)
  allocate(deb%je)
  allocate(deb%ec)
  allocate(deb%mr_mic) ; deb%mr_mic  = decompkf_eca%mr_mic
  allocate(deb%yld_enz) ; deb%yld_enz  = this%yld_enz
  allocate(deb%gmax_mic) ; deb%gmax_mic  = this%gmax_mic
  allocate(deb%yld_mic) ; deb%yld_mic  = this%yld_mic
  allocate(deb%pmax_enz) ; deb%pmax_enz  = this%pmax_enz

  associate(   &
   gB => deb%gB                       , &
   pE => deb%pE                       , &
   je => deb%je                       , &
   ec => deb%ec                       , &
   mr_mic => deb%mr_mic               , &
   yld_enz => deb%yld_enz             , &
   gmax_mic => deb%gmax_mic           , &
   yld_mic => deb%yld_mic             , &
   pmax_enz => deb%pmax_enz           , &
   nelms => summsbgc_index%nelms                       , &
   c_loc => summsbgc_index%c_loc                       , &
   vmax_mic         => decompkf_eca%vmax_mic           , &
   vmax_enz         => decompkf_eca%vmax_enz           , &
   kaff_mono_mic    => decompkf_eca%kaff_mono_mic      , &
   kaff_enz_poly    => decompkf_eca%kaff_enz_poly      , &
   kappa_mic        => decompkf_eca%kappa_mic          , &
   kaff_mono_msurf  => decompkf_eca%kaff_mono_msurf    , &
   kaff_enz_msurf   => decompkf_eca%kaff_enz_msurf     , &
   lit1           => summsbgc_index%lit1               , & !
   lit2           => summsbgc_index%lit2               , & !
   lit3           => summsbgc_index%lit3               , & !
   mic            => summsbgc_index%mic                , & !
   res            => summsbgc_index%res                , & !
   enz            => summsbgc_index%enz                , & !
   mono           => summsbgc_index%mono               , & !
   poly           => summsbgc_index%poly               , & !
   cwd            => summsbgc_index%cwd                , & !
   lwd            => summsbgc_index%lwd                , & !
   fwd            => summsbgc_index%fwd                , & !
   yld_res        => this%yld_res                      , & !
   fenz2poly      => this%fenz2poly                    , & !
   minsite        => this%minsite                      , & !
   mic_transp     => this%mic_transp                   , & !
   decay_mic0     => this%decay_mic0                   , & !
   decay_mic1     => this%decay_mic1                   , & !
   decay_enz      => this%decay_enz                    , & !
   actgB          => this%actgB                        , &
   actpE          => this%actpE                        , &
   actmr          => this%actmr                        , &
   decay_mic      => this%decay_mic                      & 
  )
  call bstatus%reset()


      y_poly = ystates((poly-1) * nelms + c_loc)
      y_mono = ystates((mono-1) * nelms + c_loc)
      y_mic=ystates((mic-1) * nelms + c_loc)
      y_res=ystates((res-1) * nelms + c_loc)
      y_enz=ystates((enz-1) * nelms + c_loc)

  !define  gB, kappa, pE, mr, residual, decay_mic
  !ec = safe_div(y_res,y_mic)
  ec = safe_div(y_res*y_mic , y_mic**(2._r8+1.e-20_r8))
  je = kappa_mic*ec

        ! Determine if there is enough C for growth after maintenance
        dc=je-mr_mic
        
        isgrw=-1
        if (dc>0._r8) then
            ! There is carbon to support growth activity
            ! Compute the actual carbon flux to support growth, the yield rate < 1, 
            ! indicating the fraction of carbon being turned into
            ! the required structure after taking off the overhead
                   ! Maximum carbon export
                   scal_c=safe_div( dc-gmax_mic*ec , safe_div(gmax_mic,yld_mic) + safe_div(pmax_enz,yld_enz) )
                   if (scal_c>=1) then
                       ! Maximum growth
                       isgrw=1
                   else
                       ! Less than maximum growth
                       isgrw=0
                   end if
        end if


        select case (isgrw)    
            case (-1)            
                ! No growth        
                ! There is no flux for growth, but penalty for mortality        
                actgB=0._r8      
                actpE=0._r8      
                actmr=je          

            case (0)    
                ! Limited growth              
                ! Potential growth rate > 0         

                ! Solve for the growth rate
                ! Actual DEB calculation
                    gB=0._r8
                    pE=0._r8
                    
                    m0=safe_div(mr_mic,gmax_mic)
                    jeg=safe_div(je,gmax_mic)
                    g0=safe_div(pmax_enz,gmax_mic)
                    ev=ec

                    !call deb_grow ( gb0, deb, fb1 )
                    !call deb_grow ( gb1, deb, fb2 )
                    !call brent(gbtemp, gb0, gb1, fb1, fb2, macheps, tol, deb, deb_grow, bstatus)

                    aa= yld_mic*yld_enz + safe_div(1._r8,ev)*(yld_enz + yld_mic*g0)
                    bb= yld_mic*g0 - safe_div(jeg - m0 , ev) * (yld_enz + yld_mic*g0) +2._r8*safe_div(g0,ev) 
                    cc= -(jeg - m0) * safe_div(2._r8,ev) * g0
                    delta= bb*bb - 4._r8*aa*cc
                    jxx=safe_div(-bb + sqrt(delta) , 2._r8*aa)
                    gB=safe_div(jeg - m0 - jxx , ev)
                    pE=(jxx - safe_div(gB,yld_mic))*yld_enz

                ! Population growth        
                actgB=gB*gmax_mic   !deb%gB
                ! Enzyme production        
                actpE=pE*gmax_mic   !deb%pE
                ! Maintenance
                actmr=mr_mic                      

            case (1)    
                ! Maximum growth       
                actgB=gmax_mic        
                ! Maximum enzyme production        
                actpE=pmax_enz        
                ! Maintenance        
                actmr=mr_mic

        end select

        ! Compute density-dependent specific microbial mortality [mol c/ut]   
              decay_mic=decay_mic0*y_mic!**0.5_r8
              !decay_mic=safe_div(decay_mic0*y_mic , decay_mic1+y_mic )

        ! Compute residual
              !residual = (kappa_mic-actgB)*y_res-(actmr + safe_div(actpE,yld_enz) + safe_div(actgB,yld_mic) )*y_mic

  k_decay(poly) = safe_div( y_enz*vmax_enz , kaff_enz_poly+y_poly+y_enz+minsite*safe_div( kaff_enz_poly , kaff_enz_msurf ))
  k_decay(mono) = safe_div( y_mic*vmax_mic*mic_transp , kaff_mono_mic+y_mono+ mic_transp*y_mic +minsite*safe_div( kaff_mono_mic ,&
  kaff_mono_msurf ))
  k_decay(mic) = decay_mic-actgB
  k_decay(res) = kappa_mic-actgB+decay_mic
  k_decay(enz) = decay_enz
  k_decay(lit1) = this%k_decay_lit1(lay)
  k_decay(lit2) = this%k_decay_lit2(lay)
  k_decay(lit3) = this%k_decay_lit3(lay)
  k_decay(cwd) = this%k_decay_cwd
  k_decay(lwd) = this%k_decay_lwd
  k_decay(fwd) = this%k_decay_fwd

  !partitioning for cascade_matrix
      this%part_mic2poly = safe_div(decay_mic , (decay_mic-actgB))

      this%part_mic2enz = safe_div(pmax_enz , (decay_mic-actgB))

      this%part_res2mono = safe_div(decay_mic , (kappa_mic-actgB+decay_mic))

      this%rate_co2 = k_decay(mono)*y_mono*(1-yld_res) + y_res*(kappa_mic-actgB) - y_mic*(actpE+actgB)
      this%cue = 1 - safe_div(this%rate_co2,k_decay(mono)*y_mono)

  end associate
  end subroutine calc_som_decay_k
  !-------------------------------------------------------------------------------
  subroutine calc_pot_min_np_flx(this, dtime, summsbgc_index, ystates, k_decay, cascade_matrix, &
    alpha_n, alpha_p, pot_decomp, pot_nn_flx, pot_np_flx)
  use BgcSummsIndexType       , only : summsbgc_index_type
  implicit none
  class(SummsSom_type)        , intent(inout) :: this
  real(r8)                    , intent(in) :: dtime
  type(summsbgc_index_type)   , intent(in) :: summsbgc_index
  real(r8)                    , intent(in) :: ystates(1:summsbgc_index%nom_tot_elms)
  real(r8)                    , intent(in) :: k_decay(1:nsummspools)
  real(r8)                    , intent(in) :: cascade_matrix(summsbgc_index%nstvars, summsbgc_index%nreactions)
  real(r8)                    , intent(in) :: alpha_n(nsummspools)
  real(r8)                    , intent(in) :: alpha_p(nsummspools)
  real(r8)                    , intent(out) :: pot_decomp(nsummspools)
  real(r8)                    , intent(out):: pot_nn_flx
  real(r8)                    , intent(out):: pot_np_flx

  integer :: reac
  integer :: reacs(nsummspools)

  associate(                                                    & !
       nom_pools => summsbgc_index%nom_pools                , & !
       nom_tot_elms=> summsbgc_index%nom_tot_elms           , & !
       lid_nh4   => summsbgc_index%lid_nh4                  , & !
       lid_minp_soluble  => summsbgc_index%lid_minp_soluble , & !
       lit1      => summsbgc_index%lit1                     , & !
       lit2      => summsbgc_index%lit2                     , & !
       lit3      => summsbgc_index%lit3                     , & !
       mic      => summsbgc_index%mic                       , & !
       res      => summsbgc_index%res                       , & !
       enz      => summsbgc_index%enz                       , & !
       mono      => summsbgc_index%mono                     , & !
       poly      => summsbgc_index%poly                     , & !
       cwd       => summsbgc_index%cwd                      , & !
       lit1_dek_reac=> summsbgc_index%lit1_dek_reac         , & !
       lit2_dek_reac=> summsbgc_index%lit2_dek_reac         , & !
       lit3_dek_reac=> summsbgc_index%lit3_dek_reac         , & !
       poly_dek_reac=> summsbgc_index%poly_dek_reac         , & !
       mono_dek_reac=> summsbgc_index%mono_dek_reac         , & !
       mic_dek_reac=> summsbgc_index%mic_dek_reac           , & !
       enz_dek_reac=> summsbgc_index%enz_dek_reac           , & !
       res_dek_reac=> summsbgc_index%res_dek_reac           , & !
       cwd_dek_reac=> summsbgc_index%cwd_dek_reac           , & !
       lwd_dek_reac=> summsbgc_index%lwd_dek_reac           , & !
       fwd_dek_reac=> summsbgc_index%fwd_dek_reac             & !
   )

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(summsbgc_index, dtime, k_decay(1:nom_pools), &
      ystates(1:nom_tot_elms), pot_decomp)

  pot_nn_flx = 0._r8; pot_np_flx = 0._r8

  reacs=(/lit1_dek_reac, lit2_dek_reac, lit3_dek_reac, &
    cwd_dek_reac, lwd_dek_reac, fwd_dek_reac, &
    poly_dek_reac, mono_dek_reac, mic_dek_reac, enz_dek_reac, res_dek_reac/)

  do reac = 1, nom_pools
    if(alpha_n(reac)>0._r8)then
      pot_nn_flx = pot_nn_flx - cascade_matrix(lid_nh4, reacs(reac)) * pot_decomp(reac)
    endif
    if(alpha_p(reac)>0._r8)then
      pot_np_flx = pot_np_flx - cascade_matrix(lid_minp_soluble, reacs(reac)) * pot_decomp(reac)
    endif
  enddo
  end associate
  end subroutine calc_pot_min_np_flx
  !-------------------------------------------------------------------------------
end module BgcSummsSOMType
