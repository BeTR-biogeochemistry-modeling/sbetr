module cdomBGCSOMType
!
!DESCRIPTION
!module defines the century decomposition
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !totally 7 om pools
  !met, cell, lig, cwd, mic, pom, humus
  !during decomposition, the cnp ratios of mic, pom, and humus are all fixed
  !the litter pools have their cnp ratios varying with nutrient status
  !We consider mic and pom as DOM, humus as humus pool
  !
  integer :: ncentpools
  type, public :: cdomSom_type
    real(r8), pointer :: icn_ratios(:)  => null()
    real(r8), pointer :: icp_ratios(:) => null()
    real(r8), pointer :: icc14_ratios(:)=> null()
    real(r8), pointer :: icc13_ratios(:)=> null()
    real(r8), pointer :: def_cn(:)=> null()
    real(r8), pointer :: def_cp(:)=> null()
    real(r8), pointer :: def_cc13(:)=> null()
    real(r8), pointer :: def_cc14(:)=> null()

    !private parameters
    real(r8) :: rf_l1s1_bgc(2)  !respired co2 as a fraction of decomposing flux
    real(r8) :: rf_l2s1_bgc(2)
    real(r8) :: rf_l3s2_bgc
    real(r8) :: rf_s2s1_bgc
    real(r8) :: rf_s3s1_bgc
    real(r8) :: rf_s1s2a_bgc(2)
    real(r8) :: rf_s1s2b_bgc(2)
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8) :: lwd_fcel
    real(r8) :: lwd_flig
    real(r8) :: fwd_fcel
    real(r8) :: fwd_flig
    real(r8) :: lit_flig
    real(r8) :: k_decay_lmet(2)
    real(r8) :: k_decay_lcel(2)
    real(r8) :: k_decay_llig(2)
    real(r8) :: k_decay_mic(2)
    real(r8) :: k_decay_pom
    real(r8) :: k_decay_humus
    real(r8) :: k_decay_cwd  !coarse root
    real(r8) :: k_decay_lwd  !large wood
    real(r8) :: k_decay_fwd  !fine branch wood
    real(r8) :: Kaff_enz
    real(r8) :: alpha_E
    real(r8) :: alpha_T
    real(r8) :: rf_doms1_bgc
    real(r8) :: k_decay_dom
    real(r8) :: Kaff_mic_dom
    real(r8) :: Kaff_CM       !dom affinity to soil minerals
    logical  :: use_c13
    logical  :: use_c14
  contains
    procedure, public :: Init
    procedure, public :: run_decomp
    procedure, public :: calc_pot_min_np_flx
    procedure, public :: stoichiometry_fix
    procedure, private :: calc_cascade_matrix
    procedure, private :: calc_som_decay_k
    procedure, private :: calc_som_decay_r
    procedure, public  :: calc_cnp_ratios
    procedure, private :: InitAllocate
    procedure, public  :: UpdateParas
    procedure, private :: calc_potential_aerobic_hr
    procedure, private :: apply_spinupf
    procedure, private :: calc_dom_rf
  end type cdomSom_type
  real(r8), private, parameter :: tiny_val=1.e-14_r8
  real(r8), private, parameter :: tiny_ncon = 1.e-15_r8
contains

  subroutine Init(this, cdom_bgc_index, biogeo_con, bstatus)

  use cdomBGCIndexType , only : cdom_bgc_index_type
  use cdomParaType       , only : cdomPara_type
  use BetrStatusType      , only : betr_status_type
  implicit none
  class(cdomSom_type)         , intent(inout) :: this
  type(cdom_bgc_index_type) , intent(in)    :: cdom_bgc_index
  type(cdomPara_type)        , intent(in)    :: biogeo_con
  type(betr_status_type)      , intent(out)   :: bstatus

  call bstatus%reset()
  ncentpools = cdom_bgc_index%nom_pools

  call this%InitAllocate()

  end subroutine Init

!------------------------------------------
  subroutine calc_dom_rf(this, cdom_bgc_index, ystates)

  !DESCRIPTION
  !compute the respiration coefficient for the dom pool during decomposition
  use cdomBGCIndexType , only : cdom_bgc_index_type
  use tracer_varcon    , only : catomw, natomw, patomw
  implicit none
  class(cdomSom_type)  , intent(inout) :: this
  type(cdom_bgc_index_type) , intent(in) :: cdom_bgc_index
  real(r8)                 , intent(in) :: ystates(1:cdom_bgc_index%nom_tot_elms)

  integer :: kc, ke
  associate(                         &
    nelms => cdom_bgc_index%nelms   , &
    c_loc => cdom_bgc_index%c_loc   , &
    e_loc => cdom_bgc_index%e_loc   , &
    dom   => cdom_bgc_index%dom       &
  )

  kc = (dom-1)*nelms + c_loc
  ke = (dom-1)*nelms + e_loc
  if(abs(ystates(kc))<tiny_val)then
    this%rf_doms1_bgc = 0._r8
  else
    this%rf_doms1_bgc = ystates(ke)/ystates(kc)
  endif
  end associate
  end subroutine calc_dom_rf

!------------------------------------------
  subroutine InitAllocate (this)
  !
  !
  implicit none
  class(cdomSom_type), intent(inout) :: this

  allocate(this%icn_ratios(ncentpools));
  allocate(this%icp_ratios(ncentpools));
  allocate(this%icc14_ratios(ncentpools)); this%icc14_ratios(:) = 0._r8
  allocate(this%icc13_ratios(ncentpools)); this%icc13_ratios(:) = 0._r8
  allocate(this%def_cn(ncentpools));
  allocate(this%def_cp(ncentpools));
  allocate(this%def_cc13(ncentpools));this%def_cc13(:) = 0._r8
  allocate(this%def_cc14(ncentpools));this%def_cc14(:) = 0._r8
  end subroutine InitAllocate
!------------------------------------------
  subroutine UpdateParas(this,cdom_bgc_index,  biogeo_con)
  !
  ! intialize model parameters
  use cdomParaType , only : cdomPara_type
  use cdomBGCIndexType , only : cdom_bgc_index_type
  use tracer_varcon    , only : catomw, natomw, patomw
  implicit none
  class(cdomSom_type)  , intent(inout) :: this
  type(cdom_bgc_index_type) , intent(in) :: cdom_bgc_index
  type(cdomPara_type) , intent(in)    :: biogeo_con

  this%k_decay_dom    = biogeo_con%k_decay_dom
  this%Kaff_mic_dom   = biogeo_con%Kaff_mic_dom
  this%Kaff_CM        = biogeo_con%Kaff_CM
  this%rf_l1s1_bgc    = biogeo_con%rf_l1s1_bgc
  this%rf_l2s1_bgc    = biogeo_con%rf_l2s1_bgc
  this%rf_l3s2_bgc    = biogeo_con%rf_l3s2_bgc
  this%rf_s2s1_bgc    = biogeo_con%rf_s2s1_bgc
  this%rf_s3s1_bgc    = biogeo_con%rf_s3s1_bgc
  this%rf_s1s2a_bgc   = biogeo_con%rf_s1s2a_bgc
  this%rf_s1s2b_bgc   = biogeo_con%rf_s1s2b_bgc
  this%cwd_fcel   = biogeo_con%cwd_fcel_bgc
  this%cwd_flig   = biogeo_con%cwd_flig_bgc
  this%lwd_fcel   = biogeo_con%lwd_fcel_bgc
  this%lwd_flig   = biogeo_con%lwd_flig_bgc
  this%fwd_fcel   = biogeo_con%fwd_fcel_bgc
  this%fwd_flig   = biogeo_con%fwd_flig_bgc

  this%k_decay_lmet   =  biogeo_con%k_decay_lmet
  this%k_decay_lcel   =  biogeo_con%k_decay_lcel
  this%k_decay_llig   =  biogeo_con%k_decay_llig
  this%k_decay_mic    =  biogeo_con%k_decay_mic
  this%k_decay_pom    =  biogeo_con%k_decay_pom
  this%k_decay_humus  =  biogeo_con%k_decay_humus
  this%k_decay_cwd    =  biogeo_con%k_decay_cwd
  this%k_decay_lwd    =  biogeo_con%k_decay_lwd
  this%k_decay_fwd    =  biogeo_con%k_decay_fwd
  this%alpha_E        =  biogeo_con%alpha_E
  this%alpha_T        =  biogeo_con%alpha_T
  this%Kaff_enz       =  biogeo_con%Kaff_enz
  this%def_cn(cdom_bgc_index%lmet) = biogeo_con%init_cn_met * natomw/catomw
  this%def_cn(cdom_bgc_index%lcel) = biogeo_con%init_cn_cel * natomw/catomw
  this%def_cn(cdom_bgc_index%llig) = biogeo_con%init_cn_lig * natomw/catomw
  this%def_cn(cdom_bgc_index%cwd)  = biogeo_con%init_cn_cwd * natomw/catomw
  this%def_cn(cdom_bgc_index%lwd)  = biogeo_con%init_cn_lwd * natomw/catomw
  this%def_cn(cdom_bgc_index%fwd)  = biogeo_con%init_cn_fwd * natomw/catomw
  this%def_cn(cdom_bgc_index%dom)  = biogeo_con%init_cn_dom * natomw/catomw

  this%def_cn(cdom_bgc_index%mic) = biogeo_con%init_cn_mic * natomw/catomw
  this%def_cn(cdom_bgc_index%pom) = biogeo_con%init_cn_pom * natomw/catomw
  this%def_cn(cdom_bgc_index%humus) = biogeo_con%init_cn_humus * natomw/catomw

  this%def_cp(cdom_bgc_index%lmet) = biogeo_con%init_cp_met * patomw/catomw
  this%def_cp(cdom_bgc_index%lcel) = biogeo_con%init_cp_cel * patomw/catomw
  this%def_cp(cdom_bgc_index%llig) = biogeo_con%init_cp_lig * patomw/catomw
  this%def_cp(cdom_bgc_index%cwd)  = biogeo_con%init_cp_cwd * patomw/catomw
  this%def_cp(cdom_bgc_index%lwd)  = biogeo_con%init_cp_lwd * patomw/catomw
  this%def_cp(cdom_bgc_index%fwd)  = biogeo_con%init_cp_fwd * patomw/catomw

  this%def_cp(cdom_bgc_index%dom)  = biogeo_con%init_cp_dom * patomw/catomw
  this%def_cp(cdom_bgc_index%mic) = biogeo_con%init_cp_mic * patomw/catomw
  this%def_cp(cdom_bgc_index%pom) = biogeo_con%init_cp_pom * patomw/catomw
  this%def_cp(cdom_bgc_index%humus) = biogeo_con%init_cp_humus * patomw/catomw

  this%use_c13=biogeo_con%use_c13
  this%use_c14=biogeo_con%use_c14

  if(this%use_c13)then
    this%def_cc13(cdom_bgc_index%lmet) = biogeo_con%init_cc13_met
    this%def_cc13(cdom_bgc_index%lcel) = biogeo_con%init_cc13_cel
    this%def_cc13(cdom_bgc_index%llig) = biogeo_con%init_cc13_lig
    this%def_cc13(cdom_bgc_index%cwd)  = biogeo_con%init_cc13_cwd
    this%def_cc13(cdom_bgc_index%lwd)  = biogeo_con%init_cc13_lwd
    this%def_cc13(cdom_bgc_index%fwd)  = biogeo_con%init_cc13_fwd
    this%def_cc13(cdom_bgc_index%mic) = biogeo_con%init_cc13_mic
    this%def_cc13(cdom_bgc_index%pom) = biogeo_con%init_cc13_pom
    this%def_cc13(cdom_bgc_index%humus) = biogeo_con%init_cc13_humus
  endif

  if(this%use_c14)then
    this%def_cc14(cdom_bgc_index%lmet) = biogeo_con%init_cc14_met
    this%def_cc14(cdom_bgc_index%lcel) = biogeo_con%init_cc14_cel
    this%def_cc14(cdom_bgc_index%llig) = biogeo_con%init_cc14_lig
    this%def_cc14(cdom_bgc_index%cwd)  = biogeo_con%init_cc14_cwd
    this%def_cc14(cdom_bgc_index%lwd)  = biogeo_con%init_cc14_lwd
    this%def_cc14(cdom_bgc_index%fwd)  = biogeo_con%init_cc14_fwd
    this%def_cc14(cdom_bgc_index%mic) = biogeo_con%init_cc14_mic
    this%def_cc14(cdom_bgc_index%pom) = biogeo_con%init_cc14_pom
    this%def_cc14(cdom_bgc_index%humus)= biogeo_con%init_cc14_humus
  endif

  end subroutine UpdateParas
!------------------------------------------

  subroutine run_decomp(this, is_surflit, cdom_bgc_index, dtime, ystates,&
      decompkf_eca, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix, &
      k_decay, pot_co2_hr, bstatus)
  !
  !DESCRIPTION
  !
  use cdomBGCIndexType , only : cdom_bgc_index_type
  use cdomBGCDecompType, only : Decompcdom_type
  use BetrStatusType      , only : betr_status_type
  use betr_ctrl           , only : betr_spinup_state
  implicit none
  class(cdomSom_type)         , intent(inout) :: this
  type(cdom_bgc_index_type) , intent(in) :: cdom_bgc_index
  real(r8)                    , intent(in) :: dtime
  real(r8)                    , intent(inout) :: ystates(1:cdom_bgc_index%nom_tot_elms)
  type(Decompcdom_type)       , intent(in) :: decompkf_eca
  logical                     , intent(in) :: is_surflit
  real(r8)                    , intent(in) :: pct_sand
  real(r8)                    , intent(in) :: pct_clay
  real(r8)                    , intent(inout) :: cascade_matrix(cdom_bgc_index%nstvars, cdom_bgc_index%nreactions)
  real(r8)                    , intent(out) :: k_decay(1:ncentpools)
  real(r8)                    , intent(out) :: pot_co2_hr
  real(r8)                    , intent(out) :: alpha_n(1:ncentpools)
  real(r8)                    , intent(out) :: alpha_p(1:ncentpools)
  type(betr_status_type)      , intent(out) :: bstatus

  !local variables
  real(r8) :: pot_om_decay_rates(1:ncentpools)
  integer :: kc, jj, lay, mic_c, dom_c

  associate(                                   &
    nelms => cdom_bgc_index%nelms,              &
    nom_tot_elms=> cdom_bgc_index%nom_tot_elms, &
    micbiom_beg=> cdom_bgc_index%micbiom_beg  , &
    dom_beg    => cdom_bgc_index%dom_beg      , &
    c_loc => cdom_bgc_index%c_loc               &
  )
  call bstatus%reset()
  mic_c = micbiom_beg - 1 + c_loc
  dom_c = dom_beg - 1 + c_loc
  if(is_surflit)then
    lay=1
  else
    lay=2
  endif
  call this%calc_cnp_ratios(cdom_bgc_index, ystates, bstatus)
  if (bstatus%check_status())return

  !calculate potential decay coefficient (1/s)
  call this%calc_som_decay_k(lay, cdom_bgc_index, decompkf_eca, ystates(mic_c), ystates(dom_c), k_decay)

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(cdom_bgc_index, dtime, k_decay(1:ncentpools), &
      ystates(1:nom_tot_elms), pot_om_decay_rates)

  do jj = 1, ncentpools
    kc = (jj-1) * nelms + c_loc
    !the following avoids over-estimation of potential hr which is used for nitri-denit estimation
    pot_om_decay_rates(jj) = min(pot_om_decay_rates(jj), ystates(kc)/dtime)
  enddo

  call this%calc_cascade_matrix(lay, cdom_bgc_index, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix)

  !calculate potential respiration rates by summarizing all om decomposition pathways
  call this%calc_potential_aerobic_hr(cdom_bgc_index, pot_om_decay_rates, &
    cascade_matrix, pot_co2_hr)

  end associate
  end subroutine run_decomp
!------------------------------------------

  subroutine calc_cascade_matrix(this, lay, cdom_bgc_index, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix)
  !
  ! DESCRIPTION
  ! calculate cascade matrix for decomposition
  ! in all the reactions, the nominal carbon oxidation status is assumed as zero, which is apparently not correct.
  ! It is also assumed the recycling of nitrogen and phosphorus during decomposition is 100%, which is likely
  ! not quite right as well.
  use cdomBGCIndexType , only : cdom_bgc_index_type
  use MathfuncMod         , only : safe_div, fpmax
  implicit none
  class(cdomSom_type),           intent(inout) :: this
  type(cdom_bgc_index_type)  , intent(in)    :: cdom_bgc_index
  integer                      , intent(in)    :: lay
  real(r8)                     , intent(in)    :: pct_sand
  real(r8)                     , intent(in)    :: pct_clay
  real(r8)                     , intent(out)   :: alpha_n(ncentpools) !indicating factor for nitrogen limitation
  real(r8)                     , intent(out)   :: alpha_p(ncentpools) !indicating factor for phosphorus limitation
  real(r8)                     , intent(inout) :: cascade_matrix(cdom_bgc_index%nstvars, cdom_bgc_index%nreactions)

  integer  :: reac,jj
  real(r8) :: f1, f2, rf_s1

  associate(                                                &
    lmet      => cdom_bgc_index%lmet                       , & !
    lcel      => cdom_bgc_index%lcel                       , & !
    llig      => cdom_bgc_index%llig                       , & !
    mic       => cdom_bgc_index%mic                        , & !
    pom       => cdom_bgc_index%pom                        , & !
    humus     => cdom_bgc_index%humus                      , & !
    dom       => cdom_bgc_index%dom                        , & !
    cwd       => cdom_bgc_index%cwd                        , & !
    lwd       => cdom_bgc_index%lwd                        , & !
    fwd       => cdom_bgc_index%fwd                        , & !
    c_loc     => cdom_bgc_index%c_loc                      , & !
    n_loc     => cdom_bgc_index%n_loc                      , & !
    p_loc     => cdom_bgc_index%p_loc                      , & !
    e_loc     => cdom_bgc_index%e_loc                      , & !
    c13_loc   => cdom_bgc_index%c13_loc                    , & !
    c14_loc   => cdom_bgc_index%c14_loc                    , & !
    nelms     => cdom_bgc_index%nelms                      , & !
    lid_o2    => cdom_bgc_index%lid_o2                     , & !
    lid_co2   => cdom_bgc_index%lid_co2                    , & !
    lid_nh4   => cdom_bgc_index%lid_nh4                    , & !
    lid_c14_co2=> cdom_bgc_index%lid_c14_co2               , & !
    lid_c13_co2=> cdom_bgc_index%lid_c13_co2               , & !
    lid_co2_hr => cdom_bgc_index%lid_co2_hr                , &
    lid_minn_nh4_immob=> cdom_bgc_index%lid_minn_nh4_immob , &
    lid_minp_immob => cdom_bgc_index%lid_minp_immob        , &
    lid_minp_soluble=> cdom_bgc_index%lid_minp_soluble     , &
    lmet_dek_reac => cdom_bgc_index%lmet_dek_reac          , &
    lcel_dek_reac => cdom_bgc_index%lcel_dek_reac          , &
    llig_dek_reac => cdom_bgc_index%llig_dek_reac          , &
    mic_dek_reac => cdom_bgc_index%mic_dek_reac            , &
    pom_dek_reac => cdom_bgc_index%pom_dek_reac            , &
    humus_dek_reac => cdom_bgc_index%humus_dek_reac        , &
    dom_dek_reac => cdom_bgc_index%dom_dek_reac            , &
    cwd_dek_reac => cdom_bgc_index%cwd_dek_reac            , &
    lwd_dek_reac => cdom_bgc_index%lwd_dek_reac            , &
    fwd_dek_reac => cdom_bgc_index%fwd_dek_reac            , &
    cwd_flig     => this%cwd_flig                         , &
    lwd_flig     => this%lwd_flig                         , &
    fwd_flig     => this%fwd_flig                         , &
    rf_l2s1_bgc  => this%rf_l2s1_bgc                      , &
    rf_l3s2_bgc  => this%rf_l3s2_bgc                      , &
    rf_s2s1_bgc  => this%rf_s2s1_bgc                      , &
    rf_s3s1_bgc  => this%rf_s3s1_bgc                      , &
    rf_l1s1_bgc  => this%rf_l1s1_bgc                      , &
    rf_s1s2a_bgc => this%rf_s1s2a_bgc                     , &
    rf_s1s2b_bgc => this%rf_s1s2b_bgc                     , &
    rf_doms1_bgc => this%rf_doms1_bgc                     , &
    debug        => cdom_bgc_index%debug                     &
  )

    alpha_n = 0._r8; alpha_p = 0._r8
    !---------------------------------------------------------------------------------
    !reactions produce dom
    !reaction1: lmet -> dom_c + dom_n + dom_p + rf_l1s1_bgc(lay) dom_e
    reac=lmet_dek_reac
    call som_to_dom(lmet, reac, rf_l1s1_bgc(lay))

    !reaction 2: lcel -> dom_c + dom_n + dom_p + rf_l2s1_bgc(lay) dom_e
    reac=lcel_dek_reac
    call som_to_dom(lcel, reac, rf_l2s1_bgc(lay))

    !reaction 3: humus -> dom_c + dom_n + dom_p + dom_e
    reac = humus_dek_reac
    call som_to_dom(humus, reac, rf_s3s1_bgc)

    !reaction 4: pom -> (1-f1)[dom_c + dom_n + dom_p] + dom_e + f1 humus
    reac = pom_dek_reac
    f1 = 0.003_r8+0.00009_r8*pct_clay
    cascade_matrix((pom-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((pom-1)*nelms+n_loc   ,reac)   = -this%icn_ratios(pom)
    cascade_matrix((pom-1)*nelms+p_loc   ,reac)   = -this%icp_ratios(pom)

    cascade_matrix((humus-1)*nelms+c_loc   ,reac) =  f1
    cascade_matrix((humus-1)*nelms+n_loc   ,reac) =  f1*this%icn_ratios(pom)
    cascade_matrix((humus-1)*nelms+p_loc   ,reac) =  f1*this%icp_ratios(pom)

    !the following may sometimes make the nutrient flux into dom negative.
    cascade_matrix((dom-1)*nelms+c_loc   ,reac)   =   1._r8-f1
    cascade_matrix((dom-1)*nelms+n_loc   ,reac)   =  (1._r8-f1)*this%icn_ratios(pom)
    cascade_matrix((dom-1)*nelms+p_loc   ,reac)   =  (1._r8-f1)*this%icp_ratios(pom)
    cascade_matrix((dom-1)*nelms+e_loc   ,reac)   =  rf_s2s1_bgc*(1._r8-f1)
    !---------------------------------------------------------------------------
    !reaction produces microbes
    !reaction 5: dom -> mic + co2
    reac = dom_dek_reac
    cascade_matrix((dom-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((dom-1)*nelms+e_loc   ,reac)   = -this%rf_doms1_bgc
    cascade_matrix((dom-1)*nelms+n_loc   ,reac)   = -this%icn_ratios(dom)
    cascade_matrix((dom-1)*nelms+p_loc   ,reac)   = -this%icp_ratios(dom)
    cascade_matrix(lid_co2               ,reac)   =  this%rf_doms1_bgc
    cascade_matrix((mic-1)*nelms+c_loc   ,reac)   = 1._r8-cascade_matrix(lid_co2  ,reac)
    cascade_matrix((mic-1)*nelms+n_loc   ,reac)   = cascade_matrix((mic-1)*nelms+c_loc,reac)* this%icn_ratios(mic)
    cascade_matrix((mic-1)*nelms+p_loc   ,reac)   = cascade_matrix((mic-1)*nelms+c_loc,reac)* this%icp_ratios(mic)

    cascade_matrix(lid_o2                 ,reac)  = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac)  = -cascade_matrix((dom-1)*nelms+n_loc   ,reac) - &
        cascade_matrix((mic-1)*nelms+n_loc,reac)
    cascade_matrix(lid_minp_soluble       ,reac)  = -cascade_matrix((dom-1)*nelms+p_loc   ,reac) - &
        cascade_matrix((mic-1)*nelms+p_loc,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)  = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_co2_hr             ,reac)  = cascade_matrix(lid_co2          ,reac)
    cascade_matrix(lid_minp_immob         ,reac)  = -cascade_matrix(lid_minp_soluble,reac)

    if(this%use_c14)then
      cascade_matrix((dom-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(dom)
      cascade_matrix(lid_c14_co2             , reac) = rf_doms1_bgc*this%icc14_ratios(dom)
      cascade_matrix((mic-1)*nelms+c14_loc   , reac) = cascade_matrix((mic-1)*nelms+c_loc,reac)*this%icc14_ratios(dom)
    endif

    if(this%use_c13)then
      cascade_matrix((lmet-1)*nelms+c13_loc  , reac) = -this%icc13_ratios(lmet)
      cascade_matrix(lid_c13_co2             , reac) = rf_doms1_bgc*this%icc13_ratios(dom)
      cascade_matrix((mic-1)*nelms+c13_loc   , reac) = cascade_matrix((mic-1)*nelms+c_loc,reac)*this%icc13_ratios(dom)
    endif

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    !---------------------------------------------------------------------------------
    !reaction 6, llig -> pom
    !because ligin decomposition is energy costly, an amount CO2 is emitted.
    !This makes the cn ratio or lignin floating
    reac = llig_dek_reac

    cascade_matrix((llig-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((llig-1)*nelms+n_loc   ,reac) = -this%icn_ratios(llig)
    cascade_matrix((llig-1)*nelms+p_loc   ,reac) = -this%icp_ratios(llig)

    cascade_matrix((pom-1)*nelms+c_loc   ,reac) =  1._r8-rf_l3s2_bgc
    cascade_matrix((pom-1)*nelms+n_loc   ,reac) =  this%icn_ratios(llig)
    cascade_matrix((pom-1)*nelms+p_loc   ,reac) =  this%icp_ratios(llig)

    cascade_matrix(lid_co2               ,reac) = rf_l3s2_bgc
    cascade_matrix(lid_o2                ,reac) = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_co2_hr            ,reac) = cascade_matrix(lid_co2    ,reac)

    if(this%use_c14)then
      cascade_matrix((llig-1)*nelms+c14_loc  , reac) = -this%icc14_ratios(llig)
      cascade_matrix(lid_c14_co2             , reac) =  rf_l3s2_bgc*this%icc14_ratios(llig)
      cascade_matrix((pom-1)*nelms+c14_loc   , reac) =  cascade_matrix((pom-1)*nelms+c_loc,reac)*this%icc14_ratios(llig)
    endif

    if(this%use_c13)then
      cascade_matrix((llig-1)*nelms+c13_loc  , reac) = -this%icc13_ratios(llig)
      cascade_matrix(lid_c13_co2             , reac) =  rf_l3s2_bgc*this%icc13_ratios(llig)
      cascade_matrix((pom-1)*nelms+c13_loc   , reac) =  cascade_matrix((pom-1)*nelms+c_loc,reac)*this%icc13_ratios(llig)
    endif

    !---------------------------------------------------------------------------------
    !reaction 7: the partition mic into pom and humus is soil texture dependent
    !SOM1 -> f1*SOM2 + f2*SOm3 + rf_s1*CO2 + (1/cn_ratios(SOM1)-f1/cn_ratios(SOM2)-f2/cn_ratios(SOM3))*min_n
    ! +(1/cp_ratios(SOM1)-f1/cp_ratios(SOM2)-f2/cp_ratios(SOM3))*min_p
    reac = mic_dek_reac
    f2=0.003_r8 + 0.00032_r8*pct_clay
    rf_s1 = rf_s1s2a_bgc(lay) + rf_s1s2b_bgc(lay) * pct_sand * 0.01_r8
    f1 = 1._r8 - rf_s1 - f2

    cascade_matrix((mic-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((mic-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(mic)
    cascade_matrix((mic-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(mic)

    cascade_matrix((humus-1)*nelms+c_loc   ,reac)  = f2
    cascade_matrix((humus-1)*nelms+n_loc   ,reac)  = f2*this%icn_ratios(mic)
    cascade_matrix((humus-1)*nelms+p_loc   ,reac)  = f2*this%icp_ratios(mic)

    cascade_matrix((pom-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((pom-1)*nelms+n_loc   ,reac) = this%icn_ratios(mic)-cascade_matrix((humus-1)*nelms+n_loc,reac)
    cascade_matrix((pom-1)*nelms+p_loc   ,reac) = this%icp_ratios(mic)-cascade_matrix((humus-1)*nelms+p_loc,reac)

    cascade_matrix(lid_co2, reac)               = rf_s1
    cascade_matrix(lid_o2                ,reac) = -cascade_matrix(lid_co2, reac)
    cascade_matrix(lid_co2_hr            ,reac) =  cascade_matrix(lid_co2 ,reac)

    if(this%use_c14)then
      cascade_matrix((mic-1)*nelms+c14_loc    , reac) = -this%icc14_ratios(mic)
      cascade_matrix(lid_c14_co2              , reac) = rf_s1*this%icc14_ratios(mic)
      cascade_matrix((pom-1)*nelms+c14_loc    , reac) = f1*this%icc14_ratios(mic)
      cascade_matrix((humus-1)*nelms+c14_loc  , reac) = f2*this%icc14_ratios(mic)
    endif

    if(this%use_c13)then
      cascade_matrix((mic-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(mic)
      cascade_matrix(lid_c13_co2             , reac) = rf_s1*this%icc13_ratios(mic)
      cascade_matrix((pom-1)*nelms+c13_loc   , reac) = f1*this%icc13_ratios(mic)
      cascade_matrix((humus-1)*nelms+c13_loc , reac) = f2*this%icc13_ratios(mic)
    endif

    !---------------------------------------------------------------------------------
    !reaction 8, the partition cwd into mic and pom
    reac = cwd_dek_reac
    !cwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(mic)-f2/cn_ratios(pom))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(mic)-f2/cp_ratios(pom))

    call wood_decomp_cascade(cwd, reac, cwd_flig)

    !---------------------------------------------------------------------------------
    !reaction 9, the partition lwd into cellulose and lignin
    reac = lwd_dek_reac
    !lwd -> cel + lignin

    call wood_decomp_cascade(lwd, reac, lwd_flig)

    !---------------------------------------------------------------------------------
    !reaction 9, the partition fwd into mic and pom
    reac = fwd_dek_reac
    !fwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(mic)-f2/cn_ratios(pom))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(mic)-f2/cp_ratios(pom))
    call wood_decomp_cascade(fwd, reac, fwd_flig)

  end associate
  contains
!------------------------------------------------------------------------------
    subroutine som_to_dom(lsom, reac, rf_som)

    implicit none
    integer , intent(in) :: lsom
    integer , intent(in) :: reac
    real(r8), intent(in) :: rf_som
    associate(                                 &
      c_loc     => cdom_bgc_index%c_loc       , & !
      n_loc     => cdom_bgc_index%n_loc       , & !
      p_loc     => cdom_bgc_index%p_loc       , & !
      e_loc     => cdom_bgc_index%e_loc       , & !
      nelms     => cdom_bgc_index%nelms       , & !
      c13_loc   => cdom_bgc_index%c13_loc     , & !
      c14_loc   => cdom_bgc_index%c14_loc     , & !
      dom       => cdom_bgc_index%dom           &
    )
    cascade_matrix((lsom-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((lsom-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(lsom)
    cascade_matrix((lsom-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(lsom)

    cascade_matrix((dom-1)*nelms+c_loc   ,reac)  = 1._r8
    cascade_matrix((dom-1)*nelms+n_loc   ,reac)  = this%icn_ratios(lsom)
    cascade_matrix((dom-1)*nelms+p_loc   ,reac)  = this%icp_ratios(lsom)
    cascade_matrix((dom-1)*nelms+e_loc   ,reac)  = rf_som

    if(this%use_c13)then
      cascade_matrix((lsom-1)*nelms+c13_loc,reac)= -this%icc13_ratios(lsom)
      cascade_matrix((dom-1)*nelms +c13_loc,reac)= this%icc13_ratios(lsom)
    endif

    if(this%use_c14)then
      cascade_matrix((lsom-1)*nelms+c14_loc,reac)= -this%icc14_ratios(lsom)
      cascade_matrix((dom-1)*nelms +c14_loc,reac)=  this%icc14_ratios(lsom)
    endif
    end associate
    end subroutine som_to_dom
!------------------------------------------------------------------------------
    subroutine wood_decomp_cascade(iwd, reac,  flig)

    implicit none
    integer , intent(in) :: iwd, reac
    real(r8), intent(in)::  flig

    real(r8) :: fcel
    associate(                                        &
      c_loc     => cdom_bgc_index%c_loc              , & !
      n_loc     => cdom_bgc_index%n_loc              , & !
      p_loc     => cdom_bgc_index%p_loc              , & !
      c13_loc   => cdom_bgc_index%c13_loc            , & !
      c14_loc   => cdom_bgc_index%c14_loc            , & !
      llig      => cdom_bgc_index%llig               , & !
      lcel      => cdom_bgc_index%lcel               , & !
      nelms     => cdom_bgc_index%nelms                & !
    )

    fcel = 1._r8-flig
    cascade_matrix((iwd-1)*nelms+c_loc    ,reac) = -1._r8
    cascade_matrix((iwd-1)*nelms+n_loc    ,reac) = -this%icn_ratios(iwd)
    cascade_matrix((iwd-1)*nelms+p_loc    ,reac) = -this%icp_ratios(iwd)

    cascade_matrix((lcel-1)*nelms+c_loc   ,reac) = fcel
    cascade_matrix((lcel-1)*nelms+n_loc   ,reac) = fcel*this%icn_ratios(iwd)
    cascade_matrix((lcel-1)*nelms+p_loc   ,reac) = fcel*this%icp_ratios(iwd)

    cascade_matrix((llig-1)*nelms+c_loc   ,reac) = flig
    cascade_matrix((llig-1)*nelms+n_loc   ,reac) = flig*this%icn_ratios(iwd)
    cascade_matrix((llig-1)*nelms+p_loc   ,reac) = flig*this%icp_ratios(iwd)

    if(this%use_c14)then
      cascade_matrix((iwd-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(iwd)
      cascade_matrix((lcel-1)*nelms+c14_loc  , reac) = fcel*this%icc14_ratios(iwd)
      cascade_matrix((llig-1)*nelms+c14_loc  , reac) = flig*this%icc14_ratios(iwd)
    endif

    if(this%use_c13)then
      cascade_matrix((iwd-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(iwd)
      cascade_matrix((lcel-1)*nelms+c13_loc  , reac) =  fcel*this%icc13_ratios(iwd)
      cascade_matrix((llig-1)*nelms+c13_loc  , reac) =  flig*this%icc13_ratios(iwd)
    endif
    end associate
    end subroutine wood_decomp_cascade
  end subroutine calc_cascade_matrix

  !-----------------------------------------------------------------------
  subroutine calc_potential_aerobic_hr(this, cdom_bgc_index, pot_decay_rates, &
    cascade_matrix, pot_co2_hr)
    !
    ! DESCRIPTION:
    ! calculate potential aerobic heteorotrophic respiration, and potential oxygen consumption based on cascade_matrix
    ! !USES:
    use MathfuncMod         , only : dot_sum
    use MathfuncMod         , only : safe_div
    use cdomBGCIndexType       , only : cdom_bgc_index_type
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(cdomSom_type), intent(inout) :: this
    type(cdom_bgc_index_type)   , intent(in) :: cdom_bgc_index
    real(r8)                , intent(in) :: pot_decay_rates(ncentpools)
    real(r8)                , intent(in) :: cascade_matrix(cdom_bgc_index%nstvars, cdom_bgc_index%nreactions)
    real(r8)                , intent(out):: pot_co2_hr
    ! !LOCAL VARIABLES:

    associate(                                         & !
         lid_co2_hr=> cdom_bgc_index%lid_co2_hr       , & !
         llig      => cdom_bgc_index%llig             , & !
         mic      => cdom_bgc_index%mic               , & !
         dom      => cdom_bgc_index%dom               , & !
         llig_dek_reac=> cdom_bgc_index%llig_dek_reac , & !
         dom_dek_reac => cdom_bgc_index%dom_dek_reac  , &
         mic_dek_reac => cdom_bgc_index%mic_dek_reac    &
    )

    pot_co2_hr = cascade_matrix(lid_co2_hr,dom_dek_reac) * pot_decay_rates(dom)  + &
          cascade_matrix(lid_co2_hr,llig_dek_reac) * pot_decay_rates(llig)       + &
          cascade_matrix(lid_co2_hr,mic_dek_reac) * pot_decay_rates(mic)         !mol CO2/m3/s
    !print*,'pot_co2_hr',pot_co2_hr,pot_decay_rates(dom),pot_decay_rates(llig),pot_decay_rates(mic)
    end associate
  end subroutine calc_potential_aerobic_hr

  !-----------------------------------------------------------------------
  subroutine calc_cnp_ratios(this, cdom_bgc_index, ystates, bstatus)
  !
  ! DESCRIPTION
  ! compute the cnp ratios for the om pools
  use BetrStatusType      , only : betr_status_type
  use MathfuncMod         , only : safe_div
  use cdomBGCIndexType       , only : cdom_bgc_index_type
  implicit none
  class(cdomSom_type), intent(inout) :: this
  type(cdom_bgc_index_type)   , intent(in) :: cdom_bgc_index
  real(r8)                      , intent(inout) :: ystates(cdom_bgc_index%nstvars)
  type(betr_status_type)      , intent(out) :: bstatus
  integer :: jj
  integer :: kc, kn, kp, kc13, kc14, kc1, kc2
  real(r8):: rat
  real(r8) :: difn
  real(r8) :: stoibal_ncon
  character(len=255) :: msg

  associate(                         &
    nelms => cdom_bgc_index%nelms, &
    c_loc => cdom_bgc_index%c_loc, &
    n_loc => cdom_bgc_index%n_loc, &
    p_loc => cdom_bgc_index%p_loc, &
    c13_loc => cdom_bgc_index%c13_loc, &
    c14_loc => cdom_bgc_index%c14_loc, &
    lcel  => cdom_bgc_index%lcel , &
    llig  => cdom_bgc_index%llig , &
    is_ompool_som => cdom_bgc_index%is_ompool_som, &
    is_dom_pool => cdom_bgc_index%is_dom_pool, &
    ompoolnames => cdom_bgc_index%ompoolnames &
  )
  !compute the respiration efficiency for dom
  call calc_dom_rf(this, cdom_bgc_index, ystates)

  !for om pools
  do jj = 1, ncentpools
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
    if(cdom_bgc_index%debug)then
       write(*,'(A,X,I2,5(X,E20.10))')'cnp',jj,ystates(kc),ystates(kn),ystates(kp),1._r8/this%icn_ratios(jj),1._r8/this%icp_ratios(jj)
    endif
    if(is_ompool_som(jj) .and. ystates(kc)>tiny_val .and. .not. is_dom_pool(jj))then
      stoibal_ncon = ystates(kc)*this%icn_ratios(jj)
      difn=ystates(kn)-stoibal_ncon
      if(difn<-tiny_ncon)then
        ystates(kn)=stoibal_ncon
      endif
!      write(msg,*)'phosphorus weirdo',jj,trim(ompoolnames(jj)),ystates(kc),ystates(kn),ystates(kp), rat, this%def_cn(jj),this%def_cp(jj),&
!         1._r8/this%icn_ratios(jj),1._r8/this%icp_ratios(jj)
!      print*,msg
 !     call bstatus%set_msg(msg,err=-1)
 !     return
    endif
    if(this%use_c14)then
      kc14 = (jj-1) * nelms + c14_loc
      this%icc14_ratios(jj) = 1._r8/this%def_cc14(jj)*(1._r8-rat)+ystates(kc14)/ystates(kc)
      if(cdom_bgc_index%debug)then
        write(*,'(A,X,I4,2(X,E20.10))') 'c14rrr som jj',jj,1._r8/this%def_cc14(jj),this%icc14_ratios(jj)
      endif
    endif
    if(this%use_c13)then
      kc13 = (jj-1) * nelms + c13_loc
      kc13 = (jj-1) * nelms + c13_loc
      this%icc13_ratios(jj) = 1._r8/this%def_cc13(jj)*(1._r8-rat)+ystates(kc13)/ystates(kc)
    endif
  enddo
  kc1 = (lcel-1)*nelms+c_loc
  kc2 = (llig-1)*nelms+c_loc
  !lignin fraction of the structural carbon
  this%lit_flig = safe_div(ystates(kc2),ystates(kc1)+ystates(kc2))

  end associate
  end subroutine calc_cnp_ratios

  !-------------------------------------------------------------------------------
  subroutine stoichiometry_fix(this, cdom_bgc_index,ystates)

  !
  ! DESCRIPTION
  ! this fixes the stoichiometric drift due to limite precision of
  ! double precision.
  use cdomBGCIndexType       , only : cdom_bgc_index_type
  implicit none
  class(cdomSom_type)           , intent(inout) :: this
  type(cdom_bgc_index_type)   , intent(in) :: cdom_bgc_index
  real(r8)                      , intent(inout) :: ystates(cdom_bgc_index%nstvars)

  associate(                             &
    c13_loc => cdom_bgc_index%c13_loc, &
    c14_loc => cdom_bgc_index%c14_loc, &
    mic  => cdom_bgc_index%mic , &
    pom  => cdom_bgc_index%pom , &
    humus  => cdom_bgc_index%humus , &
    is_ompool_som => cdom_bgc_index%is_ompool_som, &
    ompoolnames => cdom_bgc_index%ompoolnames &
  )


  !for om pools
  call stoi_fix(mic)

  call stoi_fix(pom)

  call stoi_fix(humus)

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
      nelms => cdom_bgc_index%nelms, &
      c_loc => cdom_bgc_index%c_loc, &
      n_loc => cdom_bgc_index%n_loc  &
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
  subroutine calc_som_decay_r(this, cdom_bgc_index, dtime, om_k_decay, om_pools, om_decay_rates)
    !
    ! !DESCRIPTION:
    ! calculate degradation for all different pools
    !
    ! !USES:
    use cdomBGCIndexType       , only : cdom_bgc_index_type
   implicit none
   class(cdomSom_type)     , intent(inout) :: this
   type(cdom_bgc_index_type) , intent(in)    :: cdom_bgc_index
    real(r8)  , intent(in)    :: dtime
    real(r8)  , intent(in)    :: om_k_decay(ncentpools)
    real(r8)  , intent(in)    :: om_pools(cdom_bgc_index%nom_tot_elms)
    real(r8)  , intent(out)   :: om_decay_rates(ncentpools)

    ! !LOCAL VARIABLES:
    integer :: jj, fc, c, j
    integer :: kc, kn
    associate(                                      &
         nelms => cdom_bgc_index%nelms            , &
         mic   => cdom_bgc_index%mic              , &
         nom_pools => cdom_bgc_index%nom_pools    , &
         c_loc => cdom_bgc_index%c_loc              &
    )

    !for om pools
    do jj = 1, nom_pools
      kc = (jj-1) * nelms + c_loc
      om_decay_rates(jj) = om_pools(kc) * om_k_decay(jj)
    enddo
    end associate
  end subroutine calc_som_decay_r
  !-------------------------------------------------------------------------------

  subroutine apply_spinupf(this, cdom_bgc_index, decompkf_eca, k_decay, spinup_scalar, spinup_flg)
  use cdomBGCIndexType       , only : cdom_bgc_index_type
  use cdomBGCDecompType      , only : Decompcdom_type
  use betr_varcon               , only : kyr_spinup
  implicit none
  class(cdomSom_type)     , intent(inout) :: this
  type(Decompcdom_type), intent(in) :: decompkf_eca
  type(cdom_bgc_index_type)   , intent(in)    :: cdom_bgc_index
  real(r8)                      , intent(inout) :: k_decay(ncentpools)
  real(r8)                      , intent(inout) :: spinup_scalar
  integer                       , intent(in)    :: spinup_flg
  integer :: jj

  associate(   &
   t_scalar       => decompkf_eca%t_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil temperature scalar for decomp
   w_scalar       => decompkf_eca%w_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil water scalar for decomp
   o_scalar       => decompkf_eca%o_scalar        , & ! Intput: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
   mic           => cdom_bgc_index%mic        , & !
   pom           => cdom_bgc_index%pom        , & !
   humus           => cdom_bgc_index%humus          & !
  )

  if(spinup_flg==2)then
    !cumulated more than 2 * kyr_spinup
    do jj = 1, ncentpools
      k_decay(jj) = k_decay(jj)/spinup_scalar
    enddo
  elseif(spinup_flg==1)then
    !cumulated more than kyr_spinup but less than 2 * kyr_spinup
    spinup_scalar = spinup_scalar + t_scalar * w_scalar * o_scalar / (365._r8 * 86400._r8 * kyr_spinup)
  endif

  end associate

  end subroutine apply_spinupf
  !-------------------------------------------------------------------------------
  subroutine calc_som_decay_k(this, lay, cdom_bgc_index, decompkf_eca, micb, domc, k_decay)

  use cdomBGCIndexType       , only : cdom_bgc_index_type
  use cdomBGCDecompType      , only : Decompcdom_type
  implicit none
  class(cdomSom_type)        , intent(inout) :: this
  integer                    , intent(in)  :: lay
  type(Decompcdom_type)      , intent(in)  :: decompkf_eca
  type(cdom_bgc_index_type)   , intent(in)  :: cdom_bgc_index
  real(r8)                   , intent(in)  :: micb
  real(r8)                   , intent(in)  :: domc
  real(r8)                   , intent(out) :: k_decay(ncentpools)
  integer :: jj
  real(r8), parameter :: rc=1.e-6_r8
  real(r8), parameter :: rm=3.e-6_r8
  real(r8) :: wamp, fmic
  associate(   &
   tfnr           => decompkf_eca%tfnr            , & ! Intput: [real(r8) (:,:)   ]  soil temperature scalar for decomp
   tfng           => decompkf_eca%tfng            , & ! Intput: [real(r8) (:,:)   ]  soil water scalar for decomp
   o_scalar       => decompkf_eca%o_scalar        , & ! Intput: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
   Minsurf        => decompkf_eca%Minsurf         , & !
   KM_OM          => decompkf_eca%KM_OM           , & !
   enz_modifier   => decompkf_eca%enz_modifier    , &
   filmt          => decompkf_eca%filmt           , & !
   lmet           => cdom_bgc_index%lmet          , & !
   lcel           => cdom_bgc_index%lcel          , & !
   llig           => cdom_bgc_index%llig          , & !
   mic            => cdom_bgc_index%mic           , & !
   pom            => cdom_bgc_index%pom           , & !
   humus          => cdom_bgc_index%humus         , & !
   cwd            => cdom_bgc_index%cwd           , & !
   lwd            => cdom_bgc_index%lwd           , & !
   fwd            => cdom_bgc_index%fwd           , & !
   dom            => cdom_bgc_index%dom           , & !
   alpha_E        => this%alpha_E                 , & !
   alpha_T        => this%alpha_T                 , & !
   Kaff_CM        => this%Kaff_CM                 , & !
   Kaff_mic_dom   => this%Kaff_mic_dom            , & !
   Kaff_enz       => this%Kaff_enz                  & !
  )

  fmic = alpha_E * micb / (Kaff_enz * enz_modifier+ alpha_E * micb)
  k_decay(lmet) = this%k_decay_lmet(lay) * tfng * o_scalar * fmic
  k_decay(lcel) = this%k_decay_lcel(lay) * tfng * o_scalar * fmic
  k_decay(llig) = this%k_decay_llig(lay) * tfng * o_scalar * fmic
  k_decay(mic)  = this%k_decay_mic(lay) * tfnr * o_scalar
  k_decay(pom)  = this%k_decay_pom * tfng * o_scalar * fmic
  k_decay(humus) = this%k_decay_humus * tfng * o_scalar * fmic
  k_decay(cwd)  = this%k_decay_cwd  * tfng * o_scalar * fmic
  k_decay(lwd)  = this%k_decay_lwd  * tfng * o_scalar * fmic
  k_decay(fwd)  = this%k_decay_fwd  * tfng * o_scalar * fmic

  !impose the ligin effect
  k_decay(cwd)  = k_decay(cwd) * exp(-3._r8*this%cwd_flig)
  k_decay(lwd)  = k_decay(lwd) * exp(-3._r8*this%lwd_flig)
  k_decay(fwd)  = k_decay(fwd) * exp(-3._r8*this%fwd_flig)
  k_decay(lcel) = k_decay(lcel)* exp(-3._r8*this%lit_flig)
  k_decay(llig) = k_decay(llig)* exp(-3._r8*this%lit_flig)

  !compute the decay rate of dom
  k_decay(dom) = this%k_decay_dom * alpha_T * micb / (Kaff_mic_dom * enz_modifier + &
   domc + alpha_T * micb + Minsurf * Kaff_mic_dom/(Kaff_CM*KM_OM)) * tfng * o_scalar

  !print*,'lmet',this%k_decay_lmet(lay)
  !print*,'lcel',this%k_decay_lcel(lay)
  !print*,'llig',this%k_decay_llig(lay)
  !print*,'mic',k_decay(mic)*micb
  !print*,'dom',k_decay(dom)*domc*(1._r8-this%rf_doms1_bgc)
  !print*,'humus',this%k_decay_humus
  !print*,'cwd',this%k_decay_cwd
  !print*,'lwd',this%k_decay_lwd
  !print*,'fwd',this%k_decay_fwd
  !print*,'micb',micb, fmic, o_scalar,tfng
  !print*,'k_decay',(k_decay(jj),jj=1,ncentpools)

  end associate
  end subroutine calc_som_decay_k
  !-------------------------------------------------------------------------------
  subroutine calc_pot_min_np_flx(this, dtime, cdom_bgc_index, ystates, k_decay, cascade_matrix, &
    alpha_n, alpha_p, pot_decomp, pot_nn_flx, pot_np_flx)
  use cdomBGCIndexType       , only : cdom_bgc_index_type
  implicit none
  class(cdomSom_type)         , intent(inout) :: this
  real(r8)                    , intent(in) :: dtime
  type(cdom_bgc_index_type) , intent(in) :: cdom_bgc_index
  real(r8)                    , intent(in) :: ystates(1:cdom_bgc_index%nom_tot_elms)
  real(r8)                    , intent(in) :: k_decay(1:ncentpools)
  real(r8)                    , intent(in) :: cascade_matrix(cdom_bgc_index%nstvars, cdom_bgc_index%nreactions)
  real(r8)                    , intent(in) :: alpha_n(ncentpools)
  real(r8)                    , intent(in) :: alpha_p(ncentpools)
  real(r8)                    , intent(out) :: pot_decomp(ncentpools)
  real(r8)                    , intent(out):: pot_nn_flx
  real(r8)                    , intent(out):: pot_np_flx

  integer :: reac
  integer :: reacs(ncentpools)

  associate(                                                 & !
       nom_pools => cdom_bgc_index%nom_pools                , & !
       nom_tot_elms=> cdom_bgc_index%nom_tot_elms           , & !
       lid_nh4   => cdom_bgc_index%lid_nh4                  , & !
       lid_minp_soluble  => cdom_bgc_index%lid_minp_soluble , & !
       lmet      => cdom_bgc_index%lmet                     , & !
       lcel      => cdom_bgc_index%lcel                     , & !
       llig      => cdom_bgc_index%llig                     , & !
       mic      => cdom_bgc_index%mic                       , & !
       pom      => cdom_bgc_index%pom                       , & !
       humus      => cdom_bgc_index%humus                   , & !
       cwd       => cdom_bgc_index%cwd                      , & !
       dom_dek_reac => cdom_bgc_index%dom_dek_reac          , &
       lmet_dek_reac=> cdom_bgc_index%lmet_dek_reac         , & !
       lcel_dek_reac=> cdom_bgc_index%lcel_dek_reac         , & !
       llig_dek_reac=> cdom_bgc_index%llig_dek_reac         , & !
       mic_dek_reac=> cdom_bgc_index%mic_dek_reac           , & !
       pom_dek_reac=> cdom_bgc_index%pom_dek_reac           , & !
       humus_dek_reac=> cdom_bgc_index%humus_dek_reac       , & !
       cwd_dek_reac=> cdom_bgc_index%cwd_dek_reac           , & !
       lwd_dek_reac=> cdom_bgc_index%lwd_dek_reac           , & !
       fwd_dek_reac=> cdom_bgc_index%fwd_dek_reac             & !
   )

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(cdom_bgc_index, dtime, k_decay(1:nom_pools), &
      ystates(1:nom_tot_elms), pot_decomp)

  pot_nn_flx = 0._r8; pot_np_flx = 0._r8

  !only dom decay could be nutrient limited
  reac = dom_dek_reac
  if(alpha_n(reac)>0._r8)then
      pot_nn_flx = pot_nn_flx - cascade_matrix(lid_nh4, reac) * pot_decomp(reac)
  endif
  if(alpha_p(reac)>0._r8)then
    pot_np_flx = pot_np_flx - cascade_matrix(lid_minp_soluble, reac) * pot_decomp(reac)
  endif

  end associate
  end subroutine calc_pot_min_np_flx

end module cdomBGCSOMType
