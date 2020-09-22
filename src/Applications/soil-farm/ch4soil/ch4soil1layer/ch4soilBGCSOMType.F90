module ch4soilBGCSOMType
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
  !met, cell, lig, cwd, som1, som2, som3
  !during decomposition, the cnp ratios of som1, som2, and som3 are all fixed
  !the litter pools have their cnp ratios varying with nutrient status
  !We consider som1 and som2 as DOM, som3 as humus pool
  !
  integer :: ncentpools
  type, public :: CentSom_type
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
    real(r8) :: k_decay_lit1(2)
    real(r8) :: k_decay_lit2(2)
    real(r8) :: k_decay_lit3(2)
    real(r8) :: k_decay_som1(2)
    real(r8) :: k_decay_som2
    real(r8) :: k_decay_som3
    real(r8) :: k_decay_cwd  !coarse root
    real(r8) :: k_decay_lwd  !large wood
    real(r8) :: k_decay_fwd  !fine branch wood

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

  end type CentSom_type
contains

  subroutine Init(this, ch4soil_bgc_index, biogeo_con, bstatus)

  use ch4soilBGCIndexType , only : ch4soil_bgc_index_type
  use ch4soilParaType       , only : ch4soil_para_type
  use BetrStatusType      , only : betr_status_type
  implicit none
  class(CentSom_type)         , intent(inout) :: this
  type(ch4soil_bgc_index_type) , intent(in)    :: ch4soil_bgc_index
  type(ch4soil_para_type)        , intent(in)    :: biogeo_con
  type(betr_status_type)      , intent(out)   :: bstatus

  call bstatus%reset()
  ncentpools = ch4soil_bgc_index%nom_pools

  call this%InitAllocate()

  end subroutine Init
!------------------------------------------
  subroutine InitAllocate (this)
  !
  !
  implicit none
  class(CentSom_type), intent(inout) :: this

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
  subroutine UpdateParas(this,ch4soil_bgc_index,  biogeo_con)
  !
  ! intialize model parameters
  use ch4soilParaType , only : ch4soil_para_type
  use ch4soilBGCIndexType , only : ch4soil_bgc_index_type
  use tracer_varcon    , only : catomw, natomw, patomw
  implicit none
  class(CentSom_type)  , intent(inout) :: this
  type(ch4soil_bgc_index_type) , intent(in) :: ch4soil_bgc_index
  type(ch4soil_para_type) , intent(in)    :: biogeo_con

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

  this%k_decay_lit1   =  biogeo_con%k_decay_lit1
  this%k_decay_lit2   =  biogeo_con%k_decay_lit2
  this%k_decay_lit3   =  biogeo_con%k_decay_lit3
  this%k_decay_som1   =  biogeo_con%k_decay_som1
  this%k_decay_som2   =  biogeo_con%k_decay_som2
  this%k_decay_som3   =  biogeo_con%k_decay_som3
  this%k_decay_cwd    =  biogeo_con%k_decay_cwd
  this%k_decay_lwd    =  biogeo_con%k_decay_lwd
  this%k_decay_fwd    =  biogeo_con%k_decay_fwd

  this%def_cn(ch4soil_bgc_index%lit1) = biogeo_con%init_cn_met * natomw/catomw;
  this%def_cn(ch4soil_bgc_index%lit2) = biogeo_con%init_cn_cel * natomw/catomw
  this%def_cn(ch4soil_bgc_index%lit3) = biogeo_con%init_cn_lig * natomw/catomw
  this%def_cn(ch4soil_bgc_index%cwd)  = biogeo_con%init_cn_cwd * natomw/catomw
  this%def_cn(ch4soil_bgc_index%lwd)  = biogeo_con%init_cn_lwd * natomw/catomw
  this%def_cn(ch4soil_bgc_index%fwd)  = biogeo_con%init_cn_fwd * natomw/catomw

  this%def_cn(ch4soil_bgc_index%som1) = biogeo_con%init_cn_som1 * natomw/catomw
  this%def_cn(ch4soil_bgc_index%som2) = biogeo_con%init_cn_som2 * natomw/catomw
  this%def_cn(ch4soil_bgc_index%som3) = biogeo_con%init_cn_som3 * natomw/catomw

  this%def_cp(ch4soil_bgc_index%lit1) = biogeo_con%init_cp_met * patomw/catomw
  this%def_cp(ch4soil_bgc_index%lit2) = biogeo_con%init_cp_cel * patomw/catomw
  this%def_cp(ch4soil_bgc_index%lit3) = biogeo_con%init_cp_lig * patomw/catomw
  this%def_cp(ch4soil_bgc_index%cwd)  = biogeo_con%init_cp_cwd * patomw/catomw
  this%def_cp(ch4soil_bgc_index%lwd)  = biogeo_con%init_cp_lwd * patomw/catomw
  this%def_cp(ch4soil_bgc_index%fwd)  = biogeo_con%init_cp_fwd * patomw/catomw

  this%def_cp(ch4soil_bgc_index%som1) = biogeo_con%init_cp_som1 * patomw/catomw
  this%def_cp(ch4soil_bgc_index%som2) = biogeo_con%init_cp_som2 * patomw/catomw
  this%def_cp(ch4soil_bgc_index%som3) = biogeo_con%init_cp_som3 * patomw/catomw

  this%use_c13=biogeo_con%use_c13
  this%use_c14=biogeo_con%use_c14

  if(this%use_c13)then
    this%def_cc13(ch4soil_bgc_index%lit1) = biogeo_con%init_cc13_met
    this%def_cc13(ch4soil_bgc_index%lit2) = biogeo_con%init_cc13_cel
    this%def_cc13(ch4soil_bgc_index%lit3) = biogeo_con%init_cc13_lig
    this%def_cc13(ch4soil_bgc_index%cwd)  = biogeo_con%init_cc13_cwd
    this%def_cc13(ch4soil_bgc_index%lwd)  = biogeo_con%init_cc13_lwd
    this%def_cc13(ch4soil_bgc_index%fwd)  = biogeo_con%init_cc13_fwd
    this%def_cc13(ch4soil_bgc_index%som1) = biogeo_con%init_cc13_som1
    this%def_cc13(ch4soil_bgc_index%som2) = biogeo_con%init_cc13_som2
    this%def_cc13(ch4soil_bgc_index%som3) = biogeo_con%init_cc13_som3
  endif

  if(this%use_c14)then
    this%def_cc14(ch4soil_bgc_index%lit1) = biogeo_con%init_cc14_met
    this%def_cc14(ch4soil_bgc_index%lit2) = biogeo_con%init_cc14_cel
    this%def_cc14(ch4soil_bgc_index%lit3) = biogeo_con%init_cc14_lig
    this%def_cc14(ch4soil_bgc_index%cwd)  = biogeo_con%init_cc14_cwd
    this%def_cc14(ch4soil_bgc_index%lwd)  = biogeo_con%init_cc14_lwd
    this%def_cc14(ch4soil_bgc_index%fwd)  = biogeo_con%init_cc14_fwd
    this%def_cc14(ch4soil_bgc_index%som1) = biogeo_con%init_cc14_som1
    this%def_cc14(ch4soil_bgc_index%som2) = biogeo_con%init_cc14_som2
    this%def_cc14(ch4soil_bgc_index%som3) = biogeo_con%init_cc14_som3
  endif

  end subroutine UpdateParas
!------------------------------------------

  subroutine run_decomp(this, is_surflit, ch4soil_bgc_index, dtime, ystates,&
      decompkf_eca, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix, &
      k_decay, pot_co2_hr, bstatus)
  !
  !DESCRIPTION
  !
  use ch4soilBGCIndexType , only : ch4soil_bgc_index_type
  use ch4soilBGCDecompType, only : DecompCent_type
  use BetrStatusType      , only : betr_status_type
  use betr_ctrl           , only : betr_spinup_state
  implicit none
  class(CentSom_type)         , intent(inout) :: this
  type(ch4soil_bgc_index_type) , intent(in) :: ch4soil_bgc_index
  real(r8)                    , intent(in) :: dtime
  real(r8)                    , intent(inout) :: ystates(1:ch4soil_bgc_index%nom_tot_elms)
  type(DecompCent_type)       , intent(in) :: decompkf_eca
  logical                     , intent(in) :: is_surflit
  real(r8)                    , intent(in) :: pct_sand
  real(r8)                    , intent(in) :: pct_clay
  real(r8)                    , intent(inout) :: cascade_matrix(ch4soil_bgc_index%nstvars, ch4soil_bgc_index%nreactions)
  real(r8)                    , intent(out) :: k_decay(1:ncentpools)
  real(r8)                    , intent(out) :: pot_co2_hr
  real(r8)                    , intent(out) :: alpha_n(1:ncentpools)
  real(r8)                    , intent(out) :: alpha_p(1:ncentpools)
  type(betr_status_type)      , intent(out) :: bstatus

  !local variables
  real(r8) :: pot_om_decay_rates(1:ncentpools)
  integer :: kc, jj, lay

  associate(                                      &
    nelms => ch4soil_bgc_index%nelms,              &
    nom_tot_elms=> ch4soil_bgc_index%nom_tot_elms, &
    c_loc => ch4soil_bgc_index%c_loc               &
  )
  call bstatus%reset()
  if(is_surflit)then
    lay=1
  else
    lay=2
  endif
  call this%calc_cnp_ratios(ch4soil_bgc_index, ystates, bstatus)
  if (bstatus%check_status())return
  call this%calc_som_decay_k(lay, ch4soil_bgc_index, decompkf_eca, k_decay)

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(ch4soil_bgc_index, dtime, k_decay(1:ncentpools), &
      ystates(1:nom_tot_elms), pot_om_decay_rates)

  do jj = 1, ncentpools
    kc = (jj-1) * nelms + c_loc
    !the following avoids over-estimation of potential hr which is used for nitri-denit estimation
    pot_om_decay_rates(jj) = min(pot_om_decay_rates(jj), ystates(kc)/dtime)
  enddo

  call this%calc_cascade_matrix(lay, ch4soil_bgc_index, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix)

  !calculate potential respiration rates by summarizing all om decomposition pathways
  call this%calc_potential_aerobic_hr(ch4soil_bgc_index, pot_om_decay_rates, &
    cascade_matrix, pot_co2_hr, bstatus)
  end associate
  end subroutine run_decomp
!------------------------------------------

  subroutine calc_cascade_matrix(this, lay, ch4soil_bgc_index, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix)
  !
  ! DESCRIPTION
  ! calculate cascade matrix for decomposition
  ! in all the reactions, the nominal carbon oxidation status is assumed as zero, which is apparently not correct.
  ! It is also assumed the recycling of nitrogen and phosphorus during decomposition is 100%, which is likely
  ! not quite right as well.
  use ch4soilBGCIndexType , only : ch4soil_bgc_index_type
  use MathfuncMod         , only : safe_div, fpmax
  implicit none
  class(CentSom_type),           intent(inout) :: this
  type(ch4soil_bgc_index_type)  , intent(in)    :: ch4soil_bgc_index
  integer                      , intent(in)    :: lay
  real(r8)                     , intent(in)    :: pct_sand
  real(r8)                     , intent(in)    :: pct_clay
  real(r8)                     , intent(out)   :: alpha_n(ncentpools) !indicating factor for nitrogen limitation
  real(r8)                     , intent(out)   :: alpha_p(ncentpools) !indicating factor for phosphorus limitation
  real(r8)                     , intent(inout) :: cascade_matrix(ch4soil_bgc_index%nstvars, ch4soil_bgc_index%nreactions)

  integer  :: reac,jj
  real(r8) :: f1, f2, rf_s1

  associate(                                                   &
    lit1      => ch4soil_bgc_index%lit1                       , & !
    lit2      => ch4soil_bgc_index%lit2                       , & !
    lit3      => ch4soil_bgc_index%lit3                       , & !
    som1      => ch4soil_bgc_index%som1                       , & !
    som2      => ch4soil_bgc_index%som2                       , & !
    som3      => ch4soil_bgc_index%som3                       , & !
    cwd       => ch4soil_bgc_index%cwd                        , & !
    lwd       => ch4soil_bgc_index%lwd                        , & !
    fwd       => ch4soil_bgc_index%fwd                        , & !
    c_loc     => ch4soil_bgc_index%c_loc                      , & !
    n_loc     => ch4soil_bgc_index%n_loc                      , & !
    p_loc     => ch4soil_bgc_index%p_loc                      , & !
    c13_loc   => ch4soil_bgc_index%c13_loc                    , & !
    c14_loc   => ch4soil_bgc_index%c14_loc                    , & !
    nelms     => ch4soil_bgc_index%nelms                      , & !
    lid_o2    => ch4soil_bgc_index%lid_o2                     , & !
    lid_co2   => ch4soil_bgc_index%lid_co2                    , & !
    lid_nh4   => ch4soil_bgc_index%lid_nh4                    , & !
    lid_c14_co2=> ch4soil_bgc_index%lid_c14_co2               , & !
    lid_c13_co2=> ch4soil_bgc_index%lid_c13_co2               , & !
    lid_co2_hr => ch4soil_bgc_index%lid_co2_hr                , &
    lid_minn_nh4_immob=> ch4soil_bgc_index%lid_minn_nh4_immob , &
    lid_minp_immob => ch4soil_bgc_index%lid_minp_immob        , &
    lid_minp_soluble=> ch4soil_bgc_index%lid_minp_soluble     , &
    lit1_dek_reac => ch4soil_bgc_index%lit1_dek_reac          , &
    lit2_dek_reac => ch4soil_bgc_index%lit2_dek_reac          , &
    lit3_dek_reac => ch4soil_bgc_index%lit3_dek_reac          , &
    som1_dek_reac => ch4soil_bgc_index%som1_dek_reac          , &
    som2_dek_reac => ch4soil_bgc_index%som2_dek_reac          , &
    som3_dek_reac => ch4soil_bgc_index%som3_dek_reac          , &
    cwd_dek_reac => ch4soil_bgc_index%cwd_dek_reac            , &
    lwd_dek_reac => ch4soil_bgc_index%lwd_dek_reac            , &
    fwd_dek_reac => ch4soil_bgc_index%fwd_dek_reac            , &
    cwd_fcel     => this%cwd_fcel                            , &
    cwd_flig     => this%cwd_flig                            , &
    lwd_fcel     => this%lwd_fcel                            , &
    lwd_flig     => this%lwd_flig                            , &
    fwd_fcel     => this%fwd_fcel                            , &
    fwd_flig     => this%fwd_flig                            , &
    rf_l2s1_bgc  => this%rf_l2s1_bgc                         , &
    rf_l3s2_bgc  => this%rf_l3s2_bgc                         , &
    rf_s2s1_bgc  => this%rf_s2s1_bgc                         , &
    rf_s3s1_bgc  => this%rf_s3s1_bgc                         , &
    rf_l1s1_bgc  => this%rf_l1s1_bgc                         , &
    rf_s1s2a_bgc => this%rf_s1s2a_bgc                        , &
    rf_s1s2b_bgc => this%rf_s1s2b_bgc                        , &
    debug        => ch4soil_bgc_index%debug                     &
  )

    alpha_n = 0._r8; alpha_p = 0._r8
    !---------------------------------------------------------------------------------
    !reaction1, lit1 -> s1
    reac=lit1_dek_reac
    !lit1 + 0.55*o2 -> 0.45 som1 + 0.55co2 + (1/cn_ratios(lit1) - 0.45/cn_ratios(som1))min_n+ (1/cp_ratios(lit1)-0.45/cp_ratios(som1))min_p
    cascade_matrix((lit1-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((lit1-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(lit1)
    cascade_matrix((lit1-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(lit1)

    cascade_matrix((som1-1)*nelms+c_loc   ,reac)  = 1._r8-rf_l1s1_bgc(lay)
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)  = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icn_ratios(som1)
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)  = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icp_ratios(som1)

    cascade_matrix(lid_co2                ,reac)  = rf_l1s1_bgc(lay)

    cascade_matrix(lid_o2                 ,reac)  = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac)  = -cascade_matrix((lit1-1)*nelms+n_loc   ,reac) - &
        cascade_matrix((som1-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble       ,reac)  = -cascade_matrix((lit1-1)*nelms+p_loc   ,reac) - &
        cascade_matrix((som1-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)  = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_co2_hr             ,reac)  = cascade_matrix(lid_co2           ,reac)
    cascade_matrix(lid_minp_immob         ,reac)  = -cascade_matrix(lid_minp_soluble  ,reac)

    if(this%use_c14)then
      cascade_matrix((lit1-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(lit1)
      cascade_matrix(lid_c14_co2              , reac) = rf_l1s1_bgc(lay)*this%icc14_ratios(lit1)
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc14_ratios(lit1)
    endif

    if(this%use_c13)then
      cascade_matrix((lit1-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(lit1)
      cascade_matrix(lid_c13_co2              , reac) = rf_l1s1_bgc(lay)*this%icc13_ratios(lit1)
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc13_ratios(lit1)
    endif

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    !---------------------------------------------------------------------------------
    !reaction 2, lit2 -> s1
    reac = lit2_dek_reac
    !lit2 + 0.5 o2  -> 0.5 som1 + 0.5 co2 + (1/cn_ratios(lit2)-0.5/cn_ratios(som1))min_n +(1/cp_ratios(lit2)-0.5/cp_ratios(som1))min_p
    cascade_matrix((lit2-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((lit2-1)*nelms+n_loc   ,reac)   = -this%icn_ratios(lit2)
    cascade_matrix((lit2-1)*nelms+p_loc   ,reac)   = -this%icp_ratios(lit2)

    cascade_matrix((som1-1)*nelms+c_loc   ,reac)   =  1._r8-rf_l2s1_bgc(lay)
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)   =  cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icn_ratios(som1)
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)   =  cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icp_ratios(som1)

    cascade_matrix(lid_co2                ,reac)   =  rf_l2s1_bgc(lay)
    cascade_matrix(lid_o2                 ,reac)   = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac)   = -cascade_matrix((lit2-1)*nelms+n_loc   ,reac) - &
                                                      cascade_matrix((som1-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((lit2-1)*nelms+p_loc   ,reac) - &
                                                       cascade_matrix((som1-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)   = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac)   = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac)   = cascade_matrix(lid_co2        ,reac)


    if(cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if(cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((lit2-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(lit2)
      cascade_matrix(lid_c14_co2              , reac) = rf_l2s1_bgc(lay)*this%icc14_ratios(lit2)
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc14_ratios(lit2)
    endif
    if(this%use_c13)then
      cascade_matrix((lit2-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(lit2)
      cascade_matrix(lid_c13_co2              , reac) =  rf_l2s1_bgc(lay)*this%icc13_ratios(lit2)
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) =  cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc13_ratios(lit2)
    endif

    !---------------------------------------------------------------------------------
    !reaction 3, lit3->s2
    reac = lit3_dek_reac
    !lit3 + 0.5 o2 -> 0.5 som2 + 0.5 co2 + (1/cn_ratios(lit3) - 0.5/cn_ratios(som2))min_n + (1/cp_ratios(lit3)-0.5_r8/cp_ratios(som2))minp
    cascade_matrix((lit3-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((lit3-1)*nelms+n_loc   ,reac) = -this%icn_ratios(lit3)
    cascade_matrix((lit3-1)*nelms+p_loc   ,reac) = -this%icp_ratios(lit3)

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) =  1._r8-rf_l3s2_bgc
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) =  cascade_matrix((som2-1)*nelms+c_loc,reac)*this%icn_ratios(som2)
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) =  cascade_matrix((som2-1)*nelms+c_loc,reac)*this%icp_ratios(som2)

    cascade_matrix(lid_co2                ,reac) = rf_l3s2_bgc
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((lit3-1)*nelms+n_loc   ,reac) - &
                                                    cascade_matrix((som2-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((lit3-1)*nelms+p_loc   ,reac) - &
                                                   cascade_matrix((som2-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((lit3-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(lit3)
      cascade_matrix(lid_c14_co2              , reac) =  rf_l3s2_bgc*this%icc14_ratios(lit3)
      cascade_matrix((som2-1)*nelms+c14_loc   , reac) =  cascade_matrix((som2-1)*nelms+c_loc,reac)*this%icc14_ratios(lit3)
    endif

    if(this%use_c13)then
      cascade_matrix((lit3-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(lit3)
      cascade_matrix(lid_c13_co2              , reac) =  rf_l3s2_bgc*this%icc13_ratios(lit3)
      cascade_matrix((som2-1)*nelms+c13_loc   , reac) =  cascade_matrix((som2-1)*nelms+c_loc,reac)*this%icc13_ratios(lit3)
    endif

    !---------------------------------------------------------------------------------
    !double check those stoichiometry parameters
    !reaction 4, the partition into som2 and som3 is soil texture dependent
    !SOM1 -> f1*SOM2 + f2*SOm3 + rf_s1*CO2 + (1/cn_ratios(SOM1)-f1/cn_ratios(SOM2)-f2/cn_ratios(SOM3))*min_n
    ! +(1/cp_ratios(SOM1)-f1/cp_ratios(SOM2)-f2/cp_ratios(SOM3))*min_p
    reac = som1_dek_reac
    f2=0.003_r8 + 0.00032_r8*pct_clay
    rf_s1 = rf_s1s2a_bgc(lay) + rf_s1s2b_bgc(lay) * pct_sand * 0.01_r8

    f1 = 1._r8 - rf_s1 - f2

    cascade_matrix((som1-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)  = -this%icn_ratios(som1)
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)  = -this%icp_ratios(som1)

    cascade_matrix((som3-1)*nelms+c_loc   ,reac)  = f2
    cascade_matrix((som3-1)*nelms+n_loc   ,reac)  = f2*this%icn_ratios(som3)
    cascade_matrix((som3-1)*nelms+p_loc   ,reac)  = f2*this%icp_ratios(som3)

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) = f1*this%icn_ratios(som2)
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) = f1*this%icp_ratios(som2)

    cascade_matrix(lid_co2, reac)                = rf_s1

    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2, reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((som1-1)*nelms+n_loc   ,reac)  &
                                                   -cascade_matrix((som2-1)*nelms+n_loc   ,reac)  &
                                                   -cascade_matrix((som3-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((som1-1)*nelms+p_loc   ,reac)  &
                                                   -cascade_matrix((som2-1)*nelms+p_loc   ,reac)  &
                                                   -cascade_matrix((som3-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)


    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(som1)
      cascade_matrix(lid_c14_co2              , reac) =  rf_s1*this%icc14_ratios(som1)
      cascade_matrix((som2-1)*nelms+c14_loc   , reac) =  f1*this%icc14_ratios(som1)
      cascade_matrix((som3-1)*nelms+c14_loc   , reac) =  f2*this%icc14_ratios(som1)
    endif

    if(this%use_c13)then
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(som1)
      cascade_matrix(lid_c13_co2              , reac) = rf_s1*this%icc13_ratios(som1)
      cascade_matrix((som2-1)*nelms+c13_loc   , reac) = f1*this%icc13_ratios(som1)
      cascade_matrix((som3-1)*nelms+c13_loc   , reac) = f2*this%icc13_ratios(som1)
    endif

    !---------------------------------------------------------------------------------
    !reaction 5, som2->som1, som3
    reac = som2_dek_reac
    !som2 + 0.55 o2 -> (0.45-f1) som1 + f1*som3 + 0.55co2 + (1/cn_ratios(som2)-0.42/cn_ratios(som1)-0.03/cn_ratios(som3)) + (1/cp_raitos(som2)-0.42/cp_ratios(som1)-0.03/cp_ratios(som3))
    f1 = 0.003_r8+0.00009_r8*pct_clay
    cascade_matrix((som2-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((som2-1)*nelms+n_loc   ,reac)   = -this%icn_ratios(som2)
    cascade_matrix((som2-1)*nelms+p_loc   ,reac)   = -this%icp_ratios(som2)


    cascade_matrix((som1-1)*nelms+c_loc   ,reac)   =  1._r8-rf_s2s1_bgc-f1
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)   =  cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icn_ratios(som1)
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)   =  cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icp_ratios(som1)

    cascade_matrix((som3-1)*nelms+c_loc   ,reac)   =  f1
    cascade_matrix((som3-1)*nelms+n_loc   ,reac)   =  f1*this%icn_ratios(som3)
    cascade_matrix((som3-1)*nelms+p_loc   ,reac)   =  f1*this%icp_ratios(som3)

    cascade_matrix(lid_co2                ,reac)   = rf_s2s1_bgc
    cascade_matrix(lid_o2                 ,reac)   = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac)   = -cascade_matrix((som2-1)*nelms+n_loc   ,reac) &
                                                     -cascade_matrix((som1-1)*nelms+n_loc   ,reac) &
                                                     -cascade_matrix((som3-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((som2-1)*nelms+p_loc   ,reac)  &
                                                     -cascade_matrix((som1-1)*nelms+p_loc   ,reac)  &
                                                     -cascade_matrix((som3-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)   = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac)   = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac)   = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((som2-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(som2)
      cascade_matrix(lid_c14_co2              , reac) = rf_s2s1_bgc*this%icc14_ratios(som2)
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc14_ratios(som2)
      cascade_matrix((som3-1)*nelms+c14_loc   , reac) = f1*this%icc14_ratios(som2)
    endif

    if(this%use_c13)then
      cascade_matrix((som2-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(som2)
      cascade_matrix(lid_c13_co2              , reac) = rf_s2s1_bgc*this%icc13_ratios(som2)
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc13_ratios(som2)
      cascade_matrix((som3-1)*nelms+c13_loc   , reac) = f1*this%icc13_ratios(som2)
    endif

    !---------------------------------------------------------------------------------
    !reaction 6, s3-> s1
    reac = som3_dek_reac
    !som3 + 0.55 o2 -> 0.45*som1 + 0.55co2 + (1/cn_ratios(som3)-0.45/cn_ratios(som1)) + (1/cp_ratios(som3)-0.45/cp_ratios(som1))
    cascade_matrix((som3-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((som3-1)*nelms+n_loc   ,reac) = -this%icn_ratios(som3)
    cascade_matrix((som3-1)*nelms+p_loc   ,reac) = -this%icp_ratios(som3)

    cascade_matrix((som1-1)*nelms+c_loc   ,reac) = 1._r8-rf_s3s1_bgc
    cascade_matrix((som1-1)*nelms+n_loc   ,reac) = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icn_ratios(som1)
    cascade_matrix((som1-1)*nelms+p_loc   ,reac) = cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icp_ratios(som1)

    cascade_matrix(lid_co2                ,reac) = rf_s3s1_bgc
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((som3-1)*nelms+n_loc   ,reac)  &
                                                   -cascade_matrix((som1-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((som3-1)*nelms+p_loc   ,reac)  &
                                                   -cascade_matrix((som1-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)


    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((som3-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(som3)
      cascade_matrix(lid_c14_co2              , reac) =  rf_s3s1_bgc*this%icc14_ratios(som3)
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) =  cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc14_ratios(som3)
    endif

    if(this%use_c13)then
      cascade_matrix((som3-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(som3)
      cascade_matrix(lid_c13_co2              , reac) =  rf_s3s1_bgc*this%icc13_ratios(som3)
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) =  cascade_matrix((som1-1)*nelms+c_loc,reac)*this%icc13_ratios(som3)
    endif

    !---------------------------------------------------------------------------------
    !reaction 7, the partition cwd into som1 and som2
    reac = cwd_dek_reac
    !cwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = cwd_fcel*(1._r8-rf_l2s1_bgc(lay))
    f2 = (1._r8-cwd_fcel)*(1._r8-rf_l3s2_bgc)

    call wood_decomp_cascade(cwd, reac, f1, f2)

    !---------------------------------------------------------------------------------
    !reaction 8, the partition lwd into som1 and som2
    reac = lwd_dek_reac
    !lwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = lwd_fcel*(1-rf_l2s1_bgc(lay))
    f2 = (1._r8-lwd_fcel)*(1-rf_l3s2_bgc)

    call wood_decomp_cascade(lwd, reac, f1, f2)

    !---------------------------------------------------------------------------------
    !reaction 9, the partition fwd into som1 and som2
    reac = fwd_dek_reac
    !fwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = fwd_fcel*(1-rf_l2s1_bgc(lay))
    f2 = (1._r8-fwd_fcel)*(1-rf_l3s2_bgc)

    call wood_decomp_cascade(fwd, reac, f1, f2)

  end associate
  contains

    subroutine wood_decomp_cascade(iwd, reac, f1, f2)

    implicit none
    integer , intent(in) :: iwd, reac
    real(r8), intent(in):: f1, f2
    associate(                                                   &
      c_loc     => ch4soil_bgc_index%c_loc                      , & !
      n_loc     => ch4soil_bgc_index%n_loc                      , & !
      p_loc     => ch4soil_bgc_index%p_loc                      , & !
      c13_loc   => ch4soil_bgc_index%c13_loc                    , & !
      c14_loc   => ch4soil_bgc_index%c14_loc                    , & !
      nelms     => ch4soil_bgc_index%nelms                      , & !
      lid_o2    => ch4soil_bgc_index%lid_o2                     , & !
      lid_co2   => ch4soil_bgc_index%lid_co2                    , & !
      lid_nh4   => ch4soil_bgc_index%lid_nh4                    , & !
      lid_c14_co2=> ch4soil_bgc_index%lid_c14_co2               , & !
      lid_c13_co2=> ch4soil_bgc_index%lid_c13_co2               , & !
      lid_co2_hr => ch4soil_bgc_index%lid_co2_hr                , &
      lid_minn_nh4_immob=> ch4soil_bgc_index%lid_minn_nh4_immob , &
      lid_minp_immob => ch4soil_bgc_index%lid_minp_immob        , &
      lid_minp_soluble=> ch4soil_bgc_index%lid_minp_soluble     , &
      som1      => ch4soil_bgc_index%som1                       , & !
      som2      => ch4soil_bgc_index%som2                         & !
    )
    cascade_matrix((iwd-1)*nelms+c_loc    ,reac) = -1._r8
    cascade_matrix((iwd-1)*nelms+n_loc    ,reac) = -this%icn_ratios(iwd)
    cascade_matrix((iwd-1)*nelms+p_loc    ,reac) = -this%icp_ratios(iwd)

    cascade_matrix((som1-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((som1-1)*nelms+n_loc   ,reac) = f1*this%icn_ratios(som1)
    cascade_matrix((som1-1)*nelms+p_loc   ,reac) = f1*this%icp_ratios(som1)

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) = f2
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) = f2*this%icn_ratios(som2)
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) = f2*this%icp_ratios(som2)

    cascade_matrix(lid_co2                ,reac) = 1._r8-f1-f2

    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((iwd-1)*nelms+n_loc    ,reac)  &
                                                   -cascade_matrix((som1-1)*nelms+n_loc   ,reac)  &
                                                   -cascade_matrix((som2-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((iwd-1)*nelms+p_loc    ,reac)  &
                                                   -cascade_matrix((som1-1)*nelms+p_loc   ,reac)  &
                                                   -cascade_matrix((som2-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((iwd-1)*nelms+c14_loc   , reac) = -this%icc14_ratios(iwd)
      cascade_matrix((som1-1)*nelms+c14_loc  , reac) =  f1*this%icc14_ratios(iwd)
      cascade_matrix((som2-1)*nelms+c14_loc  , reac) =  f2*this%icc14_ratios(iwd)
    endif

    if(this%use_c13)then
      cascade_matrix((iwd-1)*nelms+c13_loc   , reac) = -this%icc13_ratios(iwd)
      cascade_matrix((som1-1)*nelms+c13_loc  , reac) =  f1*this%icc13_ratios(iwd)
      cascade_matrix((som2-1)*nelms+c13_loc  , reac) =  f2*this%icc13_ratios(iwd)
    endif
    end associate
    end subroutine wood_decomp_cascade
  end subroutine calc_cascade_matrix

  !-----------------------------------------------------------------------
  subroutine calc_potential_aerobic_hr(this, ch4soil_bgc_index, pot_decay_rates, &
    cascade_matrix, pot_co2_hr, bstatus)
    !
    ! DESCRIPTION:
    ! calculate potential aerobic heteorotrophic respiration, and potential oxygen consumption based on cascade_matrix
    ! !USES:
    use MathfuncMod         , only : dot_sum
    use MathfuncMod         , only : safe_div
    use ch4soilBGCIndexType       , only : ch4soil_bgc_index_type
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(CentSom_type), intent(inout) :: this
    type(ch4soil_bgc_index_type)   , intent(in) :: ch4soil_bgc_index
    real(r8)                , intent(in) :: pot_decay_rates(ncentpools)
    real(r8)                , intent(in) :: cascade_matrix(ch4soil_bgc_index%nstvars, ch4soil_bgc_index%nreactions)
    real(r8)                , intent(out):: pot_co2_hr
    type(betr_status_type)  , intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    real(r8) :: cascade_matrix_hr(ncentpools)
    integer  :: reac

    associate(                                           & !
         nom_pools => ch4soil_bgc_index%nom_pools        , & !
         lid_co2_hr=> ch4soil_bgc_index%lid_co2_hr       , & !
         lit1      => ch4soil_bgc_index%lit1             , & !
         lit2      => ch4soil_bgc_index%lit2             , & !
         lit3      => ch4soil_bgc_index%lit3             , & !
         som1      => ch4soil_bgc_index%som1             , & !
         som2      => ch4soil_bgc_index%som2             , & !
         som3      => ch4soil_bgc_index%som3             , & !
         cwd       => ch4soil_bgc_index%cwd              , & !
         lwd       => ch4soil_bgc_index%lwd              , & !
         fwd       => ch4soil_bgc_index%fwd              , & !
         lit1_dek_reac=> ch4soil_bgc_index%lit1_dek_reac , & !
         lit2_dek_reac=> ch4soil_bgc_index%lit2_dek_reac , & !
         lit3_dek_reac=> ch4soil_bgc_index%lit3_dek_reac , & !
         som1_dek_reac=> ch4soil_bgc_index%som1_dek_reac , & !
         som2_dek_reac=> ch4soil_bgc_index%som2_dek_reac , & !
         som3_dek_reac=> ch4soil_bgc_index%som3_dek_reac , & !
         cwd_dek_reac=> ch4soil_bgc_index%cwd_dek_reac   , & !
         lwd_dek_reac=> ch4soil_bgc_index%lwd_dek_reac   , & !
         fwd_dek_reac=> ch4soil_bgc_index%fwd_dek_reac     & !
         )

    cascade_matrix_hr = 0._r8
    reac=lit1_dek_reac; cascade_matrix_hr(lit1)=cascade_matrix(lid_co2_hr,reac)
    reac=lit2_dek_reac; cascade_matrix_hr(lit2)=cascade_matrix(lid_co2_hr,reac)
    reac=lit3_dek_reac; cascade_matrix_hr(lit3)=cascade_matrix(lid_co2_hr,reac)
    reac=cwd_dek_reac ; cascade_matrix_hr(cwd) =cascade_matrix(lid_co2_hr,reac)
    reac=lwd_dek_reac ; cascade_matrix_hr(lwd) =cascade_matrix(lid_co2_hr,reac)
    reac=fwd_dek_reac ; cascade_matrix_hr(fwd) =cascade_matrix(lid_co2_hr,reac)
    reac=som1_dek_reac; cascade_matrix_hr(som1)=cascade_matrix(lid_co2_hr,reac)
    reac=som2_dek_reac; cascade_matrix_hr(som2)=cascade_matrix(lid_co2_hr,reac)
    reac=som3_dek_reac; cascade_matrix_hr(som3)=cascade_matrix(lid_co2_hr,reac)

    pot_co2_hr = dot_sum(cascade_matrix_hr, pot_decay_rates, bstatus)  !mol CO2/m3/s
    end associate
  end subroutine calc_potential_aerobic_hr

  !-----------------------------------------------------------------------
  subroutine calc_cnp_ratios(this, ch4soil_bgc_index, ystates, bstatus)
  !
  ! DESCRIPTION
  ! compute the cnp ratios for the om pools
  use BetrStatusType      , only : betr_status_type
  use MathfuncMod         , only : safe_div
  use ch4soilBGCIndexType       , only : ch4soil_bgc_index_type
  implicit none
  class(CentSom_type), intent(inout) :: this
  type(ch4soil_bgc_index_type)   , intent(in) :: ch4soil_bgc_index
  real(r8)                      , intent(inout) :: ystates(ch4soil_bgc_index%nstvars)
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
    nelms => ch4soil_bgc_index%nelms, &
    c_loc => ch4soil_bgc_index%c_loc, &
    n_loc => ch4soil_bgc_index%n_loc, &
    p_loc => ch4soil_bgc_index%p_loc, &
    c13_loc => ch4soil_bgc_index%c13_loc, &
    c14_loc => ch4soil_bgc_index%c14_loc, &
    lit2  => ch4soil_bgc_index%lit2 , &
    lit3  => ch4soil_bgc_index%lit3 , &
    is_cenpool_som => ch4soil_bgc_index%is_cenpool_som, &
    ompoolnames => ch4soil_bgc_index%ompoolnames &
  )
  call bstatus%reset()
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
    if(ch4soil_bgc_index%debug)then
       write(*,'(A,X,I2,5(X,E20.10))')'cnp',jj,ystates(kc),ystates(kn),ystates(kp),1._r8/this%icn_ratios(jj),1._r8/this%icp_ratios(jj)
    endif
    if(is_cenpool_som(jj) .and. ystates(kc)>tiny_val)then
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
      if(ch4soil_bgc_index%debug)then
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
  subroutine stoichiometry_fix(this, ch4soil_bgc_index,ystates)

  !
  ! DESCRIPTION
  ! this fixes the stoichiometric drift due to limite precision of
  ! double precision.
  use ch4soilBGCIndexType       , only : ch4soil_bgc_index_type
  implicit none
  class(CentSom_type)           , intent(inout) :: this
  type(ch4soil_bgc_index_type)   , intent(in) :: ch4soil_bgc_index
  real(r8)                      , intent(inout) :: ystates(ch4soil_bgc_index%nstvars)

  associate(                             &
    c13_loc => ch4soil_bgc_index%c13_loc, &
    c14_loc => ch4soil_bgc_index%c14_loc, &
    som1  => ch4soil_bgc_index%som1 , &
    som2  => ch4soil_bgc_index%som2 , &
    som3  => ch4soil_bgc_index%som3 , &
    is_cenpool_som => ch4soil_bgc_index%is_cenpool_som, &
    ompoolnames => ch4soil_bgc_index%ompoolnames &
  )


  !for om pools
  call stoi_fix(som1)

  call stoi_fix(som2)

  call stoi_fix(som3)

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
      nelms => ch4soil_bgc_index%nelms, &
      c_loc => ch4soil_bgc_index%c_loc, &
      n_loc => ch4soil_bgc_index%n_loc  &
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
  subroutine calc_som_decay_r(this, ch4soil_bgc_index, dtime, om_k_decay, om_pools, om_decay_rates)
    !
    ! !DESCRIPTION:
    ! calculate degradation for all different pools
    !
    ! !USES:
    use ch4soilBGCIndexType       , only : ch4soil_bgc_index_type
   implicit none
   class(CentSom_type)     , intent(inout) :: this
   type(ch4soil_bgc_index_type) , intent(in)    :: ch4soil_bgc_index
    real(r8)  , intent(in)    :: dtime
    real(r8)  , intent(in)    :: om_k_decay(ncentpools)
    real(r8)  , intent(in)    :: om_pools(ch4soil_bgc_index%nom_tot_elms)
    real(r8)  , intent(out)   :: om_decay_rates(ncentpools)

    ! !LOCAL VARIABLES:
    integer :: jj, fc, c, j
    integer :: kc, kn
    associate(                                        &
         nelms => ch4soil_bgc_index%nelms            , &
         nom_pools => ch4soil_bgc_index%nom_pools    , &
         c_loc => ch4soil_bgc_index%c_loc              &
    )

    !for om pools
    do jj = 1, nom_pools
      kc = (jj-1) * nelms + c_loc
      om_decay_rates(jj) = om_pools(kc) * om_k_decay(jj)
    enddo
    end associate
  end subroutine calc_som_decay_r
  !-------------------------------------------------------------------------------

  subroutine apply_spinupf(this, ch4soil_bgc_index, decompkf_eca, k_decay, spinup_scalar, spinup_flg)
  use ch4soilBGCIndexType       , only : ch4soil_bgc_index_type
  use ch4soilBGCDecompType      , only : DecompCent_type
  use betr_varcon               , only : kyr_spinup
  implicit none
  class(CentSom_type)     , intent(inout) :: this
  type(DecompCent_type), intent(in) :: decompkf_eca
  type(ch4soil_bgc_index_type)   , intent(in)    :: ch4soil_bgc_index
  real(r8)                      , intent(inout) :: k_decay(ncentpools)
  real(r8)                      , intent(inout) :: spinup_scalar
  integer                       , intent(in)    :: spinup_flg
  integer :: jj

  associate(   &
   t_scalar       => decompkf_eca%t_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil temperature scalar for decomp
   w_scalar       => decompkf_eca%w_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil water scalar for decomp
   o_scalar       => decompkf_eca%o_scalar        , & ! Intput: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
   som1           => ch4soil_bgc_index%som1        , & !
   som2           => ch4soil_bgc_index%som2        , & !
   som3           => ch4soil_bgc_index%som3          & !
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
  subroutine calc_som_decay_k(this, lay, ch4soil_bgc_index, decompkf_eca, k_decay)

  use ch4soilBGCIndexType       , only : ch4soil_bgc_index_type
  use ch4soilBGCDecompType      , only : DecompCent_type
  implicit none
  class(CentSom_type)     , intent(inout) :: this
  integer, intent(in) :: lay
  type(DecompCent_type), intent(in) :: decompkf_eca
  type(ch4soil_bgc_index_type)   , intent(in)    :: ch4soil_bgc_index
  real(r8)                      , intent(out)   :: k_decay(ncentpools)
  integer :: jj

  associate(   &
   t_scalar       => decompkf_eca%t_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil temperature scalar for decomp
   w_scalar       => decompkf_eca%w_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil water scalar for decomp
   o_scalar       => decompkf_eca%o_scalar        , & ! Intput: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
   depth_scalar   => decompkf_eca%depth_scalar    , & ! Intput: [real(r8) (:,:)   ]  rate constant for decomposition (1./sec)
   latacc         => decompkf_eca%latacc          , & !
   lit1           => ch4soil_bgc_index%lit1               , & !
   lit2           => ch4soil_bgc_index%lit2               , & !
   lit3           => ch4soil_bgc_index%lit3               , & !
   som1           => ch4soil_bgc_index%som1               , & !
   som2           => ch4soil_bgc_index%som2               , & !
   som3           => ch4soil_bgc_index%som3               , & !
   cwd            => ch4soil_bgc_index%cwd                , & !
   lwd            => ch4soil_bgc_index%lwd                , & !
   fwd            => ch4soil_bgc_index%fwd                  & !
  )
  k_decay(lit1) = this%k_decay_lit1(lay) * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(lit2) = this%k_decay_lit2(lay) * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(lit3) = this%k_decay_lit3(lay) * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(som1) = this%k_decay_som1(lay) * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(som2) = this%k_decay_som2 * t_scalar * w_scalar * o_scalar * depth_scalar * latacc
  k_decay(som3) = this%k_decay_som3 * t_scalar * w_scalar * o_scalar * depth_scalar * latacc
  k_decay(cwd)  = this%k_decay_cwd  * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(lwd)  = this%k_decay_lwd  * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(fwd)  = this%k_decay_fwd  * t_scalar * w_scalar * o_scalar * depth_scalar

  !impose the ligin effect
  k_decay(cwd)  = k_decay(cwd) * exp(-3._r8*this%cwd_flig)
  k_decay(lwd)  = k_decay(lwd) * exp(-3._r8*this%lwd_flig)
  k_decay(fwd)  = k_decay(fwd) * exp(-3._r8*this%fwd_flig)
  k_decay(lit2) = k_decay(lit2)* exp(-3._r8*this%lit_flig)
  k_decay(lit3) = k_decay(lit3)* exp(-3._r8*this%lit_flig)
!  print*,'lay',lay
!  print*,'decy',this%k_decay_lit1(lay),this%k_decay_lit2(lay),this%k_decay_lit3(lay),this%k_decay_som1,this%k_decay_som2,this%k_decay_som3
!  print*,'k_decay',(k_decay(jj),jj=1,ncentpools)
!  print*,'scalar',t_scalar, w_scalar, o_scalar, depth_scalar
  end associate
  end subroutine calc_som_decay_k
  !-------------------------------------------------------------------------------
  subroutine calc_pot_min_np_flx(this, dtime, ch4soil_bgc_index, ystates, k_decay, cascade_matrix, &
    alpha_n, alpha_p, pot_decomp, pot_nn_flx, pot_np_flx)
  use ch4soilBGCIndexType       , only : ch4soil_bgc_index_type
  implicit none
  class(CentSom_type)         , intent(inout) :: this
  real(r8)                    , intent(in) :: dtime
  type(ch4soil_bgc_index_type) , intent(in) :: ch4soil_bgc_index
  real(r8)                    , intent(in) :: ystates(1:ch4soil_bgc_index%nom_tot_elms)
  real(r8)                    , intent(in) :: k_decay(1:ncentpools)
  real(r8)                    , intent(in) :: cascade_matrix(ch4soil_bgc_index%nstvars, ch4soil_bgc_index%nreactions)
  real(r8)                    , intent(in) :: alpha_n(ncentpools)
  real(r8)                    , intent(in) :: alpha_p(ncentpools)
  real(r8)                    , intent(out) :: pot_decomp(ncentpools)
  real(r8)                    , intent(out):: pot_nn_flx
  real(r8)                    , intent(out):: pot_np_flx

  integer :: reac
  integer :: reacs(ncentpools)

  associate(                                                    & !
       nom_pools => ch4soil_bgc_index%nom_pools                , & !
       nom_tot_elms=> ch4soil_bgc_index%nom_tot_elms           , & !
       lid_nh4   => ch4soil_bgc_index%lid_nh4                  , & !
       lid_minp_soluble  => ch4soil_bgc_index%lid_minp_soluble , & !
       lit1      => ch4soil_bgc_index%lit1                     , & !
       lit2      => ch4soil_bgc_index%lit2                     , & !
       lit3      => ch4soil_bgc_index%lit3                     , & !
       som1      => ch4soil_bgc_index%som1                     , & !
       som2      => ch4soil_bgc_index%som2                     , & !
       som3      => ch4soil_bgc_index%som3                     , & !
       cwd       => ch4soil_bgc_index%cwd                      , & !
       lit1_dek_reac=> ch4soil_bgc_index%lit1_dek_reac         , & !
       lit2_dek_reac=> ch4soil_bgc_index%lit2_dek_reac         , & !
       lit3_dek_reac=> ch4soil_bgc_index%lit3_dek_reac         , & !
       som1_dek_reac=> ch4soil_bgc_index%som1_dek_reac         , & !
       som2_dek_reac=> ch4soil_bgc_index%som2_dek_reac         , & !
       som3_dek_reac=> ch4soil_bgc_index%som3_dek_reac         , & !
       cwd_dek_reac=> ch4soil_bgc_index%cwd_dek_reac           , & !
       lwd_dek_reac=> ch4soil_bgc_index%lwd_dek_reac           , & !
       fwd_dek_reac=> ch4soil_bgc_index%fwd_dek_reac             & !
   )

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(ch4soil_bgc_index, dtime, k_decay(1:nom_pools), &
      ystates(1:nom_tot_elms), pot_decomp)

  pot_nn_flx = 0._r8; pot_np_flx = 0._r8

  reacs=(/lit1_dek_reac, lit2_dek_reac, lit3_dek_reac, &
    cwd_dek_reac, lwd_dek_reac, fwd_dek_reac, &
    som1_dek_reac, som2_dek_reac, som3_dek_reac/)

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

end module ch4soilBGCSOMType
