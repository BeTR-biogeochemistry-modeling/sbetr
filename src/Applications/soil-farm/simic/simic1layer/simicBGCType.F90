module simicBGCType
#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  ! subroutines for stoichiometric configuration of the century bgc
  ! !History, created by Jinyun Tang, Dec, 2014.
  ! !USES:
  use bshr_assert_mod, only : shr_assert
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use bshr_assert_mod, only : shr_assert_any
  use bshr_kind_mod             , only : r8 => shr_kind_r8
  use bshr_log_mod              , only : errMsg => shr_log_errMsg
  use betr_varcon               , only : spval => bspval
  use betr_ctrl                 , only : spinup_state => betr_spinup_state
  use gbetrType                 , only : gbetr_type
  use BiogeoConType             , only : BiogeoCon_type
  use simicParaType             , only : simic_para_type
  use BetrStatusType            , only : betr_status_type
  use BeTRJarModel              , only : jar_model_type
  use simicBGCIndexType            , only : simic_bgc_index_type
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !Note:
  !Keeping simic_bgc_index as a private member is a workaround to call the ode solver
  !it increase the memory for each instance of simic_bgc_type, but enables
  !the ode function to be called by the ode solver

  type, extends(jar_model_type), public :: simic_bgc_type
    type(simic_bgc_index_type),     private :: simic_bgc_index
    real(r8), pointer                    :: ystates0(:)
    real(r8), pointer                    :: ystates1(:)
    real(r8), pointer                    :: scal_f(:)
    real(r8), pointer                    :: conv_f(:)
    real(r8), pointer                    :: conc_f(:)
    logical , private                    :: use_c13
    logical , private                    :: use_c14

    !decomposition
    real(r8) :: Kaff_EP_LIT1  !enzyme affinity for litter POM depolymerization, mol C/m3
    real(r8) :: Kaff_EP_LIT2  !enzyme affinity for litter POM depolymerization, mol C/m3
    real(r8) :: Kaff_EP_LIT3  !enzyme affinity for litter POM depolymerization, mol C/m3
    real(r8) :: Kaff_EP_CWD  !enzyme affinity for litter POM depolymerization, mol C/m3
    real(r8) :: Kaff_EP_POM  !enzyme affinity for litter POM depolymerization, mol C/m3
    real(r8) :: Kaff_BC  !microbial affinity for DOC, mol C/m3
    real(r8) :: Kaff_CM  !DOC affinity for mineral surface, mol C /m3
    real(r8) :: Kaff_ED  !enzyme affinity for dead microbial cell material, mol C/m3
    real(r8) :: cue_met  !assimilation efficiency of DOC from metabolic polymer
    real(r8) :: cue_cel  !assimilation efficiency of DOC from cellulose polymer
    real(r8) :: cue_lig  !assimilation efficiency of DOC from lignin polymer
    real(r8) :: cue_cwd  !assimilation efficiency of DOC from cwd polymer
    real(r8) :: cue_bm   !assimilation efficiency of DOC from microbial biomass
    real(r8) :: Rm0_spmic  !specific microbial maintenance respiration, 1/s
    real(r8) :: Mrt_spmic  !specific microbial mortality, 1/s
    real(r8) :: f_mic2C    !fraction of dead microbial biomass into DOC upon death
    real(r8) :: f_mic2D    !fraction of dead microbial biomass as dead cell
    real(r8) :: vmax_EP_LIT1    !maximum depolymerization rate
    real(r8) :: vmax_EP_LIT2    !maximum depolymerization rate
    real(r8) :: vmax_EP_LIT3    !maximum depolymerization rate
    real(r8) :: vmax_EP_CWD    !maximum depolymerization rate
    real(r8) :: vmax_EP_POM    !maximum depolymerization rate
    real(r8) :: vmax_EP_MD
    real(r8) :: vmax_BC    !maximum DOC assimilation rate
    real(r8) :: alpha_B2E  !scaling parameter from microbial biomass to hydrolysis enzyme
    real(r8) :: alpha_B2T !scaling parameter from microbial biomass to transporter enzyme
    real(r8) :: Kaff_EM  !enzyme affinity for mineral surface, mol C /m3
    real(r8) :: Minsurf  !mineral surface area for DOC/enzyme/microbial cell wall material adsorption
    real(r8) :: Kaff_o2
    real(r8) :: Kmort_MB
    real(r8) :: Kaff_CM_ref
    real(r8) :: Kaff_EM_ref
    real(r8) :: fpom_vmax
    real(r8) :: fpom_desorb
    real(r8), pointer :: cascade_matrix(:,:)
    real(r8), pointer :: cascade_matrixd(:,:)
    real(r8), pointer :: cascade_matrixp(:,:)
    real(r8) :: tfng
    real(r8) :: tfnr
    real(r8) :: o2_w2b
    logical,private  :: batch_mode
  contains
    procedure, public  :: init          => init_simic
    procedure, public  :: runbgc        => runbgc_simic
    procedure, public  :: UpdateParas   => UpdateParas_simic
    procedure, public  :: getvarllen    => getvarllen_simic
    procedure, public  :: getvarlist    => getvarlist_simic
    procedure, public  :: init_cold     => init_cold_simic
    procedure, private :: arenchyma_gas_transport
    procedure, private :: init_states
    procedure, private :: add_ext_input
    procedure, private :: InitAllocate
    procedure, private :: simic_rrates
    procedure, private :: calc_cascade_matrix
    procedure, private :: bgc_integrate
    procedure, private :: ode_adapt_ebbks1
    procedure, private :: sum_tot_store
  end type simic_bgc_type

  public :: create_jarmodel_simicbgc
contains

  function create_jarmodel_simicbgc()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(simic_bgc_type), pointer :: create_jarmodel_simicbgc
    class(simic_bgc_type), pointer :: bgc

    allocate(bgc)
    create_jarmodel_simicbgc => bgc

  end function create_jarmodel_simicbgc

  !-------------------------------------------------------------------------------
  function getvarllen_simic(this)result(ans)

  implicit none
  class(simic_bgc_type) , intent(inout) :: this
  integer :: ans

  ans =  this%simic_bgc_index%nstvars

  end function getvarllen_simic
  !-------------------------------------------------------------------------------
  subroutine getvarlist_simic(this, nstvars, varnames, varunits, vartypes)
  implicit none
  class(simic_bgc_type) , intent(inout) :: this
  integer, intent(in) :: nstvars
  character(len=*), intent(out) :: varnames(1:nstvars)
  character(len=*), intent(out) :: varunits(1:nstvars)
  integer         , intent(out) :: vartypes(1:nstvars)
  integer :: n

  do n = 1, nstvars
    vartypes(n) = this%simic_bgc_index%vartypes(n)
    write(varnames(n),'(A)')trim(this%simic_bgc_index%varnames(n))
    write(varunits(n),'(A)')trim(this%simic_bgc_index%varunits(n))
  enddo
  end subroutine getvarlist_simic


  !-------------------------------------------------------------------------------
  subroutine UpdateParas_simic(this,  biogeo_con, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(simic_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)     , intent(out)   :: bstatus
  integer :: sr
  character(len=256) :: msg

  call bstatus%reset()

  select type(biogeo_con)
  type is(simic_para_type)
    !decomposition
    this%Kaff_EP_LIT1  = biogeo_con%Kaff_EP_LIT
    this%Kaff_EP_LIT2  = this%Kaff_EP_LIT1 * 4._r8
    this%Kaff_EP_LIT3  = this%Kaff_EP_LIT1 * 4._r8
    this%Kaff_EP_CWD  = this%Kaff_EP_LIT1 * 10._r8
    this%Kaff_EP_POM  = biogeo_con%Kaff_EP_POM
    this%Kaff_BC  = biogeo_con%Kaff_BC
    this%Kaff_CM_ref  = biogeo_con%Kaff_CM
    this%Kaff_ED  = biogeo_con%Kaff_ED
    this%cue_met  = biogeo_con%cue_met
    this%cue_cel  = biogeo_con%cue_cel
    this%cue_lig  = biogeo_con%cue_lig
    this%cue_cwd  = biogeo_con%cue_cwd
    this%cue_bm  = biogeo_con%cue_bm
    this%Rm0_spmic= biogeo_con%Rm0_spmic
    this%Mrt_spmic= biogeo_con%Mrt_spmic
    this%f_mic2C  = biogeo_con%f_mic2C
    this%f_mic2D  = biogeo_con%f_mic2D
    this%vmax_EP_LIT1  = biogeo_con%vmax_EP_L
    this%vmax_EP_LIT2  = this%vmax_EP_LIT1 * 0.25_r8
    this%vmax_EP_LIT3  = this%vmax_EP_LIT1 * 0.25_r8
    this%vmax_EP_CWD  = this%vmax_EP_LIT1 * 0.1_r8
    this%vmax_EP_POM  = this%vmax_EP_LIT1 * 0.1_r8
    this%vmax_EP_MD  = this%vmax_EP_LIT1 * 0.1_r8
    this%vmax_BC  = biogeo_con%vmax_BC
    this%alpha_B2E= biogeo_con%alpha_B2E
    this%alpha_B2T= biogeo_con%alpha_B2T
    this%Kaff_EM_ref  = biogeo_con%Kaff_EM
    this%Kaff_o2  = biogeo_con%Kaff_o2
    this%Kmort_MB = biogeo_con%Kmort_MB
    this%fpom_vmax= biogeo_con%fpom_vmax
    this%fpom_desorb=biogeo_con%fpom_desorb
  class default
    write(msg,'(A)')'Wrong parameter type passed in for UpdateParas in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine UpdateParas_simic
  !-------------------------------------------------------------------------------
  subroutine init_simic(this,  biogeo_con, batch_mode,  bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(simic_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  logical                , intent(in) :: batch_mode
  type(betr_status_type)    , intent(out) :: bstatus

  character(len=256) :: msg
  write(this%jarname, '(A)')'simic'

  this%batch_mode=batch_mode
  select type(biogeo_con)
  type is(simic_para_type)
    call bstatus%reset()
    call this%simic_bgc_index%Init(biogeo_con%use_c13, biogeo_con%use_c14, &
     biogeo_con%non_limit, biogeo_con%nop_limit, betr_maxpatch_pft, this%batch_mode)
    call this%UpdateParas(biogeo_con, bstatus)
    if(bstatus%check_status())return

    this%use_c13 = biogeo_con%use_c13
    this%use_c14 = biogeo_con%use_c14
  class default
    call bstatus%reset()
    write(msg,'(A)')'Wrong parameter type passed in for init_simic in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select

  call this%InitAllocate(this%simic_bgc_index)
  end subroutine init_simic

  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, simic_bgc_index)

  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(simic_bgc_type)   , intent(inout) :: this
  type(simic_bgc_index_type)  , intent(in) :: simic_bgc_index

  associate(                                &
    nstvars => simic_bgc_index%nstvars        , &
    nprimvars=> simic_bgc_index%nprimvars     , &
    nreactions => simic_bgc_index%nreactions    &
  )

  allocate(this%ystates0(nstvars)); this%ystates0(:) = 0._r8
  allocate(this%ystates1(nstvars)); this%ystates1(:) = 0._r8
  allocate(this%scal_f(nprimvars));  this%scal_f(:) = 0._r8
  allocate(this%conv_f(nprimvars));  this%conv_f(:) = 0._r8
  allocate(this%conc_f(nprimvars));  this%conc_f(:) = 0._r8
  allocate(this%cascade_matrix(1:nstvars, 1:nreactions)); this%cascade_matrix(:,:) = 0._r8
  allocate(this%cascade_matrixd(1:nprimvars, 1:nreactions)); this%cascade_matrixd(:,:) = 0._r8
  allocate(this%cascade_matrixp(1:nprimvars, 1:nreactions)); this%cascade_matrixp(:,:) = 0._r8
  end associate
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------
  subroutine runbgc_simic(this,  is_surflit, dtime, bgc_forc, nstates, ystates0, ystatesf, bstatus)

  !DESCRIPTION
  !do bgc model integration for one step
  use JarBgcForcType        , only : JarBGC_forc_type
  use MathfuncMod           , only : pd_decomp
  use BetrStatusType        , only : betr_status_type
  use MathfuncMod           , only : safe_div
  use tracer_varcon         , only : catomw, natomw, patomw, input_only
  use EcosysMicDynParamMod  , only : calc_tm_factor
  implicit none
  class(simic_bgc_type)  , intent(inout) :: this
  logical                    , intent(in)    :: is_surflit
  real(r8)                   , intent(in)    :: dtime
  type(JarBGC_forc_type)     , intent(in)    :: bgc_forc
  integer                    , intent(in)    :: nstates
  real(r8)                   , intent(out)   :: ystates0(nstates)
  real(r8)                   , intent(out)   :: ystatesf(nstates)
  type(betr_status_type)     , intent(out)   :: bstatus

  real(r8)               :: time = 0._r8
  real(r8) :: yf(this%simic_bgc_index%nstvars)
  !local variables
  character(len=*),parameter :: subname = 'runbgc_simic'
  integer :: jj, id
  associate(                                         &
    ystates1       => this%ystates1                , &
    nstvars        => this%simic_bgc_index%nstvars     , &
    nreactions     => this%simic_bgc_index%nreactions  , &
    nprimvars      => this%simic_bgc_index%nprimvars   , &
    lid_doc        => this%simic_bgc_index%lid_doc , &
    lid_doc_e      => this%simic_bgc_index%lid_doc_e , &
    cascade_matrix => this%cascade_matrix          , &
    cascade_matrixp=> this%cascade_matrixp         , &
    cascade_matrixd=> this%cascade_matrixd           &
  )
  call bstatus%reset()
  this%o2_w2b        = bgc_forc%o2_w2b
  cascade_matrix(:,:) = 0._r8

  !initialize state variables
  call this%init_states(this%simic_bgc_index, bgc_forc)

  ystates0(:) = this%ystates0(:)

  call this%add_ext_input(dtime, this%simic_bgc_index, bgc_forc)

  if(.not. input_only)then
    this%Minsurf = bgc_forc%Msurf_OM
    this%Kaff_CM   = this%Kaff_CM_ref*bgc_forc%KM_OM_ref
    this%Kaff_EM   = this%Kaff_EM_ref*bgc_forc%KM_OM_ref

    call this%arenchyma_gas_transport(this%simic_bgc_index, dtime)

    call calc_tm_factor(bgc_forc%temp, bgc_forc%soilpsi, this%tfng, this%tfnr)

    call this%calc_cascade_matrix(this%simic_bgc_index,  cascade_matrix)

    call pd_decomp(nprimvars, nreactions, cascade_matrix(1:nprimvars, 1:nreactions), &
       cascade_matrixp, cascade_matrixd, bstatus)
    if(bstatus%check_status())return

    time = 0._r8
    yf(:) = ystates1(:)
    call this%ode_adapt_ebbks1(yf, nprimvars, nstvars, time, dtime, ystates1)
  endif
  if(this%batch_mode)call this%sum_tot_store(nstvars, ystates1)
  ystatesf(:) = ystates1(:)

  end associate
  end subroutine runbgc_simic
  !-------------------------------------------------------------------------------
  subroutine simic_rrates(this, simic_bgc_index, dtime,  nstates, ystates1, doc_cue, rrates)
  !
  !DESCRIPTION
  !the core of simic model
  use KineticsMod, only : ecacomplex
  use BetrStatusType  , only : betr_status_type
  implicit none
  class(simic_bgc_type)  , intent(inout) :: this
  real(r8)               , intent(in)    :: dtime
  integer                , intent(in)    :: nstates
  real(r8)               , intent(in)    :: ystates1(1:nstates)
  type(simic_bgc_index_type) , intent(in):: simic_bgc_index
  real(r8)               , intent(in)    :: doc_cue
  real(r8)               , intent(out)   :: rrates(1:simic_bgc_index%nreactions)

  real(r8) :: depolymer_l1,depolymer_l2,depolymer_l3,depolymer_cwd
  real(r8) :: depolymer_norm, depolymer_md
  real(r8) :: doc_uptake, Rh_pot, Rh_gpot, Rm_pot, Rmx
  real(r8) :: o2w, fo2, mort
  real(r8) :: vmax_EP_POM_f, vmax_EP_L1_f,vmax_EP_L2_f,vmax_EP_L3_f, vmax_EP_CWD_f
  real(r8) :: vmax_EP_MD_f
  real(r8) :: Minsurf, denorm, depolymer_pom1, depolymer_pom2
  real(r8) :: doc_sorb
  real(r8) :: kd71(1:7),ss7(1:7),ee1,siej71(1:7)
  real(r8) :: kd23(1:2,1:3),ss2(1:2),ee3(1:3),siej23(1:2,1:3)
  integer  :: jj
  type(betr_status_type)  :: bstatus
  associate(                                            &
    lit1             => simic_bgc_index%lit1              , &
    lit2             => simic_bgc_index%lit2              , &
    lit3             => simic_bgc_index%lit3              , &
    cwd              => simic_bgc_index%cwd               , &
    lid_doc          => simic_bgc_index%lid_doc           , &
    lid_pom          => simic_bgc_index%lid_pom           , &
    lid_micbl        => simic_bgc_index%lid_micbl         , &
    lid_micbd        => simic_bgc_index%lid_micbd         , &
    lid_doc_e          => simic_bgc_index%lid_doc_e           , &
    lid_o2           => simic_bgc_index%lid_o2            , &
    lid_co2          => simic_bgc_index%lid_co2           , &
    lit1_depoly_reac => simic_bgc_index%lit1_depoly_reac  , &
    lit2_depoly_reac => simic_bgc_index%lit2_depoly_reac  , &
    lit3_depoly_reac => simic_bgc_index%lit3_depoly_reac  , &
    cwd_depoly_reac  => simic_bgc_index%cwd_depoly_reac   , &
    micbd_depoly_reac=> simic_bgc_index%micbd_depoly_reac , &
    micbl_mort_reac  => simic_bgc_index%micbl_mort_reac   , &
    doc_uptake_reac  => simic_bgc_index%doc_uptake_reac   , &
    o2_resp_reac     => simic_bgc_index%o2_resp_reac      , &
    pom_desorb_reac  => simic_bgc_index%pom_desorb_reac   , &
    doc_sorb_reac    => simic_bgc_index%doc_sorb_reac     , &
    alpha_B2E        => this%alpha_B2E                , &
    alpha_B2T        => this%alpha_B2T                , &
    Kaff_o2          => this%Kaff_o2                  , &
    Kaff_EP_LIT1     => this%Kaff_EP_LIT1              , &
    Kaff_EP_LIT2     => this%Kaff_EP_LIT2              , &
    Kaff_EP_LIT3     => this%Kaff_EP_LIT3              , &
    Kaff_EP_CWD     => this%Kaff_EP_CWD              , &
    Kaff_EP_POM       => this%Kaff_EP_POM            , &
    Kaff_ED          => this%Kaff_ED                  , &
    Kaff_EM          => this%Kaff_EM                  , &
    Minsurf0         => this%Minsurf                  , &
    vmax_EP_LIT1        => this%vmax_EP_LIT1          , &
    vmax_EP_LIT2        => this%vmax_EP_LIT2          , &
    vmax_EP_LIT3        => this%vmax_EP_LIT3          , &
    vmax_EP_CWD        => this%vmax_EP_CWD             , &
    vmax_EP_POM        => this%vmax_EP_POM             , &
    vmax_EP_MD        => this%vmax_EP_MD             , &
    Mrt_spmic        => this%Mrt_spmic                , &
    Kmort_MB         => this%Kmort_MB                 , &
    f_mic2C          => this%f_mic2C                  , &
    f_mic2D          => this%f_mic2D                  , &
    Kaff_BC          => this%Kaff_BC                  , &
    Kaff_CM          => this%Kaff_CM                  , &
    Rm0_spmic        => this%Rm0_spmic                , &
    vmax_BC          => this%vmax_BC                  , &
    fpom_desorb      => this%fpom_desorb              , &
    fpom_vmax        => this%fpom_vmax                , &
    tfng             => this%tfng                     , &
    tfnr             => this%tfnr                       &
  )
  o2w = ystates1(lid_o2) / this%o2_w2b
  fo2= o2w/(Kaff_o2+o2w+alpha_B2T*ystates1(lid_micbl))
  Minsurf = max(Minsurf0 -ystates1(lid_pom),0._r8)

  kd71=(/Kaff_EP_LIT1,Kaff_EP_LIT2,Kaff_EP_LIT3,Kaff_EP_CWD,Kaff_ED,Kaff_EP_POM,Kaff_EM/)
  ss7=(/ystates1(lit1),ystates1(lit2),ystates1(lit3),ystates1(cwd),ystates1(lid_micbd),&
    ystates1(lid_pom),Minsurf/)
  ee1=ystates1(lid_micbl) * alpha_B2E

  call ecacomplex(kd71,ss7,ee1,siej71)

  depolymer_l1   = siej71(1) * vmax_EP_LIT1 * tfng
  depolymer_l2   = siej71(2) * vmax_EP_LIT2 * tfng
  depolymer_l3   = siej71(3) * vmax_EP_LIT3 * tfng
  depolymer_cwd  = siej71(4) * vmax_EP_CWD * tfng
  depolymer_md   = siej71(5) * vmax_EP_MD * tfng
  depolymer_pom1 = siej71(6) * vmax_EP_POM * tfng

  !potential respiration
  kd23(1,1:3)=(/Kaff_BC, Kaff_BC, Kaff_CM/)
  kd23(2,1:3)=(/0._r8,Kaff_EM,0._r8/)
  ee3(1:3)=(/ystates1(lid_micbl) * alpha_B2T, ystates1(lid_micbd) * alpha_B2T, Minsurf/)
  ss2(1:2)= (/ystates1(lid_doc),ystates1(lid_micbl) * alpha_B2E/)
  call ecacomplex(kd23,ss2,ee3,siej23,bstatus)

  doc_uptake = vmax_BC * siej23(1,1) * fo2 * tfng
  doc_sorb = fpom_vmax * siej23(1,2)

  !maintenance respiration
  Rm_pot = Rm0_spmic * ystates1(lid_micbl) * tfnr

  !potential growth respiration
  if(doc_uptake > Rm_pot)then
    Rh_pot = Rm_pot + (doc_uptake-Rm_pot)*(1._r8-doc_cue)
  else
    Rh_pot = doc_uptake
  endif
  Rmx = Rm_pot - doc_uptake

  !compute mortality
  mort = Mrt_spmic * ystates1(lid_micbl)
  if(Rmx > 0._r8)then
    mort = mort + Rmx
  endif
  mort = mort * ystates1(lid_micbl) / (ystates1(lid_micbl) + Kmort_MB)

  depolymer_pom2 = ystates1(lid_pom) * fpom_desorb
  !assemble the derivatives
  rrates(:) = 0._r8

  rrates(lit1_depoly_reac)     = depolymer_l1
  rrates(lit2_depoly_reac)     = depolymer_l2
  rrates(lit3_depoly_reac)     = depolymer_l3
  rrates(cwd_depoly_reac)      = depolymer_cwd
  rrates(micbd_depoly_reac)    = depolymer_md
  rrates(micbl_mort_reac)      = mort
  rrates(doc_uptake_reac)      = doc_uptake
  rrates(o2_resp_reac)         = Rh_pot
  rrates(pom_desorb_reac)      = depolymer_pom1+depolymer_pom2
  rrates(doc_sorb_reac)        = doc_sorb

  end associate
  end subroutine simic_rrates
  !-------------------------------------------------------------------------------
  subroutine init_states(this, simic_bgc_index, bgc_forc)

  use simicBGCIndexType            , only : simic_bgc_index_type
  use JarBgcForcType            , only : JarBGC_forc_type
  implicit none
  class(simic_bgc_type)  , intent(inout) :: this
  type(simic_bgc_index_type)  , intent(in) :: simic_bgc_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  associate(                               &
    lid_co2_hr  => simic_bgc_index%lid_co2_hr, &
    lid_n2 => simic_bgc_index%lid_n2,   &
    lid_o2 => simic_bgc_index%lid_o2,   &
    lid_co2 => simic_bgc_index%lid_co2, &
    lid_c13_co2 => simic_bgc_index%lid_c13_co2, &
    lid_c14_co2 => simic_bgc_index%lid_c14_co2, &
    lid_ar => simic_bgc_index%lid_ar,   &
    lid_o2_paere => simic_bgc_index%lid_o2_paere, &
    lid_n2_paere => simic_bgc_index%lid_n2_paere, &
    lid_ar_paere => simic_bgc_index%lid_ar_paere, &
    lid_co2_paere => simic_bgc_index%lid_co2_paere, &
    lid_ch4_paere => simic_bgc_index%lid_ch4_paere, &
    lid_c13_co2_paere => simic_bgc_index%lid_c13_co2_paere, &
    lid_c14_co2_paere => simic_bgc_index%lid_c14_co2_paere, &
    lid_ch4 => simic_bgc_index%lid_ch4  &
  )

  this%ystates0(:) = bgc_forc%ystates(:)
  this%ystates0(lid_co2_hr) = 0._r8
  this%ystates1(:) = this%ystates0(:)

  !set conversion parameters for arenchyma transport
  this%scal_f(lid_n2) = bgc_forc%aren_cond_n2
  this%conc_f(lid_n2) = bgc_forc%conc_atm_n2
  this%conv_f(lid_n2) = 1._r8/bgc_forc%n2_g2b

  this%scal_f(lid_o2) = bgc_forc%aren_cond_o2
  this%conc_f(lid_o2) = bgc_forc%conc_atm_o2
  this%conv_f(lid_o2) = 1._r8/bgc_forc%o2_g2b

  this%scal_f(lid_ar) = bgc_forc%aren_cond_ar
  this%conc_f(lid_ar) = bgc_forc%conc_atm_ar
  this%conv_f(lid_ar) = 1._r8/bgc_forc%ar_g2b

  this%scal_f(lid_co2) = bgc_forc%aren_cond_co2
  this%conc_f(lid_co2) = bgc_forc%conc_atm_co2
  this%conv_f(lid_co2) = 1._r8/bgc_forc%co2_g2b

  if(this%use_c13)then
    this%scal_f(lid_c13_co2) = bgc_forc%aren_cond_co2_c13
    this%conc_f(lid_c13_co2) = bgc_forc%conc_atm_co2_c13
    this%conv_f(lid_c13_co2) = 1._r8/bgc_forc%co2_g2b
  endif

  if(this%use_c14)then
    this%scal_f(lid_c14_co2) = bgc_forc%aren_cond_co2_c14
    this%conc_f(lid_c14_co2) = bgc_forc%conc_atm_co2_c14
  endif

  this%scal_f(lid_ch4) = bgc_forc%aren_cond_ch4
  this%conc_f(lid_ch4) = bgc_forc%conc_atm_ch4
  this%conv_f(lid_ch4) = 1._r8/bgc_forc%ch4_g2b

  end associate
  end subroutine init_states

  !--------------------------------------------------------------------
  subroutine add_ext_input(this, dtime, simic_bgc_index, bgc_forc)

  use simicBGCIndexType            , only : simic_bgc_index_type
  use JarBgcForcType        , only : JarBGC_forc_type
  use tracer_varcon             , only : catomw, natomw, patomw,c13atomw,c14atomw
  use MathfuncMod               , only : safe_div
  implicit none
  class(simic_bgc_type)  , intent(inout) :: this
  real(r8), intent(in) :: dtime
  type(simic_bgc_index_type)  , intent(in) :: simic_bgc_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  associate(                        &
    lit1 =>  simic_bgc_index%lit1, &
    lit2 =>  simic_bgc_index%lit2, &
    lit3 =>  simic_bgc_index%lit3, &
    lid_totinput =>  simic_bgc_index%lid_totinput, &
    cwd =>   simic_bgc_index%cwd   &
  )

  this%ystates1(lit1)=this%ystates0(lit1) +  dtime *  bgc_forc%cflx_input_litr_met/catomw

  this%ystates1(lit2)=this%ystates0(lit2) +  dtime *  bgc_forc%cflx_input_litr_cel/catomw

  this%ystates1(lit3)=this%ystates0(lit3) +  dtime *  bgc_forc%cflx_input_litr_lig/catomw

  this%ystates1(cwd)=this%ystates0(cwd) +  dtime *  bgc_forc%cflx_input_litr_cwd/catomw

  if(this%batch_mode)then
    this%ystates1(lid_totinput) = this%ystates0(lid_totinput) + dtime* &
      (bgc_forc%cflx_input_litr_met + bgc_forc%cflx_input_litr_cel + &
      bgc_forc%cflx_input_litr_lig + bgc_forc%cflx_input_litr_cwd)/catomw
  endif

  end associate
  end subroutine add_ext_input


  !--------------------------------------------------------------------
  subroutine arenchyma_gas_transport(this, simic_bgc_index, dtime)
  use simicBGCIndexType       , only : simic_bgc_index_type
  implicit none
  class(simic_bgc_type)     , intent(inout) :: this
  type(simic_bgc_index_type)  , intent(in) :: simic_bgc_index
  real(r8), intent(in) :: dtime

  integer :: j
  real(r8) :: y0
  associate(                                            &
    lid_n2            => simic_bgc_index%lid_n2           , &
    lid_o2            => simic_bgc_index%lid_o2           , &
    lid_co2           => simic_bgc_index%lid_co2          , &
    lid_c13_co2       => simic_bgc_index%lid_c13_co2      , &
    lid_c14_co2       => simic_bgc_index%lid_c14_co2      , &
    lid_ar            => simic_bgc_index%lid_ar           , &
    lid_o2_paere      => simic_bgc_index%lid_o2_paere     , &
    lid_n2_paere      => simic_bgc_index%lid_n2_paere     , &
    lid_ar_paere      => simic_bgc_index%lid_ar_paere     , &
    lid_co2_paere     => simic_bgc_index%lid_co2_paere    , &
    lid_ch4_paere     => simic_bgc_index%lid_ch4_paere    , &
    lid_c13_co2_paere => simic_bgc_index%lid_c13_co2_paere, &
    lid_c14_co2_paere => simic_bgc_index%lid_c14_co2_paere, &
    lid_ch4           => simic_bgc_index%lid_ch4            &
  )

  j = lid_o2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(simic_bgc_index%lid_o2_paere) = this%ystates1(j)-y0

  if( spinup_state == 0)then
    j = lid_n2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_bgc_index%lid_n2_paere) = this%ystates1(j)-y0

    j = lid_ar; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_bgc_index%lid_ar_paere) = this%ystates1(j)-y0

    j = lid_ch4; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_bgc_index%lid_ch4_paere) = this%ystates1(j)-y0

    j = lid_co2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_bgc_index%lid_co2_paere) = this%ystates1(j)-y0

    if(this%use_c13)then
      j = lid_c13_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(simic_bgc_index%lid_c13_co2_paere) = this%ystates1(j)-y0

    endif

    if(this%use_c14)then
      j = lid_c14_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(simic_bgc_index%lid_c14_co2_paere) = this%ystates1(j)-y0

    endif
  endif
  end associate
  contains
    subroutine exp_ode_int(dt, c1, c2, c3, y0)
    !
    ! DESCRIPTION
    ! solve dy/dt=-c1*(c2*y-c3) using analytic solution
    implicit none
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: c1
    real(r8), intent(in) :: c2
    real(r8), intent(in) :: c3
    real(r8), intent(inout) :: y0

    if(c1>0._r8)then
      y0 = c3/c2+(y0-c3/c2)*exp(-c1/c2*dtime)
    endif

    end subroutine exp_ode_int
  end subroutine arenchyma_gas_transport


  !--------------------------------------------------------------------
  subroutine calc_cascade_matrix(this, simic_bgc_index,  cascade_matrix)

  implicit none
  class(simic_bgc_type)     , intent(inout) :: this
  type(simic_bgc_index_type)    , intent(in)    :: simic_bgc_index
  real(r8)                  , intent(inout) :: cascade_matrix(1:simic_bgc_index%nstvars, 1:simic_bgc_index%nreactions)
  integer :: reac

  associate(                                &
    lit1       => simic_bgc_index%lit1        , &
    lit2       => simic_bgc_index%lit2        , &
    lit3       => simic_bgc_index%lit3        , &
    cwd        => simic_bgc_index%cwd         , &
    lid_doc    => simic_bgc_index%lid_doc     , &
    lid_micbl  => simic_bgc_index%lid_micbl   , &
    lid_micbd  => simic_bgc_index%lid_micbd   , &
    lid_doc_e  => simic_bgc_index%lid_doc_e   , &
    lid_o2     => simic_bgc_index%lid_o2      , &
    lid_co2    => simic_bgc_index%lid_co2     , &
    lid_co2_hr => simic_bgc_index%lid_co2_hr  , &
    lid_pom    => simic_bgc_index%lid_pom     , &
    lid_pom_e    => simic_bgc_index%lid_pom_e , &
    cue_met    => this%cue_met            , &
    cue_cel    => this%cue_cel            , &
    cue_lig    => this%cue_lig            , &
    cue_cwd    => this%cue_cwd            , &
    cue_bm     => this%cue_bm             , &
    f_mic2d    => this%f_mic2d            , &
    f_mic2c    => this%f_mic2c              &
  )

  reac = simic_bgc_index%lit1_depoly_reac
  cascade_matrix(lit1, reac)  = -1._r8
  cascade_matrix(lid_doc,reac)=  1._r8
  cascade_matrix(lid_doc_e,reac)=cue_met

  reac = simic_bgc_index%lit2_depoly_reac
  cascade_matrix(lit2, reac)  = -1._r8
  cascade_matrix(lid_doc,reac)=  1._r8
  cascade_matrix(lid_doc_e,reac)= cue_cel

  reac = simic_bgc_index%lit3_depoly_reac
  cascade_matrix(lit3, reac)  = -1._r8
  cascade_matrix(lid_doc,reac)=  1._r8
  cascade_matrix(lid_doc_e,reac)=cue_lig

  reac = simic_bgc_index%cwd_depoly_reac
  cascade_matrix(cwd, reac)   = -1._r8
  cascade_matrix(lid_doc,reac)=  1._r8
  cascade_matrix(lid_doc_e,reac)=cue_cwd

  reac = simic_bgc_index%micbd_depoly_reac
  cascade_matrix(lid_micbd, reac)  = -1._r8
  cascade_matrix(lid_doc,reac)     =  1._r8
  cascade_matrix(lid_doc_e,reac)   =  cue_bm

  reac = simic_bgc_index%micbl_mort_reac
  cascade_matrix(lid_micbl, reac)  = -1._r8
  cascade_matrix(lid_micbd,reac)   =  f_mic2D
  cascade_matrix(lid_doc,reac)     =  f_mic2C
  cascade_matrix(lid_doc_e,reac)   =  f_mic2C*cue_bm

  reac = simic_bgc_index%o2_resp_reac
  cascade_matrix(lid_o2, reac)     = -1._r8
  cascade_matrix(lid_co2_hr,reac)  =  1._r8
  cascade_matrix(lid_co2,reac)     =  1._r8

  reac = simic_bgc_index%doc_sorb_reac
  cascade_matrix(lid_doc, reac)   = -1._r8
  cascade_matrix(lid_pom, reac)   =  1._r8

  reac = simic_bgc_index%doc_uptake_reac
  cascade_matrix(lid_doc, reac)   = -1._r8
  cascade_matrix(lid_micbl, reac) =  1._r8

  reac = simic_bgc_index%pom_desorb_reac
  cascade_matrix(lid_pom, reac)   = -1._r8
  cascade_matrix(lid_doc, reac)   =  1._r8

  end associate
  end subroutine calc_cascade_matrix

  !-------------------------------------------------------------------------------
  subroutine ode_adapt_ebbks1(me, y0, nprimeq, neq, t, dt, y)
    ! !DESCRIPTION:
    !first order implicit bkks ode integration with the adaptive time stepping
    !This could be used as an example for the implementation of time-adaptive
    !mbbks1.
    ! !NOTE:
    ! this code should only be used for mass positive ODE integration
    use ODEMOD, only : ebbks, get_rerr, get_tscal
    implicit none
    ! !ARGUMENTS:
    class(simic_bgc_type)     , intent(inout) :: me
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: t        ! time stamp
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
    real(r8), intent(out) :: y(neq)   ! updated state variable

    ! !LOCAL VARIABLES:
    real(r8) :: yc(neq)    !coarse time stepping solution
    real(r8) :: yf(neq)    !fine time stepping solution
    real(r8) :: ycp(neq)   !temporary variable
    real(r8) :: f(neq)   ! derivative
    real(r8) :: dt2
    real(r8) :: dtr
    real(r8) :: dt05
    real(r8) :: dtmin
    real(r8) :: tt,tt2     !temporary variables
    logical  :: acc
    real(r8) :: rerr, dt_scal, pscal
    integer  :: n, nJ

    dt2=dt
    dtmin=dt/64._r8
    dtr=dt
    tt=0._r8
    !make a copy of the solution at the current time step
    y=y0
    do
       if(dt2<=dtmin)then
          call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)
          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)
          dtr=dtr-dt2
          tt=tt+dt2
          y=yc
       else
          !get coarse grid solution'
          call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)

          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)

          !get fine grid solution'
          dt05=dt2*0.5_r8
          call ebbks(y,f,nprimeq, neq,dt05, yf, pscal)
          tt2=tt+dt05
          ycp=yf
          call me%bgc_integrate(ycp, dt05, tt, nprimeq, neq, f)
          call ebbks(ycp,f,nprimeq, neq,dt05,yf,pscal)

          !determine the relative error'
          rerr=get_rerr(yc,yf, neq)*exp(1._r8-1._r8/(pscal+1.e-20))

          !determine time scalar factor'
          call get_tscal(rerr,dt_scal,acc)

          if(acc)then
             dtr=dtr-dt2
             tt=tt+dt2
             y=yf
          endif
          dt2=dt2*dt_scal
          dt2=min(dt2,dtr)
       endif
       if(abs(dtr/dt)<1.e-4_r8)exit
    enddo

  end subroutine ode_adapt_ebbks1

  !--------------------------------------------------------------------
  subroutine bgc_integrate(this, ystate, dtime, time, nprimvars, nstvars, dydt)

  use SOMStateVarUpdateMod , only : calc_dtrend_som_bgc
  use MathfuncMod          , only : lom_type, safe_div, dot_sum
  implicit none
  class(simic_bgc_type)     , intent(inout) :: this
  integer                   , intent(in) :: nstvars
  integer                   , intent(in) :: nprimvars
  real(r8)                  , intent(in) :: dtime
  real(r8)                  , intent(in) :: time
  real(r8)                  , intent(in) :: ystate(nstvars)
  real(r8)                  , intent(out) :: dydt(nstvars)

  integer  :: jj, it
  integer, parameter  :: itmax = 10
  type(lom_type) :: lom
  type(betr_status_type) :: bstatus
  logical :: lneg
  real(r8) :: rscal(1:this%simic_bgc_index%nreactions)
  real(r8) :: rrates(1:this%simic_bgc_index%nreactions)
  real(r8) :: p_dt(1:nprimvars)
  real(r8) :: d_dt(1:nprimvars)
  real(r8) :: pscal(1:nprimvars)
  real(r8) :: doc_cue, pom_cue
  integer  :: reac
  associate(                                                  &
    lid_doc_e        => this%simic_bgc_index%lid_doc_e      , &
    lid_doc          => this%simic_bgc_index%lid_doc        , &
    lid_pom_e        => this%simic_bgc_index%lid_pom_e      , &
    lid_pom          => this%simic_bgc_index%lid_pom        , &
    doc_uptake_reac  => this%simic_bgc_index%doc_uptake_reac, &
    lid_cum_closs        => this%simic_bgc_index%lid_cum_closs      , &
    lid_co2_hr       => this%simic_bgc_index%lid_co2_hr     , &
    nreactions       => this%simic_bgc_index%nreactions       &
  )
  doc_cue = safe_div(ystate(lid_doc_e), ystate(lid_doc))
  pom_cue = safe_div(ystate(lid_pom_e), ystate(lid_pom))
  !print*,'dcue',doc_cue,ystate(lid_doc_e), ystate(lid_doc)
  !print*,'pcue',pom_cue,ystate(lid_pom_e), ystate(lid_pom)
  if(pom_cue>1._r8)then
    print*,'error in pom_cue',pom_cue,ystate(lid_pom_e), ystate(lid_pom)
    stop
  endif
  call correct_cascade_matrix_doc(doc_cue, pom_cue)
  call this%simic_rrates(this%simic_bgc_index, dtime, nstvars, ystate, doc_cue, rrates)

  it=0
  rscal=0._r8
  do
    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixp(1:nprimvars, 1:nreactions), rrates, p_dt)

    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixd(1:nprimvars, 1:nreactions), rrates, d_dt)

    call lom%calc_state_pscal(nprimvars, dtime, ystate(1:nprimvars), p_dt,  d_dt, pscal, lneg, bstatus)

    if(lneg .and. it<=itmax)then
      call lom%calc_reaction_rscal(nprimvars, nreactions,  pscal, &
        this%cascade_matrixd(1:nprimvars, 1:nreactions),rscal, bstatus)
      call lom%apply_reaction_rscal(nreactions, rscal(1:nreactions), rrates(1:nreactions))
    else
      call calc_dtrend_som_bgc(nstvars, nreactions, this%cascade_matrix(1:nstvars, 1:nreactions), &
         rrates(1:nreactions), dydt)
      exit
    endif
    it = it + 1
  enddo

  if(this%batch_mode)dydt(lid_cum_closs)=dydt(lid_co2_hr)

  end associate
  contains
    subroutine correct_cascade_matrix_doc(doc_cue, pom_cue)
    implicit none
    real(r8), intent(in) :: doc_cue
    real(r8), intent(in) :: pom_cue
    integer :: reac
    associate(                                                &
      lid_doc_e        => this%simic_bgc_index%lid_doc_e    , &
      lid_pom_e        => this%simic_bgc_index%lid_pom_e      &
     )
    reac = this%simic_bgc_index%doc_sorb_reac
    this%cascade_matrix(lid_doc_e, reac)  = -doc_cue
    this%cascade_matrix(lid_pom_e, reac)  =  doc_cue
    this%cascade_matrixd(lid_doc_e, reac) = -doc_cue
    this%cascade_matrixp(lid_pom_e, reac) =  doc_cue

    reac = this%simic_bgc_index%doc_uptake_reac
    this%cascade_matrix(lid_doc_e, reac) = -doc_cue
    this%cascade_matrixd(lid_doc_e, reac)= -doc_cue

    reac = this%simic_bgc_index%pom_desorb_reac
    this%cascade_matrix(lid_pom_e, reac) = -pom_cue
    this%cascade_matrix(lid_doc_e, reac) =  pom_cue
    this%cascade_matrixd(lid_pom_e, reac)= -pom_cue
    this%cascade_matrixp(lid_doc_e, reac)=  pom_cue

    end associate
    end subroutine correct_cascade_matrix_doc
  end subroutine bgc_integrate

  !-------------------------------------------------------------------------------
  subroutine init_cold_simic(this, nstvars, ystates)

  implicit none
  class(simic_bgc_type)     , intent(inout) :: this
  integer                   , intent(in)    :: nstvars
  real(r8)                  , intent(inout) :: ystates(nstvars)

  associate(                                &
    lid_micbl => this%simic_bgc_index%lid_micbl &
  )
  ystates(lid_micbl) = 0.01_r8

  end associate
  end subroutine init_cold_simic

  !-------------------------------------------------------------------------------
  subroutine sum_tot_store(this, nstvars, ystates)
  implicit none
  class(simic_bgc_type)     , intent(inout) :: this
  integer                   , intent(in)    :: nstvars
  real(r8)                  , intent(inout) :: ystates(nstvars)
  associate(                  &
    lid_totstore => this%simic_bgc_index%lid_totstore       , &
    lid_doc          => this%simic_bgc_index%lid_doc        , &
    lid_pom          => this%simic_bgc_index%lid_pom        , &
    lit1       => this%simic_bgc_index%lit1        , &
    lit2       => this%simic_bgc_index%lit2        , &
    lit3       => this%simic_bgc_index%lit3        , &
    cwd        => this%simic_bgc_index%cwd         , &
    lid_micbl  => this%simic_bgc_index%lid_micbl   , &
    lid_micbd  => this%simic_bgc_index%lid_micbd     &

  )
  ystates(lid_totstore) =ystates(lit1) + ystates(lit2)  + ystates(lit3) + ystates(cwd) +&
    ystates(lid_doc) + ystates(lid_pom) + ystates(lid_micbl) + ystates(lid_micbd)
  end associate
  end subroutine sum_tot_store
end module simicBGCType
