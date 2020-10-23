module ecosysBGCType
#include "bshr_assert.h"
  !
  ! !DESCRIPTION:

  ! !USES:
  use bshr_kind_mod             , only : r8 => shr_kind_r8
  use bshr_log_mod              , only : errMsg => shr_log_errMsg
  use betr_varcon               , only : spval => bspval
  use betr_ctrl                 , only : spinup_state => betr_spinup_state
  use gbetrType                 , only : gbetr_type
  use BiogeoConType             , only : BiogeoCon_type
  use ecosysParaType             , only : ecosys_para_type
  use BetrStatusType            , only : betr_status_type
  use BeTRJarModel              , only : jar_model_type
  use ecosysBGCIndexType            , only : ecosys_bgc_index_type
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, extends(jar_model_type), public :: ecosys_bgc_type
  type(ecosys_bgc_index_type),     private :: ecosys_index
    real(r8), pointer                    :: ystates0(:)
    real(r8), pointer                    :: ystates1(:)
    real(r8), pointer                    :: scal_f(:)
    real(r8), pointer                    :: conv_f(:)
    real(r8), pointer                    :: conc_f(:)
    logical , private                    :: use_c13
    logical , private                    :: use_c14
    logical , private                    :: batch_mode

    !declare parameters below

    real(r8), pointer :: cascade_matrix(:,:)
    real(r8), pointer :: cascade_matrixd(:,:)
    real(r8), pointer :: cascade_matrixp(:,:)
    real(r8) :: o2_w2b
  contains
    procedure, public  :: init          => init_ecosys
    procedure, public  :: runbgc        => runbgc_ecosys
    procedure, public  :: UpdateParas   => UpdateParas_ecosys
    procedure, public  :: getvarllen    => getvarllen_ecosys
    procedure, public  :: getvarlist    => getvarlist_ecosys
    procedure, public  :: init_cold     => init_cold_ecosys
    procedure, private :: arenchyma_gas_transport
    procedure, private :: init_states
    procedure, private :: add_ext_input
    procedure, private :: InitAllocate
    procedure, private :: ecosys_rrates
    procedure, private :: calc_cascade_matrix
    procedure, private :: bgc_integrate
    procedure, private :: ode_adapt_ebbks1
  end type ecosys_bgc_type
  public :: create_jarmodel_ecosysbgc
contains

  function create_jarmodel_ecosysbgc()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(ecosys_bgc_type), pointer :: create_jarmodel_ecosysbgc
    class(ecosys_bgc_type), pointer :: bgc

    allocate(bgc)
    create_jarmodel_ecosysbgc => bgc

  end function create_jarmodel_ecosysbgc
  !-------------------------------------------------------------------------------
  function getvarllen_ecosys(this)result(ans)

  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  integer :: ans
  ans =  this%ecosys_index%nstvars

  end function getvarllen_ecosys
  !-------------------------------------------------------------------------------
  subroutine getvarlist_ecosys(this, nstvars, varnames, varunits, vartypes)
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  integer, intent(in) :: nstvars
  character(len=*), intent(out) :: varnames(1:nstvars)
  character(len=*), intent(out) :: varunits(1:nstvars)
  integer         , intent(out) :: vartypes(1:nstvars)
  !local variables
  integer :: n

  do n = 1, nstvars
    vartypes(n) = this%ecosys_index%vartypes(n)
    write(varnames(n),'(A)')trim(this%ecosys_index%varnames(n))
    write(varunits(n),'(A)')trim(this%ecosys_index%varunits(n))
  enddo
  end subroutine getvarlist_ecosys

  !-------------------------------------------------------------------------------
  subroutine UpdateParas_ecosys(this,  biogeo_con, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)      , intent(out)   :: bstatus
  !local variables
  integer :: sr
  character(len=256) :: msg
  call bstatus%reset()

  select type(biogeo_con)
  type is(ecosys_para_type)
    !pass in parameter values
  class default
    write(msg,'(A)')'Wrong parameter type passed in for UpdateParas in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine UpdateParas_ecosys
  !-------------------------------------------------------------------------------

  subroutine init_ecosys(this,  biogeo_con, batch_mode, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  logical                   , intent(in) :: batch_mode
  type(betr_status_type)      , intent(out) :: bstatus
  !local variables
  character(len=256) :: msg
  write(this%jarname, '(A)')'ecosys'

  this%batch_mode=batch_mode
  select type(biogeo_con)
  type is(ecosys_para_type)
    call bstatus%reset()
    call this%ecosys_index%Init(biogeo_con%use_c13, biogeo_con%use_c14, &
     biogeo_con%non_limit, biogeo_con%nop_limit, betr_maxpatch_pft, this%batch_mode)
    call this%UpdateParas(biogeo_con, bstatus)
    if(bstatus%check_status())return

    this%use_c13 = biogeo_con%use_c13
    this%use_c14 = biogeo_con%use_c14
  class default
    call bstatus%reset()
    write(msg,'(A)')'Wrong parameter type passed in for init_ecosys in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select

  call this%InitAllocate(this%ecosys_index)
  end subroutine init_ecosys
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, ecosys_index)

  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(ecosys_bgc_type)   , intent(inout) :: this
  type(ecosys_bgc_index_type)  , intent(in) :: ecosys_index

  associate(                                &
    nstvars => ecosys_index%nstvars        , &
    nprimvars=> ecosys_index%nprimvars     , &
    nreactions => ecosys_index%nreactions    &
  )

  allocate(this%ystates0(nstvars)); this%ystates0(:) = 0._r8
  allocate(this%ystates1(nstvars)); this%ystates1(:) = 0._r8
  allocate(this%scal_f(nprimvars));  this%scal_f(:) = 0._r8
  allocate(this%conv_f(nprimvars));  this%conv_f(:) = 0._r8
  allocate(this%conc_f(nprimvars));  this%conc_f(:) = 0._r8
  allocate(this%cascade_matrix(nstvars, nreactions)); this%cascade_matrix(:,:) = 0._r8
  allocate(this%cascade_matrixd(1:nprimvars, 1:nreactions)); this%cascade_matrixd(:,:) = 0._r8
  allocate(this%cascade_matrixp(1:nprimvars, 1:nreactions)); this%cascade_matrixp(:,:) = 0._r8
  end associate
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------
  subroutine runbgc_ecosys(this,  is_surflit, dtime, bgc_forc, nstates, ystates0, ystatesf, bstatus)

  !DESCRIPTION
  !do bgc model integration for one step
  use JarBgcForcType        , only : JarBGC_forc_type
  use MathfuncMod           , only : pd_decomp
  use BetrStatusType        , only : betr_status_type
  use MathfuncMod           , only : safe_div
  use tracer_varcon         , only : catomw, natomw, patomw
  use MathfuncMod           , only : pd_decomp
  implicit none
  class(ecosys_bgc_type)  , intent(inout) :: this
  logical                    , intent(in)    :: is_surflit
  real(r8)                   , intent(in)    :: dtime
  type(JarBGC_forc_type)     , intent(in)    :: bgc_forc
  integer                    , intent(in)    :: nstates
  real(r8)                   , intent(out)   :: ystates0(nstates)
  real(r8)                   , intent(out)   :: ystatesf(nstates)
  type(betr_status_type)     , intent(out)   :: bstatus

  !local variables
  real(r8)               :: time = 0._r8
  real(r8) :: yf(this%ecosys_index%nstvars)
  character(len=*),parameter :: subname = 'runbgc_ecosys'

  associate(                                         &
    ystates1       => this%ystates1                , &
    nstvars        => this%ecosys_index%nstvars     , &
    nreactions     => this%ecosys_index%nreactions  , &
    nprimvars      => this%ecosys_index%nprimvars   , &
    cascade_matrix => this%cascade_matrix          , &
    cascade_matrixp=> this%cascade_matrixp         , &
    cascade_matrixd=> this%cascade_matrixd           &
  )

  call bstatus%reset()
  !initialize state variables
  call this%init_states(this%ecosys_index, bgc_forc)

  ystates0(:) = this%ystates0(:)
  call this%add_ext_input(dtime, this%ecosys_index, bgc_forc)
  call this%arenchyma_gas_transport(this%ecosys_index, dtime)
  call this%calc_cascade_matrix(this%ecosys_index,  cascade_matrix)
  call pd_decomp(nprimvars, nreactions, cascade_matrix(1:nprimvars, 1:nreactions), &
     cascade_matrixp, cascade_matrixd, bstatus)
  if(bstatus%check_status())return

  time = 0._r8
  yf(:) = ystates1(:)
  call this%ode_adapt_ebbks1(yf, nprimvars, nstvars, time, dtime, ystates1)

  ystatesf(:) = ystates1(:)

  end associate
  end subroutine runbgc_ecosys
  !-------------------------------------------------------------------------------
  subroutine ecosys_rrates(this, ecosys_index, dtime,  nstates, ystates1, doc_cue, rrates)
  !
  !DESCRIPTION
  !calculate reaction rates, this subroutine should be customized
  implicit none
  class(ecosys_bgc_type)  , intent(inout) :: this
  real(r8)               , intent(in)    :: dtime
  integer                , intent(in)    :: nstates
  real(r8)               , intent(in)    :: ystates1(nstates)
  type(ecosys_bgc_index_type) , intent(in)    :: ecosys_index
  real(r8)               , intent(in)    :: doc_cue
  real(r8)               , intent(out)   :: rrates(ecosys_index%nreactions)
  !local variables
  integer  :: jj
  associate(                                            &
    lit1             => ecosys_index%lit1              , &
    lit2             => ecosys_index%lit2              , &
    lit3             => ecosys_index%lit3              , &
    cwd              => ecosys_index%cwd               , &
    lid_dom          => ecosys_index%lid_dom           , &
    lid_o2           => ecosys_index%lid_o2            , &
    lid_co2          => ecosys_index%lid_co2           , &
    lit1_depoly_reac => ecosys_index%lit1_depoly_reac  , &
    lit2_depoly_reac => ecosys_index%lit2_depoly_reac  , &
    lit3_depoly_reac => ecosys_index%lit3_depoly_reac  , &
    cwd_depoly_reac  => ecosys_index%cwd_depoly_reac     &
 )

  !assemble the derivatives
  rrates(:) = 0._r8
  end associate
  end subroutine ecosys_rrates
  !-------------------------------------------------------------------------------
  subroutine init_states(this, ecosys_index, bgc_forc)

  use ecosysBGCIndexType            , only : ecosys_bgc_index_type
  use JarBgcForcType            , only : JarBGC_forc_type
  implicit none
  class(ecosys_bgc_type)  , intent(inout) :: this
  type(ecosys_bgc_index_type)  , intent(in) :: ecosys_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  associate(                               &
    lid_co2_hr  => ecosys_index%lid_co2_hr, &
    lid_n2 => ecosys_index%lid_n2,   &
    lid_o2 => ecosys_index%lid_o2,   &
    lid_co2 => ecosys_index%lid_co2, &
    lid_c13_co2 => ecosys_index%lid_c13_co2, &
    lid_c14_co2 => ecosys_index%lid_c14_co2, &
    lid_ar => ecosys_index%lid_ar,   &
    lid_o2_paere => ecosys_index%lid_o2_paere, &
    lid_n2_paere => ecosys_index%lid_n2_paere, &
    lid_ar_paere => ecosys_index%lid_ar_paere, &
    lid_co2_paere => ecosys_index%lid_co2_paere, &
    lid_ch4_paere => ecosys_index%lid_ch4_paere, &
    lid_c13_co2_paere => ecosys_index%lid_c13_co2_paere, &
    lid_c14_co2_paere => ecosys_index%lid_c14_co2_paere, &
    lid_ch4 => ecosys_index%lid_ch4  &
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
  subroutine add_ext_input(this, dtime, ecosys_index, bgc_forc)
  !
  !DESCRIPTION
  use ecosysBGCIndexType            , only : ecosys_bgc_index_type
  use JarBgcForcType        , only : JarBGC_forc_type
  use tracer_varcon             , only : catomw, natomw, patomw,c13atomw,c14atomw
  use MathfuncMod               , only : safe_div
  implicit none
  class(ecosys_bgc_type)  , intent(inout) :: this
  real(r8), intent(in) :: dtime
  type(ecosys_bgc_index_type)  , intent(in) :: ecosys_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  associate(                        &
    lit1 =>  ecosys_index%lit1, &
    lit2 =>  ecosys_index%lit2, &
    lit3 =>  ecosys_index%lit3, &
    cwd =>   ecosys_index%cwd   &
   )

  this%ystates1(lit1)=this%ystates0(lit1) +  dtime *  bgc_forc%cflx_input_litr_met
  this%ystates1(lit2)=this%ystates0(lit2) +  dtime *  bgc_forc%cflx_input_litr_cel
  this%ystates1(lit3)=this%ystates0(lit3) +  dtime *  bgc_forc%cflx_input_litr_lig
  this%ystates1(cwd)=this%ystates0(cwd) +  dtime *  bgc_forc%cflx_input_litr_cwd

  end associate
  end subroutine add_ext_input

  !--------------------------------------------------------------------
  subroutine arenchyma_gas_transport(this, ecosys_index, dtime)
  use ecosysBGCIndexType       , only : ecosys_bgc_index_type
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  type(ecosys_bgc_index_type)  , intent(in) :: ecosys_index
  real(r8), intent(in) :: dtime

  !local variables
  integer :: j
  real(r8) :: y0
  associate(                                            &
    lid_n2            => ecosys_index%lid_n2           , &
    lid_o2            => ecosys_index%lid_o2           , &
    lid_co2           => ecosys_index%lid_co2          , &
    lid_c13_co2       => ecosys_index%lid_c13_co2      , &
    lid_c14_co2       => ecosys_index%lid_c14_co2      , &
    lid_ar            => ecosys_index%lid_ar           , &
    lid_o2_paere      => ecosys_index%lid_o2_paere     , &
    lid_n2_paere      => ecosys_index%lid_n2_paere     , &
    lid_ar_paere      => ecosys_index%lid_ar_paere     , &
    lid_co2_paere     => ecosys_index%lid_co2_paere    , &
    lid_ch4_paere     => ecosys_index%lid_ch4_paere    , &
    lid_c13_co2_paere => ecosys_index%lid_c13_co2_paere, &
    lid_c14_co2_paere => ecosys_index%lid_c14_co2_paere, &
    lid_ch4           => ecosys_index%lid_ch4            &
  )

  j = lid_o2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(ecosys_index%lid_o2_paere) = this%ystates1(j)-y0

  if( spinup_state == 0)then
    j = lid_n2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))

    this%ystates1(ecosys_index%lid_n2_paere) = this%ystates1(j)-y0
    j = lid_ar; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(ecosys_index%lid_ar_paere) = this%ystates1(j)-y0

    j = lid_ch4; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(ecosys_index%lid_ch4_paere) = this%ystates1(j)-y0

    j = lid_co2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(ecosys_index%lid_co2_paere) = this%ystates1(j)-y0

    if(this%use_c13)then
      j = lid_c13_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(ecosys_index%lid_c13_co2_paere) = this%ystates1(j)-y0
    endif

    if(this%use_c14)then
      j = lid_c14_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(ecosys_index%lid_c14_co2_paere) = this%ystates1(j)-y0
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
  subroutine calc_cascade_matrix(this, ecosys_index,  cascade_matrix)
  !
  !DESCRPTION
  !compute the cascade matrix
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  type(ecosys_bgc_index_type)    , intent(in)    :: ecosys_index
  real(r8)                  , intent(inout) :: cascade_matrix(ecosys_index%nstvars, ecosys_index%nreactions)
  !local variables
  integer :: reac

  associate(                                &
    lit1       => ecosys_index%lit1        , &
    lit2       => ecosys_index%lit2        , &
    lit3       => ecosys_index%lit3        , &
    cwd        => ecosys_index%cwd         , &
    lid_dom    => ecosys_index%lid_dom     , &
    lid_micbl  => ecosys_index%lid_micbl   , &
    lid_micbd  => ecosys_index%lid_micbd   , &
    lid_cue    => ecosys_index%lid_cue     , &
    lid_o2     => ecosys_index%lid_o2      , &
    lid_co2    => ecosys_index%lid_co2     , &
    lid_co2_hr => ecosys_index%lid_co2_hr    &
   ! cue_met    => this%cue_met            , &
   ! cue_cel    => this%cue_cel            , &
   ! cue_lig    => this%cue_lig            , &
   ! cue_cwd    => this%cue_cwd              &
   ! f_mic2d    => this%f_mic2d            , &
   ! f_mic2c    => this%f_mic2c              &
  )

  reac = ecosys_index%lit1_depoly_reac
  cascade_matrix(lit1, reac)  = -1._r8
  cascade_matrix(lid_dom,reac)=1._r8
 ! cascade_matrix(lid_cue,reac)=cue_met

  reac = ecosys_index%lit2_depoly_reac
  cascade_matrix(lit2, reac)  = -1._r8
  cascade_matrix(lid_dom,reac)=  1._r8
 ! cascade_matrix(lid_cue,reac)= cue_cel

  reac = ecosys_index%lit3_depoly_reac
  cascade_matrix(lit3, reac)  = -1._r8
  cascade_matrix(lid_dom,reac)=1._r8
 ! cascade_matrix(lid_cue,reac)=cue_lig

 ! reac = ecosys_index%cwd_depoly_reac
 ! cascade_matrix(cwd, reac)   = -1._r8
 ! cascade_matrix(lid_dom,reac)=1._r8
 ! cascade_matrix(lid_cue,reac)=cue_cwd

 ! reac = ecosys_index%micbd_depoly_reac
 ! cascade_matrix(lid_micbd, reac)  = -1._r8
 ! cascade_matrix(lid_dom,reac)     = 1._r8
 ! cascade_matrix(lid_cue,reac)     = 0.5_r8

 ! reac = ecosys_index%micbl_mort_reac
 ! cascade_matrix(lid_micbl, reac)  = -1._r8
 ! cascade_matrix(lid_micbd,reac)   =  f_mic2D
 ! cascade_matrix(lid_dom,reac)     =  f_mic2C

  reac = ecosys_index%o2_resp_reac
  cascade_matrix(lid_o2, reac)     = -1._r8
  cascade_matrix(lid_co2_hr,reac)  =  1._r8
  cascade_matrix(lid_co2,reac)  =  1._r8

 ! reac = ecosys_index%doc_uptake_reac
 ! cascade_matrix(lid_doc, reac)   = -1._r8
 ! cascade_matrix(lid_micbl, reac) = 1._r8
  end associate
  end subroutine calc_cascade_matrix

  !-------------------------------------------------------------------------------
  subroutine ode_adapt_ebbks1(me, y0, nprimeq, neq, t, dt, y)
    !
    !DESCRIPTION:
    !first order implicit bkks ode integration with the adaptive time stepping
    !This could be used as an example for the implementation of time-adaptive
    !mbbks1.
    ! !NOTE:
    ! this code should only be used for mass positive ODE integration
    use ODEMOD, only : ebbks, get_rerr, get_tscal
    implicit none
    ! !ARGUMENTS:
    class(ecosys_bgc_type)     , intent(inout) :: me
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: t        ! time stamp
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  ! number of primary variables subject to positivity constraint
    real(r8), intent(out) :: y(neq)   ! updated state variable

    ! !LOCAL VARIABLES:
    real(r8) :: yc(neq)    !coarse time stepping solution
    real(r8) :: yf(neq)    !fine time stepping solution
    real(r8) :: ycp(neq)   !temporary variable
    real(r8) :: f(neq)     ! derivative
    real(r8) :: dt2
    real(r8) :: dtr
    real(r8) :: dt05
    real(r8) :: dtmin
    real(r8) :: tt,tt2     !temporary variables
    logical  :: acc
    real(r8) :: rerr, dt_scal, pscal
    integer  :: n, nJ
    real(r8), parameter :: maxtdiv=64._r8

    dt2=dt
    dtmin=dt/maxtdiv
    dtr=dt
    tt=0._r8
    !make a copy of the solution at the current time step
    y(:)=y0(:)
    do
       if(dt2<=dtmin)then
         call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)
         call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)
         dtr=dtr-dt2
         tt=tt+dt2
         y=yc
       else
         !get coarse grid solution
         call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)
         call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)
         !get fine grid solution
         dt05=dt2*0.5_r8
         call ebbks(y,f,nprimeq, neq,dt05, yf, pscal)
         tt2=tt+dt05
         ycp=yf
         call me%bgc_integrate(ycp, dt05, tt, nprimeq, neq, f)
         call ebbks(ycp,f,nprimeq, neq,dt05,yf,pscal)
         !determine the relative error
         rerr=get_rerr(yc,yf, neq)*exp(1._r8-1._r8/(pscal+1.e-20))
         !determine time scalar factor
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
  !
  !DESCRPTION
  !core code for integrating the bgc system
  use SOMStateVarUpdateMod , only : calc_dtrend_som_bgc
  use MathfuncMod          , only : lom_type, safe_div
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  integer                   , intent(in) :: nstvars
  integer                   , intent(in) :: nprimvars
  real(r8)                  , intent(in) :: dtime
  real(r8)                  , intent(in) :: time
  real(r8)                  , intent(in) :: ystate(nstvars)
  real(r8)                  , intent(out) :: dydt(nstvars)
  !local variables
  integer  :: jj, it
  integer, parameter  :: itmax = 10
  type(lom_type) :: lom
  type(betr_status_type) :: bstatus
  logical :: lneg
  real(r8) :: rscal(1:this%ecosys_index%nreactions)
  real(r8) :: rrates(1:nstvars)
  real(r8) :: p_dt(1:nprimvars)
  real(r8) :: d_dt(1:nprimvars)
  real(r8) :: pscal(1:nprimvars)
  !real(r8) :: doc_cue

  associate(                                            &
    !lid_cue        => this%ecosys_index%lid_cue        , &
    lid_dom        => this%ecosys_index%lid_dom        , &
    doc_uptake_reac=> this%ecosys_index%doc_uptake_reac, &
    nreactions => this%ecosys_index%nreactions           &
  )

  !customize the following lines
  !doc_cue = safe_div(ystate(lid_cue), ystate(lid_doc))
  !call this%ecosys_rrates(this%ecosys_index, dtime, nstvars, ystate, doc_cue, rrates)
  !this%cascade_matrixd(lid_cue, doc_uptake_reac) =  doc_cue
  !this%cascade_matrix(lid_cue, doc_uptake_reac)  = -doc_cue

  it=0
  rscal=0._r8
  do
    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixp(1:nprimvars, 1:nreactions), rrates, p_dt)
    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixd(1:nprimvars, 1:nreactions), rrates, d_dt)
    !update the state variables
    call lom%calc_state_pscal(nprimvars, dtime, ystate(1:nprimvars), p_dt(1:nprimvars),  d_dt(1:nprimvars), &
        pscal(1:nprimvars), lneg, bstatus)
    if(lneg .and. it<=itmax)then
      call lom%calc_reaction_rscal(nprimvars, nreactions,  pscal(1:nprimvars), &
        this%cascade_matrixd(1:nprimvars, 1:nreactions),rscal, bstatus)
      call lom%apply_reaction_rscal(nreactions, rscal(1:nreactions), rrates(1:nreactions))
    else
      call calc_dtrend_som_bgc(nstvars, nreactions, this%cascade_matrix(1:nstvars, 1:nreactions), &
         rrates(1:nreactions), dydt)
      exit
    endif
    it = it + 1
  enddo
  end associate
  end subroutine bgc_integrate
  !-------------------------------------------------------------------------------
  subroutine init_cold_ecosys(this, nstvars, ystates)
  !
  !DESCRPTION
  !do a cold state initialization for batch mode simulation
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  integer                   , intent(in)    :: nstvars
  real(r8)                  , intent(inout) :: ystates(nstvars)

  !Initialize necessary state variables below
  end subroutine init_cold_ecosys
  end module ecosysBGCType
