module ecosysBGCType
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
  use BetrStatusType            , only : betr_status_type
  use BeTRJarModel              , only : jar_model_type
  use ecosysBGCIndexType       , only : ecosys_bgc_index_type
  use ecosysParaType           , only : ecosys_para_type
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !Note:
  !Keeping ecosys_bgc_index as a private member is a workaround to call the ode solver
  !it increase the memory for each instance of ecosys_bgc_type, but enables
  !the ode function to be called by the ode solver

  type, extends(jar_model_type), public :: ecosys_bgc_type
    type(ecosys_bgc_index_type), private :: ecosys_bgc_index
    real(r8), pointer                    :: ystates0(:)
    real(r8), pointer                    :: ystates1(:)
    real(r8), pointer                    :: k_decay(:)
    real(r8), pointer                    :: cascade_matrix(:,:)
    real(r8), pointer                    :: cascade_matrixd(:,:)
    real(r8), pointer                    :: cascade_matrixp(:,:)
    real(r8), pointer                    :: alpha_n(:)
    real(r8), pointer                    :: alpha_p(:)
    real(r8)                             :: pot_f_nit
    real(r8)                             :: pot_f_denit
    real(r8)                             :: rt_ar
    real(r8), pointer                    :: frac_p_sec_to_sol(:)
    real(r8), pointer                    :: minp_secondary_decay(:)
    real(r8), pointer                    :: mumax_minp_soluble_to_secondary(:)
    integer                              :: plant_ntypes
    logical                              :: nop_limit
    logical                              :: non_limit
    real(r8), pointer                    :: scal_f(:)
    real(r8), pointer                    :: conv_f(:)
    real(r8), pointer                    :: conc_f(:)
    real(r8), pointer                    :: cascade_matnh4(:)
    integer                              :: soilorder
    real(r8)                             :: msurf_nh4
    real(r8)                             :: msurf_minp
    real(r8), private                    :: c14decay_const
    real(r8), private                    :: c14decay_som_const
    real(r8), private                    :: c14decay_dom_const
    real(r8), private                    :: c14decay_pom_const
    real(r8), private                    :: c14decay_Bm_const
    logical , private                    :: use_c13
    logical , private                    :: use_c14
    real(r8), private                    :: beg_c_mass, beg_c13_mass, beg_c14_mass
    real(r8), private                    :: beg_n_mass
    real(r8), private                    :: beg_p_mass
    real(r8), private                    :: c_inflx,n_inflx, p_inflx
    logical                              :: bgc_on
    logical , private                    :: batch_mode
  contains
    procedure, public  :: init          => init_ecosys
    procedure, public  :: initVL        => initVL_ecosys
    procedure, public  :: runbgc        => runbgc_ecosys
    procedure, public  :: UpdateParas   => UpdateParas_ecosys
    procedure, public  :: sumup_cnp_msflx => sumup_cnp_msflx_ecosys
    procedure, public  :: getvarllen    => getvarllen_ecosys
    procedure, public  :: getvarlist    => getvarlist_ecosys
    procedure, public  :: init_cold     => init_cold_ecosys
    procedure, private :: calc_cascade_matrix
    procedure, private :: init_states
    procedure, private :: add_ext_input
    procedure, private :: InitAllocate
    procedure, private :: arenchyma_gas_transport
    procedure, private :: sumup_cnp_mass
    procedure, private :: bgc_integrate
    procedure, private :: c14decay
    procedure, private :: checksum_cascade
    procedure, private :: begin_massbal_check
    procedure, private :: end_massbal_check
    procedure, private :: sum_tot_store
    procedure, public  :: display_index
  end type ecosys_bgc_type

  public :: create_jarmodel_ecosys
contains

  function create_jarmodel_ecosys()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(ecosys_bgc_type), pointer :: create_jarmodel_ecosys
    class(ecosys_bgc_type), pointer :: bgc

    allocate(bgc)
    create_jarmodel_ecosys => bgc

  end function create_jarmodel_ecosys

  !-------------------------------------------------------------------------------
  function getvarllen_ecosys(this)result(ans)

  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  integer :: ans

  ans = this%ecosys_bgc_index%nstvars

  end function getvarllen_ecosys
  !-------------------------------------------------------------------------------
  subroutine getvarlist_ecosys(this, nstvars, varnames, varunits, vartypes)
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  integer, intent(in) :: nstvars
  character(len=*), intent(out) :: varnames(1:nstvars)
  character(len=*), intent(out) :: varunits(1:nstvars)
  integer         , intent(out) :: vartypes(1:nstvars)
  integer :: n

  do n = 1, nstvars
    vartypes(n) = this%ecosys_bgc_index%vartypes(n)
    write(varnames(n),'(A)')trim(this%ecosys_bgc_index%varnames(n))
    write(varunits(n),'(A)')trim(this%ecosys_bgc_index%varunits(n))
  enddo
  end subroutine getvarlist_ecosys

  !-------------------------------------------------------------------------------
  subroutine begin_massbal_check(this)

  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  real(r8) :: c_mass0, n_mass0, p_mass0
  print*,'begin_massbal_check'
  call this%sumup_cnp_msflx(this%ystates1, this%beg_c_mass, &
    this%beg_n_mass, this%beg_p_mass)

  print*,this%beg_c_mass,this%beg_n_mass,this%beg_p_mass

  call this%sumup_cnp_msflx(this%ystates0, c_mass0, &
    n_mass0, p_mass0)
  print*,c_mass0,n_mass0,p_mass0
  end subroutine begin_massbal_check
  !-------------------------------------------------------------------------------

  subroutine end_massbal_check(this, header)

  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  character(len=*), intent(in) :: header
  real(r8) :: c_mass, n_mass, p_mass
  real(r8) :: c_flx,n_flx,p_flx
  type(betr_status_type) :: bstatus
  real(r8) :: dmass_c, dmass_n, dmass_p
  real(r8) :: emass_c, emass_n, emass_p
  real(r8) :: remass_c, remass_n, remass_p

  real(r8), parameter :: tiny_val=1.e-10_r8
  print*,trim(header)
  call this%sumup_cnp_msflx(this%ystates1, c_mass, n_mass, p_mass,c_flx,n_flx,p_flx, bstatus)

  dmass_c=c_mass - this%beg_c_mass
  dmass_n=n_mass - this%beg_n_mass
  dmass_p=p_mass - this%beg_p_mass

  emass_c=dmass_c + c_flx - this%c_inflx
  emass_n=dmass_n + n_flx - this%n_inflx
  emass_p=dmass_p + p_flx - this%p_inflx

  remass_c=emass_c/max(abs(dmass_c),tiny_val)
  remass_n=emass_n/max(abs(dmass_n),tiny_val)
  remass_p=emass_p/max(abs(dmass_p),tiny_val)

  write(*,'(A)')'------------------------------------------------------------------------------------------------------'
  write(*,'(A)')'type             beg_mass             end_mass            dmass                inflx            outflx'
  write(*, '(A,5(X,E20.10))')'c_mass bal=',this%beg_c_mass, c_mass, emass_c, this%c_inflx, c_flx
  write(*, '(A,5(X,E20.10))')'n_mass bal=',this%beg_n_mass, n_mass, emass_n, this%n_inflx, n_flx
  write(*, '(A,5(X,E20.10))')'p_mass bal=',this%beg_p_mass, p_mass, emass_p, this%p_inflx, p_flx
  if(maxval((/abs(remass_c),abs(remass_n),abs(remass_p)/))>1.e-3_r8)stop
  end subroutine end_massbal_check

  !-------------------------------------------------------------------------------
  subroutine UpdateParas_ecosys(this,  biogeo_con, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)     , intent(out)   :: bstatus
  integer :: sr
  character(len=256) :: msg

  call bstatus%reset()

  select type(biogeo_con)
  type is(ecosys_para_type)

    if(this%use_c14)then
      this%c14decay_const=biogeo_con%c14decay_const
      this%c14decay_som_const=biogeo_con%c14decay_som_const
      this%c14decay_dom_const=biogeo_con%c14decay_dom_const
      this%c14decay_pom_const=biogeo_con%c14decay_pom_const
      this%c14decay_Bm_const=biogeo_con%c14decay_Bm_const
    endif


  class default
    write(msg,'(A)')'Wrong parameter type passed in for UpdateParas in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine UpdateParas_ecosys
  !-------------------------------------------------------------------------------
  subroutine init_ecosys(this,  biogeo_con,  batch_mode, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  logical                   , intent(in) :: batch_mode
  type(betr_status_type)    , intent(out) :: bstatus

  character(len=256) :: msg
  write(this%jarname, '(A)')'ecosys'

  this%bgc_on=.true.
  this%batch_mode=batch_mode
  select type(biogeo_con)
  type is(ecosys_para_type)
    call bstatus%reset()
    call this%ecosys_bgc_index%Init(biogeo_con%use_c13, biogeo_con%use_c14, &
     biogeo_con%non_limit, biogeo_con%nop_limit, betr_maxpatch_pft, this%batch_mode)

    this%nop_limit=biogeo_con%nop_limit
    this%non_limit=biogeo_con%non_limit

    call this%InitAllocate(this%ecosys_bgc_index)

    this%use_c13 = biogeo_con%use_c13

    this%use_c14 = biogeo_con%use_c14

    call this%UpdateParas(biogeo_con, bstatus)
    if(bstatus%check_status())return
  class default
    call bstatus%reset()
    write(msg,'(A)')'Wrong parameter type passed in for init_ecosys in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine init_ecosys
  !-------------------------------------------------------------------------------
  subroutine initVL_ecosys(this,  biogeo_con,  batch_mode, ecosys_bgc_index, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  logical                   , intent(in) :: batch_mode
  class(ecosys_bgc_index_type), intent(in) :: ecosys_bgc_index
  type(betr_status_type)    , intent(out) :: bstatus

  character(len=256) :: msg
  write(this%jarname, '(A)')'ecosys'

  this%bgc_on=.true.
  this%batch_mode=batch_mode
  select type(biogeo_con)
  type is(ecosys_para_type)
    call bstatus%reset()
    call this%ecosys_bgc_index%hcopy(ecosys_bgc_index)

    this%nop_limit=biogeo_con%nop_limit
    this%non_limit=biogeo_con%non_limit

    call this%InitAllocate(this%ecosys_bgc_index)

    if(bstatus%check_status())return

    this%use_c13 = biogeo_con%use_c13

    this%use_c14 = biogeo_con%use_c14

    call this%UpdateParas(biogeo_con, bstatus)
    if(bstatus%check_status())return
  class default
    call bstatus%reset()
    write(msg,'(A)')'Wrong parameter type passed in for init_ecosys in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine initVL_ecosys
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, ecosys_bgc_index)
  use ecosysBGCIndexType , only : ecosys_bgc_index_type
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(ecosys_bgc_type)   , intent(inout) :: this
  type(ecosys_bgc_index_type) , intent(in):: ecosys_bgc_index

  associate(                                   &
    nstvars => ecosys_bgc_index%nstvars  ,     &
    nreactions => ecosys_bgc_index%nreactions, &
    nprimvars => ecosys_bgc_index%nprimvars    &
  )

  allocate(this%ystates0(nstvars)); this%ystates0(:) = 0._r8
  allocate(this%ystates1(nstvars)); this%ystates1(:) = 0._r8
  allocate(this%scal_f(nprimvars));  this%scal_f(:) = 0._r8
  allocate(this%conv_f(nprimvars));  this%conv_f(:) = 0._r8
  allocate(this%conc_f(nprimvars));  this%conc_f(:) = 0._r8

  end associate
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------


  subroutine checksum_cascade(this, ecosys_bgc_index)

  use ecosysBGCIndexType       , only : ecosys_bgc_index_type

  implicit none
  ! !ARGUMENTS:
  class(ecosys_bgc_type)     , intent(in) :: this
  type(ecosys_bgc_index_type)   , intent(in) :: ecosys_bgc_index

  real(r8) :: resc, resn, resp
  integer  :: reac,jj

  print*,'checksum'
  print*,'som1c ','som1n ','som1p'

  end subroutine checksum_cascade

  !-------------------------------------------------------------------------------
  subroutine runbgc_ecosys(this,  is_surflit, dtime, bgc_forc, nstates, ystates0, ystatesf, bstatus)

  !DESCRIPTION
  !do bgc model integration for one step
  use JarBgcForcType        , only : JarBGC_forc_type
  use MathfuncMod           , only : pd_decomp
  use BetrStatusType        , only : betr_status_type
  use MathfuncMod           , only : safe_div
  use tracer_varcon         , only : catomw, natomw, patomw
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
  real(r8) :: pot_co2_hr
  real(r8) :: pot_f_nit_mol_per_sec
  real(r8) :: n2_n2o_ratio_denit
  real(r8) :: yf(this%ecosys_bgc_index%nstvars)
  real(r8) :: o2_decomp_depth
  real(r8) :: time
  real(r8) :: frc_c13, frc_c14
  real(r8) :: c_mass1, n_mass1, p_mass1
  real(r8) :: c_mass2, n_mass2, p_mass2
  real(r8) :: c_flx, n_flx, p_flx
  integer :: jj
  character(len=*),parameter :: subname = 'runbgc_ecosys'
  associate(                                            &
    pctsand   => bgc_forc%pct_sand             , &  !sand in %
    rt_ar     => bgc_forc%rt_ar                , &  !root autotrophic respiration
    rt_ar_c13 => bgc_forc%rt_ar_c13            , &  !root autotrophic respiration
    rt_ar_c14 => bgc_forc%rt_ar_c14            , &  !root autotrophic respiration
    lid_nh3x   => this%ecosys_bgc_index%lid_nh3x        , &  !position id of nh4
    lid_hno3   => this%ecosys_bgc_index%lid_hno3        , &  !
    lid_o2    => this%ecosys_bgc_index%lid_o2         , &  !
    nstvars   => this%ecosys_bgc_index%nstvars        , &
    nprimvars => this%ecosys_bgc_index%nprimvars      , &
    nreactions=> this%ecosys_bgc_index%nreactions     , &
    lid_plant_minn_nh4  => this%ecosys_bgc_index%lid_plant_minn_nh4       , &
    lid_plant_minn_no3  => this%ecosys_bgc_index%lid_plant_minn_no3       , &
    ystates1 => this%ystates1                           &
  )
  this%ecosys_bgc_index%debug = bgc_forc%debug
  if(this%ecosys_bgc_index%debug)print*,'enter runbgc_ecosys'
  this%rt_ar = rt_ar
  frc_c13 = safe_div(rt_ar_c13,rt_ar); frc_c14 = safe_div(rt_ar_c14,rt_ar)
  call bstatus%reset()

  !initialize state variables
  call this%init_states(this%ecosys_bgc_index, bgc_forc)

  ystates0(:) = this%ystates0(:)

  !add all external input
   call this%add_ext_input(dtime, this%ecosys_bgc_index, bgc_forc, &
      this%c_inflx, this%n_inflx, this%p_inflx)


  if(this%batch_mode)call this%sum_tot_store(nstvars, ystates1)
  ystatesf(:) = ystates1(:)
  end associate
  end subroutine runbgc_ecosys
  !-------------------------------------------------------------------------------
  subroutine c14decay(this, ecosys_bgc_index, dtime, ystates1)

  !apply c14 decay to om pools
  use ecosysBGCIndexType       , only : ecosys_bgc_index_type

  implicit none
  ! !ARGUMENTS:
  class(ecosys_bgc_type)     , intent(in) :: this
  type(ecosys_bgc_index_type)   , intent(in) :: ecosys_bgc_index
  real(r8), intent(in) :: dtime
  real(r8), intent(inout) :: ystates1(:)

  integer :: jj
  integer :: kc14


  contains
    subroutine somc14_decay(ebeg, eend, nelms, c14_loc, decay_const)
    implicit none
    integer, intent(in) :: ebeg, eend, nelms
    integer, intent(in) :: c14_loc
    real(r8),intent(in) :: decay_const
    integer :: jj, kc14
    do jj = ebeg, eend, nelms
      kc14=jj-1+c14_loc
      ystates1(kc14) = ystates1(kc14)*exp(- decay_const * dtime)
    enddo

    end subroutine somc14_decay
  end subroutine c14decay
  !-------------------------------------------------------------------------------
  subroutine calc_cascade_matrix(this,ecosys_bgc_index, cascade_matrix, frc_c13, frc_c14)
    !
    ! !DESCRIPTION:
    ! calculate cascade matrix for the decomposition model
    !
    ! !USES:
    use MathfuncMod               , only : safe_div
    use ecosysBGCIndexType       , only : ecosys_bgc_index_type
    use betr_ctrl                 , only : spinup_state => betr_spinup_state
    implicit none
    ! !ARGUMENTS:
    class(ecosys_bgc_type)     , intent(in) :: this
    type(ecosys_bgc_index_type)   , intent(in) :: ecosys_bgc_index
    real(r8)                      , intent(inout)   :: cascade_matrix(ecosys_bgc_index%nstvars, ecosys_bgc_index%nreactions)
    real(r8)                      , intent(in) :: frc_c13, frc_c14
    ! !LOCAL VARIABLES:
    real(r8) :: ftxt, f1, f2
    integer :: k, reac

  end subroutine calc_cascade_matrix
  !--------------------------------------------------------------------
  subroutine init_states(this, ecosys_bgc_index, bgc_forc)

  use ecosysBGCIndexType       , only : ecosys_bgc_index_type
  use JarBgcForcType            , only : JarBGC_forc_type
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  type(ecosys_bgc_index_type)  , intent(in) :: ecosys_bgc_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  integer :: j
  associate(                           &
    lid_n2 => ecosys_bgc_index%lid_n2, &
    lid_o2 => ecosys_bgc_index%lid_o2, &
    lid_co2x => ecosys_bgc_index%lid_co2x, &
    lid_c13_co2x => ecosys_bgc_index%lid_c13_co2x, &
    lid_c14_co2x => ecosys_bgc_index%lid_c14_co2x, &
    lid_n2o => ecosys_bgc_index%lid_n2o, &
    lid_ar => ecosys_bgc_index%lid_ar, &
    lid_ch4 => ecosys_bgc_index%lid_ch4,  &
    lid_co2_hr => ecosys_bgc_index%lid_co2_hr &
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

  this%scal_f(lid_co2x) = bgc_forc%aren_cond_co2
  this%conc_f(lid_co2x) = bgc_forc%conc_atm_co2
  this%conv_f(lid_co2x) = 1._r8/bgc_forc%co2_g2b

  if(this%use_c13)then
    this%scal_f(lid_c13_co2x) = bgc_forc%aren_cond_co2_c13
    this%conc_f(lid_c13_co2x) = bgc_forc%conc_atm_co2_c13
    this%conv_f(lid_c13_co2x) = 1._r8/bgc_forc%co2_g2b
  endif

  if(this%use_c14)then
    this%scal_f(lid_c14_co2x) = bgc_forc%aren_cond_co2_c14
    this%conc_f(lid_c14_co2x) = bgc_forc%conc_atm_co2_c14
  endif

  this%scal_f(lid_ch4) = bgc_forc%aren_cond_ch4
  this%conc_f(lid_ch4) = bgc_forc%conc_atm_ch4
  this%conv_f(lid_ch4) = 1._r8/bgc_forc%ch4_g2b

  this%scal_f(lid_n2o) = bgc_forc%aren_cond_n2o
  this%conc_f(lid_n2o) = bgc_forc%conc_atm_n2o
  this%conv_f(lid_n2o) = 1._r8/bgc_forc%n2o_g2b

  this%plant_ntypes = bgc_forc%plant_ntypes


  end associate
  end subroutine init_states
  !--------------------------------------------------------------------
  subroutine add_ext_input(this, dtime, ecosys_bgc_index, bgc_forc, c_inf, n_inf, p_inf)
  use ecosysBGCIndexType       , only : ecosys_bgc_index_type
  use JarBgcForcType        , only : JarBGC_forc_type
  use tracer_varcon             , only : catomw, natomw, patomw,c13atomw,c14atomw
  use MathfuncMod               , only : safe_div
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  real(r8), intent(in) :: dtime
  type(ecosys_bgc_index_type)  , intent(in) :: ecosys_bgc_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc
  real(r8), optional, intent(out) :: c_inf, n_inf, p_inf
  integer :: kc, kn, kp,kc13,kc14
  integer :: jj
  real(r8):: totp, pmin_frac, pmin_cleave
  real(r8):: totc, totn
  real(r8):: c_inf_loc, n_inf_loc, p_inf_loc
  associate(                        &
    lid_nh3x=> ecosys_bgc_index%lid_nh3x, &
    lid_hno3=> ecosys_bgc_index%lid_hno3  &
  )
  c_inf_loc=0._r8; n_inf_loc=0._r8; p_inf_loc=0._r8

  if(present(c_inf))c_inf=c_inf*catomw
  if(present(n_inf))n_inf=n_inf*natomw
  if(present(p_inf))p_inf=p_inf*patomw
  end associate

  end subroutine add_ext_input


  !--------------------------------------------------------------------
  subroutine bgc_integrate(this, ystate, dtime, time, nprimvars, nstvars, dydt)
  !
  !DESCRIPTION
  !calculate the reaction rates
  !In current implementation, no active adsorption of NH4 is involved.
  !The inorganic phosphorus does the transition from soluble->secondary->occlude
  use SOMStateVarUpdateMod , only : calc_dtrend_som_bgc
  use MathfuncMod          , only : lom_type, safe_div
  implicit none
  class(ecosys_bgc_type) , intent(inout) :: this
  integer                   , intent(in) :: nstvars
  integer                   , intent(in) :: nprimvars
  real(r8)                  , intent(in) :: dtime
  real(r8)                  , intent(in) :: time
  real(r8)                  , intent(in) :: ystate(nstvars)
  real(r8)                  , intent(out) :: dydt(nstvars)

  !local variables
  real(r8) :: mic_pot_nn_flx  !potential nitrogen uptake to support decomposition
  real(r8) :: mic_pot_np_flx  !potential phosphorus uptake to support decomposition
  real(r8) :: rrates(this%ecosys_bgc_index%nreactions)
  real(r8) :: p_dt(1:this%ecosys_bgc_index%nprimvars)
  real(r8) :: d_dt(1:this%ecosys_bgc_index%nprimvars)
  real(r8) :: dydt1(nstvars)
  real(r8) :: pscal(1:nprimvars)
  real(r8) :: rscal(1:this%ecosys_bgc_index%nreactions)
  real(r8) :: dminn
  integer  :: jj, it
  integer, parameter  :: itmax = 10
  type(lom_type) :: lom
  type(betr_status_type) :: bstatus
  logical :: lneg
  real(r8) :: scal

  associate(                                                                      &
    lid_nh3x => this%ecosys_bgc_index%lid_nh3x                                  , &
    lid_hno3 => this%ecosys_bgc_index%lid_hno3                                  , &
    lid_co2_hr => this%ecosys_bgc_index%lid_co2_hr                              , &
    lid_cum_closs=> this%ecosys_bgc_index%lid_cum_closs                         , &
    lid_plant_minn_no3_pft=> this%ecosys_bgc_index%lid_plant_minn_no3_pft       , &
    lid_plant_minn_nh4_pft=> this%ecosys_bgc_index%lid_plant_minn_nh4_pft       , &
    lid_plant_minp_pft=> this%ecosys_bgc_index%lid_plant_minp_pft               , &
    lid_plant_minp    => this%ecosys_bgc_index%lid_plant_minp                   , &
    lid_plant_minn_nh4 => this%ecosys_bgc_index%lid_plant_minn_nh4              , &
    lid_plant_minn_no3 => this%ecosys_bgc_index%lid_plant_minn_no3              , &
    lid_supp_minp => this%ecosys_bgc_index%lid_supp_minp                        , &
    lid_supp_minn => this%ecosys_bgc_index%lid_supp_minn                          &
  )

  dydt(:) = 0._r8
  rrates(:) = 0._r8

  if(lid_supp_minp>0)then
    !check for mineral phosphorous
  endif
  if(lid_supp_minn>0)then

  endif
  if(this%batch_mode)dydt(lid_cum_closs)=dydt(lid_co2_hr)
  end associate
  end subroutine bgc_integrate
  !--------------------------------------------------------------------
  subroutine arenchyma_gas_transport(this, ecosys_bgc_index, dtime)
  use ecosysBGCIndexType       , only : ecosys_bgc_index_type
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  type(ecosys_bgc_index_type)  , intent(in) :: ecosys_bgc_index
  real(r8), intent(in) :: dtime

  integer :: j
  real(r8) :: y0
  associate(                             &
    lid_n2 => ecosys_bgc_index%lid_n2,   &
    lid_o2 => ecosys_bgc_index%lid_o2,   &
    lid_co2x => ecosys_bgc_index%lid_co2x, &
    lid_c13_co2x => ecosys_bgc_index%lid_c13_co2x, &
    lid_c14_co2x => ecosys_bgc_index%lid_c14_co2x, &
    lid_n2o => ecosys_bgc_index%lid_n2o, &
    lid_ar => ecosys_bgc_index%lid_ar,   &
    lid_ch4 => ecosys_bgc_index%lid_ch4  &
  )

  j = lid_o2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(ecosys_bgc_index%lid_o2_paere) = this%ystates1(j)-y0

  j = lid_n2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(ecosys_bgc_index%lid_n2_paere) = this%ystates1(j)-y0

    j = lid_ar; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(ecosys_bgc_index%lid_ar_paere) = this%ystates1(j)-y0

    j = lid_ch4; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(ecosys_bgc_index%lid_ch4_paere) = this%ystates1(j)-y0

    j = lid_co2x; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(ecosys_bgc_index%lid_co2x_paere) = this%ystates1(j)-y0

    if(this%use_c13)then
      j = lid_c13_co2x; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(ecosys_bgc_index%lid_c13_co2x_paere) = this%ystates1(j)-y0
    endif

    if(this%use_c14)then
      j = lid_c14_co2x; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(ecosys_bgc_index%lid_c14_co2x_paere) = this%ystates1(j)-y0

    endif
    j = lid_n2o; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(ecosys_bgc_index%lid_n2o_paere) = this%ystates1(j)-y0
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
  subroutine sumup_cnp_msflx_ecosys(this, ystates1, c_mass, n_mass, p_mass,c_flx,n_flx,p_flx, bstatus)

  use tracer_varcon, only : catomw, natomw, patomw
  implicit none
  class(ecosys_bgc_type)     , intent(in) :: this
  real(r8), intent(in)  :: ystates1(:)
  real(r8), intent(out) :: c_mass, n_mass, p_mass
  real(r8), optional, intent(out) :: c_flx, n_flx, p_flx
  type(betr_status_type), optional, intent(out)   :: bstatus
  !local variables

  integer :: kc, kn, kp, jj
  associate(                                        &
    lid_co2_hr => this%ecosys_bgc_index%lid_co2_hr  &
  )
  if(present(bstatus)) &
    SHR_ASSERT_ALL((size(ystates1) == this%ecosys_bgc_index%nstvars), errMsg(mod_filename,__LINE__),bstatus)

  c_mass = 0._r8; n_mass = 0._r8; p_mass = 0._r8;
  if(present(c_flx))c_flx=0._r8
  if(present(n_flx))n_flx=0._r8
  if(present(p_flx))p_flx=0._r8


  c_mass = c_mass * catomw
  n_mass = n_mass * natomw
  p_mass = p_mass * patomw

  if(present(c_flx))then
    c_flx = c_flx + ystates1(lid_co2_hr)
    c_flx = c_flx * catomw
  endif

  if(present(n_flx))then
    n_flx = n_flx * natomw
  endif

  if(present(p_flx))then
    p_flx = p_flx * patomw
  endif

  end associate

  end subroutine sumup_cnp_msflx_ecosys

  !--------------------------------------------------------------------
  subroutine sumup_cnp_mass(this, header, c_mass, n_mass, p_mass)
  use tracer_varcon         , only : catomw, natomw, patomw
  implicit none
  class(ecosys_bgc_type)     , intent(in) :: this
  character(len=*), intent(in) :: header
  real(r8), intent(out) :: c_mass, n_mass, p_mass
  !local variables

  integer :: kc, kn, kp, jj
  associate(                                   &
    lid_nh3x => this%ecosys_bgc_index%lid_nh3x, &
    lid_hno3 => this%ecosys_bgc_index%lid_hno3, &
    lid_hno2 => this%ecosys_bgc_index%lid_hno2, &
    lid_h0po4=> this%ecosys_bgc_index%lid_h0po4,&
    lid_h1po4=> this%ecosys_bgc_index%lid_h1po4,&
    lid_h2po4=> this%ecosys_bgc_index%lid_h2po4,&
    lid_h3po4=> this%ecosys_bgc_index%lid_h3po4,&
    ystates1 => this%ystates1                   &
  )
  print*,header
  c_mass = 0._r8; n_mass = 0._r8; p_mass = 0._r8;

  n_mass = n_mass + ystates1(lid_nh3x) + ystates1(lid_hno3)+ ystates1(lid_hno2)

  c_mass = c_mass * catomw
  n_mass = n_mass * natomw
  p_mass = p_mass * patomw
  if(p_mass>1.e10_r8)then
     print*,'sum cnp bad mass',p_mass
     stop
  endif
  end associate
  end subroutine sumup_cnp_mass

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
    class(ecosys_bgc_type),  intent(inout)  :: me
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
    yc=y0
    call me%bgc_integrate(yc, dt, tt, nprimeq, neq, f)
    call ebbks(yc, f, nprimeq, neq, dt, y, pscal)
    return
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

  !-------------------------------------------------------------------------------
  subroutine init_cold_ecosys(this, nstvars, ystates)
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  integer , intent(in) :: nstvars
  real(r8), intent(inout) :: ystates(nstvars)
  ystates(:)=0._r8
  if(nstvars>=0)continue
  end subroutine init_cold_ecosys

  !-------------------------------------------------------------------------------
  subroutine sum_tot_store(this, nstvars, ystates)
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  integer                   , intent(in)    :: nstvars
  real(r8)                  , intent(inout) :: ystates(nstvars)

  associate(                  &
    lid_totstore => this%ecosys_bgc_index%lid_totstore   &
  )

  end associate
  end subroutine sum_tot_store
  !-------------------------------------------------------------------------------
  subroutine display_index(this)
  implicit none
  class(ecosys_bgc_type)     , intent(inout) :: this
  print*,'no3',this%ecosys_bgc_index%lid_plant_minn_no3_pft
  print*,'nh4', this%ecosys_bgc_index%lid_plant_minn_nh4_pft

  end subroutine display_index
end module ecosysBGCType
