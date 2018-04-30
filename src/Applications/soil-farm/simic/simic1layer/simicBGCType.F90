module simicBGCType
#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  ! subroutines for stoichiometric configuration of the century bgc
  ! !History, created by Jinyun Tang, Dec, 2014.
  ! !USES:
  use bshr_kind_mod             , only : r8 => shr_kind_r8
  use bshr_log_mod              , only : errMsg => shr_log_errMsg
  use betr_varcon               , only : spval => bspval
  use betr_ctrl                 , only : spinup_state => betr_spinup_state
  use gbetrType                 , only : gbetr_type
  use BiogeoConType             , only : BiogeoCon_type
  use simicParaType             , only : simic_para_type
  use BetrStatusType            , only : betr_status_type
  use BeTRJarModel              , only : jar_model_type
  use simicBGCIndexType            , only : simic_index_type
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !Note:
  !Keeping simic_index as a private member is a workaround to call the ode solver
  !it increase the memory for each instance of simic_bgc_type, but enables
  !the ode function to be called by the ode solver

  type, extends(jar_model_type), public :: simic_bgc_type
    type(simic_index_type),     private :: simic_index
    real(r8), pointer                    :: ystates0(:)
    real(r8), pointer                    :: ystates1(:)
    real(r8), pointer                    :: scal_f(:)
    real(r8), pointer                    :: conv_f(:)
    real(r8), pointer                    :: conc_f(:)
    logical , private                    :: use_c13
    logical , private                    :: use_c14
  contains
    procedure, public  :: init          => init_simic
    procedure, public  :: runbgc        => runbgc_simic
    procedure, public  :: UpdateParas   => UpdateParas_simic
    procedure, public  :: getvarllen    => getvarllen_simic
    procedure, public  :: getvarlist    => getvarlist_simic
    procedure, private :: arenchyma_gas_transport
    procedure, private :: init_states
    procedure, private :: add_ext_input
    procedure, private :: InitAllocate
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
  print*,'nvar',this%simic_index%nstvars
  ans =  this%simic_index%nstvars

  end function getvarllen_simic
  !-------------------------------------------------------------------------------
  subroutine getvarlist_simic(this, nstvars, varnames, varunits)
  implicit none
  class(simic_bgc_type) , intent(inout) :: this
  integer, intent(in) :: nstvars
  character(len=*), intent(out) :: varnames(1:nstvars)
  character(len=*), intent(out) :: varunits(1:nstvars)
  integer :: n

  do n = 1, nstvars
    write(varnames(n),'(A)')trim(this%simic_index%varnames(n))
    write(varunits(n),'(A)')trim(this%simic_index%varunits(n))
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

  class default
    write(msg,'(A)')'Wrong parameter type passed in for UpdateParas in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine UpdateParas_simic
  !-------------------------------------------------------------------------------
  subroutine init_simic(this,  biogeo_con,  bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(simic_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)    , intent(out) :: bstatus

  character(len=256) :: msg
  write(this%jarname, '(A)')'simic'

  select type(biogeo_con)
  type is(simic_para_type)
    call bstatus%reset()
    call this%simic_index%Init(biogeo_con%use_c13, biogeo_con%use_c14, &
     biogeo_con%non_limit, biogeo_con%nop_limit, betr_maxpatch_pft)
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

  call this%InitAllocate(this%simic_index)
  end subroutine init_simic

  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, simic_index)

  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(simic_bgc_type)   , intent(inout) :: this
  type(simic_index_type)  , intent(in) :: simic_index

  associate(                              &
    nstvars => simic_index%nstvars      , &
    nprimvars=> simic_index%nprimvars     &
  )

  allocate(this%ystates0(nstvars)); this%ystates0(:) = 0._r8
  allocate(this%ystates1(nstvars)); this%ystates1(:) = 0._r8
  allocate(this%scal_f(nprimvars));  this%scal_f(:) = 0._r8
  allocate(this%conv_f(nprimvars));  this%conv_f(:) = 0._r8
  allocate(this%conc_f(nprimvars));  this%conc_f(:) = 0._r8
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
  use tracer_varcon         , only : catomw, natomw, patomw
  implicit none
  class(simic_bgc_type)  , intent(inout) :: this
  logical                    , intent(in)    :: is_surflit
  real(r8)                   , intent(in)    :: dtime
  type(JarBGC_forc_type)     , intent(in)    :: bgc_forc
  integer                    , intent(in)    :: nstates
  real(r8)                   , intent(out)   :: ystates0(nstates)
  real(r8)                   , intent(out)   :: ystatesf(nstates)
  type(betr_status_type)     , intent(out)   :: bstatus

  !local variables
  character(len=*),parameter :: subname = 'runbgc_simic'

  associate(                                            &
    ystates1 => this%ystates1                           &
  )
  call bstatus%reset()

  !initialize state variables
  call this%init_states(this%simic_index, bgc_forc)

  ystates0(:) = this%ystates0(:)

  call this%add_ext_input(dtime, this%simic_index, bgc_forc)

  call this%arenchyma_gas_transport(this%simic_index, dtime)
  ystatesf(:) = ystates1(:)
  end associate
  end subroutine runbgc_simic

  !-------------------------------------------------------------------------------
  subroutine init_states(this, simic_index, bgc_forc)

  use simicBGCIndexType            , only : simic_index_type
  use JarBgcForcType            , only : JarBGC_forc_type
  implicit none
  class(simic_bgc_type)  , intent(inout) :: this
  type(simic_index_type)  , intent(in) :: simic_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  associate(                               &
    lid_co2_hr  => simic_index%lid_co2_hr, &
    lid_n2 => simic_index%lid_n2,   &
    lid_o2 => simic_index%lid_o2,   &
    lid_co2 => simic_index%lid_co2, &
    lid_c13_co2 => simic_index%lid_c13_co2, &
    lid_c14_co2 => simic_index%lid_c14_co2, &
    lid_ar => simic_index%lid_ar,   &
    lid_o2_paere => simic_index%lid_o2_paere, &
    lid_n2_paere => simic_index%lid_n2_paere, &
    lid_ar_paere => simic_index%lid_ar_paere, &
    lid_co2_paere => simic_index%lid_co2_paere, &
    lid_ch4_paere => simic_index%lid_ch4_paere, &
    lid_c13_co2_paere => simic_index%lid_c13_co2_paere, &
    lid_c14_co2_paere => simic_index%lid_c14_co2_paere, &
    lid_ch4 => simic_index%lid_ch4  &
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
  subroutine add_ext_input(this, dtime, simic_index, bgc_forc)

  use simicBGCIndexType            , only : simic_index_type
  use JarBgcForcType        , only : JarBGC_forc_type
  use tracer_varcon             , only : catomw, natomw, patomw,c13atomw,c14atomw
  use MathfuncMod               , only : safe_div
  implicit none
  class(simic_bgc_type)  , intent(inout) :: this
  real(r8), intent(in) :: dtime
  type(simic_index_type)  , intent(in) :: simic_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  associate(                        &
    lit1 =>  simic_index%lit1, &
    lit2 =>  simic_index%lit2, &
    lit3 =>  simic_index%lit3, &
    cwd =>   simic_index%cwd   &
  )

  this%ystates1(lit1)=this%ystates0(lit1) +  dtime *  bgc_forc%cflx_input_litr_met

  this%ystates1(lit2)=this%ystates0(lit2) +  dtime *  bgc_forc%cflx_input_litr_cel

  this%ystates1(lit3)=this%ystates0(lit3) +  dtime *  bgc_forc%cflx_input_litr_lig

  this%ystates1(cwd)=this%ystates0(cwd) +  dtime *  bgc_forc%cflx_input_litr_cwd


  end associate
  end subroutine add_ext_input


  !--------------------------------------------------------------------
  subroutine arenchyma_gas_transport(this, simic_index, dtime)
  use simicBGCIndexType       , only : simic_index_type
  implicit none
  class(simic_bgc_type)     , intent(inout) :: this
  type(simic_index_type)  , intent(in) :: simic_index
  real(r8), intent(in) :: dtime

  integer :: j
  real(r8) :: y0
  associate(                             &
    lid_n2 => simic_index%lid_n2,   &
    lid_o2 => simic_index%lid_o2,   &
    lid_co2 => simic_index%lid_co2, &
    lid_c13_co2 => simic_index%lid_c13_co2, &
    lid_c14_co2 => simic_index%lid_c14_co2, &
    lid_ar => simic_index%lid_ar,   &
    lid_o2_paere => simic_index%lid_o2_paere, &
    lid_n2_paere => simic_index%lid_n2_paere, &
    lid_ar_paere => simic_index%lid_ar_paere, &
    lid_co2_paere => simic_index%lid_co2_paere, &
    lid_ch4_paere => simic_index%lid_ch4_paere, &
    lid_c13_co2_paere => simic_index%lid_c13_co2_paere, &
    lid_c14_co2_paere => simic_index%lid_c14_co2_paere, &
    lid_ch4 => simic_index%lid_ch4  &
  )

  j = lid_o2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(simic_index%lid_o2_paere) = this%ystates1(j)-y0

  if( spinup_state == 0)then
    j = lid_n2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_index%lid_n2_paere) = this%ystates1(j)-y0

    j = lid_ar; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_index%lid_ar_paere) = this%ystates1(j)-y0

    j = lid_ch4; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_index%lid_ch4_paere) = this%ystates1(j)-y0

    j = lid_co2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(simic_index%lid_co2_paere) = this%ystates1(j)-y0

    if(this%use_c13)then
      j = lid_c13_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(simic_index%lid_c13_co2_paere) = this%ystates1(j)-y0

    endif

    if(this%use_c14)then
      j = lid_c14_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(simic_index%lid_c14_co2_paere) = this%ystates1(j)-y0

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

end module simicBGCType
