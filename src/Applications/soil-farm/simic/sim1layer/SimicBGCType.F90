module SimicBGCType
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
  use SimParaType               , only : SimPara_type
  use BetrStatusType            , only : betr_status_type
  use BeTRJarModel              , only : jar_model_type
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !Note:
  !Keeping centurybgc_index as a private member is a workaround to call the ode solver
  !it increase the memory for each instance of simicbgc_type, but enables
  !the ode function to be called by the ode solver

  type, extends(jar_model_type), public :: simicbgc_type
    real(r8), pointer                    :: ystates0(:)
    real(r8), pointer                    :: ystates1(:)

  contains
    procedure, public  :: init          => init_simic
    procedure, public  :: runbgc        => runbgc_simic
    procedure, public  :: UpdateParas   => UpdateParas_simic
    procedure, public  :: getvarllen    => getvarllen_simic
    procedure, public  :: getvarlist    => getvarlist_simic
  end type simicbgc_type

  public :: create_jarmodel_simicbgc
contains

  function create_jarmodel_simicbgc()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(simicbgc_type), pointer :: create_jarmodel_simicbgc
    class(simicbgc_type), pointer :: bgc

    allocate(bgc)
    create_jarmodel_simicbgc => bgc

  end function create_jarmodel_simicbgc

  !-------------------------------------------------------------------------------
  function getvarllen_simic(this)result(ans)

  implicit none
  class(simicbgc_type) , intent(inout) :: this
  integer :: ans

  ans = 1

  end function getvarllen_simic
  !-------------------------------------------------------------------------------
  subroutine getvarlist_simic(this, nstvars, varnames, varunits)
  implicit none
  class(simicbgc_type) , intent(inout) :: this
  integer, intent(in) :: nstvars
  character(len=*), intent(out) :: varnames(1:nstvars)
  character(len=*), intent(out) :: varunits(1:nstvars)
  integer :: n

  do n = 1, nstvars
!    write(varnames(n),'(A)')trim(this%centurybgc_index%varnames(n))
!    write(varunits(n),'(A)')trim(this%centurybgc_index%varunits(n))
  enddo
  end subroutine getvarlist_simic


  !-------------------------------------------------------------------------------
  subroutine UpdateParas_simic(this,  biogeo_con, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(simicbgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)     , intent(out)   :: bstatus
  integer :: sr
  character(len=256) :: msg

  call bstatus%reset()

  select type(biogeo_con)
  type is(SimPara_type)


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
  class(simicbgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)    , intent(out) :: bstatus

  character(len=256) :: msg
  write(this%jarname, '(A)')'simic'

  select type(biogeo_con)
  type is(SimPara_type)
    call bstatus%reset()

    call this%UpdateParas(biogeo_con, bstatus)
    if(bstatus%check_status())return
  class default
    call bstatus%reset()
    write(msg,'(A)')'Wrong parameter type passed in for init_simic in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine init_simic

  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this)

  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(simicbgc_type)   , intent(inout) :: this


!  allocate(this%ystates0(nstvars)); this%ystates0(:) = 0._r8
!  allocate(this%ystates1(nstvars)); this%ystates1(:) = 0._r8


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
  class(simicbgc_type)  , intent(inout) :: this
  logical                    , intent(in)    :: is_surflit
  real(r8)                   , intent(in)    :: dtime
  type(JarBGC_forc_type)     , intent(in)    :: bgc_forc
  integer                    , intent(in)    :: nstates
  real(r8)                   , intent(out)   :: ystates0(nstates)
  real(r8)                   , intent(out)   :: ystatesf(nstates)
  type(betr_status_type)     , intent(out)   :: bstatus

  !local variables
  character(len=*),parameter :: subname = 'runbgc_simic'



  end subroutine runbgc_simic

  !-------------------------------------------------------------------------------

end module SimicBGCType
