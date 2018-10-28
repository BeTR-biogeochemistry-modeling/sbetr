module fakeParaType
use bshr_kind_mod   , only : r8 => shr_kind_r8
use BiogeoContype   , only : BiogeoCon_type
implicit none

  private
  character(len=*), private, parameter :: filename = &
       __FILE__
  type, public :: fake_para_type
  !declare variables here

  contains
    procedure, public  :: Init     => fake_para_Init
    procedure, public  :: readPars => fake_para_readPars
    procedure, public  :: printPars=> fake_para_printPars
    procedure, private :: fake_InitAllocate
    procedure, private :: set_defpar_default
  end type fake_para_type

  type(fake_para_type), public :: fake_para
  public :: create_jarpars_fake
contains

  function create_jarpars_fake()
  ! DESCRIPTION
  ! constructor
  implicit none
  class(fake_para_type), pointer :: create_jarpars_fake
  class(fake_para_type), pointer :: bgc

  allocate(bgc)
  create_jarpars_fake => bgc

  end function create_jarpars_fake
 !--------------------------------------------------------------------
  subroutine fake_para_Init(this, bstatus)
  !
  !DESCRIPTION
  !initialize default parameters

  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  implicit none
  class(fake_para_type), intent(inout) :: this
  type(betr_status_type)       , intent(out) :: bstatus


  call this%fake_InitAllocate()
  call this%set_defpar_default()

  end subroutine fake_para_Init
 !--------------------------------------------------------------------
  subroutine fake_InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(fake_para_type), intent(inout) :: this
  !allocate memory for necessary varaibles

  end subroutine fake_InitAllocate
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)
  !
  !DESCRIPTION
  !set default value for relevant parameters
  use tracer_varcon      , only : natomw,patomw
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  implicit none
  class(fake_para_type), intent(inout) :: this


  end subroutine set_defpar_default
 !--------------------------------------------------------------------
  subroutine fake_para_readPars(this, ncid, bstatus)
  !
  !DESCRIPTION
  !read parameters from input file
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use ncdio_pio       , only : file_desc_t, ncd_io
  use BetrStatusType  , only : betr_status_type
  use betr_ctrl       , only : betr_spinup_state
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  use tracer_varcon   , only : natomw,patomw
  implicit none
  class(fake_para_type), intent(inout) :: this
  type(file_desc_t)            , intent(inout) :: ncid  ! pio netCDF file id
  type(betr_status_type)       , intent(out)   :: bstatus

  !local variables
  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr(1) ! temporary to read in constant
  real(r8)           :: temparr(1:2,1:1)
  character(len=100) :: tString ! temp. var for reading

  call bstatus%reset()
  return

  end subroutine fake_para_readPars
!--------------------------------------------------------------------
  subroutine fake_para_printPars(this)

  implicit none
  class(fake_para_type), intent(inout) :: this

  end subroutine fake_para_printPars
end module fakeParaType
