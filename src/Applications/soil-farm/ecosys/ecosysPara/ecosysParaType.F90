module ecosysParaType
use bshr_kind_mod   , only : r8 => shr_kind_r8
use BiogeoContype   , only : BiogeoCon_type
implicit none

  private
  character(len=*), private, parameter :: filename = &
       __FILE__
  type, public, extends(BiogeoCon_type) :: ecosys_para_type
  !declare variables here

  real(r8) :: ORAD    !
  real(r8) :: BIOS    !
  real(r8) :: BIOA    !
  real(r8) :: DCKI    !
  real(r8) :: RCCX    !
  real(r8) :: RCCQ    !
  real(r8) :: RCCZ    !
  real(r8) :: RCCY    !
  real(r8) :: FPRIM   !
  real(r8) :: FPRIMM  !
  real(r8) :: OMGR    !
  real(r8) :: OQKI    !
  real(r8) :: H2KI    !
  real(r8) :: OAKI    !
  real(r8) :: COMKI   !
  real(r8) :: COMKM   !
  real(r8) :: CKC     !
  real(r8) :: FOSCZ0  !
  real(r8) :: FOSCZL  !
  real(r8) :: FMN     !
  real(r8) :: DCKM0   !
  real(r8) :: DCKML   !

  real(r8) :: VMXO    !
  real(r8) :: VMXF    !
  real(r8) :: VMXM    !
  real(r8) :: VMXH    !
  real(r8) :: VMXN    !
  real(r8) :: VMX4    !
  real(r8) :: VMXC    !
  real(r8) :: OQKM    !
  real(r8) :: OQKA    !
  real(r8) :: OQKAM   !
  real(r8) :: CCKM    !
  real(r8) :: CCK4    !
  real(r8) :: ZHKM    !
  real(r8) :: ZHKI    !
  real(r8) :: ZNKM    !
  real(r8) :: Z3KM    !
  real(r8) :: Z2KM    !
  real(r8) :: Z1KM    !
  real(r8) :: Z4MX    !
  real(r8) :: Z4KU    !
  real(r8) :: Z4MN    !
  real(r8) :: ZOMX    !
  real(r8) :: ZOKU    !
  real(r8) :: ZOMN    !
  real(r8) :: HPMX    !
  real(r8) :: HPKU    !
  real(r8) :: HPMN    !
  real(r8) :: ZFKM    !
  real(r8) :: H2KM    !
  real(r8) :: ECNH    !
  real(r8) :: ECNO    !
  real(r8) :: ECN3    !
  real(r8) :: ECN2    !
  real(r8) :: ECN1    !
  real(r8) :: RNFNI   !
  real(r8) :: ECHO    !
  real(r8) :: VMKI    !
  real(r8) :: VHKI    !
  real(r8) :: OXKA    !
  real(r8) :: EDNH    !
  real(r8) :: EDNA    !

  real(r8) :: EOMC    !
  real(r8) :: EOMD    !
  real(r8) :: EOMG    !
  real(r8) :: EOMF    !
  real(r8) :: EOMH    !
  real(r8) :: EOMN    !
  real(r8) :: GO2X    !
  real(r8) :: GH4X    !
  real(r8) :: GCHX    !
  real(r8) :: GO2A    !
  real(r8) :: GC4X    !
  real(r8) :: GCOX    !
  real(r8) :: GNOX    !
  real(r8) :: GN2X    !
  real(r8) :: EN2X    !
  real(r8) :: EN2Y    !
  real(r8) :: EO2X    !
  real(r8) :: EH4X    !
  real(r8) :: EO2G    !
  real(r8) :: EO2D    !
  real(r8) :: ENFX    !
  real(r8) :: ENOX    !
  real(r8) :: EO2A    !

  real(r8) :: TSORP   !
  real(r8) :: HSORP   !
  real(r8) :: DOSA    !
  real(r8) :: DOSX    !

  real(r8) :: SPOHC   !
  real(r8) :: SPOHA   !
  real(r8) :: RMOM    !

  real(r8) :: SPOSC(4,0:4)
  real(r8) :: SPORC(2)
  real(r8) :: SPOMC(2)
  real(r8) :: DCKX(0:4)
  real(r8) :: EN2F(7)
  contains
    procedure, public  :: Init     => ecosys_para_Init
    procedure, public  :: readPars => ecosys_para_readPars
    procedure, public  :: printPars=> ecosys_para_printPars
    procedure, private :: ecosys_InitAllocate
    procedure, private :: set_defpar_default
  end type ecosys_para_type

  type(ecosys_para_type), public :: ecosys_para
  public :: create_jarpars_ecosys
contains

  function create_jarpars_ecosys()
  ! DESCRIPTION
  ! constructor
  implicit none
  class(ecosys_para_type), pointer :: create_jarpars_ecosys
  class(ecosys_para_type), pointer :: bgc

  allocate(bgc)
  create_jarpars_ecosys => bgc

  end function create_jarpars_ecosys
 !--------------------------------------------------------------------
  subroutine ecosys_para_Init(this, bstatus)
  !
  !DESCRIPTION
  !initialize default parameters

  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  implicit none
  class(ecosys_para_type), intent(inout) :: this
  type(betr_status_type)       , intent(out) :: bstatus

  call this%bcon_Init(bstatus)
  if(bstatus%check_status())return

  call this%ecosys_InitAllocate()
  call this%set_defpar_default()

  end subroutine ecosys_para_Init
 !--------------------------------------------------------------------
  subroutine ecosys_InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(ecosys_para_type), intent(inout) :: this
  !allocate memory for necessary varaibles

  end subroutine ecosys_InitAllocate
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)
  !
  !DESCRIPTION
  !set default value for relevant parameters
  use tracer_varcon      , only : natomw,patomw
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  implicit none
  class(ecosys_para_type), intent(inout) :: this

  !SUBSTRATE DECOMPOSITION BY MICROBIAL POPULATIONS
  this%ORAD=1.0E-06_r8
  this%BIOS=1.0E-06_r8/(4.19_r8*this%ORAD**3._r8)
  this%BIOA=this%BIOS*12.57_r8*this%ORAD**2._r8
  this%DCKI=2.5_r8
  this%RCCX=0.833_r8
  this%RCCQ=0.833_r8
  this%RCCZ=0.167_r8
  this%RCCY=0.833_r8
  this%FPRIM=5.0E-02_r8
  this%FPRIMM=1.0E-06_r8
  this%OMGR=2.5E-01_r8
  this%OQKI=1.2E+03_r8
  this%H2KI=1.0_r8      !
  this%OAKI=12.0_r8
  this%COMKI=1.0E-03_r8
  this%COMKM=1.0E-04_r8  !MM constant
  this%CKC=1.0E-03_r8
  this%FOSCZ0=2.0E-02_r8
  this%FOSCZL=5.0E-06_r8
  this%FMN=1.0E-03_r8    !minimum value for fraction of substrate demanded by each microbial or plant population
  this%DCKM0=5.0E+03_r8
  this%DCKML=1.0E+03_r8
!     SPECIFIC RESPIRATION RATES, M-M UPTAKE CONSTANTS,
!     STOICHIOMETRIC CONSTANTS FOR MICROBIAL REDOX REACTIONS

  this%VMXO=0.125_r8  !maximum specific oxidation rate of nonstructural C by symbiotic N2 fixers
  this%VMXF=0.125_r8
  this%VMXM=0.125_r8
  this%VMXH=0.375_r8  !maximum specific oxidation rate of NH3 by nitrifiers	h-1
  this%VMXN=0.25_r8
  this%VMX4=0.375_r8
  this%VMXC=0.125_r8
  this%OQKM=1.2E+01_r8   ! MM constant for DOC uptake, gC/m3
  this%OQKA=1.2E+01_r8   ! MM constant for acetate uptake, gC/m3
  this%OQKAM=1.2E+01_r8
  this%CCKM=0.15_r8
  this%CCK4=1.2E-03_r8
  this%ZHKM=1.4_r8       ! MM constant for nh4 uptake, gN/m3
  this%ZHKI=7.0E+03_r8
  this%ZNKM=1.4_r8       !mm constatnt for no2 uptake, gN/m3
  this%Z3KM=1.4_r8       !Km for NO3 reduction by denitrifiers	g m-3
  this%Z2KM=1.4_r8       !Km for NO2 reduction by denitrifiers	g m-3
  this%Z1KM=0.14_r8      !Km for N2O reduction by denitrifiers	g m-3
  this%Z4MX=5.0E-03_r8
  this%Z4KU=0.40_r8
  this%Z4MN=0.0125_r8
  this%ZOMX=5.0E-03_r8
  this%ZOKU=0.35_r8
  this%ZOMN=0.03_r8
  this%HPMX=1.0E-03_r8
  this%HPKU=0.075_r8
  this%HPMN=0.002_r8
  this%ZFKM=0.14_r8
  this%H2KM=0.01_r8
  this%ECNH=0.30_r8   !ratio of CO2-C fixed per NH3-N oxidized	g g-1
  this%ECNO=0.10_r8   !ratio of CO2-C fixed per NO2-N oxidized	g g-1
  this%ECN3=0.857_r8  !ratio of DOC oxidized to NO3-N reduced by denitrifiers	g g-1
  this%ECN2=0.857_r8  !ratio of DOC oxidized to NO2-N reduced by denitrifiers	g g-1
  this%ECN1=0.429_r8  !ratio of DOC oxidized to N2O-N reduced by denitrifiers	g g-1
  this%RNFNI=2.0E-04_r8   !nitrification inhibition effect
  this%ECHO=0.75_r8
  this%VMKI=0.25_r8       !inhibition constant for transfer of electrons from O2 to NO3 during denitrification	g m-3 h-1
  this%VHKI=15.0_r8
  this%OXKA=0.16_r8
  this%EDNH=1.00_r8
  this%EDNA=1.00_r8

!     ENERGY REQUIREMENTS FOR MICROBIAL GROWTH AND
!     ENERGY YIELDS FROM REDUCTION OF O2, OC, CH4, NO3, N2
  this%EOMC=25.0_r8
  this%EOMD=37.5_r8
  this%EOMG=37.5_r8
  this%EOMF=37.5_r8
  this%EOMH=25.0_r8
  this%EOMN=75.0_r8
  this%GO2X=37.5_r8
  this%GH4X=66.5_r8
  this%GCHX=4.50_r8
  this%GO2A=this%GO2X-this%GCHX
  this%GC4X=3.00_r8
  this%GCOX=11.00_r8  !energy yield of H2 oxidation and CO2 reduction	kJ g C-1
  this%GNOX=10.0_r8
  this%GN2X=187.5_r8
  this%EN2X=  this%GO2X/this%GN2X  !ratio of N2 fixed per DOC oxidized by aerobic diazotrophs, g/g
  this%EN2Y=  this%GCHX/ this%GN2X !ratio of N2 fixed per DOC oxidized by anaerobic diazotrophs, g/g
  this%EO2X= 1.0_r8/(1.0_r8+this%GO2X/this%EOMC)  !ratio of DOC oxidization to DOC uptake by aerobic micobial biomass, g/g
  this%EH4X= 1.0_r8/(1.0_r8+this%GH4X/this%EOMC)
  this%EO2G= 1.0_r8/(1.0_r8+this%GO2X/this%EOMG)  !ratio of DOC oxidization to DOC uptake by fungal biomass, g/g
  this%EO2D= 1.0_r8/(1.0_r8+this%GO2X/this%EOMD)  !ratio of DOC oxidization to DOC uptake by denitrifier biomass reducing O2, g/g
  this%ENFX= 1.0_r8/(1.0_r8+this%GO2X/this%EOMN)  !ratio of DOC oxidization to DOC uptake by diazotrophic biomass, g/g
  this%ENOX= 1.0_r8/(1.0_r8+this%GNOX/this%EOMC)  !ratio of DOC oxidization to DOC uptake by denitrifier biomass reducing NOx, g/g
  this%EO2A= 1.0_r8/(1.0_r8+this%GO2A/this%EOMC)  !ratio of acetate oxidization to acetate uptake by aerobic micobial biomass, g/g

!     SORPTION RATE CONSTANTS
  this%TSORP=0.5_r8
  this%HSORP=1.0_r8
  this%DOSA=0.5_r8
  this%DOSX=5.0E-03_r8

!     SPECIFIC DECOMPOSITION RATES
  this%SPOHC=0.25_r8
  this%SPOHA=0.25_r8
  this%RMOM=0.010_r8

  this%SPOSC(:,0)=(/7.50_r8,7.50_r8,1.50_r8,0.50_r8/)
  this%SPOSC(:,1)=(/7.50_r8,7.50_r8,1.50_r8,0.50_r8/)
  this%SPOSC(:,2)=(/7.50_r8,7.50_r8,1.50_r8,0.50_r8/)
  this%SPOSC(:,3)=(/0.05_r8,0.00_r8,0.00_r8,0.00_r8/)
  this%SPOSC(:,4)=(/0.05_r8,0.0167_r8,0.00_r8,0.00_r8/)
  this%SPORC=(/7.5_r8,1.5_r8/)
  this%SPOMC=(/1.e-02_r8,5.0e-4_r8/)
  this%DCKX(0:4)=(/0.50_r8,0.50_r8,0.50_r8,0.00_r8,0.00_r8/)

! MICROBIAL C:N:P RATIOS DURING HUMIFICATION
  this%EN2F=(/0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,this%EN2X,this%EN2Y/)
  end subroutine set_defpar_default
 !--------------------------------------------------------------------
  subroutine ecosys_para_readPars(this, ncid, bstatus)
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
  class(ecosys_para_type), intent(inout) :: this
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
  !read shared parameters
  call this%readPars_bgc(ncid, bstatus)

  end subroutine ecosys_para_readPars
!--------------------------------------------------------------------
  subroutine ecosys_para_printPars(this)

  implicit none
  class(ecosys_para_type), intent(inout) :: this
  call this%prtPars_bgc()
  end subroutine ecosys_para_printPars
end module ecosysParaType
