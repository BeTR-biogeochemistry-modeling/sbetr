module simicParaType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use BiogeoContype       , only : BiogeoCon_type
implicit none

 private
  character(len=*), private, parameter :: filename = &
       __FILE__

 type, public, extends(BiogeoCon_type) :: simic_para_type

  !decomposition
  real(r8) :: Kaff_EP_LIT  !enzyme affinity for litter POM depolymerization, mol C/m3
  real(r8) :: Kaff_EP_POM  !enzyme affinity for sorbed POM depolymerization, mol C/m3
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
  real(r8) :: vmax_EP_L    !maximum depolymerization rate for litter
  real(r8) :: vmax_BC    !maximum DOC assimilation rate
  real(r8) :: alpha_B2E  !scaling parameter from microbial biomass to enzyme
  real(r8) :: alpha_B2T
  real(r8) :: Kaff_EM  !enzyme affinity for mineral surface, mol C /m3
  real(r8) :: Minsurf  !mineral surface area for DOC/enzyme/microbial cell wall material adsorption
  real(r8) :: Kaff_o2  !oxygen affinity
  real(r8) :: Kmort_MB
  real(r8) :: fpom_vmax
  real(r8) :: fpom_desorb
 contains
   procedure, public  :: Init => simic_para_init
   procedure, public  :: readPars => simic_para_readPars
   procedure, public  :: printPars=> simic_para_printPars
   procedure, private :: InitAllocate_simpara
   procedure, private :: set_defpar_default
   procedure, public  :: apply_spinup_factor
   procedure, public  :: set_spinup_factor
 end type simic_para_type

 type(simic_para_type), public :: simic_para
 public :: create_jarpars_simic
contains

  function create_jarpars_simic()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(simic_para_type), pointer :: create_jarpars_simic
    class(simic_para_type), pointer :: bgc

    allocate(bgc)
    create_jarpars_simic => bgc

  end function create_jarpars_simic

  !--------------------------------------------------------------------
  subroutine simic_para_init(this, bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  implicit none
  class(simic_para_type), intent(inout) :: this
  type(betr_status_type)   , intent(out) :: bstatus

  call this%bcon_Init(bstatus)
  if(bstatus%check_status())return

  call this%InitAllocate_simpara()

  call this%set_defpar_default()

  if(betr_spinup_state/=0)then
    call this%apply_spinup_factor()
  endif
  end subroutine simic_para_init
  !--------------------------------------------------------------------
  subroutine InitAllocate_simpara(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(simic_para_type), intent(inout) :: this


  !the following will be actually calculated from CNP bgc
  end subroutine InitAllocate_simpara
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)

  use tracer_varcon      , only : natomw,patomw
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  implicit none
  class(simic_para_type), intent(inout) :: this
  real(r8) :: half_life

  !the following note is from daycent
  !3.90000           'DEC1(1)'  surface litter, structural, 1/y
  !4.90000           'DEC1(2)'  soil, structural, 1/y
  !14.80000          'DEC2(1)'  surface litter, metabolic, 1/y
  !18.50000          'DEC2(2)'  soil, metabolic, 1/y
  !6.00000           'DEC3(1)'  surface litter, som1, 1/y
  !7.30000           'DEC3(2)'  soil, som1, 1/y
  !0.00450           'DEC4'     soil, som3, 1/y
  !0.20000           'DEC5'     soil, som2, 1/y
  !1.50000           'DECW1'    dead fine branch
  !0.50000           'DECW2'    dead large wood component
  !0.60000           'DECW3'    dead coarse root component
  !0.45000           'PS1CO2(1)', surface
  !0.55000           'PS1CO2(2)', soil
  !0.60000           'P1CO2A(1)'
  !0.17000           'P1CO2A(2)'
  !0.00000           'P1CO2B(1)'
  !0.68000           'P1CO2B(2)'
  !p1co2 = p1co2a + p1co2 * sand, som1->som2,som3, co2 resp frac, on surface, som1->som2
  !0.45000           'PS1CO2(1)', surface litter, struct to som1 and som2
  !0.55000           'PS1CO2(2)', soil, struct to som1 and som2
  !0.30000           'RSPLIG', fraction of lignin flow lost to respiration
  !0.55000           'PMCO2(1)' , surface litter, meta to som1, co2 resp frac
  !0.55000           'PMCO2(2)' , soil meta to som1, co2 resp frac
  !0.55000           'P2CO2', som2->som1,som3, co2 resp frac
  !0.55000           'P3CO2', som3->som1, co2 resp frac

  half_life = 5568._r8 ! yr
  half_life = half_life * year_sec

  !decomposition
  this%Kaff_EP_LIT  = 1.e-3_r8
  this%Kaff_EP_POM  = 4.e-2_r8
  this%Kaff_ED  = 1.e-3_r8
  this%Kaff_BC  = 1.e-5_r8
  this%cue_met  = 0.5_r8
  this%cue_cel  = 0.4_r8
  this%cue_lig  = 0.3_r8
  this%cue_cwd  = 0.3_r8
  this%cue_bm   = 0.5_r8
  this%Rm0_spmic= 1.e-7_r8  !s-1 /
  this%Mrt_spmic= 4.e-8_r8
  this%f_mic2C = 0.3_r8
  this%f_mic2D = 0.7_r8
  this%vmax_EP_L = 4.e-5_r8
  this%vmax_BC = 3.e-5_r8
  this%alpha_B2E = 0.05_r8
  this%alpha_B2T = 0.05_r8
  this%Kaff_CM  = 1.e-1_r8
  this%Kaff_EM  = 1.e-1_r8
  this%Minsurf  = 100._r8
  this%Kaff_o2  = 0.22_r8
  this%Kmort_MB = 0.01_r8
  this%fpom_vmax= 2.e-6_r8              !s-1
  this%fpom_desorb=this%fpom_vmax*3.e-3_r8
  end subroutine set_defpar_default


  !--------------------------------------------------------------------
  subroutine apply_spinup_factor(this)
  use betr_ctrl, only : betr_spinup_state
  implicit none
  class(simic_para_type), intent(inout) :: this

  call this%set_spinup_factor()

  if(betr_spinup_state==1)then

  endif
  end subroutine apply_spinup_factor

  !--------------------------------------------------------------------

  subroutine simic_para_readPars(this, ncid, bstatus)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use ncdio_pio       , only : file_desc_t, ncd_io
  use BetrStatusType  , only : betr_status_type
  use betr_ctrl       , only : betr_spinup_state
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  use tracer_varcon   , only : natomw,patomw
  implicit none
  class(simic_para_type), intent(inout) :: this
  type(file_desc_t)    , intent(inout)  :: ncid  ! pio netCDF file id
  type(betr_status_type) , intent(out) :: bstatus

  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr(1) ! temporary to read in constant
  real(r8)           :: temparr(1:2,1:1)
  character(len=100) :: tString ! temp. var for reading
  call bstatus%reset()

  return

  if(betr_spinup_state/=0)then
    call this%apply_spinup_factor()
  endif

  end subroutine simic_para_readPars

!--------------------------------------------------------------------
  subroutine set_spinup_factor(this)

  implicit none
  class(simic_para_type), intent(inout) :: this
  real(r8) :: k_decay_ref



  end subroutine set_spinup_factor

!--------------------------------------------------------------------
  subroutine simic_para_printPars(this)

  implicit none
  class(simic_para_type), intent(inout) :: this

  call this%prtPars_bgc()


  end subroutine simic_para_printPars

end module simicParaType
