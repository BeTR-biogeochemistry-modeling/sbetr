module ch4soilParaType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use BiogeoContype       , only : BiogeoCon_type
implicit none

 private
  character(len=*), private, parameter :: filename = &
       __FILE__

 type, public, extends(BiogeoCon_type) :: ch4soil_para_type

  !decomposition
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding

  real(r8) :: rf_l1s1_bgc(2)
  real(r8) :: rf_l2s1_bgc(2)
  real(r8) :: rf_l3s2_bgc
  real(r8) :: rf_s2s1_bgc
  real(r8) :: rf_s3s1_bgc
  real(r8) :: rf_s1s2a_bgc(2)
  real(r8) :: rf_s1s2b_bgc(2)

  real(r8) :: k_decay_lit1(2)
  real(r8) :: k_decay_lit2(2)
  real(r8) :: k_decay_lit3(2)
  real(r8) :: k_decay_som1(2)
  real(r8) :: k_decay_som2
  real(r8) :: k_decay_som3
  real(r8) :: k_decay_cwd
  real(r8) :: k_decay_lwd
  real(r8) :: k_decay_fwd
  real(r8) :: k_m_o2_bgc  !MM parameter for O2 consumption

  !nitrification-denitrification
  real(r8) :: nitrif_n2o_loss_frac
  real(r8) :: organic_max
  real(r8) :: rij_kro_a
  real(r8) :: rij_kro_alpha
  real(r8) :: rij_kro_beta
  real(r8) :: rij_kro_gamma
  real(r8) :: rij_kro_delta
  real(r8) :: surface_tension_water
  real(r8) :: minpsi_bgc

  real(r8) :: c14decay_som_const
  real(r8) :: c14decay_dom_const
  real(r8) :: c14decay_pom_const
  real(r8) :: c14decay_bm_const
  real(r8) :: k_nitr_max

  real(r8) :: init_cn_som1
  real(r8) :: init_cn_som2
  real(r8) :: init_cn_som3

  real(r8) :: init_cp_som1
  real(r8) :: init_cp_som2
  real(r8) :: init_cp_som3

  real(r8) :: init_cc14_som1
  real(r8) :: init_cc14_som2
  real(r8) :: init_cc14_som3

  real(r8) :: init_cc13_som1
  real(r8) :: init_cc13_som2
  real(r8) :: init_cc13_som3
  real(r8) :: c14decay_const
  real(r8) :: km_decomp_no3
  real(r8) :: km_decomp_nh4
  real(r8) :: km_decomp_p
  real(r8) :: km_den
  real(r8) :: km_nit
  real(r8) :: vmax_decomp_n
  real(r8) :: vmax_decomp_p
  real(r8) :: vmax_den
  real(r8) :: vmax_nit
  real(r8), pointer :: spinup_factor(:)

 contains
   procedure, public  :: Init => centpara_init
   procedure, public  :: readPars => centpara_readPars
   procedure, public  :: printPars=> centpara_printPars
   procedure, private :: InitAllocate_centpara
   procedure, private :: set_defpar_default
   procedure, public  :: apply_spinup_factor
   procedure, public  :: set_spinup_factor
 end type ch4soil_para_type

 type(ch4soil_para_type), public :: ch4soil_para
 public :: create_jarpars_ch4soil
contains

  function create_jarpars_ch4soil()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(ch4soil_para_type), pointer :: create_jarpars_ch4soil
    class(ch4soil_para_type), pointer :: bgc

    allocate(bgc)
    create_jarpars_ch4soil => bgc

  end function create_jarpars_ch4soil

  !--------------------------------------------------------------------
  subroutine centpara_init(this, bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  implicit none
  class(ch4soil_para_type), intent(inout) :: this
  type(betr_status_type)   , intent(out) :: bstatus

  call this%bcon_Init(bstatus)
  if(bstatus%check_status())return

  call this%InitAllocate_centpara()

  call this%set_defpar_default()

  if(betr_spinup_state/=0)then
    call this%apply_spinup_factor()
  endif
  end subroutine centpara_init
  !--------------------------------------------------------------------
  subroutine InitAllocate_centpara(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(ch4soil_para_type), intent(inout) :: this


  allocate(this%spinup_factor(9))
  !the following will be actually calculated from CNP bgc
  end subroutine InitAllocate_centpara
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)

  use tracer_varcon      , only : natomw,patomw
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  implicit none
  class(ch4soil_para_type), intent(inout) :: this
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
  this%c14decay_const = - log(0.5_r8) / half_life
  this%c14decay_som_const  =this%c14decay_const
  this%c14decay_dom_const  =this%c14decay_const
  this%c14decay_pom_const  =this%c14decay_const
  this%c14decay_Bm_const  =this%c14decay_const

  !decomposition
  this%Q10                   = 2._r8
  this%froz_q10              = 10._r8
  this%decomp_depth_efolding = 1._r8

  !following is based on century parameterization
  this%rf_l1s1_bgc           = (/0.55_r8, 0.55_r8/)
  this%rf_l2s1_bgc           = (/0.45_r8, 0.55_r8/)
  this%rf_l3s2_bgc           = 0.3_r8
  this%rf_s2s1_bgc           = 0.55_r8
  this%rf_s3s1_bgc           = 0.55_r8
  this%rf_s1s2a_bgc          = (/0.60_r8,0.17_r8/)
  this%rf_s1s2b_bgc          = (/0._r8, 0.68_r8/)

  !following is based on century parameterization
  this%k_decay_lit1          = (/14.8_r8,18.5_r8/)/year_sec    !1/second
  this%k_decay_lit2          = (/3.9_r8 ,4.9_r8/) /year_sec    !1/second
  this%k_decay_lit3          = (/3.9_r8 ,4.9_r8/) /year_sec    !1/second
  this%k_decay_som1          = (/6.7_r8, 7.3_r8/) /year_sec    !1/second
  this%k_decay_som2          = 0.2_r8/year_sec                 !1/second
  this%k_decay_som3          = 0.0045_r8/year_sec              !1/second
  this%k_decay_cwd           = 0.6_r8/year_sec                 !1/second
  this%k_decay_fwd           = 1.5_r8/year_sec                 !1/second
  this%k_decay_lwd           = 0.5_r8/year_sec                 !1/second

  !nitrification-denitrification
  this%nitrif_n2o_loss_frac  = 1.e-4_r8   !Arah and Vinten, 1995
  this%organic_max           = 160._r8    !organic matter content (kg/m3) where soil is assumed to act like peat
  this%rij_kro_a             = 1.5e-10_r8 ! Arah and Vinten, 1995
  this%rij_kro_alpha         = 1.26_r8    ! Arah and Vinten, 1995
  this%rij_kro_beta          = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_gamma         = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_delta         = 0.85_r8    ! Arah and Vinten, 1995
  this%surface_tension_water = 73.e-3_r8  ! (J/m^2), Arah and Vinten, 1995
  this%minpsi_bgc            = -10._r8    ! M Pa
  this%k_m_o2_bgc            = 0.22_r8
  this%k_nitr_max            = 1.1574074e-06_r8
  this%init_cp_som1 = 110._r8 !mass based
  this%init_cp_som2 = 320._r8 !mass based
  this%init_cp_som3 = 114._r8 !mass based
  this%init_cn_som1 = 8._r8   !mass based
  this%init_cn_som2 = 11._r8  !mass based
  this%init_cn_som3 = 11._r8  !mass based

  this%init_cc14_som1= 0._r8
  this%init_cc14_som2= 0._r8
  this%init_cc14_som3= 0._r8

  this%init_cc13_som1= 0._r8
  this%init_cc13_som2= 0._r8
  this%init_cc13_som3= 0._r8

  this%km_decomp_no3 = 0.41_r8/natomw
  this%km_decomp_nh4 = 0.18_r8/natomw
  this%km_decomp_p   = 0.2_r8/patomw
  this%km_den        = 0.11_r8/natomw
  this%km_nit        = 0.76_r8/natomw

  this%vmax_decomp_n= 5.e-3_r8/3600._r8  ! 1/s
  this%vmax_decomp_p  = 1.e-3_r8/3600._r8  ! 1/s
  this%vmax_den = 1.8_r8 /86400._r8     ! 1/s
  this%vmax_nit = 0.67_r8/86400._r8     ! 1/s

  end subroutine set_defpar_default


  !--------------------------------------------------------------------
  subroutine apply_spinup_factor(this)
  use betr_ctrl, only : betr_spinup_state
  implicit none
  class(ch4soil_para_type), intent(inout) :: this


  call this%set_spinup_factor()

  if(betr_spinup_state==1)then
    this%k_decay_lit1 = this%k_decay_lit1 * this%spinup_factor(1)
    this%k_decay_lit2 = this%k_decay_lit2 * this%spinup_factor(2)
    this%k_decay_lit3 = this%k_decay_lit3 * this%spinup_factor(3)
    this%k_decay_cwd = this%k_decay_cwd * this%spinup_factor(4)
    this%k_decay_lwd = this%k_decay_lwd * this%spinup_factor(5)
    this%k_decay_fwd = this%k_decay_fwd * this%spinup_factor(6)
    this%k_decay_som1 = this%k_decay_som1 * this%spinup_factor(7)
    this%k_decay_som3 = this%k_decay_som3 * this%spinup_factor(8)
    this%k_decay_som2 = this%k_decay_som2 * this%spinup_factor(9)

    this%c14decay_Bm_const   =this%c14decay_Bm_const * this%spinup_factor(7)
    this%c14decay_som_const  =this%c14decay_som_const * this%spinup_factor(8)
    this%c14decay_dom_const  =this%c14decay_dom_const * this%spinup_factor(9)

  endif
  end subroutine apply_spinup_factor

  !--------------------------------------------------------------------

  subroutine centpara_readPars(this, ncid, bstatus)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use ncdio_pio       , only : file_desc_t, ncd_io
  use BetrStatusType  , only : betr_status_type
  use betr_ctrl       , only : betr_spinup_state
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  use tracer_varcon   , only : natomw,patomw
  implicit none
  class(ch4soil_para_type), intent(inout) :: this
  type(file_desc_t)    , intent(inout)  :: ncid  ! pio netCDF file id
  type(betr_status_type) , intent(out) :: bstatus

  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr(1) ! temporary to read in constant
  real(r8)           :: temparr(1:2,1:1)
  character(len=100) :: tString ! temp. var for reading
  call bstatus%reset()

  call this%readPars_bgc(ncid, bstatus)
  if(bstatus%check_status())return
  tString='surface_tension_water'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__),err=-1)
  if(bstatus%check_status())return
  this%surface_tension_water=tempr(1)

  tString='rij_kro_a'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__),err=-1)
  if(bstatus%check_status())return
  this%rij_kro_a=tempr(1)

  tString='rij_kro_alpha'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_alpha=tempr(1)

  tString='rij_kro_beta'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_beta=tempr(1)

  tString='rij_kro_gamma'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_gamma=tempr(1)

  tString='rij_kro_delta'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_delta=tempr(1)

  tString='rf_l1s1_bgc'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_l1s1_bgc(1:2)=temparr(1:2,1)

  tString='rf_l2s1_bgc'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_l2s1_bgc(1:2)=temparr(1:2,1)

  tString='rf_l3s2_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_l3s2_bgc=tempr(1)

  tString='rf_s2s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_s2s1_bgc=tempr(1)

  tString='rf_s3s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_s3s1_bgc=tempr(1)

  tString='rf_s1s2a_bgc'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_s1s2a_bgc(1:2)=temparr(1:2,1)

  tString='rf_s1s2b_bgc'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_s1s2b_bgc(1:2)=temparr(1:2,1)

  tString='k_decay_cwd'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_cwd=tempr(1)

  tString='k_decay_lwd'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_lwd=tempr(1)

  tString='k_decay_fwd'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_fwd=tempr(1)

  tString='k_decay_lit1'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_lit1(1:2) = temparr(1:2,1)

  tString='k_decay_lit2'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_lit2(1:2) = temparr(1:2,1)

  tString='k_decay_lit3'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_lit3(1:2) = temparr(1:2,1)

  tString='k_decay_som1'
  call ncd_io(trim(tString),temparr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_som1(1:2)=temparr(1:2,1)

  tString='k_decay_som2'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_som2=tempr(1)

  tString='k_decay_som3'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_decay_som3=tempr(1)

  tString='froz_q10'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%froz_q10=tempr(1)

  tString='decomp_depth_efolding'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%decomp_depth_efolding=tempr(1)

  tString='Q10'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%Q10=tempr(1)

  tString='minpsi_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%minpsi_bgc=tempr(1)

  tString='k_m_o2_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_m_o2_bgc=tempr(1)

  tString='organic_max'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%organic_max=tempr(1)

  tString='k_nitr_max'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_nitr_max=tempr(1)

  tString='km_den'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%km_den=tempr(1)/natomw

  tString='km_nit'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%km_nit=tempr(1)/natomw

  tString='km_decomp_nh4'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%km_decomp_nh4=tempr(1)/natomw

  tString='km_decomp_no3'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%km_decomp_no3=tempr(1)/natomw

  tString='km_decomp_p'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%km_decomp_p=tempr(1)/patomw

  tString='vmax_decomp_n'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%vmax_decomp_n=tempr(1)/3600._r8

  tString='vmax_decomp_p'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%vmax_decomp_p=tempr(1)/3600._r8

  tString='vmax_den'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%vmax_den=tempr(1)

  tString='vmax_nit'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%vmax_nit=tempr(1)

!   the following are purposely commented out to use default parameters from century, Jinyun Tang, Feb 21, 2018
  this%k_decay_lit1(1:2)  = this%k_decay_lit1(1:2) /year_sec    !1/second
  this%k_decay_lit2(1:2)  = this%k_decay_lit2(1:2) /year_sec    !1/second
  this%k_decay_lit3(1:2)  = this%k_decay_lit3(1:2) /year_sec    !1/second
  this%k_decay_som1(1:2)  = this%k_decay_som1(1:2) /year_sec    !1/second
  this%k_decay_som2       = this%k_decay_som2/year_sec          !1/second
  this%k_decay_som3       = this%k_decay_som3/year_sec          !1/second
  this%k_decay_cwd        = this%k_decay_cwd /year_sec          !1/second
  this%k_decay_fwd        = this%k_decay_fwd /year_sec          !1/second
  this%k_decay_lwd        = this%k_decay_lwd /year_sec          !1/second

  if(betr_spinup_state/=0)then
    call this%apply_spinup_factor()
  endif

  end subroutine centpara_readPars

!--------------------------------------------------------------------
  subroutine set_spinup_factor(this)

  implicit none
  class(ch4soil_para_type), intent(inout) :: this
  real(r8) :: k_decay_ref

  !the order is, lit1, lit2, lit3, cwd, lwd, fwd, som1, som3, som2
  this%spinup_factor(1) = 1._r8
  this%spinup_factor(2) = 1._r8
  this%spinup_factor(3) = 1._r8

  this%spinup_factor(4) = 1._r8
  this%spinup_factor(5) = 1._r8
  this%spinup_factor(6) = 1._r8

  k_decay_ref=this%k_decay_som1(1)
  this%spinup_factor(7) = k_decay_ref/this%k_decay_som1(1)
  this%spinup_factor(8) = k_decay_ref/this%k_decay_som3
  this%spinup_factor(9) = k_decay_ref/this%k_decay_som2

  end subroutine set_spinup_factor

!--------------------------------------------------------------------
  subroutine centpara_printPars(this)

  implicit none
  class(ch4soil_para_type), intent(inout) :: this

  call this%prtPars_bgc()

  !decomposition
  print*,'Q10=',this%Q10
  print*,'froz_q10=',this%froz_q10
  print*,'decomp_depth_efolding=',this%decomp_depth_efolding

  !following is based on century parameterization
  print*,'rf_l1s1_bgc =',this%rf_l1s1_bgc
  print*,'rf_l2s1_bgc =',this%rf_l2s1_bgc
  print*,'rf_l3s2_bgc =',this%rf_l3s2_bgc
  print*,'rf_s2s1_bgc =',this%rf_s2s1_bgc
  print*,'rf_s3s1_bgc =',this%rf_s3s1_bgc
  print*,'rf_s1s2a_bgc=',this%rf_s1s2a_bgc
  print*,'rf_s1s2b_bgc=',this%rf_s1s2b_bgc


  !following is based on century parameterization
  print*,'k_decay_lit1=',this%k_decay_lit1
  print*,'k_decay_lit2=',this%k_decay_lit2
  print*,'k_decay_lit3=',this%k_decay_lit3
  print*,'k_decay_som1=',this%k_decay_som1
  print*,'k_decay_som2=',this%k_decay_som2
  print*,'k_decay_som3=',this%k_decay_som3
  print*,'k_decay_cwd=',this%k_decay_cwd
  print*,'k_decay_fwd=',this%k_decay_fwd
  print*,'k_decay_lwd=',this%k_decay_lwd

  !nitrification-denitrification
  print*,'nitrif_n2o_loss_frac=',this%nitrif_n2o_loss_frac
  print*,'organic_max=',this%organic_max
  print*,'rij_kro_a =',this%rij_kro_a
  print*,'rij_kro_alpha=',this%rij_kro_alpha
  print*,'rij_kro_beta=',this%rij_kro_beta
  print*,'rij_kro_gamma =',this%rij_kro_gamma
  print*,'rij_kro_delta=',this%rij_kro_delta
  print*,'surface_tension_water=',this%surface_tension_water

  print*,'init_cp_som1=',this%init_cp_som1
  print*,'init_cp_som2=',this%init_cp_som2
  print*,'init_cp_som3=',this%init_cp_som3
  print*,'init_cn_som1=',this%init_cn_som1
  print*,'init_cn_som2=',this%init_cn_som2
  print*,'init_cn_som3=',this%init_cn_som3

  print*,'km_den = ', this%km_den
  print*,'km_nit = ', this%km_nit
  print*,'km_decomp_nh4=', this%km_decomp_nh4
  print*,'km_decomp_no3=', this%km_decomp_no3
  print*,'km_decomp_p  =', this%km_decomp_p
  print*,'vmax_decomp_p=', this%vmax_decomp_p
  print*,'vmax_decomp_n=', this%vmax_decomp_n
  print*,'vmax_den=', this%vmax_den
  print*,'vmax_nit=', this%vmax_nit
  end subroutine centpara_printPars

end module ch4soilParaType
