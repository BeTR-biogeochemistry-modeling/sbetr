module v1ecaParaType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use BiogeoContype       , only : BiogeoCon_type
implicit none

 private
  character(len=*), private, parameter :: filename = &
       __FILE__

 type, public, extends(BiogeoCon_type) :: v1eca_para_type

  !decomposition
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding

  real(r8) :: rf_l1s1_bgc
  real(r8) :: rf_l2s1_bgc
  real(r8) :: rf_l3s2_bgc
  real(r8) :: rf_s2s1_bgc
  real(r8) :: rf_s3s1_bgc
  real(r8) :: rf_s2s3_bgc
  real(r8) :: rf_cwdl2_bgc
  real(r8) :: rf_cwdl3_bgc
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
   procedure, public  :: Init => v1eca_init
   procedure, public  :: readPars => v1eca_readPars
   procedure, public  :: printPars=> v1eca_printPars
   procedure, private :: InitAllocate_v1eca
   procedure, private :: set_defpar_default
   procedure, public  :: apply_spinup_factor
   procedure, public  :: set_spinup_factor
 end type v1eca_para_type

 type(v1eca_para_type), public :: v1eca_para
 public :: create_jarpars_v1eca
contains

  function create_jarpars_v1eca()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(v1eca_para_type), pointer :: create_jarpars_v1eca
    class(v1eca_para_type), pointer :: bgc

    allocate(bgc)
    create_jarpars_v1eca => bgc

  end function create_jarpars_v1eca

  !--------------------------------------------------------------------
  subroutine v1eca_init(this, bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  implicit none
  class(v1eca_para_type), intent(inout) :: this
  type(betr_status_type)   , intent(out) :: bstatus

  call this%bcon_Init(bstatus)
  if(bstatus%check_status())return

  call this%InitAllocate_v1eca()

  call this%set_defpar_default()

  if(betr_spinup_state/=0)then
    call this%apply_spinup_factor()
  endif
  end subroutine v1eca_init
  !--------------------------------------------------------------------
  subroutine InitAllocate_v1eca(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(v1eca_para_type), intent(inout) :: this


  allocate(this%spinup_factor(9))
  !the following will be actually calculated from CNP bgc
  end subroutine InitAllocate_v1eca
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)

  use tracer_varcon      , only : natomw,patomw
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  implicit none
  class(v1eca_para_type), intent(inout) :: this
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
  this%rf_l1s1_bgc           = 0.55_r8
  this%rf_l2s1_bgc           = 0.45_r8
  this%rf_l3s2_bgc           = 0.3_r8
  this%rf_s2s1_bgc           = 0.55_r8
  this%rf_s3s1_bgc           = 0.55_r8


  !nitrification-denitrification
  this%nitrif_n2o_loss_frac  = 6.e-4_r8   !Arah and Vinten, 1995
  this%organic_max           = 160._r8    !organic matter content (kg/m3) where soil is assumed to act like peat
  this%rij_kro_a             = 1.5e-10_r8 ! Arah and Vinten, 1995
  this%rij_kro_alpha         = 1.26_r8    ! Arah and Vinten, 1995
  this%rij_kro_beta          = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_gamma         = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_delta         = 0.85_r8    ! Arah and Vinten, 1995
  this%surface_tension_water = 73.e-3_r8  ! (J/m^2), Arah and Vinten, 1995
  this%minpsi_bgc            = -10._r8    ! M Pa
  this%k_m_o2_bgc            = 0.0022_r8
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
  class(v1eca_para_type), intent(inout) :: this

  return
  call this%set_spinup_factor()


  end subroutine apply_spinup_factor

  !--------------------------------------------------------------------

  subroutine v1eca_readPars(this, ncid, bstatus)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use ncdio_pio       , only : file_desc_t, ncd_io
  use BetrStatusType  , only : betr_status_type
  use betr_ctrl       , only : betr_spinup_state
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  use tracer_varcon   , only : natomw,patomw
  implicit none
  class(v1eca_para_type), intent(inout) :: this
  type(file_desc_t)    , intent(inout)  :: ncid  ! pio netCDF file id
  type(betr_status_type) , intent(out) :: bstatus

  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr ! temporary to read in constant
  real(r8)           :: temparr(1:2,1:1)
  character(len=100) :: tString ! temp. var for reading
  call bstatus%reset()

  call this%readPars_bgc(ncid, bstatus)
  if(bstatus%check_status())return

  tString='k_nitr_max'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_nitr_max=tempr

  tString='surface_tension_water'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%surface_tension_water=tempr

  tString='rij_kro_a'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_a=tempr

  tString='rij_kro_alpha'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_alpha=tempr

  tString='rij_kro_beta'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_beta=tempr

  tString='rij_kro_gamma'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_gamma=tempr

  tString='rij_kro_delta'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rij_kro_delta=tempr


  tString='cn_s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%init_cn_som1=tempr

  tString='cn_s2_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%init_cn_som2=tempr

  tString='cn_s3_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%init_cn_som3=tempr

!!! read in phosphorus variables - note that these NP ratio parameters for BGC  will have
!!! to be added in the parameter file

  tString='np_s1_new'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%init_cp_som1=this%init_cn_som1*tempr

  tString='np_s2_new'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%init_cp_som2=this%init_cn_som2*tempr

  tString='np_s3_new'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%init_cp_som3=this%init_cn_som3*tempr

  tString='rf_l1s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_l1s1_bgc=tempr

  tString='rf_l2s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_l2s1_bgc=tempr

  tString='rf_l3s2_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_l3s2_bgc=tempr

  tString='rf_s2s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_s2s1_bgc=tempr

  tString='rf_s2s3_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  this%rf_s2s3_bgc=tempr

  tString='rf_s3s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  this%rf_s3s1_bgc=tempr

  tString='rf_cwdl2_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  this%rf_cwdl2_bgc=tempr

  tString='rf_cwdl3_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  this%rf_cwdl3_bgc=tempr

  tString='organic_max'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  this%organic_max=tempr


  end subroutine v1eca_readPars

!--------------------------------------------------------------------
  subroutine set_spinup_factor(this)

  implicit none
  class(v1eca_para_type), intent(inout) :: this
  real(r8) :: k_decay_ref

  end subroutine set_spinup_factor

!--------------------------------------------------------------------
  subroutine v1eca_printPars(this)

  implicit none
  class(v1eca_para_type), intent(inout) :: this

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
  end subroutine v1eca_printPars

end module v1ecaParaType
