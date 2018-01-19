module BiogeoConType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
implicit none

 private
  character(len=*), private, parameter :: filename = &
       __FILE__
  real(r8), private, parameter :: year_sec=86400._r8*365._r8
 type, public :: BiogeoCon_type

  !decomposition
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding

  real(r8) :: rf_l1s1_bgc
  real(r8) :: rf_l2s1_bgc
  real(r8) :: rf_l3s2_bgc
  real(r8) :: rf_s2s1_bgc
  real(r8) :: rf_s3s1_bgc
  real(r8) :: cwd_fcel_bgc
  real(r8) :: cwd_flig_bgc
  real(r8) :: lwd_fcel_bgc
  real(r8) :: lwd_flig_bgc
  real(r8) :: fwd_fcel_bgc
  real(r8) :: fwd_flig_bgc
  real(r8) :: k_decay_lit1
  real(r8) :: k_decay_lit2
  real(r8) :: k_decay_lit3
  real(r8) :: k_decay_som1
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
  real(r8) :: init_cn_met
  real(r8) :: init_cn_cel
  real(r8) :: init_cn_lig
  real(r8) :: init_cn_cwd
  real(r8) :: init_cn_lwd
  real(r8) :: init_cn_fwd
  real(r8) :: init_cn_som1
  real(r8) :: init_cn_som2
  real(r8) :: init_cn_som3

  real(r8) :: init_cp_met
  real(r8) :: init_cp_cel
  real(r8) :: init_cp_lig
  real(r8) :: init_cp_cwd
  real(r8) :: init_cp_lwd
  real(r8) :: init_cp_fwd
  real(r8) :: init_cp_som1
  real(r8) :: init_cp_som2
  real(r8) :: init_cp_som3

  real(r8) :: init_cc13_met
  real(r8) :: init_cc13_cel
  real(r8) :: init_cc13_lig
  real(r8) :: init_cc13_cwd
  real(r8) :: init_cc13_lwd
  real(r8) :: init_cc13_fwd
  real(r8) :: init_cc13_som1
  real(r8) :: init_cc13_som2
  real(r8) :: init_cc13_som3

  real(r8) :: init_cc14_met
  real(r8) :: init_cc14_cel
  real(r8) :: init_cc14_lig
  real(r8) :: init_cc14_cwd
  real(r8) :: init_cc14_lwd
  real(r8) :: init_cc14_fwd
  real(r8) :: init_cc14_som1
  real(r8) :: init_cc14_som2
  real(r8) :: init_cc14_som3
  real(r8) :: c14decay_const
  real(r8) :: c14decay_som_const
  real(r8) :: c14decay_dom_const
  real(r8) :: c14decay_bm_const
  real(r8) :: k_nitr_max
  logical :: use_c13
  logical :: use_c14
  logical :: nop_limit                                              !switch for P limitation
  logical :: non_limit
  !ECA nutrient competition
  real(r8), pointer :: vmax_minp_soluble_to_secondary(:)  => null() !maximum conversion rate of soluble P into secondary P

  !inorganic phosphorus cycling
  real(r8), pointer :: frac_p_sec_to_sol(:)    => null()   !fraction of released secondary phosphorus that goes into soluble form
  real(r8), pointer :: minp_secondary_decay(:) => null()   !decay rate of secondary phosphorus
  real(r8), pointer :: spinup_factor(:)
 contains
   procedure, public  :: Init
   procedure, private :: InitAllocate
   procedure, private :: set_defpar_default
   procedure, public  :: apply_spinup_factor
   procedure, private :: ReadNamelist
   procedure, public  :: readPars
   procedure, public  :: set_spinup_factor
 end type BiogeoCon_type

 type(BiogeoCon_type), public :: bgc_con_eca
contains
  !--------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  character(len=betr_namelist_buffer_size_ext) , intent(in)    :: namelist_buffer
  type(betr_status_type)                   , intent(out) :: bstatus

  call bstatus%reset()

  call this%InitAllocate()

  call this%set_defpar_default()

  !update parameter from namelist
  call this%ReadNamelist(namelist_buffer, bstatus)

  if(betr_spinup_state/=0)then
    call this%apply_spinup_factor()
  endif

  end subroutine Init
  !--------------------------------------------------------------------
  subroutine InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(BiogeoCon_type), intent(inout) :: this


  allocate(this%minp_secondary_decay(0:betr_max_soilorder))
  allocate(this%vmax_minp_soluble_to_secondary(0:betr_max_soilorder))
  allocate(this%frac_p_sec_to_sol(0:betr_max_soilorder))
  allocate(this%spinup_factor(9))
  !the following will be actually calculated from CNP bgc
  end subroutine InitAllocate
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  real(r8) :: half_life

  half_life = 5568._r8 ! yr
  half_life = half_life * 86400._r8 * 365._r8
  this%c14decay_const = - log(0.5_r8) / half_life
  this%c14decay_som_const  =this%c14decay_const
  this%c14decay_dom_const  =this%c14decay_const
  this%c14decay_Bm_const  =this%c14decay_const
  !decomposition
  this%Q10                   = 2._r8
  this%froz_q10              = 10._r8
  this%decomp_depth_efolding = 1._r8

  !following is based on Table 15.4 in CLM4.5 tech note
  this%rf_l1s1_bgc           = 0.55_r8
  this%rf_l2s1_bgc           = 0.5_r8
  this%rf_l3s2_bgc           = 0.5_r8
  this%rf_s2s1_bgc           = 0.55_r8
  this%rf_s3s1_bgc           = 0.55_r8
  this%cwd_fcel_bgc          = 0.76_r8
  this%cwd_flig_bgc          = 0.24_r8

  !following is based on Table 15.3 in CLM4.5 tech note
  this%k_decay_lit1          = 1._r8/(0.066_r8*86400._r8*365._r8)   !1/second
  this%k_decay_lit2          = 1._r8/(0.25_r8*86400._r8*365._r8)    !1/second
  this%k_decay_lit3          = 1._r8/(0.25_r8*86400._r8*365._r8)    !1/second
  this%k_decay_som1          = 1._r8/(0.17_r8*86400._r8*365._r8)    !1/second
  this%k_decay_som2          = 1._r8/(6.1_r8*86400._r8*365._r8)     !1/second
  this%k_decay_som3          = 1._r8/(270._r8*86400._r8*365._r8)    !1/second
  this%k_decay_cwd           = 1._r8/(4.1_r8*86400._r8*365._r8)     !1/second
  this%k_decay_fwd           = 1._r8/(4.1_r8*86400._r8*365._r8)     !1/second
  this%k_decay_lwd           = 1._r8/(4.1_r8*86400._r8*365._r8)     !1/second

  this%init_cn_met  = 90._r8  !mass based
  this%init_cn_cel  = 90._r8  !mass based
  this%init_cn_lig  = 90._r8  !mass based
  this%init_cn_cwd  = 90._r8  !mass based
  this%init_cn_fwd  = 90._r8  !mass based
  this%init_cn_lwd  = 90._r8  !mass based
  this%init_cn_som1 = 8._r8   !mass based
  this%init_cn_som2 = 11._r8  !mass based
  this%init_cn_som3 = 11._r8  !mass based

  this%init_cp_met  = 1600._r8
  this%init_cp_cel  = 2000._r8
  this%init_cp_lig  = 2500._r8
  this%init_cp_cwd  = 4500._r8
  this%init_cp_lwd  = 4500._r8!mass based
  this%init_cp_fwd  = 4500._r8!mass based
  this%init_cp_som1 = 110._r8 !mass based
  this%init_cp_som2 = 320._r8 !mass based
  this%init_cp_som3 = 114._r8 !mass based

  !nitrification-denitrification
  this%nitrif_n2o_loss_frac  = 1.e-4_r8   !Arah and Vinten, 1995
  this%organic_max           = 160._r8    !organic matter content (kg/m3) where soil is assumed to act like peat
  this%rij_kro_a             = 1.5e-10_r8 ! Arah and Vinten, 1995
  this%rij_kro_alpha         = 1.26_r8    ! Arah and Vinten, 1995
  this%rij_kro_beta          = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_gamma         = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_delta         = 0.85_r8    ! Arah and Vinten, 1995
  this%surface_tension_water = 73.e-3_r8  ! (J/m^2), Arah and Vinten, 1995

  !ECA nutrient competition
  this%vmax_minp_soluble_to_secondary(:) = 1.e-9_r8  !1/s on the order of 1.e-9 1/s
  !inorganic phosphorus cycling
  !Note: (1._r8-frac_p_sec_to_sol)*minp_secondary_decay = occlusion rate
  this%frac_p_sec_to_sol(:)              = 0.95_r8    !fraction of released secondary phosphorus that goes into soluble form, occ rate is on the order of  1.e-13 1/s
  this%minp_secondary_decay(:)           = 1.e-11_r8  !decay rate of secondary phosphorus, 1/s, this is on the order of 1.e-11 1/s

  this%use_c13 = .false.
  this%use_c14 = .false.
  this%nop_limit=.false.
  this%non_limit=.false.
  this%init_cc13_met = 0._r8
  this%init_cc13_cel = 0._r8
  this%init_cc13_lig = 0._r8
  this%init_cc13_cwd = 0._r8
  this%init_cc13_som1= 0._r8
  this%init_cc13_som2= 0._r8
  this%init_cc13_som3= 0._r8

  this%init_cc14_met = 0._r8
  this%init_cc14_cel = 0._r8
  this%init_cc14_lig = 0._r8
  this%init_cc14_cwd = 0._r8
  this%init_cc14_som1= 0._r8
  this%init_cc14_som2= 0._r8
  this%init_cc14_som3= 0._r8
  end subroutine set_defpar_default

  !--------------------------------------------------------------------
  subroutine ReadNamelist(this, namelist_buffer, bstatus)
  !
  ! DESCRIPTION
  ! reading bgc parameters
  ! will be updated later
  use betr_constants , only : stdout, betr_string_length_long, betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : iulog => biulog
  use betr_ctrl      , only : betr_spinup_state
  use bshr_log_mod   , only : errMsg => shr_log_errMsg
  use tracer_varcon  , only : use_c13_betr, use_c14_betr, is_nitrogen_active, is_phosphorus_active
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  character(len=betr_namelist_buffer_size_ext) , intent(in)    :: namelist_buffer
  type(betr_status_type), intent(out) :: bstatus

  !
  ! !LOCAL VARIABLES:
  integer                                :: nml_error
  character(len=betr_string_length_long) :: ioerror_msg
  real(r8) :: tau_decay_lit1
  real(r8) :: tau_decay_lit2
  real(r8) :: tau_decay_lit3
  real(r8) :: tau_decay_som1
  real(r8) :: tau_decay_som2
  real(r8) :: tau_decay_som3
  real(r8) :: tau_decay_cwd
  real(r8) :: tau_decay_fwd
  real(r8) :: tau_decay_lwd

  call bstatus%reset()

  !years
  tau_decay_lit1          = 0.066_r8
  tau_decay_lit2          = 0.25_r8
  tau_decay_lit3          = 0.25_r8
  tau_decay_som1          = 0.17_r8
  tau_decay_som2          = 6.1_r8
  tau_decay_som3          = 270._r8
  tau_decay_cwd           = 4.1_r8
  tau_decay_fwd           = 4.1_r8
  tau_decay_lwd           = 4.1_r8

  this%use_c13 = use_c13_betr
  this%use_c14 = use_c14_betr
  this%nop_limit=.not. is_phosphorus_active
  this%non_limit=.not. is_nitrogen_active
  this%k_decay_lit1          = 1._r8/(tau_decay_lit1*year_sec)    !1/second
  this%k_decay_lit2          = 1._r8/(tau_decay_lit2*year_sec)    !1/second
  this%k_decay_lit3          = 1._r8/(tau_decay_lit3*year_sec)    !1/second
  this%k_decay_som1          = 1._r8/(tau_decay_som1*year_sec)    !1/second
  this%k_decay_som2          = 1._r8/(tau_decay_som2*year_sec)    !1/second
  this%k_decay_som3          = 1._r8/(tau_decay_som3*year_sec)    !1/second
  this%k_decay_cwd           = 1._r8/(tau_decay_cwd*year_sec)     !1/second
  this%k_decay_fwd           = 1._r8/(tau_decay_fwd*year_sec)     !1/second
  this%k_decay_lwd           = 1._r8/(tau_decay_lwd*year_sec)     !1/second


  end subroutine ReadNamelist

  !--------------------------------------------------------------------
  subroutine apply_spinup_factor(this)
  use betr_ctrl, only : betr_spinup_state
  implicit none
  class(BiogeoCon_type), intent(inout) :: this


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

  subroutine readPars(this, ncid, bstatus)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use ncdio_pio       , only : file_desc_t, ncd_io
  use BetrStatusType  , only : betr_status_type
  use betr_ctrl       , only : betr_spinup_state
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  type(file_desc_t)    , intent(inout)  :: ncid  ! pio netCDF file id
  type(betr_status_type) , intent(out) :: bstatus

  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr ! temporary to read in constant
  character(len=100) :: tString ! temp. var for reading
  real(r8) :: tau_decay_lit1
  real(r8) :: tau_decay_lit2
  real(r8) :: tau_decay_lit3
  real(r8) :: tau_decay_som1
  real(r8) :: tau_decay_som2
  real(r8) :: tau_decay_som3
  real(r8) :: tau_decay_cwd
  real(r8) :: tau_decay_fwd
  real(r8) :: tau_decay_lwd

  call bstatus%reset()

  tString='surface_tension_water'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__),err=-1)
  if(bstatus%check_status())return
  this%surface_tension_water=tempr

  tString='rij_kro_a'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__),err=-1)
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

  tString='rf_s3s1_bgc'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%rf_s3s1_bgc=tempr

  tString='cwd_fcel'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%cwd_fcel_bgc=tempr

  tString='lwd_fcel'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%lwd_fcel_bgc=tempr

  tString='fwd_fcel'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%fwd_fcel_bgc=tempr

  tString='cwd_flig'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%cwd_flig_bgc=tempr

  tString='lwd_flig'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%lwd_flig_bgc=tempr

  tString='fwd_flig'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%fwd_flig_bgc=tempr

  tString='tau_cwd'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_cwd=tempr

  tString='tau_fwd'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_fwd=tempr

  tString='tau_lwd'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_lwd=tempr

  tString='tau_l1'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_lit1 = tempr

  tString='tau_l2'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_lit2 = tempr

  tString='tau_l3'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_lit3 = tempr

  tString='tau_s1'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_som1=tempr

  tString='tau_s2'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_som2=tempr

  tString='tau_s3'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  tau_decay_som3=tempr

  tString='froz_q10'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%froz_q10=tempr

  tString='decomp_depth_efolding'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%decomp_depth_efolding=tempr

  tString='q10_hr'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%Q10=tempr

  tString='minpsi_hr'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%minpsi_bgc=tempr

  tString='k_m_o2'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_m_o2_bgc=tempr

  tString='organic_max'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%organic_max=tempr

  tString='k_nitr_max'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%k_nitr_max=tempr

  call ncd_io('vmax_minp_soluble_to_secondary',this%vmax_minp_soluble_to_secondary, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=' ERROR: error in reading in soil order vmax_minp_soluble_to_secondary'//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return

  call ncd_io('minp_secondary_decay',this%minp_secondary_decay, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=' ERROR: error in reading in soil order minp_secondary_decay'//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return

  call ncd_io('frac_p_sec_to_sol',this%frac_p_sec_to_sol, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=' ERROR: error in reading in soil order frac_p_sec_to_sol'//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return

  this%k_decay_lit1          = 1._r8/(tau_decay_lit1*year_sec)    !1/second
  this%k_decay_lit2          = 1._r8/(tau_decay_lit2*year_sec)    !1/second
  this%k_decay_lit3          = 1._r8/(tau_decay_lit3*year_sec)    !1/second
  this%k_decay_som1          = 1._r8/(tau_decay_som1*year_sec)    !1/second
  this%k_decay_som2          = 1._r8/(tau_decay_som2*year_sec)    !1/second
  this%k_decay_som3          = 1._r8/(tau_decay_som3*year_sec)    !1/second
  this%k_decay_cwd           = 1._r8/(tau_decay_cwd*year_sec)     !1/second
  this%k_decay_fwd           = 1._r8/(tau_decay_fwd*year_sec)     !1/second
  this%k_decay_lwd           = 1._r8/(tau_decay_lwd*year_sec)     !1/second

  if(betr_spinup_state/=0)then
    call this%apply_spinup_factor()
  endif

  end subroutine readPars

!--------------------------------------------------------------------
  subroutine set_spinup_factor(this)

  implicit none
  class(BiogeoCon_type), intent(inout) :: this

  !the order is, lit1, lit2, lit3, cwd, lwd, fwd, som1, som3, som2
  this%spinup_factor(1) = 1._r8
  this%spinup_factor(2) = 1._r8
  this%spinup_factor(3) = 1._r8

  this%spinup_factor(4) = 1._r8
  this%spinup_factor(5) = 1._r8
  this%spinup_factor(6) = 1._r8

  this%spinup_factor(7) = this%k_decay_lit1/this%k_decay_som1
  this%spinup_factor(8) = this%k_decay_lit1/this%k_decay_som3
  this%spinup_factor(9) = this%k_decay_lit1/this%k_decay_som2

  end subroutine set_spinup_factor
end module BiogeoConType
