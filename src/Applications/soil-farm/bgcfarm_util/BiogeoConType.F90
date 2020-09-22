module BiogeoConType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
implicit none

 private
  character(len=*), private, parameter :: filename = &
       __FILE__
 type, public :: BiogeoCon_type


  real(r8) :: cwd_fcel_bgc
  real(r8) :: cwd_flig_bgc
  real(r8) :: lwd_fcel_bgc
  real(r8) :: lwd_flig_bgc
  real(r8) :: fwd_fcel_bgc
  real(r8) :: fwd_flig_bgc

  real(r8) :: init_cn_met
  real(r8) :: init_cn_cel
  real(r8) :: init_cn_lig
  real(r8) :: init_cn_cwd
  real(r8) :: init_cn_lwd
  real(r8) :: init_cn_fwd

  real(r8) :: init_cp_met
  real(r8) :: init_cp_cel
  real(r8) :: init_cp_lig
  real(r8) :: init_cp_cwd
  real(r8) :: init_cp_lwd
  real(r8) :: init_cp_fwd

  real(r8) :: init_cc13_met
  real(r8) :: init_cc13_cel
  real(r8) :: init_cc13_lig
  real(r8) :: init_cc13_cwd
  real(r8) :: init_cc13_lwd
  real(r8) :: init_cc13_fwd

  real(r8) :: init_cc14_met
  real(r8) :: init_cc14_cel
  real(r8) :: init_cc14_lig
  real(r8) :: init_cc14_cwd
  real(r8) :: init_cc14_lwd
  real(r8) :: init_cc14_fwd

  logical :: use_c13
  logical :: use_c14
  logical :: nop_limit                                              !switch for P limitation
  logical :: non_limit
  !ECA nutrient competition
  real(r8), pointer :: vmax_minp_soluble_to_secondary(:)  => null() !maximum conversion rate of soluble P into secondary P

  !inorganic phosphorus cycling
  real(r8), pointer :: frac_p_sec_to_sol(:)    => null()   !fraction of released secondary phosphorus that goes into soluble form
  real(r8), pointer :: minp_secondary_decay(:) => null()   !decay rate of secondary phosphorus

  real(r8), pointer :: E_weath(:)=> null()
  real(r8) :: T_ref_weath
  real(r8), pointer :: b_weath(:)=> null()
  real(r8), pointer :: f_shield(:) => null()
  real(r8), pointer :: P_weip(:) => null()

 contains
   procedure, public  :: Init => bgc_con_init
   procedure, public  :: readPars => readPars_bcon
   procedure, public  :: printPars=> printPars_bcon
   procedure, public  :: bcon_Init
   procedure, private :: InitAllocate_bcon
   procedure, private :: set_defpar_default_base
   procedure, public  :: readPars_bgc
   procedure, public  :: prtPars_bgc
   procedure, public  :: deep_copy_bgc
 end type BiogeoCon_type

contains

  !--------------------------------------------------------------------
  subroutine bgc_con_init(this,  bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  use tracer_varcon  , only : use_c13_betr, use_c14_betr, is_nitrogen_active, is_phosphorus_active
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  type(betr_status_type)                   , intent(out) :: bstatus

  call bstatus%reset()

  end subroutine bgc_con_init
  !--------------------------------------------------------------------

  subroutine readPars_bcon(this, ncid, bstatus)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use ncdio_pio       , only : file_desc_t, ncd_io
  use BetrStatusType  , only : betr_status_type
  use betr_ctrl       , only : betr_spinup_state
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  type(file_desc_t)    , intent(inout)  :: ncid  ! pio netCDF file id
  type(betr_status_type) , intent(out) :: bstatus

  call this%readPars_bgc(ncid, bstatus)
  end subroutine readPars_bcon
  !--------------------------------------------------------------------
  subroutine bcon_Init(this,  bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : betr_spinup_state
  use tracer_varcon  , only : use_c13_betr, use_c14_betr, is_nitrogen_active, is_phosphorus_active
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  type(betr_status_type)                   , intent(out) :: bstatus

  call bstatus%reset()

  call this%InitAllocate_bcon()
  call this%set_defpar_default_base()

  this%use_c13 = use_c13_betr
  this%use_c14 = use_c14_betr
  this%nop_limit=.not. is_phosphorus_active
  this%non_limit=.not. is_nitrogen_active

  this%cwd_fcel_bgc          = 0.76_r8
  this%cwd_flig_bgc          = 0.24_r8

  end subroutine bcon_Init
  !--------------------------------------------------------------------
  subroutine InitAllocate_bcon(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder

  implicit none
  class(BiogeoCon_type), intent(inout) :: this

  allocate(this%minp_secondary_decay(0:betr_max_soilorder))
  allocate(this%vmax_minp_soluble_to_secondary(0:betr_max_soilorder))
  allocate(this%frac_p_sec_to_sol(0:betr_max_soilorder))

  allocate(this%E_weath(1:16))
  allocate(this%b_weath(1:16))
  allocate(this%f_shield(1:16))
  allocate(this%P_weip(1:16))

  !the following will be actually calculated from CNP bgc
  end subroutine InitAllocate_bcon
  !--------------------------------------------------------------------
  subroutine set_defpar_default_base(this)
  use bshr_const_mod  , only : year_sec=>SHR_CONST_YEARSECS
  implicit none
  class(BiogeoCon_type), intent(inout) :: this

  this%init_cn_met  = 90._r8  !mass based
  this%init_cn_cel  = 90._r8  !mass based
  this%init_cn_lig  = 90._r8  !mass based
  this%init_cn_cwd  = 90._r8  !mass based
  this%init_cn_fwd  = 90._r8  !mass based
  this%init_cn_lwd  = 90._r8  !mass based


  this%init_cp_met  = 900._r8
  this%init_cp_cel  = 900._r8
  this%init_cp_lig  = 900._r8
  this%init_cp_cwd  = 900._r8
  this%init_cp_lwd  = 4500._r8!mass based
  this%init_cp_fwd  = 4500._r8!mass based

  !ECA nutrient competition
  this%vmax_minp_soluble_to_secondary(:) = 1.e-7_r8  !1/s on the order of 1.e-7 1/s
  !inorganic phosphorus cycling
  !Note: (1._r8-frac_p_sec_to_sol)*minp_secondary_decay = occlusion rate
  this%frac_p_sec_to_sol(:)              = 0.95_r8    !fraction of released secondary phosphorus that goes into soluble form, occ rate is on the order of  1.e-13 1/s
  this%minp_secondary_decay(:)           = 2.e-8_r8   !decay rate of secondary phosphorus, 1/s, which gives occlusion rate about 1./(25 yrs)

  this%use_c13 = .false.
  this%use_c14 = .false.
  this%nop_limit=.false.
  this%non_limit=.false.
  this%init_cc13_met = 0._r8
  this%init_cc13_cel = 0._r8
  this%init_cc13_lig = 0._r8
  this%init_cc13_cwd = 0._r8

  this%init_cc14_met = 0._r8
  this%init_cc14_cel = 0._r8
  this%init_cc14_lig = 0._r8
  this%init_cc14_cwd = 0._r8

  end subroutine set_defpar_default_base

  !--------------------------------------------------------------------

  subroutine readPars_bgc(this, ncid, bstatus)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use ncdio_pio       , only : file_desc_t, ncd_io
  use BetrStatusType  , only : betr_status_type
  use betr_ctrl       , only : betr_spinup_state, bgc_type
  use betr_varcon     , only : betr_maxpatch_pft, betr_max_soilorder

  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  type(file_desc_t)    , intent(inout)  :: ncid  ! pio netCDF file id
  type(betr_status_type) , intent(out) :: bstatus

  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr(1) ! temporary to read in constant
  real(r8)           :: temparr(1:25)
  real(r8)           :: temparrs(1:16)
  character(len=100) :: tString ! temp. var for reading

  call bstatus%reset()
  if(index(bgc_type,'type1_bgc')/=0)return
  tString='cwd_fcel'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//' '//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%cwd_fcel_bgc=tempr(1)

  tString='cwd_flig'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) call bstatus%set_msg(msg=trim(errCode)//trim(tString)//' '//errMsg(__FILE__, __LINE__), err=-1)
  if(bstatus%check_status())return
  this%cwd_flig_bgc=tempr(1)

  if(index(bgc_type,'type1_bgc')/=0)return

  tString='lwd_fcel'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) then ! call bstatus%set_msg(msg=trim(errCode)//trim(tString)//' '//errMsg(__FILE__, __LINE__), err=-1)
    this%lwd_fcel_bgc=this%cwd_fcel_bgc
  else
    this%lwd_fcel_bgc=tempr(1)
  endif

  tString='fwd_fcel'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) then !call bstatus%set_msg(msg=trim(errCode)//trim(tString)//' '//errMsg(__FILE__, __LINE__), err=-1)
    this%fwd_fcel_bgc=this%cwd_fcel_bgc
  else
    this%fwd_fcel_bgc=tempr(1)
  endif

  tString='lwd_flig'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) then !call bstatus%set_msg(msg=trim(errCode)//trim(tString)//' '//errMsg(__FILE__, __LINE__), err=-1)
    this%lwd_flig_bgc=this%cwd_flig_bgc
  else
    this%lwd_flig_bgc=tempr(1)
  endif

  tString='fwd_flig'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( .not. readv ) then !call bstatus%set_msg(msg=trim(errCode)//trim(tString)//' '//errMsg(__FILE__, __LINE__), err=-1)
    this%fwd_flig_bgc=this%cwd_flig_bgc
  else
    this%fwd_flig_bgc=tempr(1)
  endif

  call ncd_io('vmax_minp_soluble_to_secondary',temparr, 'read', ncid, readvar=readv)
  if (readv ) then !call bstatus%set_msg(msg=' ERROR: error in reading in soil order vmax_minp_soluble_to_secondary '//errMsg(__FILE__, __LINE__), err=-1)
    this%vmax_minp_soluble_to_secondary(1:betr_max_soilorder)=temparr(1:betr_max_soilorder)
  endif

  call ncd_io('minp_secondary_decay',temparr, 'read', ncid, readvar=readv)
  if ( readv ) then !call bstatus%set_msg(msg=' ERROR: error in reading in soil order minp_secondary_decay '//errMsg(__FILE__, __LINE__), err=-1)
    this%minp_secondary_decay(1:betr_max_soilorder)=temparr(1:betr_max_soilorder)
  endif

  call ncd_io('frac_p_sec_to_sol',temparr, 'read', ncid, readvar=readv)
  if ( readv ) then !call bstatus%set_msg(msg=' ERROR: error in reading in soil order frac_p_sec_to_sol '//errMsg(__FILE__, __LINE__), err=-1)
    this%frac_p_sec_to_sol(1:betr_max_soilorder)=temparr(1:betr_max_soilorder)
  endif

  call ncd_io('E_weath',temparrs, 'read', ncid, readvar=readv)
  if ( readv ) then !call bstatus%set_msg(msg=' ERROR: error in reading in soil order E_weath '//errMsg(__FILE__, __LINE__), err=-1)
    this%E_weath(1:16) = temparrs(1:16)
  endif

  tString='T_ref_weath'
  call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
  if ( readv ) then !call bstatus%set_msg(msg=trim(errCode)//trim(tString)//' '//errMsg(__FILE__, __LINE__), err=-1)
    this%T_ref_weath=tempr(1)
  endif

  call ncd_io('b_weath',temparrs, 'read', ncid, readvar=readv)
  if ( readv ) then !call bstatus%set_msg(msg=' ERROR: error in reading in soil order b_weath'//errMsg(__FILE__, __LINE__), err=-1)
    this%b_weath(1:16) = temparrs(1:16)
  endif

  call ncd_io('f_shield',temparrs, 'read', ncid, readvar=readv)
  if ( readv ) then !call bstatus%set_msg(msg=' ERROR: error in reading in soil order f_shield'//errMsg(__FILE__, __LINE__), err=-1)
    this%f_shield(1:16) = temparrs(1:16)
  endif

  call ncd_io('P_weip',temparrs, 'read', ncid, readvar=readv)
  if ( readv ) then !call bstatus%set_msg(msg=' ERROR: error in reading in soil order P_weip'//errMsg(__FILE__, __LINE__), err=-1)
    this%P_weip(1:16) = temparrs(1:16)
  endif
  end subroutine readPars_bgc

  !--------------------------------------------------------------------
  subroutine prtPars_bgc(this)
  use betr_varcon     , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(BiogeoCon_type), intent(inout) :: this

  integer :: j
  print*,'init_cn_met=',this%init_cn_met
  print*,'init_cn_cel=',this%init_cn_cel
  print*,'init_cn_lig=',this%init_cn_lig
  print*,'init_cn_cwd=',this%init_cn_cwd
  print*,'init_cn_fwd=',this%init_cn_fwd
  print*,'init_cn_lwd=',this%init_cn_lwd

  print*,'init_cp_met=',this%init_cp_met
  print*,'init_cp_cel=',this%init_cp_cel
  print*,'init_cp_lig=',this%init_cp_lig
  print*,'init_cp_cwd=',this%init_cp_cwd
  print*,'init_cp_lwd=',this%init_cp_lwd
  print*,'init_cp_fwd=',this%init_cp_fwd

  !ECA nutrient competition
  print*,'vmax_minp_soluble_to_secondary=',(this%vmax_minp_soluble_to_secondary(j),j=1,betr_max_soilorder)
  print*,'frac_p_sec_to_sol=',(this%frac_p_sec_to_sol(j),j=1,betr_max_soilorder)
  print*,'minp_secondary_decay=',(this%minp_secondary_decay(j),j=1,betr_max_soilorder)


  end subroutine prtPars_bgc

  !--------------------------------------------------------------------
  subroutine printPars_bcon(this)
  use betr_varcon     , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(BiogeoCon_type), intent(inout) :: this

  call this%prtPars_bgc()

  end subroutine printPars_bcon
  !--------------------------------------------------------------------
  subroutine deep_copy_bgc(this,mother)
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  class(BiogeoCon_type), intent(in) :: mother


  this%cwd_fcel_bgc = mother%cwd_fcel_bgc
  this%cwd_flig_bgc = mother%cwd_flig_bgc
  this%lwd_fcel_bgc = mother%lwd_fcel_bgc
  this%lwd_flig_bgc = mother%lwd_flig_bgc
  this%fwd_fcel_bgc = mother%fwd_fcel_bgc
  this%fwd_flig_bgc = mother%fwd_flig_bgc

  this%init_cn_met  = mother%init_cn_met
  this%init_cn_cel  = mother%init_cn_cel
  this%init_cn_lig  = mother%init_cn_lig
  this%init_cn_cwd  = mother%init_cn_cwd
  this%init_cn_lwd  = mother%init_cn_lwd
  this%init_cn_fwd  = mother%init_cn_fwd

  this%init_cp_met  = mother%init_cp_met
  this%init_cp_cel  = mother%init_cp_cel
  this%init_cp_lig  = mother%init_cp_cel
  this%init_cp_cwd  = mother%init_cp_cwd
  this%init_cp_lwd  = mother%init_cp_lwd
  this%init_cp_fwd  = mother%init_cp_fwd

  this%init_cc13_met= mother%init_cc13_met
  this%init_cc13_cel= mother%init_cc13_cel
  this%init_cc13_lig= mother%init_cc13_lig
  this%init_cc13_cwd= mother%init_cc13_cwd
  this%init_cc13_lwd= mother%init_cc13_lwd
  this%init_cc13_fwd= mother%init_cc13_fwd

  this%init_cc14_met= mother%init_cc14_met
  this%init_cc14_cel= mother%init_cc14_cel
  this%init_cc14_lig= mother%init_cc14_lig
  this%init_cc14_cwd= mother%init_cc14_cwd
  this%init_cc14_lwd= mother%init_cc14_lwd
  this%init_cc14_fwd= mother%init_cc14_fwd

  this%use_c13 = mother%use_c13
  this%use_c14 = mother%use_c14
  this%nop_limit = mother%nop_limit                                              !switch for P limitation
  this%non_limit = mother%non_limit
  !ECA nutrient competition
  this%vmax_minp_soluble_to_secondary(:) = mother%vmax_minp_soluble_to_secondary(:)
  !inorganic phosphorus cycling
  this%frac_p_sec_to_sol(:) = mother%frac_p_sec_to_sol(:)
  this%minp_secondary_decay(:) = mother%minp_secondary_decay(:)
  this%E_weath(:) = mother%E_weath(:)
  this%T_ref_weath = mother%T_ref_weath
  this%b_weath(:) = mother%b_weath(:)
  this%f_shield(:) = mother%f_shield(:)
  this%P_weip(:) = mother%P_weip(:)
  end subroutine deep_copy_bgc


end module BiogeoConType
