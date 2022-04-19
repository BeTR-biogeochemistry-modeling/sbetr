module ApplicationsFactory
  !
  ! !DESCRIPTION:
  !  factory to load the specific bgc reaction modules
  !
  ! History:
  !  Created by Jinyun Tang, April 28, 2016
  !
  ! !USES:
  !
  use bshr_kind_mod          , only : r8 => shr_kind_r8
  use bshr_log_mod           , only : errMsg => shr_log_errMsg
  use BGCReactionsMod        , only : bgc_reaction_type
  use PlantSoilBGCMod        , only : plant_soilbgc_type

  implicit none
  character(len=*), parameter :: mod_filename = &
       __FILE__
  private
  public :: create_betr_usr_application
  public :: AppLoadParameters
  public :: AppInitParameters
  public :: AppSetSpinup
  public :: AppCopyParas
contains

  subroutine create_betr_usr_application(bgc_reaction, plant_soilbgc, method, asoibgc, bstatus)
  !DESCRIPTION
  !create betr applications
  !
  use betr_constants  , only : betr_errmsg_len
  use BetrStatusType  , only : betr_status_type
  implicit none
  class(bgc_reaction_type),  allocatable, intent(out) :: bgc_reaction
  class(plant_soilbgc_type), allocatable, intent(out) :: plant_soilbgc
  character(len=*),                       intent(in)  :: method
  logical,                                intent(out) :: asoibgc
  type(betr_status_type), intent(out) :: bstatus

  call create_bgc_reaction_type(bgc_reaction, method,asoibgc,bstatus)

  if(bstatus%check_status())return

  call create_plant_soilbgc_type(plant_soilbgc, method,bstatus)

  end subroutine create_betr_usr_application

!-------------------------------------------------------------------------------

  subroutine create_bgc_reaction_type(bgc_reaction, method, asoibgc, bstatus)
    !
    ! !DESCRIPTION:
    ! create and return an object of bgc_reaction
    !
    ! !USES:
    use BGCReactionsMod , only : bgc_reaction_type
    use betr_ctrl       , only : iulog  => biulog, inloop_reaction, bgc_type
    use betr_constants  , only : betr_errmsg_len
    use BetrStatusType  , only : betr_status_type
    !begin_appadd
#if (defined SBETR)
    use ecacnpBGCReactionsType, only : ecacnp_bgc_reaction_type
    use ch4soilBGCReactionsType, only : ch4soil_bgc_reaction_type
    use cdomBGCReactionsType  , only : cdom_bgc_reaction_type
    use simicBGCReactionsType , only : simic_bgc_reaction_type
    use kecaBGCReactionsType  , only : keca_bgc_reaction_type
#endif
    use v1ecaBGCReactionsType, only : v1eca_bgc_reaction_type
    !end_appadd

    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_type),  allocatable, intent(inout) :: bgc_reaction
    character(len=*), intent(in)          :: method
    type(betr_status_type), intent(out)   :: bstatus
    logical,                intent(out)   :: asoibgc
    character(len=*), parameter           :: subname = 'create_bgc_reaction_type'
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    asoibgc = .false.
    select case(trim(method))
    !begin_appadd
#if (defined SBETR)
    case ("ecacnp","ecacnp_mosart")
       asoibgc=.true.;allocate(bgc_reaction, source=ecacnp_bgc_reaction_type())
       bgc_type='type2_bgc'
    case ("ch4soil")
       asoibgc=.true.;allocate(bgc_reaction, source=ch4soil_bgc_reaction_type())
    case ("cdom","cdom_mosart")
       asoibgc=.true.;allocate(bgc_reaction, source=cdom_bgc_reaction_type())
    case ("simic")
       asoibgc=.true.;allocate(bgc_reaction, source=simic_bgc_reaction_type())
    case ("keca")
       asoibgc=.true.;allocate(bgc_reaction, source=keca_bgc_reaction_type())
#endif
    case ("v1eca")
       asoibgc=.true.;allocate(bgc_reaction, source=v1eca_bgc_reaction_type())
       inloop_reaction=.false.; bgc_type='type1_bgc'
    !end_appadd
    case default
       write(msg,*)subname //' ERROR: unknown method: ', method
       msg = trim(msg)//new_line('A')//errMsg(mod_filename, __LINE__)
       call bstatus%set_msg(msg=msg,err=-1)
    end select
  end subroutine create_bgc_reaction_type
  !-------------------------------------------------------------------------------

  subroutine create_plant_soilbgc_type(plant_soilbgc, method, bstatus)

  !DESCRIPTION
  !create and return an object of plant_soilbgc
  !
  !USES
  use PlantSoilBGCMod , only : plant_soilbgc_type
  use betr_ctrl       , only : iulog  => biulog
  use betr_constants  , only : betr_errmsg_len
  use BetrStatusType  , only : betr_status_type
  !begin_appadd
#if (defined SBETR)
  use ecacnpPlantSoilBGCType, only : ecacnp_plant_soilbgc_type
  use ch4soilPlantSoilBGCType, only : ch4soil_plant_soilbgc_type
  use cdomPlantSoilBGCType  , only : cdom_plant_soilbgc_type
  use simicPlantSoilBGCType , only : simic_plant_soilbgc_type
  use kecaPlantSoilBGCType  , only : keca_plant_soilbgc_type
#endif
  use v1ecaPlantSoilBGCType, only : v1eca_plant_soilbgc_type
  !end_appadd
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_type), allocatable, intent(inout) :: plant_soilbgc
  character(len=*), intent(in)          :: method
  type(betr_status_type), intent(out)   :: bstatus

  character(len=*)          , parameter   :: subname = 'create_plant_soilbgc_type'
  character(len=betr_errmsg_len) :: msg

  call bstatus%reset()

  select case(trim(method))
  !begin_appadd
#if (defined SBETR)
  case ("ecacnp","ecacnp_mosart")
     allocate(plant_soilbgc, source=ecacnp_plant_soilbgc_type())
  case ("ch4soil")
     allocate(plant_soilbgc, source=ch4soil_plant_soilbgc_type())
  case ("cdom","cdom_mosart")
     allocate(plant_soilbgc, source=cdom_plant_soilbgc_type())
  case ("simic")
     allocate(plant_soilbgc, source=simic_plant_soilbgc_type())
  case ("keca")
     allocate(plant_soilbgc, source=keca_plant_soilbgc_type())
#endif
  case ("v1eca","v1eca_mosart")
     allocate(plant_soilbgc, source=v1eca_plant_soilbgc_type())
  !end_appadd
  case default
     write(msg, *)subname //' ERROR: unknown method: ', method
     msg = trim(msg)//new_line('A')//errMsg(mod_filename, __LINE__)
     call bstatus%set_msg(msg=msg,err=-1)
  end select

  end subroutine create_plant_soilbgc_type

  !-------------------------------------------------------------------------------
  subroutine AppLoadParameters(ncid, bstatus)
  !
  ! DESCRIPTION
  ! read in the parameters for specified bgc implementation
  !begin_appadd
#if (defined SBETR)
  use ecacnpParaType   , only : ecacnp_para
  use ch4soilParaType   , only : ch4soil_para
  use cdomParaType     , only : cdom_para
  use simicParaType    , only : simic_para
  use kecaParaType     , only : keca_para
#endif
  use v1ecaParaType   , only : v1eca_para
  !end_appadd
  use tracer_varcon    , only : reaction_method
  use ncdio_pio        , only : file_desc_t
  use BetrStatusType   , only : betr_status_type
  implicit none
  type(file_desc_t), intent(inout)  :: ncid
  type(betr_status_type) , intent(out) :: bstatus

   select case (trim(reaction_method))
  !begin_appadd
#if (defined SBETR)
   case ("ecacnp","ecacnp_mosart")
     call ecacnp_para%readPars(ncid, bstatus)
   case ("ch4soil")
     call ch4soil_para%readPars(ncid, bstatus)
   case ("cdom","cdom_mosart")
     call cdom_para%readPars(ncid, bstatus)
   case ("simic")
     call simic_para%readPars(ncid, bstatus)
   case ("keca")
     call keca_para%readPars(ncid, bstatus)
#endif
   case ("v1eca","v1eca_mosart")
     call v1eca_para%readPars(ncid, bstatus)
   !end_appadd
   case default
     !do nothing
   end select

  end subroutine  AppLoadParameters

  !-------------------------------------------------------------------------------
  subroutine AppInitParameters(reaction_method, bstatus)
  !
  ! DESCRIPTION
  ! read in the parameters for specified bgc implementation
  !begin_appadd
#if (defined SBETR)
  use ecacnpParaType   , only : ecacnp_para
  use ch4soilParaType   , only : ch4soil_para
  use cdomParaType     , only : cdom_para
  use simicParaType    , only : simic_para
  use kecaParaType     , only : keca_para
#endif
  use v1ecaParaType   , only : v1eca_para
  !end_appadd
  use betr_constants   , only : betr_namelist_buffer_size_ext
  use BetrStatusType   , only : betr_status_type
  use tracer_varcon    , only : nparcols
  implicit none
  character(len=*), intent(in) :: reaction_method
  type(betr_status_type), intent(out)   :: bstatus
  character(len=255) :: msg

   call bstatus%reset()

   select case (trim(reaction_method))
   !begin_appadd
#if (defined SBETR)
   case ("ecacnp","ecacnp_mosart")
     nparcols=0
     call ecacnp_para%Init(bstatus)
   case ("ch4soil")
     call ch4soil_para%Init(bstatus)
   case ("cdom","cdom_mosart")
     call cdom_para%Init(bstatus)
   case ("simic")
     call simic_para%Init(bstatus)
   case ("keca")
     call keca_para%Init(bstatus)
#endif
   case ("v1eca")
     call v1eca_para%Init(bstatus)
   !end_appadd
   case default
     !do nothing
   end select

  end subroutine  AppInitParameters
  !-------------------------------------------------------------------------------
  subroutine AppSetSpinup()

  ! set spinup strategies
  !begin_appadd
#if (defined SBETR)
  use ecacnpParaType  , only : ecacnp_para
  use ch4soilParaType  , only : ch4soil_para
  use cdomParaType    , only : cdom_para
  use kecaParaType    , only : keca_para
#endif
  use v1ecaParaType  , only : v1eca_para
  !end_appadd
  use tracer_varcon   , only : reaction_method
  implicit none

  select case (trim(reaction_method))
  !begin_appadd
#if (defined SBETR)
  case ("ecacnp","ecacnp_mosart")
     call  ecacnp_para%set_spinup_factor()
  case ("ch4soil")
     call  ch4soil_para%set_spinup_factor()
  case ("cdom","cdom_mosart")
     call cdom_para%set_spinup_factor()
  case ("keca")
     call keca_para%set_spinup_factor()
#endif
  case ("v1eca","v1eca_mosart")
     call  v1eca_para%set_spinup_factor()
  !end_appadd
  end select

  end subroutine AppSetSpinup

  !-------------------------------------------------------------------------------

  subroutine AppCopyParas(begc, endc, bstatus)
  !
  ! DESCRIPTION
  ! read in the parameters for specified bgc implementation
  !begin_appadd
  use tracer_varcon    , only : reaction_method, nparcols
#if (defined SBETR)
  use ecacnpParaType   , only : ecacnp_para,ecacnp_paras
  use ch4soilParaType   , only : ch4soil_para
  use cdomParaType     , only : cdom_para
  use simicParaType    , only : simic_para
  use kecaParaType     , only : keca_para
#endif
  use v1ecaParaType   , only : v1eca_para
  !end_appadd
  use betr_constants   , only : betr_namelist_buffer_size_ext
  use BetrStatusType   , only : betr_status_type
  implicit none
  integer, intent(in) :: begc, endc
  type(betr_status_type), intent(out)   :: bstatus
  integer :: fl
   call bstatus%reset()

   select case (trim(reaction_method))
   !begin_appadd
#if (defined SBETR)
   case ("ecacnp","ecacnp_mosart")
    allocate(ecacnp_paras(begc:endc))
    nparcols=endc-begc+1
    do fl = begc, endc
      call ecacnp_paras(fl)%Init(bstatus)
      if(bstatus%check_status())return
      call ecacnp_paras(fl)%deep_copy(ecacnp_para)
    enddo
!   case ("ch4soil")
!     call ch4soil_para%Init(bstatus)
!   case ("cdom","cdom_mosart")
!     call cdom_para%Init(bstatus)
!   case ("simic")
!     call simic_para%Init(bstatus)
!   case ("keca")
!     call keca_para%Init(bstatus)
#endif
!   case ("v1eca")
!     call v1eca_para%Init(bstatus)
   !end_appadd
   case default
     !do nothing
   end select

  end subroutine AppCopyParas

end module ApplicationsFactory
