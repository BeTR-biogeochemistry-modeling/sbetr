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
    use betr_ctrl       , only : iulog  => biulog
    use betr_constants  , only : betr_errmsg_len
    use BetrStatusType  , only : betr_status_type
    use ecacnpBGCReactionsType, only : ecacnp_bgc_reaction_type
    use cdomBGCReactionsType, only : cdom_bgc_reaction_type
    use SimicBGCReactionsType, only : simic_bgc_reaction_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_type),  allocatable, intent(inout) :: bgc_reaction
    character(len=*), intent(in)          :: method
    type(betr_status_type), intent(out)   :: bstatus
    logical,                intent(out)   :: asoibgc
    character(len=*), parameter           :: subname = 'create_bgc_reaction_type'
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    select case(trim(method))
    case ("ecacnp")
       asoibgc=.true.
       allocate(bgc_reaction, source=ecacnp_bgc_reaction_type())
    case ("cdom")
       asoibgc=.true.
       allocate(bgc_reaction, source=cdom_bgc_reaction_type())
    case ("simic")
       asoibgc=.true.
       allocate(bgc_reaction, source=simic_bgc_reaction_type())
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
  use ecacnpPlantSoilBGCType, only : ecacnp_plant_soilbgc_type
  use cdomPlantSoilBGCType, only : cdom_plant_soilbgc_type
  use SimicPlantSoilBGCType , only : simic_plant_soilbgc_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_type), allocatable, intent(inout) :: plant_soilbgc
  character(len=*), intent(in)          :: method
  type(betr_status_type), intent(out)   :: bstatus

  character(len=*)          , parameter   :: subname = 'create_plant_soilbgc_type'
  character(len=betr_errmsg_len) :: msg

  call bstatus%reset()

  select case(trim(method))
  case ("ecacnp")
     allocate(plant_soilbgc, source=ecacnp_plant_soilbgc_type())
  case ("cdom")
     allocate(plant_soilbgc, source=cdom_plant_soilbgc_type())
  case ("simic")
     allocate(plant_soilbgc, source=simic_plant_soilbgc_type())
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
  use ecacnpParaType   , only : ecacnp_para
  use cdomParaType   , only : cdom_para
  use simicParaType  , only : simic_para
  use tracer_varcon  , only : reaction_method
  use ncdio_pio      , only : file_desc_t
  use BetrStatusType , only : betr_status_type
  implicit none
  type(file_desc_t), intent(inout)  :: ncid
  type(betr_status_type) , intent(out) :: bstatus

   select case (trim(reaction_method))
   case ("ecacnp")
     call ecacnp_para%readPars(ncid, bstatus)
   case ("cdom")
     call cdom_para%readPars(ncid, bstatus)
   case ("simic")
     call simic_para%readPars(ncid, bstatus)
   case default
     !do nothing
   end select

  end subroutine  AppLoadParameters

  !-------------------------------------------------------------------------------
  subroutine AppInitParameters(reaction_method, bstatus)
  !
  ! DESCRIPTION
  ! read in the parameters for specified bgc implementation
  use ecacnpParaType   , only : ecacnp_para
  use cdomParaType   , only : cdom_para
  use simicParaType  , only : simic_para
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  implicit none
  character(len=*), intent(in) :: reaction_method
  type(betr_status_type), intent(out)   :: bstatus
  character(len=255) :: msg

   call bstatus%reset()

   select case (trim(reaction_method))
   case ("ecacnp")
     call  ecacnp_para%Init(bstatus)
     !do nothing
   case ("cdom")
     call cdom_para%Init(bstatus)
   case ("simic")
     call simic_para%Init(bstatus)
   case default
     !do nothing
   end select

  end subroutine  AppInitParameters
  !-------------------------------------------------------------------------------
  subroutine AppSetSpinup()

  use ecacnpParaType  , only : ecacnp_para
  use cdomParaType    , only : cdom_para
  use tracer_varcon  , only : reaction_method
  implicit none

  select case (trim(reaction_method))
  case ("ecacnp")
     call  ecacnp_para%set_spinup_factor()
  case ("cdom")
     call cdom_para%set_spinup_factor()
  end select

  end subroutine AppSetSpinup
end module ApplicationsFactory
