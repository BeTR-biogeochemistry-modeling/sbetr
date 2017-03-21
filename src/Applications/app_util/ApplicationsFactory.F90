module ApplicationsFactory
  !
  ! !DESCRIPTION:
  !  factory to load the specific bgc reaction modules
  !
  ! History:
  !  Created by Jinyun Tang, April 28, 2016
  !
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
contains

  subroutine create_betr_usr_application(bgc_reaction, plant_soilbgc, method, bstatus)
  !DESCRIPTION
  !create betr applications
  !
  use betr_constants  , only : betr_errmsg_len
  use BetrStatusType  , only : betr_status_type
  implicit none
  class(bgc_reaction_type),  allocatable, intent(out) :: bgc_reaction
  class(plant_soilbgc_type), allocatable, intent(out) :: plant_soilbgc
  character(len=*),                       intent(in)  :: method
  type(betr_status_type), intent(out) :: bstatus


  allocate(bgc_reaction, source=create_bgc_reaction_type(method,bstatus))
  if(bstatus%check_status())return

  allocate(plant_soilbgc, source=create_plant_soilbgc_type(method,bstatus))

  end subroutine create_betr_usr_application

!-------------------------------------------------------------------------------

  function create_bgc_reaction_type(method, bstatus) result(bgc_reaction)
    !
    ! !DESCRIPTION:
    ! create and return an object of bgc_reaction
    !
    ! !USES:
    use BGCReactionsMod , only : bgc_reaction_type
    use betr_ctrl       , only : iulog  => biulog
    use betr_constants  , only : betr_errmsg_len
    use BetrStatusType  , only : betr_status_type
    use BGCReactionsCentECACnpType, only : bgc_reaction_CENTURY_ECACNP_type
    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(in)          :: method
    type(betr_status_type), intent(out)   :: bstatus

    ! temporary varaibles
    class(bgc_reaction_type), allocatable :: bgc_reaction
    character(len=*), parameter           :: subname = 'create_bgc_reaction_type'
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    select case(trim(method))
    case ("eca_cnp")
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_ECACNP_type())
    case default
       write(msg,*)subname //' ERROR: unknown method: ', method
       msg = trim(msg)//new_line('A')//errMsg(mod_filename, __LINE__)
       call bstatus%set_msg(msg=msg,err=-1)
    end select
  end function create_bgc_reaction_type
  !-------------------------------------------------------------------------------

  function create_plant_soilbgc_type(method, bstatus)result(plant_soilbgc)

  !DESCRIPTION
  !create and return an object of plant_soilbgc
  !
  !USES
  use PlantSoilBGCMod , only : plant_soilbgc_type
  use betr_ctrl       , only : iulog  => biulog
  use betr_constants  , only : betr_errmsg_len
  use BetrStatusType  , only : betr_status_type
  use PlantSoilBgcCnpType, only : plant_soilbgc_cnp_type
!  use PlantSoilBgcDcnpType, only : plant_soilbgc_dcnp_type

  implicit none
  ! !ARGUMENTS:
  character(len=*), intent(in)          :: method
  type(betr_status_type), intent(out)   :: bstatus
  !temporary variables
  class(plant_soilbgc_type) , allocatable :: plant_soilbgc
  character(len=*)          , parameter   :: subname = 'create_plant_soilbgc_type'
  character(len=betr_errmsg_len) :: msg

  call bstatus%reset()

  select case(trim(method))
  case ("eca_cnp")
     allocate(plant_soilbgc, source=plant_soilbgc_cnp_type())
!  case ("eca_dcnp")
!     allocate(plant_soilbgc, source=plant_soilbgc_dcnp_type())
  case default
     write(msg, *)subname //' ERROR: unknown method: ', method
     msg = trim(msg)//new_line('A')//errMsg(mod_filename, __LINE__)
     call bstatus%set_msg(msg=msg,err=-1)
  end select

  end function create_plant_soilbgc_type
end module ApplicationsFactory
