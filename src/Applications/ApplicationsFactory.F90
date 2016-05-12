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

  subroutine create_betr_usr_application(bgc_reaction, plant_soilbgc, method)
  !DESCRIPTION
  !create betr applications
  !
  implicit none
  class(bgc_reaction_type),  allocatable, intent(out) :: bgc_reaction
  class(plant_soilbgc_type), allocatable, intent(out) :: plant_soilbgc
  character(len=*),                       intent(in)  :: method

  allocate(bgc_reaction, source=create_bgc_reaction_type(method))

  allocate(plant_soilbgc, source=create_plant_soilbgc_type(method))

  end subroutine create_betr_usr_application

!-------------------------------------------------------------------------------

  function create_bgc_reaction_type(method) result(bgc_reaction)
    !
    ! !DESCRIPTION:
    ! create and return an object of bgc_reaction
    !
    ! !USES:
    use BGCReactionsMod , only : bgc_reaction_type
    use babortutils     , only : endrun
    use betr_ctrl       , only : iulog  => biulog

    ! !ARGUMENTS:
    character(len=*), intent(in)          :: method
    ! temporary varaibles
    class(bgc_reaction_type), allocatable :: bgc_reaction
    character(len=*), parameter           :: subname = 'create_bgc_reaction_type'

    select case(trim(method))

    case default

       write(iulog,*)subname //' ERROR: unknown method: ', method
       call endrun(msg=errMsg(mod_filename, __LINE__))
    end select
  end function create_bgc_reaction_type
  !-------------------------------------------------------------------------------

  function create_plant_soilbgc_type(method)result(plant_soilbgc)

  !DESCRIPTION
  !create and return an object of plant_soilbgc
  !
  !USES
  use PlantSoilBGCMod , only : plant_soilbgc_type
  use babortutils     , only : endrun
  use betr_ctrl       , only : iulog  => biulog
  implicit none
  ! !ARGUMENTS:
  character(len=*), intent(in)          :: method
  !temporary variables
  class(plant_soilbgc_type) , allocatable :: plant_soilbgc
  character(len=*)          , parameter   :: subname = 'create_plant_soilbgc_type'

  select case(trim(method))

  case default
     write(*, *)subname //' ERROR: unknown method: ', method
     call endrun(msg=errMsg(mod_filename, __LINE__))
  end select

  end function create_plant_soilbgc_type
end module ApplicationsFactory
