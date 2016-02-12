module BGCReactionsFactoryMod
  !
  ! !DESCRIPTION:
  !  factory to load the specific bgc reaction modules
  !
  ! History:
  !  Created by Jinyun Tang, Oct 2, 2014
  !
  !
  ! !USES:
  !
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use tracer_varcon         , only : is_active_betr_bgc
  implicit none
  save
  private

  public  :: create_betr_application
  private :: create_bgc_reaction_type
  private :: create_plant_soilbgc_type
contains

  subroutine create_betr_application(bgc_reaction, plant_soilbgc, method)


  class(bgc_reaction_type),  allocatable, intent(out) :: bgc_reaction
  class(plant_soilbgc_type), allocatable, intent(out) :: plant_soilbgc
  character(len=*), intent(in)          :: method


  allocate(bgc_reaction, source=create_bgc_reaction_type(method))

  allocate(plant_soilbgc, source=create_plant_soilbgc_type(method))

  end subroutine create_betr_application
!-------------------------------------------------------------------------------

  function create_bgc_reaction_type(method) result(bgc_reaction)
    !
    ! !DESCRIPTION:
    ! create and return an object of bgc_reaction
    !
    ! !USES:
    use BGCReactionsMod             , only : bgc_reaction_type
    use BGCReactionsMockRunType     , only : bgc_reaction_mock_run_type
    use BGCReactionsCentECACNPType  , only : bgc_reaction_CENTURY_ECACNP_type
    use abortutils                  , only : endrun
    use clm_varctl                  , only : iulog

    ! !ARGUMENTS:
    class(bgc_reaction_type), allocatable :: bgc_reaction
    character(len=*), intent(in)          :: method
    character(len=*), parameter           :: subname = 'create_bgc_reaction_type'

    select case(trim(method))
    case ("mock_run")
       allocate(bgc_reaction, source=bgc_reaction_mock_run_type())
    case ("century_bgcECA")
       is_active_betr_bgc = .true.
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_ECA_type())
    case default
       write(iulog,*)subname //' ERROR: unknown method: ', method
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
  end function create_bgc_reaction_type
  !-------------------------------------------------------------------------------

  function create_plant_soilbgc_type(method)result(plant_soilbgc)

  use PlantSoilBGCMod             , only : plant_soilbgc_type
  use PlantSoilBGCMockRunType     , only : plant_soilbgc_mock_run_type
  use plant_soilbgc_ECACPL        , only : plant_soilbgc_ECA_run_type
  use abortutils                  , only : endrun
  use clm_varctl                  , only : iulog


  ! !ARGUMENTS:
  class(plant_soilbgc_type), allocatable :: plant_soilbgc
  character(len=*), intent(in)          :: method
  character(len=*), parameter           :: subname = 'create_plant_soilbgc_type'


  select case(trim(method))
  case ("mock_run")
     allocate(plant_soilbgc, source=plant_soilbgc_mock_run_type())
  case default
     write(iulog,*)subname //' ERROR: unknown method: ', method
     call endrun(msg=errMsg(__FILE__, __LINE__))
  end select

  end function create_plant_soilbgc_type
end module BGCReactionsFactoryMod
