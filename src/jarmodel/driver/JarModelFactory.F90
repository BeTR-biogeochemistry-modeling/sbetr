module JarModelFactory
  !
  ! !DESCRIPTION:
  !  factory to load the specific betr jarmodel
  !
  ! !USES:
  !
  use abortutils  , only : endrun
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use  betr_ctrl  , only : iulog => biulog
  implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  public :: create_jar_model
  public :: create_jar_pars
contains

!-------------------------------------------------------------------------------
  function create_jar_model(jarmodel_name)result(jarmodel)
    !DESCRIPTION
    ! create a betr simulation object
    !
    !USES
    use ecacnpBGCType, only : create_jarmodel_ecacnp
    use BeTRJarModel      , only : jar_model_type
    implicit none
    !ARGUMENTS
    character(len=*), intent(in)  :: jarmodel_name
    class(jar_model_type), pointer :: jarmodel

    select case(trim(jarmodel_name))
       case ("ecacnp")
          allocate(jarmodel, source=create_jarmodel_ecacnp())
       case default
          write(iulog, *) "ERROR: unknown jarmodel type '", &
               trim(jarmodel_name), "'."
          call endrun(msg=errMsg(mod_filename, __LINE__))
    end select
  end function create_jar_model

!-------------------------------------------------------------------------------
  function create_jar_pars(jarmodel_name)result(jarpars)
    !DESCRIPTION
    ! create a betr simulation object
    !
    !USES
    use ecacnpBGCType, only : create_jarmodel_ecacnp
    use BiogeoConType , only : BiogeoCon_type
    use ecacnpParaType  , only : create_jarpars_ecacnp
    implicit none
    !ARGUMENTS
    character(len=*), intent(in)  :: jarmodel_name
    class(BiogeoCon_type), pointer :: jarpars

    select case(trim(jarmodel_name))
       case ("ecacnp")
          allocate(jarpars, source=create_jarpars_ecacnp())
       case default
          write(iulog, *) "ERROR: unknown jarmodel type '", &
               trim(jarmodel_name), "'."
          call endrun(msg=errMsg(mod_filename, __LINE__))
    end select
  end function create_jar_pars
end module JarModelFactory
