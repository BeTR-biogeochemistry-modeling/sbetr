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
    !begin_appadd
    use ch4soilBGCType, only : create_jarmodel_ch4soil
    use ecacnpBGCType, only : create_jarmodel_ecacnp
    use simicBGCType,  only : create_jarmodel_simicbgc
    use cdomBGCType, only : create_jarmodel_cdom
    use v1ecaBGCType, only : create_jarmodel_v1eca
    !end_appadd
    use BeTRJarModel      , only : jar_model_type
    implicit none
    !ARGUMENTS
    character(len=*), intent(in)  :: jarmodel_name
    class(jar_model_type), pointer :: jarmodel

    select case(trim(jarmodel_name))
    !begin_appadd
       case ("ch4soil")
          allocate(jarmodel, source=create_jarmodel_ch4soil())
       case ("ecacnp")
          allocate(jarmodel, source=create_jarmodel_ecacnp())
       case ("simic")
          allocate(jarmodel, source=create_jarmodel_simicbgc())
       case ("cdom")
          allocate(jarmodel, source=create_jarmodel_cdom())
       case ("v1eca")
          allocate(jarmodel, source=create_jarmodel_v1eca())
    !end_appadd
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
    use BiogeoConType , only : BiogeoCon_type
    !begin_appadd
    use ch4soilParaType  , only : create_jarpars_ch4soil
    use ecacnpBGCType, only : create_jarmodel_ecacnp
    use ecacnpParaType  , only : create_jarpars_ecacnp
    use simicParaType   , only : create_jarpars_simic
    use cdomParaType   , only : create_jarpars_cdom
    use v1ecaParaType  , only : create_jarpars_v1eca
    !end_appadd
    implicit none
    !ARGUMENTS
    character(len=*), intent(in)  :: jarmodel_name
    class(BiogeoCon_type), pointer :: jarpars

    select case(trim(jarmodel_name))
    !begin_appadd
       case ("ch4soil")
          allocate(jarpars, source=create_jarpars_ch4soil())
       case ("ecacnp")
          allocate(jarpars, source=create_jarpars_ecacnp())
       case ("simic")
          allocate(jarpars, source=create_jarpars_simic())
       case ("cdom")
          allocate(jarpars, source=create_jarpars_cdom())
       case ("v1eca")
          allocate(jarpars, source=create_jarpars_v1eca())
    !end_appadd
       case default
          write(iulog, *) "ERROR: unknown jarmodel type '", &
               trim(jarmodel_name), "'."
          call endrun(msg=errMsg(mod_filename, __LINE__))
    end select
  end function create_jar_pars
end module JarModelFactory
