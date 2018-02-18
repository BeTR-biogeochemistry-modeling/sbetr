program main

  use BeTRJarModel   , only : jar_model_type
  use JarModelFactory, only : create_jar_model
implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  class(jar_model_type),  pointer :: jarmodel

  character(len=24) :: jarmtype

  jarmtype = 'ecacnp'
  !create the jar type

  jarmodel => create_jar_model(jarmtype)

  !initialize the parameters


end program main
