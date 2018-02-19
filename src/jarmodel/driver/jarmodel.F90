program main

  use BeTRJarModel     , only : jar_model_type
  use JarModelFactory  , only : create_jar_model
  use JarBgcForcType   , only : JarBGC_forc_type
  use BetrStatusType   , only : betr_status_type
implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  class(jar_model_type),  pointer :: jarmodel
  type(JarBGC_forc_type) :: bgc_forc
  type(betr_status_type) :: bstatus
  character(len=24) :: jarmtype

  jarmtype = 'ecacnp'
  !create the jar type

  jarmodel => create_jar_model(jarmtype)

!  call jarmodel%init(biogeo_con,  bstatus)

  !initialize the parameters
!  call jarmodel%UpdateParas(biogeo_con, bstatus)

end program main
