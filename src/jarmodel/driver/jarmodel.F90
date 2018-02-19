program main

  use BeTRJarModel     , only : jar_model_type
  use JarModelFactory  , only : create_jar_model, create_jar_pars
  use JarBgcForcType   , only : JarBGC_forc_type
  use BetrStatusType   , only : betr_status_type
  use betr_ctrl        , only : iulog => biulog
  use BiogeoConType    , only : BiogeoCon_type
  use abortutils       , only : endrun
implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  class(jar_model_type),  pointer :: jarmodel
  class(BiogeoCon_type),  pointer :: jarpars
  type(JarBGC_forc_type) :: bgc_forc
  type(betr_status_type) :: bstatus
  character(len=24) :: jarmtype
  character(len=*), parameter :: namelist_buffer="junk_data"

  jarmtype = 'ecacnp'
  !create the jar type

  jarmodel => create_jar_model(jarmtype)
  jarpars  => create_jar_pars(jarmtype)

  !initialize model parameters
  call jarpars%Init(namelist_buffer, bstatus)
  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif

  !initialize model
  call jarmodel%init(jarpars,  bstatus)
  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif

  !initialize the parameters
  call jarmodel%UpdateParas(jarpars, bstatus)
  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif


end program main
