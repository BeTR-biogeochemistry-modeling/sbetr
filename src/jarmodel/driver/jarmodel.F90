program main

  use BeTRJarModel     , only : jar_model_type
  use JarModelFactory  , only : create_jar_model, create_jar_pars
  use JarBgcForcType   , only : JarBGC_forc_type
  use BetrStatusType   , only : betr_status_type
  use betr_ctrl        , only : iulog => biulog
  use BiogeoConType    , only : BiogeoCon_type
  use abortutils       , only : endrun
  use histMod          , only : histf_type
  use histMod          , only : hist_var_str_len, hist_unit_str_len, hist_freq_str_len
  use  BeTR_TimeMod    , only : betr_time_type
  use shr_kind_mod     , only : r8 => shr_kind_r8
implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  class(jar_model_type),  pointer :: jarmodel
  class(BiogeoCon_type),  pointer :: jarpars
  type(JarBGC_forc_type) :: bgc_forc
  type(betr_status_type) :: bstatus
  type(histf_type) :: hist
  type(betr_time_type) :: timer
  character(len=24) :: jarmtype
  character(len=*), parameter :: namelist_buffer="junk_data"
  character(len=hist_var_str_len) , allocatable :: varl(:)
  character(len=hist_unit_str_len), allocatable :: unitl(:)
  character(len=hist_freq_str_len), allocatable :: freql(:)
  real(r8) :: dt
  integer :: nvars

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

  nvars=jarmodel%getvarllen()
  allocate(varl(nvars)); allocate(unitl(nvars)); allocate(freql(nvars))
  freql(:) = 'day'

  !initialize timer
  call timer%Init(namelist_buffer=namelist_buffer)

  dt=timer%get_step_size()
  !call hist%init(varl, unitl, freql, 'jarmodel')

end program main
