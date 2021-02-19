
program main

  use betr_constants   , only : stdout, betr_filename_length, betr_namelist_buffer_size
  use betr_utils       , only : remove_filename_extension, namelist_to_buffer
implicit none


  integer :: arg_count
  integer :: args
  character(len=betr_filename_length)      :: namelist_filename
  character(len=betr_namelist_buffer_size) :: namelist_buffer
  character(len=betr_filename_length)      :: base_filename

  arg_count = command_argument_count()
  if (arg_count /= 1) then
     write(stdout, '(a, i3)') 'ERROR: must pass exactly one arguement one the command line, received ', arg_count
     call usage()
     call abort()
  end if

  call get_command_argument(1, namelist_filename)
  base_filename = remove_filename_extension(namelist_filename)
  write(stdout, '(a, a)') 'Using base filename for output : ', trim(base_filename)

  write(stdout, '(a, a)') 'Reading namelist filename : ', trim(namelist_filename)
  namelist_buffer = ''

  call namelist_to_buffer(namelist_filename, namelist_buffer)

  call run_model(namelist_buffer)

end program main



! ----------------------------------------------------------------------
subroutine usage()
  !DESCRIPTION
  !display something
  use betr_constants, only : stdout
 implicit none
  write(stdout, *) 'jarmodel - standalone driver for BeTR jarmodel library.'
  write(stdout, *) 'usage: jarmodel namelist_filename'
end subroutine usage

! ----------------------------------------------------------------------
subroutine run_model(namelist_buffer)
  use BeTRJarModel     , only : jar_model_type
  use JarModelFactory  , only : create_jar_model, create_jar_pars
  use JarBgcForcType   , only : JarBGC_forc_type
  use BetrStatusType   , only : betr_status_type
  use betr_ctrl        , only : iulog => biulog
  use BiogeoConType    , only : BiogeoCon_type
  use abortutils       , only : endrun
  use histMod          , only : histf_type
  use histMod          , only : hist_var_str_len, hist_unit_str_len, hist_freq_str_len
  use BeTR_TimeMod     , only : betr_time_type
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use SetJarForcMod    , only : SetJarForc, setJarStates
  use SoilForcType     , only : soil_forc_type
  use AtmForcType      , only : atm_forc_type
  use OMForcType       , only : om_forc_type
  use NutForcType      , only : nut_forc_type
  use ForcDataType     , only : load_forc, init_forc
  use betr_constants   , only : betr_string_length_long
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  implicit none
  character(len=*), intent(in) :: namelist_buffer

  !LOCAL VARIABLES
  character(len=*), parameter :: mod_filename = &
       __FILE__
  class(jar_model_type),  pointer :: jarmodel
  class(BiogeoCon_type),  pointer :: jarpars
  type(JarBGC_forc_type) :: bgc_forc
  type(soil_forc_type)   :: soil_forc
  type(om_forc_type)     :: om_forc
  type(nut_forc_type)    :: nut_forc
  type(atm_forc_type)    :: atm_forc
  type(betr_status_type) :: bstatus
  type(histf_type) :: hist
  type(betr_time_type) :: timer
  character(len=24) :: jarmodel_name
  character(len=hist_var_str_len) , allocatable :: varl(:)
  character(len=hist_unit_str_len), allocatable :: unitl(:)
  character(len=hist_freq_str_len), allocatable :: freql(:)
  integer                         , allocatable :: vartypes(:)
  real(r8), allocatable :: ystatesf(:,:)
  real(r8), allocatable :: ystates0l(:)
  real(r8), allocatable :: ystatesfl(:)

  real(r8) :: dtime
  integer :: nvars
  character(len=14) :: yymmddhhss
  logical :: is_surflit = .false.  !logical switch for litter decomposition
  logical :: nitrogen_stress=.false.
  logical :: phosphorus_stress=.false.
  integer                                :: nml_error
  character(len=betr_string_length_long) :: ioerror_msg
  character(len=64) :: case_id
  character(len=256):: gname
  character(len=hist_freq_str_len) :: hist_freq
  integer :: ncols, jj

  namelist / jar_driver / jarmodel_name, case_id, is_surflit, &
     nitrogen_stress, phosphorus_stress, hist_freq
  logical :: batch_mode


  ncols =1
  case_id=''
  hist_freq='day'
  if ( .true. )then
     ioerror_msg=''
     read(namelist_buffer, nml=jar_driver, iostat=nml_error, iomsg=ioerror_msg)
     if (nml_error /= 0) then
        call endrun(msg="ERROR reading forcing_inparm namelist "//errmsg(mod_filename, __LINE__))
     end if
  end if
  !print*,'jarmodel_name',jarmodel_name
  !create the jar type

  jarmodel => create_jar_model(jarmodel_name)
  jarpars  => create_jar_pars(jarmodel_name)

  !initialize model parameters
  call jarpars%Init(bstatus)
  jarpars%nop_limit = .not. phosphorus_stress; jarpars%non_limit=.not. nitrogen_stress
  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif
  batch_mode=.true.
  !initialize model
  call jarmodel%init(jarpars, batch_mode, bstatus)
  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif

  !initialize the parameters
  call jarmodel%UpdateParas(jarpars, bstatus)
  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif

  !set variables
  nvars=jarmodel%getvarllen()
  allocate(varl(nvars)); allocate(unitl(nvars)); allocate(freql(nvars)); allocate(vartypes(nvars))
  call jarmodel%getvarlist(nvars, varl, unitl, vartypes)
  freql(:) = hist_freq
  allocate(ystatesf(1:ncols,1:nvars));ystatesf(:,:)=0._r8
  allocate(ystates0l(1:nvars));ystates0l(:) = 0._r8
  allocate(ystatesfl(1:nvars));ystatesfl(:) = 0._r8

  call jarmodel%init_cold(nvars, ystatesf(1,:))

  !initialize timer
  call timer%Init(namelist_buffer=namelist_buffer)

  dtime=timer%get_step_size()

  if(len(trim(case_id))==0)then
    write(gname,'(A)')'jarmodel'//'.'//trim(jarmodel_name)
  else
    write(gname,'(A)')'jarmodel'//'.'//trim(case_id)//'.'//trim(jarmodel_name)
  endif
  call hist%init(ncols,varl, unitl, vartypes, freql, gname, dtime)

  call bgc_forc%Init(nvars)
  !read in forcing
  call init_forc(namelist_buffer=namelist_buffer)

  !'load constant forcing'
  call load_forc(soil_forc, atm_forc)

  !'set constant forcing'
  call SetJarForc(bgc_forc, soil_forc)

  !print*,'run the model'
  do
    !print*,'update forcing'
    call load_forc(om_forc, nut_forc, atm_forc, soil_forc, timer%tstep)

    call SetJarForc(bgc_forc, om_forc, nut_forc, atm_forc, soil_forc, dtime, jarpars)

    call setJarStates(bgc_forc, ystatesfl)

    !print*,'runbgc'
    call jarmodel%runbgc(is_surflit, dtime, bgc_forc, nvars, &
        ystates0l, ystatesfl, bstatus)
    !print*,'update_time_stamp'
    call timer%update_time_stamp()
    !print*,'hist_wrap'
    do jj = 1, nvars
    ystatesf(1,jj)=ystatesfl(jj)
    enddo
    call hist%hist_wrap(ystatesf, timer)

    if(timer%its_a_new_year())then
      write(iulog,*)'year ',timer%get_cur_year()
      !write(iulog,*)'day ',timer%get_cur_day()
    endif
    if(timer%its_time_to_exit())exit
  enddo
  call timer%get_ymdhs(yymmddhhss)
  call hist%histrst('jarmodel', 'write', yymmddhhss)
end subroutine run_model
