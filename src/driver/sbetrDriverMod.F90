module sbetrDriverMod
!
! DESCRIPTION
! module holding subroutines to do out of clm betr application
! created by Jinyun Tang
  use shr_kind_mod , only : r8 => shr_kind_r8
  use BeTR_TimeMod , only : betr_time_type

  implicit none

  private
  save
  public :: sbetrBGC_driver

  character(len=*), parameter :: mod_filename = &
       __FILE__

contains

  subroutine sbetrBGC_driver(base_filename, namelist_buffer)
  !
  !DESCRIPTION
  !driver subroutine for sbetrBGC
  !
  !the rtm is done using the strang splitting approach (Strang, 1968)
  !
  use BeTRSimulationCLM     , only : betr_simulation_clm_type
  use BeTRSimulationStandalone, only : betr_simulation_standalone_type
  use BeTRSimulationALM     , only : betr_simulation_alm_type
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use clm_varpar            , only : nlevtrc_soil
  use decompMod             , only : bounds_type
  use ncdio_pio             , only : file_desc_t
  use clm_instMod           , only : atm2lnd_vars
  use clm_instMod           , only : canopystate_vars
  use clm_instMod           , only : carbonflux_vars
  use clm_instMod           , only : chemstate_vars
  use clm_instMod           , only : soilhydrology_vars
  use clm_instMod           , only : soilstate_vars
  use clm_instMod           , only : temperature_vars
  use clm_instMod           , only : waterflux_vars
  use clm_instMod           , only : waterstate_vars
  use ColumnType            , only : col
  use clmgridMod            , only : init_clm_vertgrid
  use clm_initializeMod     , only : initialize
  use BeTRSimulation        , only : betr_simulation_type
  use BeTRSimulationFactory , only : create_betr_simulation
  use betr_constants        , only : betr_namelist_buffer_size, betr_string_length_long, betr_filename_length
  use ForcingDataType       , only : ForcingData_type
  use BeTR_GridMod          , only : betr_grid_type
  use TracerParamsMod       , only : tracer_param_init
  use spmdMod               , only : spmd_init
  use LandunitType          , only : lun
  use PatchType             , only : pft
  use landunit_varcon       , only : istsoil
  use clm_varpar            , only : nlevsno, nlevsoi
  implicit none
  !arguments
  character(len=betr_filename_length)      , intent(in) :: base_filename
  character(len=betr_namelist_buffer_size) , intent(in) :: namelist_buffer


  !local variables
  class(betr_simulation_type), pointer :: simulation
  real(r8)                             :: dtime    !model time step
  real(r8)                             :: dtime2   !half of the model time step
  integer                              :: record

  character(len=80)                    :: subname = 'sbetrBGC_driver'

  type(bounds_type)                    :: bounds
  integer                              :: lbj, ubj
  logical :: continue_run
  type(file_desc_t)                      :: ncid
  character(len=betr_string_length_long) :: simulator_name
  character(len=betr_string_length_long) :: restfname
  class(ForcingData_type), allocatable :: forcing_data
  class(betr_grid_type), allocatable :: grid_data
  class(betr_time_type), allocatable :: time_vars
  character(len=256) :: restfile
  integer :: nstep
  !initialize parameters
  call read_name_list(namelist_buffer, simulator_name, continue_run)
  simulation => create_betr_simulation(simulator_name)

  !set up mask
  bounds%begc = 1
  bounds%endc = 1
  bounds%begp = 1
  bounds%endp = 1
  bounds%begl = 1
  bounds%endl = 1
  bounds%lbj  = 1
  bounds%ubj  = nlevtrc_soil

  !set up grid
  allocate(grid_data)
  call grid_data%Init(namelist_buffer)
  call init_clm_vertgrid(grid_data%nlevgrnd)

  call initialize(bounds)

  allocate(time_vars)
  call time_vars%Init(namelist_buffer)

  !x print*,'read forcing'
  allocate(forcing_data)
  call forcing_data%ReadData(namelist_buffer, grid_data)

  lbj = 1
  ubj = nlevtrc_soil

  lun%itype(1) = istsoil
  col%landunit(1) = 1
  col%gridcell(1) = 1
  col%npfts(1)    = 1
  col%pfti(1)     = 1
  pft%landunit(1) = 1
  pft%column(1)   = 1
  pft%itype(1)    = 2

  call spmd_init
  !print*,'set filters to load initialization data from input'
  call simulation%BeTRSetFilter()

  if(continue_run)then
    ! print*,'continue from restart file'
    call read_restinfo(restfname, nstep)
    !set current step
    call time_vars%set_nstep(nstep-1)
    ! print*,'back 1 step'
    call time_vars%print_cur_time()
  endif

  !print*,'obtain waterstate_vars for initilizations that need it'
  call forcing_data%UpdateForcing(grid_data,                                            &
       bounds, lbj, ubj, simulation%num_soilc, simulation%filter_soilc, time_vars, col, &
       pft, atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars    ,   &
       waterflux_vars, temperature_vars, chemstate_vars, simulation%jtops)

  !print*,'initial water state variable output',time_vars%tstep
  call calc_qadv(ubj, simulation%num_soilc, &
       simulation%filter_soilc, waterstate_vars)

  !print*,'base_filename:',trim(base_filename)
  call  simulation%Init(base_filename, namelist_buffer, bounds, lun, col, pft, waterstate_vars)

  select type(simulation)
  class is (betr_simulation_standalone_type)
    print*,'simulation using standalone-betr'
  class is (betr_simulation_alm_type)
    print*,'simulation using alm-betr'
  class is (betr_simulation_clm_type)
    print*,'simulation using clm-betr'
  class default
  end select

  !read initial condition from restart file is needed
  if(continue_run)then
    call simulation%BeTRRestartOpen(restfname, flag='read', ncid=ncid)
    call simulation%BeTRRestart(bounds, ncid, simulation%num_soilc, simulation%filter_soilc, flag='read')
    call simulation%BeTRRestartClose(ncid)
    !the following aligns forcing data with correct time stamp
    call time_vars%set_nstep(nstep)
    call time_vars%print_cur_time()
    record = 0
  else
    record = -1
    call time_vars%proc_initstep()
  endif


  do
    record = record + 1

    !print*,'prepare for diagnosing water flux'
    call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, waterstate_vars=waterstate_vars)

    call simulation%PreDiagSoilColWaterFlux(simulation%num_soilc,  simulation%filter_soilc)

    ! print*,'update forcing for betr'
    !set envrionmental forcing by reading foring data: temperature, moisture, atmospheric resistance
    !from either user specified file or clm history file

    call forcing_data%UpdateForcing(grid_data,                                            &
      bounds, lbj, ubj, simulation%num_soilc, simulation%filter_soilc, time_vars, col, pft, &
      atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars,                   &
      waterflux_vars, temperature_vars, chemstate_vars, simulation%jtops)

    call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, waterstate_vars=waterstate_vars, &
      waterflux_vars=waterflux_vars, soilhydrology_vars = soilhydrology_vars)

    !print*,'diagnose water flux'
    call simulation%DiagAdvWaterFlux(time_vars, simulation%num_soilc, &
      simulation%filter_soilc)

    !now assign back waterflux_vars
    call simulation%RetrieveBiogeoFlux(bounds, 1, nlevsoi, waterflux_vars=waterflux_vars)

    !no calculation in the first step
    if(record==0)cycle
    call simulation%BeginMassBalanceCheck(bounds)

    !print*,'without drainage'
    !the following call could be lsm specific, so that
    !different lsm could use different definitions of input
    !variables, e.g. clm doesn't use cnstate_vars as public variables
    call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi,  &
      carbonflux_vars=carbonflux_vars,       &
      waterstate_vars=waterstate_vars,         waterflux_vars=waterflux_vars,         &
      temperature_vars=temperature_vars,       soilhydrology_vars=soilhydrology_vars, &
      atm2lnd_vars=atm2lnd_vars,               canopystate_vars=canopystate_vars,     &
      chemstate_vars=chemstate_vars,           soilstate_vars=soilstate_vars)
    call simulation%StepWithoutDrainage(time_vars, bounds, col, pft)

    !print*,'with drainge'
    !set forcing variable for drainage
    call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi,&
       waterflux_vars=waterflux_vars )
    call simulation%StepWithDrainage(bounds, col)

    !print*,'do mass balance check'
    call simulation%MassBalanceCheck(time_vars, bounds)

    !specific for water tracer transport
    !call simulation%ConsistencyCheck(bounds, ubj, simulation%num_soilc,    &
    !  simulation%filter_soilc, waterstate_vars)

    !update time stamp
    call time_vars%update_time_stamp()

    !print*,'write output'
    call simulation%WriteOfflineHistory(bounds, record, simulation%num_soilc,  &
       simulation%filter_soilc, time_vars, waterflux_vars%qflx_adv_col)

    call time_vars%proc_nextstep()
    !print*,'write restart file? is not functionning at the moment'
    if(time_vars%its_time_to_write_restart()) then
       !set restfname
       nstep = time_vars%get_nstep()
       write(restfname,'(A,I8.8,A)')trim(base_filename)//'.',nstep,'.rst.nc'
       call write_restinfo(restfname, nstep)
       !print*, 'open restart file'
       call simulation%BeTRRestartOpen(restfname, flag='write', ncid=ncid)
       ! print*,'define restart file'
       call simulation%BeTRRestart(bounds, ncid, simulation%num_soilc, simulation%filter_soilc, flag='define')
       ! print*,'write restart file'
       call simulation%BeTRRestart(bounds, ncid, simulation%num_soilc, simulation%filter_soilc, flag='write')
       ! print*,'close file'
       call simulation%BeTRRestartClose(ncid)
    endif

    if(time_vars%its_time_to_exit()) then
       exit
    end if

  enddo

  call simulation%WriteRegressionOutput(waterflux_vars%qflx_adv_col)

end subroutine sbetrBGC_driver

! ----------------------------------------------------------------------

  subroutine read_name_list(namelist_buffer, simulator_name_arg, continue_run)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod        , only : masterproc, mpicom
    use clm_varctl     , only : iulog
    use abortutils     , only : endrun
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use betr_constants , only : stdout, betr_string_length_long, betr_namelist_buffer_size

    implicit none
    ! !ARGUMENTS:
    character(len=betr_namelist_buffer_size) , intent(in)  :: namelist_buffer
    character(len=betr_string_length_long)   , intent(out) :: simulator_name_arg
    logical, intent(out) :: continue_run
    !
    ! !LOCAL VARIABLES:
    integer                                :: nml_error
    character(len=*), parameter            :: subname = 'read_name_list'
    character(len=betr_string_length_long) :: simulator_name
    character(len=betr_string_length_long) :: ioerror_msg

    !-----------------------------------------------------------------------

    namelist / sbetr_driver / simulator_name, continue_run

    continue_run=.false.
    simulator_name = ''

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=sbetr_driver, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading sbetr_driver namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' sbetr driver :'
       write(stdout, *)
       write(stdout, *) ' sbetr_driver namelist settings :'
       write(stdout, *)
       write(stdout, sbetr_driver)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    simulator_name_arg = simulator_name

  end subroutine read_name_list



  !-------------------------------------------------------------------------------

  subroutine calc_qadv(ubj, numf, filter, waterstate_vars)

    !
    ! description
    ! calculate advective velocity between different layers
    ! this function is now not used, and will be deleted
    ! USES
    use WaterstateType  , only : waterstate_type
    use ColumnType      , only : column_type
    implicit none
    !ARGUMENTS
    integer, intent(in) :: numf
    integer, intent(in) :: filter(:)
    integer, intent(in) :: ubj
    type(waterstate_type)   , intent(inout) :: waterstate_vars

    integer :: j, fc, c

    do j = 1, ubj
       do fc = 1, numf
          c = filter(fc)
          waterstate_vars%h2osoi_liq_old(c, j) = waterstate_vars%h2osoi_liq_col(c, j)
          waterstate_vars%h2osoi_ice_old(c, j) = waterstate_vars%h2osoi_ice_col(c, j)
       end do
    end do

  end subroutine calc_qadv
  !-------------------------------------------------------------------------------
  subroutine read_restinfo(restfile, nstep)
  implicit none
  character(len=*), intent(out) :: restfile
  integer, intent(out) :: nstep

  integer :: rpt_unit, rpt_error

  rpt_unit=20
  open(unit=rpt_unit,file='rpoint.betr', status='old',&
          action='read', form='formatted', iostat=rpt_error)
  read(rpt_unit,*)restfile,nstep
  print*,'restart file: ', trim(restfile), 'step=',nstep
  close(rpt_unit)

  end subroutine read_restinfo

  !-------------------------------------------------------------------------------
  subroutine write_restinfo(fname, nstep)

  implicit none
  character(len=*), intent(in) :: fname
  integer, intent(in) :: nstep

  integer :: rpt_unit, rpt_error
  character(len=50) :: rpt_filename='rpoint.betr'

  !write rpoint.betr
  rpt_unit=20
  open(unit=rpt_unit,file=trim(rpt_filename), status='replace',&
     action='write', form='formatted', iostat=rpt_error)
  write(rpt_unit,*)trim(fname), nstep
  close(rpt_unit)
  end subroutine write_restinfo
end module sbetrDriverMod
