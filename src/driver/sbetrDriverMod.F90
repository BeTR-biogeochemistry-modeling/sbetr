module sbetrDriverMod
!
! DESCRIPTION
! module holding subroutines to do out of clm betr application
! created by Jinyun Tang
  use shr_kind_mod        , only : r8 => shr_kind_r8

  use BeTR_TimeMod, only : betr_time_type

  implicit none

  private
  save
  public :: sbetrBGC_driver

  character(len=*), parameter :: mod_filename = __FILE__

contains

  subroutine sbetrBGC_driver(base_filename, namelist_buffer)
  !
  !DESCRIPTION
  !driver subroutine for sbetrBGC
  !
  !the rtm is done using the strang splitting approach (Strang, 1968)
  !
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varpar          , only : nlevtrc_soil
  use decompMod           , only : bounds_type

  use clm_instMod, only : atm2lnd_vars
  use clm_instMod, only : canopystate_vars
  use clm_instMod, only : carbonflux_vars
  use clm_instMod, only : chemstate_vars
  use clm_instMod, only : cnstate_vars
  use clm_instMod, only : soilhydrology_vars
  use clm_instMod, only : soilstate_vars
  use clm_instMod, only : temperature_vars
  use clm_instMod, only : waterflux_vars
  use clm_instMod, only : waterstate_vars

  use ColumnType          , only : col
  use clmgridMod           , only : init_clm_vertgrid
  use clm_initializeMod    , only : initialize

  use BeTRSimulation, only : betr_simulation_type
  use BeTRSimulationFactory, only : create_betr_simulation

  use betr_constants, only : betr_namelist_buffer_size, betr_string_length_long, betr_filename_length

  use ForcingDataType, only : ForcingData_type
  use BeTR_GridMod, only : betr_grid_type

  use TracerParamsMod     , only : tracer_param_init
  use spmdMod             , only : spmd_init
  use LandunitType        , only : lun
  use PatchType           , only : pft
  use landunit_varcon     , only : istsoil

  implicit none

  character(len=betr_filename_length), intent(in) :: base_filename
  character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer


  !local variables
  class(betr_simulation_type), pointer :: simulation
  real(r8) :: dtime    !model time step
  real(r8) :: dtime2   !half of the model time step
  integer  :: record

  character(len=80) :: subname = 'sbetrBGC_driver'

  type(bounds_type) :: bounds
  integer :: lbj, ubj

  character(len=betr_string_length_long) :: simulator_name
  class(ForcingData_type), allocatable :: forcing_data
  class(betr_grid_type), allocatable :: grid_data
  class(betr_time_type), allocatable :: time_vars

  !initialize parameters
  call read_name_list(namelist_buffer, simulator_name)
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

  allocate(forcing_data)
  call forcing_data%ReadData(namelist_buffer, grid_data)

  lbj = 1
  ubj = nlevtrc_soil

  lun%itype(1) = istsoil
  col%landunit(1) = 1
  col%gridcell(1) = 1
  pft%landunit(1) = 1
  pft%column(1)   = 1
  pft%itype(1)    = 2

  call spmd_init
  !set filters to load initialization data from input
  call simulation%BeTRSetFilter()

  !obtain waterstate_vars for initilizations that need it
  call forcing_data%UpdateForcing(grid_data, &
       bounds, lbj, ubj, simulation%num_soilc, simulation%filter_soilc, time_vars, col, &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars             , &
       waterflux_vars, temperature_vars, chemstate_vars, simulation%jtops)

  !calculate advective velocity
  call calc_qadv(forcing_data, ubj, record, &
      simulation%num_soilc, simulation%filter_soilc, &
      time_vars, waterstate_vars, waterflux_vars)

  call  simulation%Init(base_filename, namelist_buffer, bounds, waterstate_vars, cnstate_vars)

  record = -1

  call time_vars%proc_initstep()
  do
    record = record + 1

    call simulation%PreDiagSoilColWaterFlux(bounds, simulation%num_soilc, &
      simulation%filter_soilc, waterstate_vars)

    !set envrionmental forcing by reading foring data: temperature, moisture, atmospheric resistance
    !from either user specified file or clm history file

    call forcing_data%UpdateForcing(grid_data, &
         bounds, lbj, ubj, simulation%num_soilc, simulation%filter_soilc, time_vars, col, &
      atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars, &
      waterflux_vars, temperature_vars, chemstate_vars, simulation%jtops)

    !calculate advective velocity
    !call calc_qadv(forcing_data, ubj, record, &
    !     simulation%num_soilc, simulation%filter_soilc, &
    !     time_vars, waterstate_vars, waterflux_vars)

    call simulation%DiagAdvWaterFlux(time_vars, bounds, &
      simulation%num_soilc, simulation%filter_soilc, &
      waterstate_vars, soilhydrology_vars, waterflux_vars)

    !no calculation in the first step
    if(record==0)cycle
    call simulation%BeginMassBalanceCheck(bounds)

    call simulation%StepWithoutDrainage(time_vars, bounds, col, &
         atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
         temperature_vars, waterflux_vars, chemstate_vars, &
         cnstate_vars, canopystate_vars, carbonflux_vars)

    call simulation%StepWithDrainage(bounds, col)

    !do mass balance check
    call simulation%MassBalanceCheck(time_vars, bounds)

    !specific for water tracer transport
    !call simulation%ConsistencyCheck(bounds, ubj, simulation%num_soilc, &
    !  simulation%filter_soilc, waterstate_vars)

    !update time stamp
    call time_vars%update_time_stamp()

    !write output
    call simulation%WriteHistory(record, lbj, ubj, time_vars, waterflux_vars%qflx_adv_col)

    !write restart file? is not functionning at the moment
    !if(its_time_to_write_restart(time_vars)) call rest_write(tracerstate_vars, tracercoeff_vars, tracerflux_vars, time_vars)
    call time_vars%proc_nextstep()
    if(time_vars%its_time_to_exit()) then
       exit
    end if

  enddo

  call simulation%WriteRegressionOutput(waterflux_vars%qflx_adv_col)

end subroutine sbetrBGC_driver

! ----------------------------------------------------------------------

  subroutine read_name_list(namelist_buffer, simulator_name_arg)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use clm_varctl   , only : iulog
    use abortutils      , only : endrun
    use shr_log_mod     , only : errMsg => shr_log_errMsg

    use betr_constants, only : stdout, betr_string_length_long, betr_namelist_buffer_size

    implicit none
    ! !ARGUMENTS:
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    character(len=betr_string_length_long), intent(out) :: simulator_name_arg
    !
    ! !LOCAL VARIABLES:
    integer :: nml_error
    character(len=*), parameter :: subname = 'read_name_list'
    character(len=betr_string_length_long) :: simulator_name
    character(len=betr_string_length_long) :: ioerror_msg


    !-----------------------------------------------------------------------

    namelist / sbetr_driver / simulator_name

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

!X!  subroutine rest_write(tracerstate_vars, tracercoeff_vars, tracerflux_vars, time_vars)
!X!  !
!X!  ! DESCRIPTION
!X!  ! write restart file
!X!
!X!  use tracerfluxType,  only : tracerflux_type
!X!  use tracerstatetype, only : tracerstate_type
!X!  use tracercoeffType, only : tracercoeff_type
!X!
!X!  implicit none
!X!  type(tracercoeff_type), intent(in) :: tracercoeff_vars
!X!  type(tracerflux_type) , intent(in) :: tracerflux_vars
!X!  type(tracerstate_type), intent(in) :: tracerstate_vars
!X!  type(betr_time_type)       , intent(in) :: time_vars
!X!  character(len=80) :: subname = 'rest_write'
!X!
!X!  end subroutine rest_write

  !-------------------------------------------------------------------------------

  subroutine calc_qadv(forcing_data, ubj, record, numf, filter, ttime, waterstate_vars, waterflux_vars)

    !
    ! description
    !calculate advective velocity between different layers
    !
    ! USES
    use WaterstateType    , only : waterstate_type
    use WaterfluxType     , only : waterflux_type
    use ColumnType        , only : column_type

    use ForcingDataType, only : ForcingData_type

    implicit none

    class(ForcingData_type), intent(inout) :: forcing_data
    integer, intent(in) :: numf
    integer, intent(in) :: filter(:)
    integer, intent(in) :: ubj

    integer                 , intent(in)    :: record
    type(betr_time_type)    , intent(in)    :: ttime
    type(waterstate_type)   , intent(inout) :: waterstate_vars
    type(waterflux_type)    , intent(inout) :: waterflux_vars

    integer :: j, fc, c
    real(r8):: dmass    !kg/m2 = mm H2O/m2

    associate( &
         tstep => ttime%tstep, &
         dtime => ttime%delta_time &
         )


      !now obtain the advective fluxes between different soil layers
      !dstorage = (h2o_new-h2o)/dt = qin-qout-qtran_dep
    if (record >= 0) then
       do fc = 1, numf
          c = filter(fc)

          do j = ubj, 1, -1
             if (j == ubj) then
                waterflux_vars%qflx_adv_col(c,j) = forcing_data%discharge(tstep) * 1.e-3_r8
             else
                dmass = (waterstate_vars%h2osoi_ice_col(c,j+1) + waterstate_vars%h2osoi_liq_col(c,j+1)) - &
                  (waterstate_vars%h2osoi_ice_old(c,j+1) + waterstate_vars%h2osoi_liq_old(c,j+1))
                waterflux_vars%qflx_adv_col(c,j) = dmass * 1.e-3_r8 / dtime + &
                     waterflux_vars%qflx_adv_col(c,j+1) + waterflux_vars%qflx_rootsoi_col(c,j+1)
             end if
          end do
          waterflux_vars%qflx_infl_col(c) = forcing_data%infiltration(tstep)

          !now correct the infiltration
          waterflux_vars%qflx_gross_infl_soil_col(c) = &
               (waterstate_vars%h2osoi_liq_col(c,1) - waterstate_vars%h2osoi_liq_old(c,1)) / dtime + &
               (waterflux_vars%qflx_rootsoi_col(c,1) + waterflux_vars%qflx_adv_col(c,1)) * 1.e3_r8

          !the following may have some problem, because it also includes contributions from
          !dew, and sublimation
          if ( waterflux_vars%qflx_gross_infl_soil_col(c) > 0._r8) then
             waterflux_vars%qflx_adv_col(c,0) = waterflux_vars%qflx_gross_infl_soil_col(c) * 1.e-3_r8
             waterflux_vars%qflx_gross_evap_soil_col(c) = 0._r8
          else
             waterflux_vars%qflx_gross_evap_soil_col(c) = waterflux_vars%qflx_gross_infl_soil_col(c)
             waterflux_vars%qflx_adv_col(c,0) = 0._r8
             waterflux_vars%qflx_gross_infl_soil_col(c) = 0._r8
          end if
       end do
    end if
    do j = 1, ubj
       do fc = 1, numf
          c = filter(fc)
          waterstate_vars%h2osoi_liq_old(c, j) = waterstate_vars%h2osoi_liq_col(c, j)
          waterstate_vars%h2osoi_ice_old(c, j) = waterstate_vars%h2osoi_ice_col(c, j)
       end do
    end do

  end associate
end subroutine calc_qadv
end module sbetrDriverMod
