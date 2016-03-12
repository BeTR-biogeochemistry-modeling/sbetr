module sbetrDriverMod
!
! DESCRIPTION
! module holding subroutines to do out of clm betr application
! created by Jinyun Tang
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use ncdio_pio
  use CLMForcType         , only : clmforc_vars

  use BeTR_TimeMod, only : betr_time_type
  
  implicit none

  private
  save
  public :: sbetrBGC_driver

  character(len=*), parameter :: mod_filename = __FILE__

  character(len=*), parameter :: histfilename="betr_output.nc"      !this will be changed
contains

  subroutine sbetrBGC_driver(namelist_buffer)
  !
  !DESCRIPTION
  !driver subroutine for sbetrBGC
  !
  !the rtm is done using the strang splitting approach (Strang, 1968)
  !
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varpar          , only : nlevtrc_soil, nlevgrnd
  use decompMod           , only : bounds_type
  use clm_instMod
  use ColumnType          , only : col
  use clmgridMod           , only : init_clm_vertgrid
  use clm_varpar           , only : nlevgrnd
  use clm_initializeMod    , only : initialize

  use BeTRSimulation, only : betr_simulation_type
  use BeTRSimulationFactory, only : create_betr_simulation

  use betr_constants, only : betr_namelist_buffer_size, betr_string_length_long

!X!  use betr_standalone_cpl , only : betr_initialize_standalone
!X!  use betr_standalone_cpl , only : run_betr_one_step_without_drainage_standalone
!X!  use betr_standalone_cpl , only : run_betr_one_step_with_drainage_standalone


  use TracerParamsMod     , only : tracer_param_init
  use spmdMod             , only : spmd_init
  use LandunitType        , only : lun
  use PatchType           , only : pft
  use landunit_varcon     , only : istsoil
  use clm_time_manager    , only : proc_initstep, proc_nextstep
  use accumulMod


  implicit none
  character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer

  !local variables
  class(betr_simulation_type), pointer :: simulation
  real(r8) :: dtime    !model time step
  real(r8) :: dtime2   !half of the model time step
  integer  :: record

  character(len=80) :: subname = 'sbetrBGC_driver'

  type(betr_time_type) :: time_vars
  type(bounds_type) :: bounds
  integer :: lbj, ubj
  integer :: num_soilp
  integer, allocatable :: filter_soilp(:)
  integer :: num_jtops 
  integer, allocatable :: jtops(:)
  integer :: num_soilc
  integer, allocatable :: filter_soilc(:)

  character(len=betr_string_length_long) :: simulator_name

  !set up mask
  bounds%begc = 1
  bounds%endc = 1
  bounds%begp = 1
  bounds%endp = 1
  bounds%begl = 1
  bounds%endl = 1
  bounds%lbj  = 1
  bounds%ubj  = nlevgrnd

  num_soilc = 1
  allocate(filter_soilc(num_soilc))
  filter_soilc(:) = 1
  num_jtops = 1
  allocate(jtops(num_jtops))
  jtops(:) = 1
  !set up grid
  call init_clm_vertgrid(nlevgrnd)

  call initialize(bounds)


  dtime = 1800._r8   !half hourly time step
  time_vars%tstep = 1
  time_vars%time  = 0._r8
  time_vars%time_end = dtime*48._r8*365._r8*2._r8
  time_vars%restart_dtime = 1800*2

  call clmforc_vars%LoadForcingData(namelist_buffer)

  lbj = 1
  ubj = nlevgrnd

  lun%itype(1) = istsoil
  col%landunit(1) = 1
  col%gridcell(1) = 1
  pft%landunit(1) = 1
  pft%column(1)   = 1
  pft%itype(1)    = 2

  num_soilp = 1
  allocate(filter_soilp(num_soilp)); filter_soilp(:) = 1

  call spmd_init

  !initialize parameters
  call read_name_list(namelist_buffer, simulator_name)
  simulation => create_betr_simulation(simulator_name)
  call  simulation%Init(namelist_buffer, bounds, waterstate_vars)
  !X!call betr_initialize_standalone(bounds, lbj, ubj)

  !create output file
  call hist_htapes_create(histfilename,nlevtrc_soil, num_soilc, simulation%betr%tracers)

  record = -1
  !in all calculations, ubj is set to nlevtrc_soil
  ubj = nlevtrc_soil
  call proc_initstep()
  do
    record = record + 1
    !set envrionmental forcing by reading foring data: temperature, moisture, atmospheric resistance
    !from either user specified file or clm history file

    call read_betrforcing(bounds, lbj, ubj, num_soilc, filter_soilc, time_vars, col, &
      atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars             , &
      waterflux_vars, temperature_vars, chemstate_vars, jtops)

    !calculate advective velocity
    call calc_qadv(ubj, record, num_soilc, filter_soilc, dtime, time_vars, col, waterstate_vars, waterflux_vars)

    !no calculation in the first step
    if(record==0)cycle
    call simulation%BeginMassBalanceCheck(bounds,  num_soilc, filter_soilc)

    call simulation%StepWithoutDrainage(bounds,  num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
         atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, temperature_vars, waterflux_vars, chemstate_vars, &
         cnstate_vars, canopystate_vars, carbonflux_vars)

    call simulation%StepWithDrainage(bounds,  num_soilc, filter_soilc, &
         jtops, waterflux_vars, col)

    !do mass balance check
    call simulation%MassBalanceCheck(bounds,  num_soilc, filter_soilc)

    !update time stamp
    call update_time_stamp(time_vars, dtime)

    !write output
    call hist_write(record, lbj, ubj, simulation%betr%tracerfluxes, &
         simulation%betr%tracerstates, time_vars, &
         simulation%betr%tracers)

    !write restart file? is not functionning at the moment
    !if(its_time_to_write_restart(time_vars)) call rest_write(tracerstate_vars, tracercoeff_vars, tracerflux_vars, time_vars)
    call proc_nextstep()
    if(its_time_to_exit(time_vars))exit

  enddo

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

  function its_time_to_write_restart(ttime)result(ans)
  !
  ! DESCRIPTION
  ! decide if to write restart file
  !

  implicit none
  type(betr_time_type), intent(in) :: ttime
  logical :: ans

  character(len=80) :: subname = 'its_time_to_write_restart'



  ans = (mod(ttime%time,ttime%restart_dtime) == 0)
  end function its_time_to_write_restart

!-------------------------------------------------------------------------------
  function its_time_to_exit(ttime)result(ans)
  !
  ! DESCRIPTION
  ! decide if to exit the loop
  !

  implicit none
  type(betr_time_type), intent(in) :: ttime
  logical :: ans


  character(len=80) :: subname = 'its_time_to_exit'

  ans= (ttime%time .eq. ttime%time_end)


  end function its_time_to_exit

!-------------------------------------------------------------------------------
  subroutine update_time_stamp(ttime, dtime)
  !
  ! DESCRIPTION
  !
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  type(betr_time_type), intent(inout) :: ttime
  real(r8),           intent(in) :: dtime

  character(len=80) :: subname='update_time_stamp'

  ttime%time=ttime%time+dtime

  ttime%tstep = ttime%tstep + 1
  if(mod(ttime%tstep, 48*365)==0)ttime%tstep = 48*365

  end subroutine update_time_stamp

!-------------------------------------------------------------------------------
  subroutine hist_write(record, lbj, ubj, tracerflux_vars, tracerstate_vars, time_vars, betrtracer_vars)
  !
  ! DESCRIPTION
  ! output hist file
  !
  use TracerFluxType,  only : tracerflux_type
  use TracerStateType, only : tracerstate_type
  use BeTRTracerType , only : BeTRTracer_Type
  implicit none
  integer                , intent(in) :: record
  integer                , intent(in) :: lbj,ubj
  type(tracerflux_type)  , intent(in) :: tracerflux_vars
  type(tracerstate_type) , intent(in) :: tracerstate_vars
  type(betr_time_type)        , intent(in) :: time_vars
  type(BeTRTracer_Type)  , intent(in) :: betrtracer_vars


  type(file_desc_t)     :: ncid
  integer :: jj
  character(len=80) :: subname='hist_write'

  associate(                                                    &
    ntracers          =>  betrtracer_vars%ntracers            , &
    ngwmobile_tracers =>  betrtracer_vars%ngwmobile_tracers   , &
    is_volatile       =>  betrtracer_vars%is_volatile         , &
    is_h2o            =>  betrtracer_vars%is_h2o              , &
    is_isotope        =>  betrtracer_vars%is_isotope          , &
    volatileid        =>  betrtracer_vars%volatileid          , &
    tracernames       =>  betrtracer_vars%tracernames           &
  )

  call ncd_pio_openfile_for_write(ncid,histfilename)

  if (mod(time_vars%time ,86400._r8)==0) then
     print*,'day', time_vars%time/86400._r8
  end if
  call ncd_putvar(ncid, "time", record, time_vars%time)

  do jj = 1, ntracers
    if(jj<= ngwmobile_tracers)then
      call ncd_putvar(ncid,trim(tracernames(jj))//'_TRACER_CONC_MOIBLE',&
        record,tracerstate_vars%tracer_conc_mobile_col(1:1,lbj:ubj,jj))

        if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then

          call ncd_putvar(ncid, trim(tracernames(jj))//'_TRACER_P_GAS_FRAC',&
            record,tracerstate_vars%tracer_P_gas_frac_col(1:1,lbj:ubj,volatileid(jj)))
        endif
        if(is_volatile(jj))then
          call ncd_putvar(ncid, trim(tracernames(jj))//'_FLX_SURFEMI', &
            record, tracerflux_vars%tracer_flx_surfemi_col(1:1, volatileid(jj)))
        endif
    else
      call ncd_putvar(ncid, trim(tracernames(jj))//'_TRACER_CONC_SOLID_PASSIVE', &
      record, tracerstate_vars%tracer_conc_solid_passive_col(1:1,lbj:ubj,jj-ngwmobile_tracers))
    endif
  enddo

  call ncd_putvar(ncid, 'TRACER_P_GAS', record, tracerstate_vars%tracer_P_gas_col(1:1,lbj:ubj))
  call ncd_pio_closefile(ncid)
  end associate
  end subroutine hist_write
!-------------------------------------------------------------------------------
  subroutine rest_write(tracerstate_vars, tracercoeff_vars, tracerflux_vars, time_vars)
  !
  ! DESCRIPTION
  ! write restart file

  use tracerfluxType,  only : tracerflux_type
  use tracerstatetype, only : tracerstate_type
  use tracercoeffType, only : tracercoeff_type

  implicit none
  type(tracercoeff_type), intent(in) :: tracercoeff_vars
  type(tracerflux_type) , intent(in) :: tracerflux_vars
  type(tracerstate_type), intent(in) :: tracerstate_vars
  type(betr_time_type)       , intent(in) :: time_vars
  character(len=80) :: subname = 'rest_write'





  end subroutine rest_write


!-------------------------------------------------------------------------------

  subroutine hist_htapes_create(histfilename,nlevtrc_soil, ncol, betrtracer_vars)
  !
  ! DESCRIPTIONS
  ! create history file and define output variables
  use netcdf
  use histFileMod, only : hist_file_create, hist_def_fld1d, hist_def_fld2d

  use BeTRTracerType, only: BeTRTracer_Type
  !
  !ARGUMENTS
  implicit none
  character(len=*)     , intent(in) :: histfilename
  integer              , intent(in) :: nlevtrc_soil, ncol
  type(BeTRTracer_Type), intent(in) :: betrtracer_vars

  integer            :: jj, kk
  type(file_desc_t)  :: ncid
  character(len=255) :: subname = 'hist_htapes_create'

  associate(                                                    &
    ntracers          =>  betrtracer_vars%ntracers            , &
    ngwmobile_tracers =>  betrtracer_vars%ngwmobile_tracers   , &
    is_volatile       =>  betrtracer_vars%is_volatile         , &
    is_h2o            =>  betrtracer_vars%is_h2o              , &
    is_isotope        =>  betrtracer_vars%is_isotope          , &
    volatileid        =>  betrtracer_vars%volatileid          , &
    tracernames       =>  betrtracer_vars%tracernames           &
  )

  call ncd_pio_createfile(ncid, histfilename)

  call hist_file_create(ncid,nlevtrc_soil, ncol)


  call  hist_def_fld2d(ncid, varname="TRACER_P_GAS", nf90_type=nf90_float,dim1name="ncol", &
      dim2name="levgrnd",long_name="total gas pressure", units="Pa")

  do jj = 1, ntracers
    if(jj<= ngwmobile_tracers)then
      call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_CONC_MOIBLE',&
        nf90_type=nf90_float, dim1name="ncol", &
        dim2name="levgrnd", long_name=trim(tracernames(jj))//"tracer concentrations", &
        units="mol m-3")

        if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then

          call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_P_GAS_FRAC',&
            nf90_type=nf90_float, dim1name="ncol", &
            dim2name="levgrnd", long_name='fraction of gas phase contributed by '//trim(tracernames(jj)), &
            units="none")

        endif

        if(is_volatile(jj))then

          call hist_def_fld1d (ncid, varname=trim(tracernames(jj))//'_FLX_SURFEMI', units='mol/m2/s', &
            nf90_type=nf90_float,  dim1name="ncol", &
            long_name='loss from surface emission for '//trim(tracernames(jj)))
        endif
    else
      call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_CONC_SOLID_PASSIVE',&
        nf90_type=nf90_float, dim1name="ncol", &
        dim2name="levgrnd", long_name=trim(tracernames(jj))//"tracer concentrations", &
        units="mol m-3")
    endif
  enddo


  call ncd_enddef(ncid)

  call ncd_pio_closefile(ncid)

  end associate
  end subroutine hist_htapes_create

!-------------------------------------------------------------------------------

  subroutine read_betrforcing(bounds, lbj, ubj, numf, filter, ttime, col, atm2lnd_vars, &
    soilhydrology_vars, soilstate_vars,waterstate_vars,waterflux_vars, &
    temperature_vars,chemstate_vars, jtops)
  !
  ! DESCRIPTIONS
  ! read environmental forcing to run betr
  ! for clm forced runs, it will forcing from history files
  !
  use TemperatureType   , only : temperature_type
  use WaterstateType    , only : waterstate_type
  use WaterfluxType     , only : waterflux_type
  use SoilStateType     , only : soilstate_type
  use ChemStateType     , only : chemstate_type
  use ColumnType        , only : column_type
  use clmgridMod        , only : dzsoi, zisoi
  use decompMod         , only : bounds_type
  use SoilHydrologyType , only : soilhydrology_type
  use atm2lndType       , only : atm2lnd_type
  implicit none
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)
  integer, intent(in) :: lbj, ubj
  type(betr_time_type),             intent(in) :: ttime
  type(chemstate_type),     intent(inout) :: chemstate_vars
  type(atm2lnd_type),       intent(inout) :: atm2lnd_vars
  type(soilstate_type),     intent(inout) :: soilstate_vars
  type(waterstate_type),    intent(inout) :: waterstate_vars
  type(waterflux_type),     intent(inout) :: waterflux_vars
  type(temperature_type),   intent(inout) :: temperature_vars
  type(column_type),        intent(inout) :: col
  type(soilhydrology_type), intent(inout) :: soilhydrology_vars
  integer,               intent(inout)    :: jtops(bounds%begc:bounds%endc)

  integer :: j, fc, c
  character(len=255) :: subname='read_betrforcing'



  associate(                   &
    tstep => ttime%tstep       &
  )
  !setup top boundary
  do fc = 1, numf
    c = filter(fc)
    jtops(c) = 1
    soilhydrology_vars%zwts_col(c) = 10._r8
    atm2lnd_vars%forc_pbot_downscaled_col(c) = clmforc_vars%pbot(tstep)             ! 1 atmos
  enddo

  !set up forcing variables
  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(j>=jtops(c))then
        waterstate_vars%h2osoi_liqvol_col(c,j) = clmforc_vars%h2osoi_liqvol(tstep,j)
        waterstate_vars%air_vol_col(c,j)       = clmforc_vars%watsat(j)-clmforc_vars%h2osoi_liqvol(tstep,j)
        soilstate_vars%eff_porosity_col(c,j)   = clmforc_vars%watsat(j)-clmforc_vars%h2osoi_icevol(tstep,j)
        soilstate_vars%bsw_col(c,j)            = clmforc_vars%bsw(j)
        temperature_vars%t_soisno_col(c,j)     = clmforc_vars%t_soi(tstep,j)
        waterflux_vars%qflx_rootsoi_col(c,j)   = clmforc_vars%qflx_tran_dep(tstep,j)  !water exchange between soil and root, mm/H2O/s
        col%dz(c,j) = dzsoi(j)
        col%zi(c,j) = zisoi(j)
        chemstate_vars%soil_pH(c,j)=7._r8
        !set drainage to zero
        !set surface runoff to zero
        waterflux_vars%qflx_surf_col(c) = 0._r8
        waterflux_vars%qflx_drain_vr_col(c,j) = 0._r8

      endif
    enddo
  enddo

  do fc = 1, numf
      c = filter(fc)
      waterflux_vars%qflx_totdrain_col(c) = 0._r8
      waterflux_vars%qflx_infl_col(c)  = clmforc_vars%qflx_infl(tstep)              !infiltration flux, mm H2O/s
      col%zi(c,0) = zisoi(0)
      waterflux_vars%qflx_gross_infl_soil_col(c) = 0._r8
  enddo


  do j = 1, ubj
    do fc = 1, numf
      c = filter(fc)
      waterstate_vars%h2osoi_liq_col(c,j)    = clmforc_vars%h2osoi_liq(tstep,j)
      waterstate_vars%h2osoi_ice_col(c,j)    = clmforc_vars%h2osoi_ice(tstep,j)
    enddo
  enddo
  end associate
  end subroutine read_betrforcing

  !-------------------------------------------------------------------------------

  subroutine calc_qadv(ubj, record, numf, filter, dtime, ttime, col, waterstate_vars, waterflux_vars)

  !
  ! description
  !calculate advective velocity between different layers
  !
  ! USES
  use WaterstateType    , only : waterstate_type
  use WaterfluxType     , only : waterflux_type
  use ColumnType        , only : column_type

  ! ARGUMENTS
  implicit none
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)
  integer, intent(in) :: ubj

  integer                 , intent(in)    :: record
  type(betr_time_type)    , intent(in)    :: ttime
  real(r8)                , intent(in)    :: dtime
  type(waterstate_type)   , intent(inout) :: waterstate_vars
  type(waterflux_type)    , intent(inout) :: waterflux_vars
  type(column_type)       , intent(inout) :: col

  integer :: j, fc, c
  real(r8):: dmass    !kg/m2 = mm H2O/m2

  associate(                   &
    tstep => ttime%tstep       &
  )


  !now obtain the advective fluxes between different soil layers
  !dstorage = (h2o_new-h2o)/dt = qin-qout-qtran_dep
  if(record > 0)then
    do fc = 1, numf
      c = filter(fc)
      waterflux_vars%qflx_adv_col(c,0) = clmforc_vars%qflx_infl(tstep)
      if(waterflux_vars%qflx_adv_col(c,0)< 0._r8)then
        waterflux_vars%qflx_gross_evap_soil_col(c) = -waterflux_vars%qflx_adv_col(c,0)
        waterflux_vars%qflx_adv_col(c,0) = 0._r8
      endif
      do j = 1,ubj

        dmass=(waterstate_vars%h2osoi_ice_col(c,j)+waterstate_vars%h2osoi_liq_col(c,j))- &
           (waterstate_vars%h2osoi_ice_old(c,j)+waterstate_vars%h2osoi_liq_old(c,j))

        waterflux_vars%qflx_adv_col(c,j)= waterflux_vars%qflx_adv_col(c,j-1) - waterflux_vars%qflx_rootsoi_col(c,j) &
          - dmass/dtime
      enddo
    enddo

    !now convert all flux unit into m/s

    do fc = 1, numf
      c = filter(fc)
      waterflux_vars%qflx_adv_col(c,0) = waterflux_vars%qflx_adv_col(c,0)*1.e-3_r8
      do j = 1,ubj
        waterflux_vars%qflx_adv_col(c,j)= waterflux_vars%qflx_adv_col(c,j)*1.e-3_r8
        waterflux_vars%qflx_rootsoi_col(c,j) = waterflux_vars%qflx_rootsoi_col(c,j) * 1.e-3_r8
      enddo
    enddo
  endif
  do j = 1, ubj
    do fc = 1, numf
      c = filter(fc)
      waterstate_vars%h2osoi_liq_old(c,j)    = waterstate_vars%h2osoi_liq_col(c,j)
      waterstate_vars%h2osoi_ice_old(c,j)    = waterstate_vars%h2osoi_ice_col(c,j)
    enddo
  enddo

  end associate
  end subroutine calc_qadv
end module sbetrDriverMod
