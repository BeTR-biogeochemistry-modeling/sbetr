module sbetrDriverMod
!
! DESCRIPTION
! module holding subroutines to do out of clm betr application
! created by Jinyun Tang
  use shr_kind_mod, only : r8 => shr_kind_r8
implicit none

  private
  save
  public :: sbetrBGC_driver


  type,public:: time_type
    real(r8) :: time_end
    real(r8) :: time
    real(r8) :: restart_dtime
  end type time_type

  character(len=255), parameter :: histfilename="betr_output.nc"      !this will be changed
contains

  subroutine sbetrBGC_driver(bounds, num_soilc, filter_soilc,time_vars)
  !
  !DESCRIPTION
  !driver subroutine for sbetrBGC
  !
  !the rtm is done using the strang splitting approach (Strang, 1968)
  !
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varpar          , only : nlevtrc_soil
  use decompMod           , only : bounds_type
  use clm_initializeMod   , only : soilstate_vars, waterstate_vars, temperature_vars, &
                                   waterflux_vars, chemstate_vars, atm2lnd_vars, soilhydrology_vars, &
                                   canopystate_vars, carbonflux_vars, carbonstate_vars, cnstate_vars, &
                                   nitrogenflux_vars, nitrogenstate_vars
  use ColumnType          , only : col
  use betr_initializeMod  , only : betr_initialize, betrtracer_vars, tracercoeff_vars,  bgc_reaction, &
                                    tracerflux_vars, tracerState_vars, tracerboundarycond_vars, plantsoilnutrientflux_vars
  use BetrBGCMod          , only : run_betr_one_step_without_drainage, betrbgc_init
  use TracerParamsMod     , only : tracer_param_init
  use ncdio_pio           , only : file_desc_t, ncd_pio_closefile
  use spmdMod             , only : spmd_init
  use accumulMod
  use TracerBalanceMod
  implicit none
  type(bounds_type),        intent(in) :: bounds                                ! bounds
  integer,                  intent(in) :: num_soilc                                  ! number of columns in column filter
  integer,                  intent(in) :: filter_soilc(:)                             ! column filter

  type(time_type),       intent(inout) :: time_vars

  !local variables

  real(r8) :: dtime    !model time step
  real(r8) :: dtime2   !half of the model time step
  integer  :: record

  character(len=80) :: subname = 'sbetrBGC_driver'

  integer :: num_soilp
  integer, allocatable :: filter_soilp(:)
  integer :: lbj, ubj
  integer :: jtops(bounds%begc:bounds%endc)


  dtime = 1800._r8


  jtops(:) = 999  !this will be replaced with nan when I figured out how to do it, Jinyun Tang, June 17, 2014

  lbj = 1
  ubj = nlevtrc_soil

  call spmd_init
  !initialize parameters
  call betr_initialize(bounds, lbj, ubj, waterstate_vars)

  !initialize the betr parameterization module
  call tracer_param_init(bounds)

  !initialize the betrBGC module
  call betrbgc_init(bounds, betrtracer_vars)

  !set up model time, in CLM, this will be clm_inparm, but one has to
  !set it different for the offline betr code

  !

  !create output file

  call hist_htapes_create(histfilename,nlevtrc_soil, num_soilc, betrtracer_vars)

  record = 0
  do
    !set envrionmental forcing by reading foring data: temperature, moisture, atmospheric resistance
    !from either user specified file or clm history file
    !and advective velocity at interfatemperature_vars
    call read_betrforcing(bounds, lbj, ubj, num_soilc, filter_soilc, time_vars, col, &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars,  &
      waterflux_vars, temperature_vars, chemstate_vars, jtops)

    call run_betr_one_step_without_drainage(bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
         atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, temperature_vars, waterflux_vars, chemstate_vars, &
         cnstate_vars, canopystate_vars, carbonflux_vars, betrtracer_vars, bgc_reaction, tracerboundarycond_vars, &
         tracercoeff_vars, tracerstate_vars, tracerflux_vars, plantsoilnutrientflux_vars)

    !do mass balance check
    !update time stamp
    call update_time_stamp(time_vars, dtime)
    record = record + 1
    !write output
    call hist_write(record, tracerflux_vars, tracerstate_vars, time_vars)

    !write restart file?
    if(its_time_to_write_restart(time_vars)) call rest_write(tracerstate_vars, tracercoeff_vars, tracerflux_vars, time_vars)

    if(its_time_to_exit(time_vars))exit

  enddo

  end subroutine sbetrBGC_driver
!-------------------------------------------------------------------------------

  function its_time_to_write_restart(ttime)result(ans)
  !
  ! DESCRIPTION
  ! decide if to write restart file
  !

  implicit none
  type(time_type), intent(in) :: ttime
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
  type(time_type), intent(in) :: ttime
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
  type(time_type), intent(inout) :: ttime
  real(r8),           intent(in) :: dtime

  character(len=80) :: subname='update_time_stamp'


  ttime%time=ttime%time+dtime

  print*,'time',ttime%time
  end subroutine update_time_stamp

!-------------------------------------------------------------------------------
  subroutine hist_write(record, tracerflux_vars, tracerstate_vars, time_vars)
  !
  ! DESCRIPTION
  ! output hist file
  !
  use TracerFluxType,  only : tracerflux_type
  use TracerStateType, only : tracerstate_type
  use ncdio_pio      , only : ncd_pio_openfile_for_write, file_desc_t
  use ncdio_pio      , only : ncd_pio_closefile, ncd_putvar
  implicit none
  integer                , intent(in) :: record
  type(tracerflux_type)  , intent(in) :: tracerflux_vars
  type(tracerstate_type) , intent(in) :: tracerstate_vars
  type(time_type)        , intent(in) :: time_vars


  type(file_desc_t)     :: ncid
  character(len=80) :: subname='hist_write'

  call ncd_pio_openfile_for_write(ncid,histfilename)

  call ncd_putvar(ncid, "time", record, time_vars%time)

  call ncd_pio_closefile(ncid)

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
  type(tracerflux_type),  intent(in) :: tracerflux_vars
  type(tracerstate_type), intent(in) :: tracerstate_vars
  type(time_type)         , intent(in) :: time_vars
  character(len=80) :: subname = 'rest_write'





  end subroutine rest_write


!-------------------------------------------------------------------------------

  subroutine hist_htapes_create(histfilename,nlevtrc_soil, ncol, betrtracer_vars)
  !
  ! DESCRIPTIONS
  ! create history file and define output variables
  use netcdf
  use histFileMod, only : hist_file_create, hist_def_fld1d, hist_def_fld2d
  use ncdio_pio,   only : file_desc_t, ncd_pio_createfile,  &
                          ncd_pio_closefile, ncd_enddef, ncd_putvar
  use BeTRTracerType, only: BeTRTracer_Type
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
    volatileid        =>  betrtracer_vars%volatileid          , &
    tracernames       =>  betrtracer_vars%tracernames           &
  )

  call ncd_pio_createfile(ncid, histfilename)

  call hist_file_create(ncid,nlevtrc_soil, ncol)


  call  hist_def_fld2d(ncid, varname="TRACER_P_GAS", nf90_type=nf90_float,dim1name="levgrnd", &
      dim2name="ncol",long_name="total gas pressure", units="Pa")
  print*,'ntracers',ntracers
  do jj = 1, ntracers
    if(jj<= ngwmobile_tracers)then
      call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_CONC_MOIBLE',&
        nf90_type=nf90_float, dim1name="levgrnd", &
        dim2name="ncol", long_name=trim(tracernames(jj))//"tracer concentrations", &
        units="mol m-3")
    else
      call hist_def_fld2d(ncid, varname=trim(tracernames(jj))//'_TRACER_CONC_SOLID_PASSIVE',&
        nf90_type=nf90_float, dim1name="levgrnd", &
        dim2name="ncol", long_name=trim(tracernames(jj))//"tracer concentrations", &
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
  type(time_type),             intent(in) :: ttime
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



  !setup top boundary
  do fc = 1, numf
    c = filter(fc)
    jtops(c) = 1
    soilhydrology_vars%zwts_col(c) = 10._r8
    atm2lnd_vars%forc_pbot_downscaled_col(c) = 1.01325e5_r8                     ! 1 atmos
  enddo
  !set up forcing variables
  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(j>=jtops(c))then
        waterstate_vars%h2osoi_liqvol_col(c,j) = 0.3_r8
        waterstate_vars%air_vol_col(c,j) = 0.2_r8
        soilstate_vars%eff_porosity_col(c,j)  = 0.5_r8
        soilstate_vars%bsw_col(c,j) = 3._r8
        col%dz(c,j) = dzsoi(j)
        col%zi(c,j) = zisoi(j)
        temperature_vars%t_soisno_col(c,j) = 298._r8
        chemstate_vars%soil_pH(c,j)=7._r8
        soilhydrology_vars%fracice_col(c,j) = 0._r8                             !no ice
        waterflux_vars%qflx_adv_col(c,j) = 1.e-6_r8                                 !advective tracer flux, m H2O/s
        waterflux_vars%qflx_rootsoi_col(c,j) = 0._r8                                !water exchange between soil and root
      endif
    enddo
  enddo
  do fc = 1, numf
      waterflux_vars%qflx_infl_col(c)  = 1.e-3_r8                                 !infiltration flux, mm H2O/s
      col%zi(c,0) = zisoi(0)
  enddo

  end subroutine read_betrforcing

!-------------------------------------------------------------------------------
  subroutine read_clmforcing(ttime, soilstate_vars, waterstate_vars, temperature_vars)
  !
  ! DESCRIPTION
  !
  ! reading offline soil biophysical forcing from CLM output
  ! by default, all data are at half ourly time step
  !
  use TemperatureType, only : temperature_type
  use WaterStateType , only : waterstate_type
  use SoilStateType  , only : soilstate_type

  implicit none
  type(soilstate_type), intent(inout) :: soilstate_vars
  type(waterstate_type), intent(inout) :: waterstate_vars
  type(temperature_type), intent(inout) :: temperature_vars
  type(time_type),             intent(in) :: ttime

  character(len=255) :: subname ='read_clmforcing'


  end subroutine read_clmforcing
end module sbetrDriverMod
