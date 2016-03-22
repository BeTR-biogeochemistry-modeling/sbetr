module BeTRSimulation
  !
  ! !DESCRIPTION:
  !  BeTR simulation base class.
  !
  !  BeTR simulation class are API definitions, mapping data
  !  structures from a specific LSM, e.g. CLM, ALM, into BeTR data
  !  structures.
  !
  use abortutils                  , only : endrun
  use clm_varctl                  , only : iulog, use_cn
  use shr_log_mod                 , only : errMsg => shr_log_errMsg
  use tracer_varcon               , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use BeTR_decompMod              , only : betr_bounds_type
  use decompMod                   , only : bounds_type

  ! !USES:
  use BetrType                  , only : betr_type
  use betr_constants            , only : betr_string_length
  use betr_constants            , only : betr_filename_length
  use betr_regression_module, only : betr_regression_type
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = __FILE__

  type, public :: betr_simulation_type
     type(betr_type), public :: betr
     character(len=betr_string_length), private :: base_filename
     character(len=betr_filename_length), private :: hist_filename

     type(betr_regression_type), private :: regression

     integer, public :: num_soilp
     integer, public, allocatable :: filter_soilp(:)
     integer, public :: num_jtops
     integer, public, allocatable :: jtops(:)
     integer, public :: num_soilc
     integer, public, allocatable :: filter_soilc(:)

     ! FIXME(bja, 201603) most of these types should be private!

     ! NOTE(bja, 201603) BeTR types only, no LSM specific types here!

   contains
     procedure, public :: BeTRInit
     procedure, public :: BeTRSetFilter
     procedure, public :: Init => BeTRSimulationInit
     procedure, public :: ReadNameList => BeTRSimulationReadNameList
     procedure, public :: RestartInit => BeTRSimulationRestartInit
     procedure, public :: ConsistencyCheck => BeTRSimulationConsistencyCheck
     procedure, public :: StepWithoutDrainage => BeTRSimulationStepWithoutDrainage
     procedure, public :: StepWithDrainage => BeTRSimulationStepWithDrainage
     procedure, public :: BeginMassBalanceCheck => BeTRSimulationBeginMassBalanceCheck
     procedure, public :: MassBalanceCheck      => BeTRSimulationMassBalanceCheck
     procedure, public :: CreateHistory => hist_htapes_create
     procedure, public :: WriteHistory => hist_write
     procedure, public :: WriteRegressionOutput
  end type betr_simulation_type

  public :: BeTRSimulationInit

contains

  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationInit(this, base_filename, namelist_buffer, &
       bounds, waterstate, cnstate)
    ! Dummy routine for inheritance purposes. don't use.

    use WaterstateType        , only : waterstate_type
    use CNStateType           , only : cnstate_type

    use betr_constants, only : betr_namelist_buffer_size, betr_filename_length


    implicit none

    class(betr_simulation_type), intent(inout) :: this
    character(len=betr_filename_length), intent(in) :: base_filename
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer

    type(bounds_type)    , intent(in) :: bounds
    type(waterstate_type), intent(inout) :: waterstate
    type(cnstate_type), intent(inout) :: cnstate

    character(len=*), parameter :: subname = 'BeTRSimulationInit'


    call endrun(msg="ERROR "//subname//" unimplemented. "//errmsg(mod_filename, __LINE__))


  end subroutine BeTRSimulationInit

!-------------------------------------------------------------------------------
  subroutine BeTRSetFilter(this)


  implicit none

    class(betr_simulation_type), intent(inout) :: this
    this%num_jtops = 1
    allocate(this%jtops(this%num_jtops))
    this%jtops(:) = 1

    this%num_soilc = 1
    allocate(this%filter_soilc(this%num_soilc))
    this%filter_soilc(:) = 1

    this%num_soilp = 1
    allocate(this%filter_soilp(this%num_soilp))
    this%filter_soilp(:) = 1

  end subroutine BeTRSetFilter
!-------------------------------------------------------------------------------

  subroutine BeTRInit(this, base_filename, namelist_buffer, &
       betr_bounds, betr_waterstate, betr_cnstate)
    !
    use BeTR_WaterStateType, only : betr_waterstate_type
    use BeTR_CNStateType, only : betr_cnstate_type

    use betr_constants, only : betr_namelist_buffer_size
    use betr_constants, only : betr_filename_length

    implicit none

    class(betr_simulation_type), intent(inout) :: this
    character(len=betr_filename_length), intent(in) :: base_filename
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer

    type(betr_bounds_type)    , intent(in) :: betr_bounds
    type(betr_waterstate_type), intent(inout) :: betr_waterstate
    type(betr_cnstate_type), intent(inout) :: betr_cnstate

    character(len=*), parameter :: subname = 'BeTRInit'

    this%base_filename = base_filename
    
    call this%betr%Init(namelist_buffer, betr_bounds, betr_waterstate, betr_cnstate)

    call this%CreateHistory(betr_nlevtrc_soil, this%num_soilc)

    call this%regression%Init(base_filename, namelist_buffer)

  end subroutine BeTRInit

  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationReadNameList(this, filename)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_type), intent(inout) :: this
    character(len=*), intent(IN) :: filename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                    ! error code
    integer :: unitn                   ! unit for namelist file
    character(len=betr_string_length) :: reaction_method

    character(len=*), parameter :: subname = 'BeTRSimulationReadNameList'

    !-----------------------------------------------------------------------

    namelist / betr_inparm / reaction_method

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in betr_inparm  namelist'
       call opnfil (filename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'betr_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, betr_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading betr_inparm namelist"//errmsg(filename, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if
    ! Broadcast namelist variables read in
    call shr_mpi_bcast(reaction_method, mpicom)

  end subroutine BeTRSimulationReadNameList

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartInit(this, bounds, ncid, flag)
    !
    !! DESCRIPTION
    ! initialize for restart run
    ! !USES:
    use ncdio_pio, only : file_desc_t
    implicit none

    class(betr_simulation_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    class(file_desc_t), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in)    :: flag ! 'read' or 'write'

    type(betr_bounds_type)     :: betr_bounds
    integer :: lbj, ubj

    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    call this%betr%tracerstates%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betr%tracers)

    call this%betr%tracerfluxes%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betr%tracers)

    call this%betr%tracercoeffs%Restart(betr_bounds, ncid, flag=flag, betrtracer_vars=this%betr%tracers)
  end subroutine BeTRSimulationRestartInit


  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithoutDrainage(this, bounds, col, &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)

    use SoilStateType, only : soilstate_type
    use WaterStateType, only : Waterstate_Type
    use TemperatureType, only : temperature_type
    use ChemStateType, only : chemstate_type
    use WaterfluxType, only : waterflux_type
    use ColumnType, only : column_type
    use atm2lndType, only : atm2lnd_type
    use SoilHydrologyType, only : soilhydrology_type
    use BeTR_CNStateType, only : betr_cnstate_type
    use BeTR_CarbonFluxType, only : betr_carbonflux_type
    use CNStateType, only : cnstate_type
    use CNCarbonFluxType, only : carbonflux_type
    use CanopyStateType, only : canopystate_type
    use BeTR_PatchType, only : betr_pft
    use PatchType, only : pft
    use pftvarcon, only : crop

    implicit none

    class(betr_simulation_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds ! bounds


    type(column_type), intent(in) :: col ! column type
    type(Waterstate_Type), intent(in) :: waterstate_vars ! water state variables
    type(soilstate_type), intent(in) :: soilstate_vars ! column physics variable
    type(temperature_type), intent(in) :: temperature_vars ! energy state variable
    type(chemstate_type), intent(in) :: chemstate_vars
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    type(cnstate_type), intent(inout) :: cnstate_vars
    type(canopystate_type), intent(in) :: canopystate_vars
    type(carbonflux_type), intent(in) :: carbonflux_vars
    type(waterflux_type), intent(inout) :: waterflux_vars !temporary variables


  end subroutine BeTRSimulationStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithDrainage(this, bounds, waterflux_vars, col)

    use ColumnType, only : column_type
    use MathfuncMod, only : safe_div
    use WaterFluxType, only : waterflux_type

    implicit none

    class(betr_simulation_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(waterflux_type)    , intent(in) :: waterflux_vars
    type(column_type), intent(in) :: col ! column type

  end subroutine BeTRSimulationStepWithDrainage


  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationBeginMassBalanceCheck(this, bounds)

    use TracerBalanceMod, only : begin_betr_tracer_massbalance

    implicit none

    class(betr_simulation_type) :: this
    type(bounds_type), intent(in) :: bounds

    type(betr_bounds_type)     :: betr_bounds
    integer  :: lbj, ubj
    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    call begin_betr_tracer_massbalance(betr_bounds, lbj, ubj, &
         this%num_soilc, this%filter_soilc, &
         this%betr%tracers, this%betr%tracerstates, &
         this%betr%tracerfluxes)

  end  subroutine BeTRSimulationBeginMassBalanceCheck
  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationMassBalanceCheck(this, bounds)

    use TracerBalanceMod, only : betr_tracer_massbalance_check

    implicit none

    class(betr_simulation_type) :: this
    type(bounds_type), intent(in) :: bounds

    integer :: lbj, ubj
    type(betr_bounds_type)     :: betr_bounds

    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    call betr_tracer_massbalance_check(betr_bounds, lbj, ubj, &
         this%num_soilc, this%filter_soilc, &
         this%betr%tracers, this%betr%tracerstates, &
         this%betr%tracerfluxes)
  end subroutine BeTRSimulationMassBalanceCheck

!-------------------------------------------------------------------------------

  subroutine hist_htapes_create(this, nlevtrc_soil, ncol)
  !
  ! DESCRIPTIONS
  ! create history file and define output variables

    use netcdf
    use ncdio_pio
    use histFileMod, only : hist_file_create, hist_def_fld1d, hist_def_fld2d

    !
    !ARGUMENTS
    implicit none

    class(betr_simulation_type), intent(inout) :: this
    integer              , intent(in) :: nlevtrc_soil, ncol

    integer            :: jj, kk
    type(file_desc_t)  :: ncid
    character(len=*), parameter :: subname = 'hist_htapes_create'

    associate( &
         ntracers => this%betr%tracers%ntracers, &
         ngwmobile_tracers => this%betr%tracers%ngwmobile_tracers, &
         is_volatile => this%betr%tracers%is_volatile, &
         is_h2o => this%betr%tracers%is_h2o, &
         is_isotope => this%betr%tracers%is_isotope, &
         volatileid => this%betr%tracers%volatileid, &
         tracernames => this%betr%tracers%tracernames &
         )

      this%hist_filename = trim(this%base_filename) // '.output.nc'
      
      call ncd_pio_createfile(ncid, this%hist_filename)

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
  subroutine hist_write(me, record, lbj, ubj, time_vars)
    !
    ! DESCRIPTION
    ! output hist file
    !
    use shr_kind_mod        , only : r8 => shr_kind_r8
    use ncdio_pio

    use BeTR_TimeMod, only : betr_time_type

    implicit none

    class(betr_simulation_type), intent(inout) :: me
    integer, intent(in) :: record
    integer, intent(in) :: lbj,ubj
    type(betr_time_type), intent(in) :: time_vars

    type(file_desc_t) :: ncid
    integer :: jj
    character(len=*), parameter :: subname='hist_write'

    associate( &
         ntracers => me%betr%tracers%ntracers, &
         ngwmobile_tracers => me%betr%tracers%ngwmobile_tracers, &
         is_volatile => me%betr%tracers%is_volatile, &
         is_h2o => me%betr%tracers%is_h2o, &
         is_isotope => me%betr%tracers%is_isotope, &
         volatileid => me%betr%tracers%volatileid, &
         tracernames => me%betr%tracers%tracernames &
         )

      call ncd_pio_openfile_for_write(ncid, me%hist_filename)

      if (mod(time_vars%time, 86400._r8)==0) then
         print*,'day', time_vars%time/86400._r8
      end if
      call ncd_putvar(ncid, "time", record, time_vars%time)

      do jj = 1, ntracers
         if(jj<= ngwmobile_tracers)then
            call ncd_putvar(ncid,trim(tracernames(jj))//'_TRACER_CONC_MOIBLE',&
                 record,me%betr%tracerstates%tracer_conc_mobile_col(1:1,lbj:ubj,jj))

            if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then

               call ncd_putvar(ncid, trim(tracernames(jj))//'_TRACER_P_GAS_FRAC',&
                    record,me%betr%tracerstates%tracer_P_gas_frac_col(1:1,lbj:ubj,volatileid(jj)))
            endif
            if(is_volatile(jj))then
               call ncd_putvar(ncid, trim(tracernames(jj))//'_FLX_SURFEMI', &
                    record, me%betr%tracerfluxes%tracer_flx_surfemi_col(1:1, volatileid(jj)))
            endif
         else
            call ncd_putvar(ncid, trim(tracernames(jj))//'_TRACER_CONC_SOLID_PASSIVE', &
                 record, me%betr%tracerstates%tracer_conc_solid_passive_col(1:1,lbj:ubj,jj-ngwmobile_tracers))
         endif
      enddo

      call ncd_putvar(ncid, 'TRACER_P_GAS', record, me%betr%tracerstates%tracer_P_gas_col(1:1,lbj:ubj))
      call ncd_pio_closefile(ncid)
    end associate
  end subroutine hist_write

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationConsistencyCheck(this, &
     bounds, ubj, num_soilc, filter_soilc, waterstate_vars)
    use WaterStateType, only : Waterstate_Type
  implicit none
    class(betr_simulation_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: ubj
    type(Waterstate_Type), intent(in) :: waterstate_vars ! water state variables

  end subroutine BeTRSimulationConsistencyCheck

  !---------------------------------------------------------------------------------

  subroutine WriteRegressionOutput(this)

    use betr_constants, only : betr_string_length
    
    implicit none

    class(betr_simulation_type), intent(inout) :: this

    integer :: jj, tt, begc, endc
    character(len=betr_string_length) :: category

    ! FIXME(bja, 201603) need a way to categorize output variables,
    ! e.g. concentration, velocity, etc. Hard coding for now.
    category = 'concentration'

    begc = 1
    endc = 1
    
    if (this%regression%write_regression_output) then
       call this%regression%OpenOutput()
       do tt = 1, this%betr%tracers%ntracers
          if (tt < this%betr%tracers%ngwmobile_tracers) then
             call this%regression%WriteData(category, &
                  this%betr%tracers%tracernames(tt), &
                  this%betr%tracerstates%tracer_conc_mobile_col(begc, :, tt))
          end if
       end do
       call this%regression%CloseOutput()
    end if
  end subroutine WriteRegressionOutput
end module BeTRSimulation
