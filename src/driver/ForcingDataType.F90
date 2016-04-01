module ForcingDataType

  use bshr_kind_mod        , only : r8 => shr_kind_r8
  use bshr_const_mod       , only : rhoice => SHR_CONST_RHOICE
  use bshr_const_mod       , only : rhoh2o => SHR_CONST_RHOFW

  use betr_constants, only : betr_filename_length
  use betr_constants, only : betr_string_length, betr_string_length_long
  use betr_constants, only : betr_namelist_buffer_size

  implicit none

  private

  character(len=*), private, parameter :: mod_filename = __FILE__

  type, public :: ForcingData_type
     character(len=betr_filename_length) :: grid_filename
     character(len=betr_filename_length) :: forcing_filename
     character(len=betr_string_length) :: forcing_format ! file format: netcdf, namelist, csv, etc.
     character(len=betr_string_length) :: forcing_type ! steady state, transient


     integer :: num_levels
     integer :: num_time
     integer :: num_columns

     real(r8), pointer :: t_soi(:,:)
     real(r8), pointer :: h2osoi_liqvol(:,:)
     real(r8), pointer :: h2osoi_liq(:,:)
     real(r8), pointer :: h2osoi_ice(:,:)
     real(r8), pointer :: qflx_infl(:)       !surface infiltration, mm/s
     real(r8), pointer :: qflx_rootsoi(:,:)  !transpiration at depth, m/s
     real(r8), pointer :: pbot(:)            !amtospheric pressure, Pa
     real(r8), pointer :: tbot(:)            !atmoshperic temperature, kelvin
     real(r8), pointer :: dzsoi(:)           !node thickness
     real(r8), pointer :: zsoi(:)            !node depth of each numerical node
     real(r8), pointer :: bsw(:)             !clap-hornberg parameter
     real(r8), pointer :: watsat(:)          !saturated volumetric water content
     real(r8), pointer :: h2osoi_icevol(:,:)
     real(r8), pointer :: qbot(:)            !water flux at bottom boundary, mm/s

   contains

     procedure, public :: Init
     procedure, public :: ReadData
     procedure, public :: ReadGridData
     procedure, public :: ReadForcingData
!     procedure, public :: GetForcing
     procedure, private :: InitAllocate
     procedure, private :: ReadNameList
  end type ForcingData_type

contains

  subroutine Init(this, diml, dimt)

    class(ForcingData_type) :: this

    integer, intent(in) :: diml, dimt

    call this%InitAllocate(diml, dimt)

  end subroutine Init


  !------------------------------------------------------------------------
  subroutine InitAllocate(this, dim_levels, dim_time)

    class(ForcingData_type) :: this
    !at this moment the variable size is fixed

    integer, intent(in) :: dim_levels, dim_time

    this%num_columns = 1
    this%num_levels=dim_levels
    this%num_time=dim_time

    allocate(this%t_soi(this%num_time, this%num_levels))
    allocate(this%h2osoi_liqvol(this%num_time, this%num_levels))
    allocate(this%h2osoi_icevol(this%num_time, this%num_levels))
    allocate(this%h2osoi_liq(this%num_time, this%num_levels))
    allocate(this%h2osoi_ice(this%num_time, this%num_levels))
    allocate(this%qflx_infl(this%num_time))
    allocate(this%qflx_rootsoi(this%num_time, this%num_levels))
    allocate(this%zsoi(this%num_levels))
    allocate(this%dzsoi(this%num_levels))
    allocate(this%pbot(this%num_time))
    allocate(this%qbot(this%num_time))
    allocate(this%tbot(this%num_time))
    allocate(this%watsat(this%num_levels))
    allocate(this%bsw(this%num_levels))

  end subroutine InitAllocate


  !------------------------------------------------------------------------
  subroutine ReadData(this, namelist_buffer)

    use ncdio_pio

    implicit none

    class(ForcingData_type) :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer

    character(len=250) :: ncf_in_filename_forc
    type(file_desc_t)  :: ncf_in_forc

    integer :: num_levels, num_time

    call this%ReadNameList(namelist_buffer)

    ! FIXME(bja, 201603) Ugh, have to crack open the netcdf file to
    ! retreive dimensions.... should we read from the namelist?
    ncf_in_filename_forc=trim(this%forcing_filename)
    call ncd_pio_openfile(ncf_in_forc, ncf_in_filename_forc, mode=ncd_nowrite)
    num_levels = get_dim_len(ncf_in_forc, 'levgrnd')
    num_time = get_dim_len(ncf_in_forc, 'time')
    call ncd_pio_closefile(ncf_in_forc)

    call this%Init(num_levels, num_time)
    call this%ReadGridData()
    call this%ReadForcingData()

  end subroutine ReadData

  !------------------------------------------------------------------------
  subroutine ReadGridData(this)

    use ncdio_pio

    implicit none

    class(ForcingData_type) :: this

    character(len=250) :: ncf_in_filename_grid
    type(file_desc_t)  :: ncf_in_grid
    real(r8), allocatable :: data(:,:)
    integer :: j2


    ncf_in_filename_grid=trim(this%grid_filename)
    call ncd_pio_openfile(ncf_in_grid, ncf_in_filename_grid, mode=ncd_nowrite)

    allocate(data(this%num_columns, this%num_levels))
    call ncd_getvar(ncf_in_grid, 'DZSOI', data)
    do j2 = 1, this%num_levels
       this%dzsoi(j2) = data(this%num_columns, j2)
    enddo

    call ncd_getvar(ncf_in_grid, 'ZSOI', data)
    do j2 = 1, this%num_levels
       this%zsoi(j2) = data(this%num_columns, j2)
    enddo

    call ncd_getvar(ncf_in_grid, 'WATSAT', data)
    do j2 = 1, this%num_levels
       this%watsat(j2) = data(this%num_columns, j2)
    enddo

    call ncd_getvar(ncf_in_grid, 'BSW', data)
    do j2 = 1, this%num_levels
       this%bsw(j2) = data(this%num_columns, j2)
    enddo

    call ncd_pio_closefile(ncf_in_grid)

    deallocate(data)

  end subroutine ReadGridData

  !------------------------------------------------------------------------
  subroutine ReadForcingData(this)

    use ncdio_pio

    implicit none

    class(ForcingData_type) :: this

    character(len=250) :: ncf_in_filename_forc
    type(file_desc_t)  :: ncf_in_forc
    real(r8), allocatable :: data_2d(:,:,:)
    real(r8), allocatable :: data_1d(:,:)
    integer :: j1, j2

    ncf_in_filename_forc=trim(this%forcing_filename)
    call ncd_pio_openfile(ncf_in_forc, ncf_in_filename_forc, mode=ncd_nowrite)

    allocate(data_2d(this%num_columns, this%num_levels, 1:this%num_time))
    allocate(data_1d(this%num_columns, this%num_time))

    call ncd_getvar(ncf_in_forc, 'TSOI', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%t_soi(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'H2OSOI', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%h2osoi_liqvol(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'SOILICE', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time

          this%h2osoi_icevol(j1, j2) = data_2d(this%num_columns, j2, j1)/rhoice/this%dzsoi(j2)
          this%h2osoi_liqvol(j1, j2) = this%h2osoi_liqvol(j1, j2) - this%h2osoi_icevol(j1, j2)
          this%h2osoi_ice(j1, j2) = data_2d(this%num_columns, j2, j1)
          this%h2osoi_liq(j1, j2) = this%h2osoi_liqvol(j1, j2)*this%dzsoi(j2)*rhoh2o
       enddo
    enddo


    call ncd_getvar(ncf_in_forc, 'QINFL', data_1d)
    do j1 =1, this%num_time
       this%qflx_infl(j1) = data_1d(this%num_columns, j1)
    enddo

    call ncd_getvar(ncf_in_forc, 'PBOT', data_1d)
    do j1 =1, this%num_time
       this%pbot(j1) = data_1d(this%num_columns, j1)
    enddo

    call ncd_getvar(ncf_in_forc, 'TBOT', data_1d)
    do j1 =1, this%num_time
       this%tbot(j1) = data_1d(this%num_columns, j1)
    enddo

    call ncd_getvar(ncf_in_forc, 'QCHARGE', data_1d)
    do j1 =1, this%num_time
       this%qbot(j1) = data_1d(this%num_columns, j1) * 1.e-3_r8  !convert int m/s
    enddo

    call ncd_getvar(ncf_in_forc, 'QFLX_ROOTSOI', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%qflx_rootsoi(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_pio_closefile(ncf_in_forc)

    deallocate(data_2d)
    deallocate(data_1d)

  end subroutine ReadForcingData

  ! ----------------------------------------------------------------------

  subroutine ReadNameList(this, namelist_buffer)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use clm_varctl   , only : iulog
    use abortutils      , only : endrun
    use shr_log_mod     , only : errMsg => shr_log_errMsg

    use betr_constants, only : stdout

    implicit none
    ! !ARGUMENTS:
    class(ForcingData_type) :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    !
    ! !LOCAL VARIABLES:
    integer :: nml_error
    character(len=*), parameter :: subname = 'ReadNameList'
    character(len=betr_filename_length) :: grid_filename
    character(len=betr_filename_length) :: forcing_format, forcing_type, forcing_filename
    character(len=betr_string_length_long) :: ioerror_msg


    !-----------------------------------------------------------------------

    namelist / forcing_inparm / grid_filename, &
         forcing_type, forcing_filename, forcing_format

    grid_filename = ''
    forcing_format = ''
    forcing_type = ''
    forcing_filename = ''

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=forcing_inparm, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading forcing_inparm namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr climate forcing :'
       write(stdout, *)
       write(stdout, *) ' forcing_inparm namelist settings :'
       write(stdout, *)
       write(stdout, forcing_inparm)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    this%grid_filename = trim(grid_filename)
    this%forcing_type = trim(forcing_type)
    this%forcing_format = trim(forcing_format)
    this%forcing_filename = trim(forcing_filename)

  end subroutine ReadNameList



end module ForcingDataType
