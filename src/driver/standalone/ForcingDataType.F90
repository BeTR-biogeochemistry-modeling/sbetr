module ForcingDataType

  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use bshr_const_mod , only : rhoice => SHR_CONST_RHOICE
  use bshr_const_mod , only : rhoh2o => SHR_CONST_RHOFW
  use betr_constants , only : betr_filename_length
  use betr_constants , only : betr_string_length, betr_string_length_long
  use betr_constants , only : betr_namelist_buffer_size

  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  character(len=*), parameter          :: steady_state_name = 'steady state'
  character(len=*), parameter          :: transient_name = 'transient'
  integer, parameter                   :: transient = 1
  integer, parameter                   :: steady_state = 2

  type, public :: ForcingData_type
     character(len=betr_filename_length) :: forcing_filename
     character(len=betr_string_length)   :: forcing_format ! file format: netcdf, namelist, csv, etc.
     character(len=betr_string_length)   :: forcing_type_name ! steady state, transient
     integer                             :: forcing_type
     integer                             :: num_levels
     integer                             :: num_time
     integer                             :: num_columns
     real(r8), allocatable               :: t_soi(:,:)
     real(r8), allocatable               :: h2osoi_liqvol(:,:)
     real(r8), allocatable               :: h2osoi_icevol(:,:)
     real(r8), allocatable               :: h2osoi_liq(:,:)
     real(r8), allocatable               :: h2osoi_ice(:,:)
     real(r8), allocatable               :: qflx_infl(:)       !surface infiltration, mm/s
     real(r8), allocatable               :: qflx_rootsoi(:,:)   !transpiration at depth, m/s
     real(r8), allocatable               :: pbot(:)            !amtospheric pressure, Pa
     real(r8), allocatable               :: tbot(:)            !atmoshperic temperature, kelvin
     real(r8), allocatable               :: qbot(:)            !water flux at bottom boundary, mm/s
     real(r8), allocatable               :: finundated(:)
     real(r8), allocatable               :: nflx_nh4_vr(:,:)
     real(r8), allocatable               :: nflx_no3_vr(:,:)
     real(r8), allocatable               :: pflx_po4_vr(:,:)
     real(r8), allocatable               :: cflx_met_vr(:,:)
     real(r8), allocatable               :: cflx_cel_vr(:,:)
     real(r8), allocatable               :: cflx_lig_vr(:,:)
     real(r8), allocatable               :: cflx_cwd_vr(:,:)
     real(r8), allocatable               :: nflx_met_vr(:,:)
     real(r8), allocatable               :: nflx_cel_vr(:,:)
     real(r8), allocatable               :: nflx_lig_vr(:,:)
     real(r8), allocatable               :: nflx_cwd_vr(:,:)
     real(r8), allocatable               :: pflx_met_vr(:,:)
     real(r8), allocatable               :: pflx_cel_vr(:,:)
     real(r8), allocatable               :: pflx_lig_vr(:,:)
     real(r8), allocatable               :: pflx_cwd_vr(:,:)
     real(r8), allocatable               :: qflx_runoff_col(:)
     real(r8), allocatable               :: rr_vr(:,:)
     logical                             :: use_rootsoit
   contains
     procedure, public :: Init
     procedure, public :: ReadData
     procedure, public :: ReadForcingData
     procedure, public :: ReadCNPData
     procedure, public :: UpdateForcing
     procedure, public :: UpdateCNPForcing
     procedure, public :: discharge
     procedure, public :: infiltration
     procedure, private :: InitAllocate
     procedure, private :: InitAllocate_CNP
     procedure, private :: ReadNameList
     procedure, public  :: Destroy
  end type ForcingData_type

contains

  subroutine Init(this, dim_levels, dim_time)
    !DESCRIPTION
    !initialize
    implicit none
    !ARGUMENTS
    class(ForcingData_type), intent(inout) :: this
    integer, intent(in)     :: dim_levels, dim_time

    this%num_columns = 1
    this%num_levels  =dim_levels
    this%num_time    =dim_time

    select case (trim(this%forcing_type_name))
       case (transient_name)
          this%forcing_type = transient
       case (steady_state_name)
          this%forcing_type = steady_state
          ! ignore the values scraped from netcdf and just use one time level.
          this%num_time = 1
       case default
          this%forcing_type = transient
          write(*, *) 'WARNING: no forcing data type specified, using transient.'
       end select

    call this%InitAllocate()

  end subroutine Init


  !------------------------------------------------------------------------
  subroutine Destroy(this)
  !DESCRIPTION
  !allocate memory
  implicit none
  class(ForcingData_type), intent(inout) :: this
  !at this moment the variable size is fixed

  if(allocated(this%t_soi))deallocate(this%t_soi)
  if(allocated(this%h2osoi_liqvol))deallocate(this%h2osoi_liqvol)
  if(allocated(this%h2osoi_icevol))deallocate(this%h2osoi_icevol)
  if(allocated(this%h2osoi_liq))deallocate(this%h2osoi_liq)
  if(allocated(this%h2osoi_ice))deallocate(this%h2osoi_ice)
  if(allocated(this%qflx_infl))deallocate(this%qflx_infl)
  if(allocated(this%qflx_rootsoi))deallocate(this%qflx_rootsoi)
  if(allocated(this%pbot))deallocate(this%pbot)
  if(allocated(this%qbot))deallocate(this%qbot)
  if(allocated(this%tbot))deallocate(this%tbot)
  if(allocated(this%finundated))deallocate(this%finundated)
  if(allocated(this%nflx_nh4_vr))deallocate(this%nflx_nh4_vr)
  if(allocated(this%nflx_no3_vr))deallocate(this%nflx_no3_vr)
  if(allocated(this%pflx_po4_vr))deallocate(this%pflx_po4_vr)
  if(allocated(this%cflx_met_vr))deallocate(this%cflx_met_vr)
  if(allocated(this%cflx_cel_vr))deallocate(this%cflx_cel_vr)
  if(allocated(this%cflx_lig_vr))deallocate(this%cflx_lig_vr)
  if(allocated(this%cflx_cwd_vr))deallocate(this%cflx_cwd_vr)
  if(allocated(this%rr_vr))deallocate(this%rr_vr)
  if(allocated(this%nflx_met_vr))deallocate(this%nflx_met_vr)
  if(allocated(this%nflx_cel_vr))deallocate(this%nflx_cel_vr)
  if(allocated(this%nflx_lig_vr))deallocate(this%nflx_lig_vr)
  if(allocated(this%nflx_cwd_vr))deallocate(this%nflx_cwd_vr)
  if(allocated(this%pflx_met_vr))deallocate(this%pflx_met_vr)
  if(allocated(this%pflx_cel_vr))deallocate(this%pflx_cel_vr)
  if(allocated(this%pflx_lig_vr))deallocate(this%pflx_lig_vr)
  if(allocated(this%pflx_cwd_vr))deallocate(this%pflx_cwd_vr)
  end subroutine Destroy

  !------------------------------------------------------------------------
  subroutine InitAllocate(this)
    !DESCRIPTION
    !allocate memory
    implicit none
    class(ForcingData_type), intent(inout) :: this
    !at this moment the variable size is fixed

    allocate(this%t_soi(this%num_time, this%num_levels))
    allocate(this%h2osoi_liqvol(this%num_time, this%num_levels))
    allocate(this%h2osoi_icevol(this%num_time, this%num_levels))
    allocate(this%h2osoi_liq(this%num_time, this%num_levels))
    allocate(this%h2osoi_ice(this%num_time, this%num_levels))
    allocate(this%qflx_infl(this%num_time))
    allocate(this%qflx_rootsoi(this%num_time, this%num_levels))
    allocate(this%pbot(this%num_time))
    allocate(this%qbot(this%num_time))
    allocate(this%tbot(this%num_time))
    allocate(this%finundated(this%num_time)); this%finundated(:)=0._r8
    allocate(this%qflx_runoff_col(this%num_time)); this%qflx_runoff_col(:)=0._r8
  end subroutine InitAllocate
  !------------------------------------------------------------------------
  subroutine InitAllocate_CNP(this)
    !DESCRIPTION
    !allocate memory
    implicit none
    class(ForcingData_type), intent(inout) :: this
    !at this moment the variable size is fixed

    allocate(this%cflx_met_vr(this%num_time, this%num_levels))
    allocate(this%cflx_cel_vr(this%num_time, this%num_levels))
    allocate(this%cflx_lig_vr(this%num_time, this%num_levels))
    allocate(this%cflx_cwd_vr(this%num_time, this%num_levels))
    allocate(this%rr_vr(this%num_time, this%num_levels)); this%rr_vr(:,:) = 0._r8
    allocate(this%nflx_met_vr(this%num_time, this%num_levels))
    allocate(this%nflx_cel_vr(this%num_time, this%num_levels))
    allocate(this%nflx_lig_vr(this%num_time, this%num_levels))
    allocate(this%nflx_cwd_vr(this%num_time, this%num_levels))
    allocate(this%pflx_met_vr(this%num_time, this%num_levels))
    allocate(this%pflx_cel_vr(this%num_time, this%num_levels))
    allocate(this%pflx_lig_vr(this%num_time, this%num_levels))
    allocate(this%pflx_cwd_vr(this%num_time, this%num_levels))
    allocate(this%nflx_nh4_vr(this%num_time, this%num_levels))
    allocate(this%nflx_no3_vr(this%num_time, this%num_levels))
    allocate(this%pflx_po4_vr(this%num_time, this%num_levels))

  end subroutine InitAllocate_CNP

  !------------------------------------------------------------------------
  subroutine ReadCNPData(this)
    !read infomration about forcing data
    !USES
    use ncdio_pio    , only : file_desc_t
    use ncdio_pio    , only : ncd_nowrite
    use ncdio_pio    , only : ncd_pio_openfile
    use ncdio_pio    , only : get_dim_len
    use ncdio_pio    , only : ncd_getvar, Var_desc_t
    use ncdio_pio    , only : ncd_pio_closefile, check_var
    use babortutils  , only : endrun
    use bshr_log_mod , only : errMsg => shr_log_errMsg
    use BeTR_GridMod , only : betr_grid_type
  implicit none
    class(ForcingData_type), intent(inout)  :: this

    character(len=250)    :: ncf_in_filename_forc
    type(file_desc_t)     :: ncf_in_forc
    real(r8), allocatable :: data_2d(:,:,:)
    integer :: j1, j2
    type(Var_desc_t)  :: vardesc
    logical :: readvar

    call this%InitAllocate_CNP()

    ncf_in_filename_forc=trim(this%forcing_filename)
    call ncd_pio_openfile(ncf_in_forc, ncf_in_filename_forc, mode=ncd_nowrite)

    allocate(data_2d(this%num_columns, this%num_levels, 1:this%num_time))

    call ncd_getvar(ncf_in_forc, 'CFLX_INPUT_LITR_MET_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%cflx_met_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'CFLX_INPUT_LITR_CEL_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%cflx_cel_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'CFLX_INPUT_LITR_LIG_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%cflx_lig_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'CFLX_INPUT_LITR_CWD_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%cflx_cwd_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call check_var(ncf_in_forc, 'RR_vr', vardesc, readvar,print_err=.false.)
    if(readvar)then
      call ncd_getvar(ncf_in_forc, 'RR_vr', data_2d)
      do j2 = 1, this%num_levels
         do j1 = 1, this%num_time
           this%rr_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
         enddo
      enddo
    endif

    call ncd_getvar(ncf_in_forc, 'NFLX_INPUT_LITR_MET_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%nflx_met_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'NFLX_INPUT_LITR_CEL_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%nflx_cel_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'NFLX_INPUT_LITR_LIG_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%nflx_lig_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'NFLX_INPUT_LITR_CWD_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%nflx_cwd_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'PFLX_INPUT_LITR_MET_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%pflx_met_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'PFLX_INPUT_LITR_CEL_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%pflx_cel_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'PFLX_INPUT_LITR_LIG_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%pflx_lig_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'PFLX_INPUT_LITR_CWD_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%pflx_cwd_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'NFLX_MINN_INPUT_NH4_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%nflx_nh4_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'NFLX_MINN_INPUT_NO3_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%nflx_no3_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_getvar(ncf_in_forc, 'PFLX_MINP_INPUT_PO4_vr', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%pflx_po4_vr(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    call ncd_pio_closefile(ncf_in_forc)

    deallocate(data_2d)
  end subroutine ReadCNPData
  !------------------------------------------------------------------------
  subroutine ReadData(this, namelist_buffer, grid)
    !DESCRIPTION
    !read infomration about forcing data
    !USES
    use ncdio_pio    , only : file_desc_t
    use ncdio_pio    , only : ncd_nowrite
    use ncdio_pio    , only : ncd_pio_openfile
    use ncdio_pio    , only : get_dim_len
    use ncdio_pio    , only : ncd_getvar
    use ncdio_pio    , only : ncd_pio_closefile
    use babortutils  , only : endrun
    use bshr_log_mod , only : errMsg => shr_log_errMsg
    use BeTR_GridMod , only : betr_grid_type
    implicit none
    !ARGUMENTS
    class(ForcingData_type), intent(inout)  :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    class(betr_grid_type), intent(in)                    :: grid
    !TEMPORARY VARIABLES
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

    if (grid%nlevgrnd /= num_levels) then
       call endrun(msg="ERROR inconsistent vertical levels between "//&
            "grid and forcing. "//errmsg(mod_filename, __LINE__))

    end if

    call this%Init(num_levels, num_time)
    call this%ReadForcingData(grid)
    !x print*,'read data tsoi',this%t_soi(1,:)
  end subroutine ReadData

  !------------------------------------------------------------------------
  subroutine ReadForcingData(this, grid)
    !DESCRIPTION
    !read forcing data
    !USES
    use ncdio_pio    , only : file_desc_t
    use ncdio_pio    , only : ncd_nowrite
    use ncdio_pio    , only : ncd_pio_openfile
    use ncdio_pio    , only : ncd_getvar, ncd_getatt
    use ncdio_pio    , only : ncd_pio_closefile
    use BeTR_GridMod , only : betr_grid_type
    implicit none
    !ARGUMENTS
    class(ForcingData_type), intent(inout)  :: this
    class(betr_grid_type), intent(in) :: grid
    !TEMPORARY VARIABLES
    character(len=250)    :: ncf_in_filename_forc
    type(file_desc_t)     :: ncf_in_forc
    real(r8), allocatable :: data_2d(:,:,:)
    real(r8), allocatable :: data_1d(:,:)
    integer               :: j1, j2
    real(r8) :: tommps
    character(len=9) :: units

    ncf_in_filename_forc=trim(this%forcing_filename)
    call ncd_pio_openfile(ncf_in_forc, ncf_in_filename_forc, mode=ncd_nowrite)

    allocate(data_2d(this%num_columns, this%num_levels, 1:this%num_time))
    allocate(data_1d(this%num_columns, this%num_time))

    !X!write(*, *) 'Reading TSOI'
    call ncd_getvar(ncf_in_forc, 'TSOI', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%t_soi(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    !X!write(*, *) 'Reading H2OSOI'
    call ncd_getvar(ncf_in_forc, 'H2OSOI', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%h2osoi_liqvol(j1, j2) = data_2d(this%num_columns, j2, j1)
       enddo
    enddo

    !X!write(*, *) 'Reading SOILICE'
    call ncd_getvar(ncf_in_forc, 'SOILICE', data_2d)
    do j2 = 1, this%num_levels
       do j1 = 1, this%num_time
          this%h2osoi_icevol(j1, j2) = data_2d(this%num_columns, j2, j1)/rhoice/grid%dzsoi(j2)
          this%h2osoi_liqvol(j1, j2) = this%h2osoi_liqvol(j1, j2) - this%h2osoi_icevol(j1, j2)
          this%h2osoi_ice(j1, j2) = data_2d(this%num_columns, j2, j1)
          this%h2osoi_liq(j1, j2) = this%h2osoi_liqvol(j1, j2)*grid%dzsoi(j2)*rhoh2o
       enddo
    enddo


    !X!write(*, *) 'Reading QINFL'
    call ncd_getvar(ncf_in_forc, 'QINFL', data_1d)
    do j1 =1, this%num_time
       this%qflx_infl(j1) = data_1d(this%num_columns, j1)
    enddo

    !X!write(*, *) 'Reading PBOT'
    call ncd_getvar(ncf_in_forc, 'PBOT', data_1d)
    do j1 =1, this%num_time
       this%pbot(j1) = data_1d(this%num_columns, j1)
    enddo

    !X!write(*, *) 'Reading TBOT'
    call ncd_getvar(ncf_in_forc, 'TBOT', data_1d)
    do j1 =1, this%num_time
       this%tbot(j1) = data_1d(this%num_columns, j1)
    enddo

    !X!write(*, *) 'Reading QCHARGE'
    call ncd_getvar(ncf_in_forc, 'QCHARGE', data_1d)
    do j1 =1, this%num_time
       this%qbot(j1) = data_1d(this%num_columns, j1)  ! mm/s
    enddo

    if(this%use_rootsoit)then
      !X!write(*, *) 'Reading QFLX_ROOTSOI'
      call ncd_getvar(ncf_in_forc, 'QFLX_ROOTSOI', data_2d)
      call ncd_getatt(ncf_in_forc,'QFLX_ROOTSOI','units',units)

      if(trim(units)=='m/s')then
        tommps=1.e3_r8
      else
        tommps=1._r8
      endif
      do j2 = 1, this%num_levels
        do j1 = 1, this%num_time
          this%qflx_rootsoi(j1, j2) = data_2d(this%num_columns, j2, j1)*1.e3_r8
        enddo
      enddo
    else
      do j2 = 1, this%num_levels
        do j1 = 1, this%num_time
          this%qflx_rootsoi(j1, j2)  = 0._r8
        enddo
      enddo
    endif
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

    use bshr_nl_mod    , only : shr_nl_find_group_name
    use betr_ctrl      , only : iulog => biulog
    use babortutils    , only : endrun
    use bshr_log_mod   , only : errMsg => shr_log_errMsg
    use betr_constants , only : stdout
    implicit none
    ! !ARGUMENTS:
    class(ForcingData_type), intent(inout)  :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    !
    ! !LOCAL VARIABLES:
    integer                                :: nml_error
    character(len=*), parameter            :: subname = 'ReadNameList'
    character(len=betr_filename_length)    :: forcing_format, forcing_type_name, forcing_filename
    character(len=betr_string_length_long) :: ioerror_msg
    logical                                :: use_rootsoit
    !-----------------------------------------------------------------------

    namelist / forcing_inparm / &
         forcing_type_name, forcing_filename, forcing_format, use_rootsoit

    forcing_format    = ''
    forcing_type_name = transient_name
    forcing_filename  = ''
    use_rootsoit =.true.
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

    this%forcing_type_name = trim(forcing_type_name)
    this%forcing_format    = trim(forcing_format)
    this%forcing_filename  = trim(forcing_filename)
    this%use_rootsoit       = use_rootsoit
  end subroutine ReadNameList

  ! ----------------------------------------------------------------------

  subroutine UpdateForcing(this, grid, bounds, lbj, ubj, numf, filter, ttime, col, pft, &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars,waterflux_vars, &
       temperature_vars,chemstate_vars, plantMicKinetics_vars, jtops)
    !
    ! DESCRIPTIONS
    ! read environmental forcing to run betr
    ! for clm forced runs, it will forcing from history files
    !
    ! USES
    use TemperatureType   , only : temperature_type
    use WaterstateType    , only : waterstate_type
    use WaterfluxType     , only : waterflux_type
    use SoilStateType     , only : soilstate_type
    use ChemStateType     , only : chemstate_type
    use ColumnType        , only : column_type
    use PatchType         , only : patch_type
    use decompMod         , only : bounds_type
    use SoilHydrologyType , only : soilhydrology_type
    use atm2lndType       , only : atm2lnd_type
    use BeTR_TimeMod      , only : betr_time_type
    use BeTR_GridMod      , only : betr_grid_type
    use betr_varcon       , only : betr_maxpatch_pft, denh2o=> bdenh2o, denice => bdenice
    use PlantMicKineticsMod, only : PlantMicKinetics_type
    implicit none
    !arguments
    class(ForcingData_type)  , intent(in)    :: this
    class(betr_grid_type)    , intent(in)    :: grid
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: numf
    integer                  , intent(in)    :: filter(:)
    integer                  , intent(in)    :: lbj, ubj
    type(betr_time_type)     , intent(in)    :: ttime
    type(chemstate_type)     , intent(inout) :: chemstate_vars
    type(atm2lnd_type)       , intent(inout) :: atm2lnd_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(column_type)        , intent(inout) :: col
    type(patch_type)         , intent(in)    :: pft
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(PlantMicKinetics_type), intent(inout) :: PlantMicKinetics_vars
    integer                  , intent(inout) :: jtops(bounds%begc:bounds%endc)

    integer            :: j, fc, c, tstep, p, pi
    character(len=255) :: subname='update_forcing'

    !X!write(*, *) 'Updating forcing data'

    if (this%forcing_type == steady_state) then
       tstep = 1
    else
       tstep = ttime%tstep
    end if
    tstep=mod(tstep,this%num_time)
    if(tstep==0)tstep=this%num_time

    !setup top boundary
    do fc = 1, numf
       c = filter(fc)
       jtops(c) = 1
       soilhydrology_vars%zwts_col(c) = 10._r8
       soilhydrology_vars%qcharge_col(c) = this%discharge(tstep)
       atm2lnd_vars%forc_pbot_downscaled_col(c) = this%pbot(tstep)  ! 1 atmos
       atm2lnd_vars%forc_t_downscaled_col(c)    = this%tbot(tstep)  ! 2 atmos temperature
    enddo

    !set up forcing variables
    do j = lbj, ubj
       do fc = 1, numf
          c = filter(fc)
          if(j>=jtops(c))then
             !take the minimum to avoid incompatible combination between uniform steady state forcing with
             !expoential grid.
             waterstate_vars%h2osoi_liqvol(c,j) = min(this%h2osoi_liqvol(tstep,j),grid%watsat(j))
             waterstate_vars%air_vol(c,j)       = grid%watsat(j)-waterstate_vars%h2osoi_liqvol(c,j)
             waterstate_vars%h2osoi_icevol(c,j) = this%h2osoi_icevol(tstep,j)
             soilstate_vars%eff_porosity_col(c,j)   = grid%watsat(j)-this%h2osoi_icevol(tstep,j)
             temperature_vars%t_soisno(c,j)     = this%t_soi(tstep,j)
             waterflux_vars%qflx_rootsoi(c,j)   = this%qflx_rootsoi(tstep,j)  !water exchange between soil and root, m/H2O/s
             do pi = 1, betr_maxpatch_pft
               if (pi <= col%npfts(c)) then
                 p = col%pfti(c) + pi - 1
                 if (pft%active(p)) then
                   waterflux_vars%qflx_rootsoi(p,j) = waterflux_vars%qflx_rootsoi(c,j)
                 endif
               endif
             enddo
             col%dz(c,j)                            = grid%dzsoi(j)
             col%zi(c,j)                            = grid%zisoi(j)
             col%z(c,j)                             = grid%zsoi(j)
             chemstate_vars%soil_pH(c,j)            = 7._r8

             !set drainage to zero
             !set surface runoff to zero
             waterflux_vars%qflx_surf(c)        = 0._r8
             waterflux_vars%qflx_drain_vr(c,j)  = 0._r8
          endif
       enddo
    enddo

    do fc = 1, numf
       c = filter(fc)
       waterflux_vars%qflx_totdrain(c)        = 0._r8
       col%zi(c,0)                            = grid%zisoi(0)
       waterflux_vars%qflx_snow2topsoi(c)     = 0._r8
       waterflux_vars%qflx_h2osfc2topsoi(c)   = 0._r8
       waterflux_vars%qflx_gross_infl_soil(c) = 0._r8
       waterflux_vars%qflx_gross_evap_soil(c) = 0._r8
       waterflux_vars%qflx_runoff(c) = this%qflx_runoff_col(c)
    enddo

    do j = 1, ubj
       do fc = 1, numf
          c = filter(fc)
          waterstate_vars%h2osoi_liq(c,j) = this%h2osoi_liq(tstep,j)
          waterstate_vars%h2osoi_ice(c,j) = this%h2osoi_ice(tstep,j)

          waterstate_vars%h2osoi_vol(c,j) = waterstate_vars%h2osoi_liqvol(c,j) + &
            waterstate_vars%h2osoi_icevol(c,j)
          waterstate_vars%h2osoi_vol(c,j) = min(waterstate_vars%h2osoi_vol(c,j), grid%watsat(j))
       enddo
    enddo

    do j = 1, ubj
      do fc = 1, numf
        PlantMicKinetics_vars%minsurf_dom_compet_vr_col(c,j)=grid%msurf_OM(j)
        PlantMicKinetics_vars%km_minsurf_dom_vr_col(c,j)=grid%KM_OM(j)
      enddo
    enddo
  end subroutine UpdateForcing

! ----------------------------------------------------------------------
  subroutine UpdateCNPForcing(this, lbj, ubj, numf, filter, ttime, &
    carbonflux_vars, c13_cflx_vars, c14_cflx_vars, nitrogenflux_vars, &
    phosphorusflux_vars, plantMicKinetics_vars)

    use CNNitrogenFluxType, only : nitrogenflux_type
    use CNCarbonFluxType  , only : carbonflux_type
    use PhosphorusFluxType, only : phosphorusflux_type
    use PlantMicKineticsMod, only : PlantMicKinetics_type
    use BeTR_TimeMod      , only : betr_time_type
  implicit none
  !arguments
    class(ForcingData_type) , intent(in) :: this
    integer                  , intent(in)    :: lbj, ubj
    integer                  , intent(in)    :: numf
    integer                  , intent(in)    :: filter(:)
    type(betr_time_type)     , intent(in)    :: ttime
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_cflx_vars
    type(carbonflux_type)    , intent(inout) :: c14_cflx_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type), intent(inout) :: phosphorusflux_vars
    type(PlantMicKinetics_type), intent(inout):: plantMicKinetics_vars

    integer :: c, fc, j, tstep

    tstep = ttime%tstep

    do j = lbj, ubj
      do fc = 1, numf
        c = filter(fc)
        carbonflux_vars%cflx_input_litr_met_vr(c,j) = this%cflx_met_vr(tstep,j)
        carbonflux_vars%cflx_input_litr_cel_vr(c,j) = this%cflx_cel_vr(tstep,j)
        carbonflux_vars%cflx_input_litr_lig_vr(c,j) = this%cflx_lig_vr(tstep,j)
        carbonflux_vars%cflx_input_litr_cwd_vr(c,j) = this%cflx_cwd_vr(tstep,j)
        carbonflux_vars%rr_vr(c,j)                  = this%rr_vr(tstep,j)

        nitrogenflux_vars%nflx_input_litr_met_vr(c,j) = this%nflx_met_vr(tstep,j)
        nitrogenflux_vars%nflx_input_litr_cel_vr(c,j) = this%nflx_cel_vr(tstep,j)
        nitrogenflux_vars%nflx_input_litr_lig_vr(c,j) = this%nflx_lig_vr(tstep,j)
        nitrogenflux_vars%nflx_input_litr_cwd_vr(c,j) = this%nflx_cwd_vr(tstep,j)

        phosphorusflux_vars%pflx_input_litr_met_vr(c,j) = this%cflx_met_vr(tstep,j)/1600._r8
        phosphorusflux_vars%pflx_input_litr_cel_vr(c,j) = this%cflx_cel_vr(tstep,j)/2000._r8
        phosphorusflux_vars%pflx_input_litr_lig_vr(c,j) = this%cflx_lig_vr(tstep,j)/2500._r8
        phosphorusflux_vars%pflx_input_litr_cwd_vr(c,j) = this%cflx_cwd_vr(tstep,j)/4500._r8

        nitrogenflux_vars%nflx_minn_input_nh4_vr(c,j) = this%nflx_nh4_vr(tstep,j)
        nitrogenflux_vars%nflx_minn_input_no3_vr(c,j) = this%nflx_no3_vr(tstep,j)
        phosphorusflux_vars%pflx_minp_input_po4_vr(c,j) = this%pflx_po4_vr(tstep,j)

      enddo
    enddo
  end subroutine UpdateCNPForcing

  ! ----------------------------------------------------------------------

  function discharge(this, tstep) result(flux)
   !DESCRIPTION
   !return discharge at bottom of the soil column
    implicit none
    !arguments
    class(ForcingData_type) , intent(in) :: this
    integer                 , intent(in) :: tstep

    real(r8) :: flux
    integer  :: index

    if (this%forcing_type == steady_state) then
       index = 1
    else
       index = tstep
    end if

    flux = this%qbot(index)

  end function discharge

  ! ----------------------------------------------------------------------

  function infiltration(this, tstep) result(flux)
   !DESCRIPTION
   !return infiltration
   implicit none
    class(ForcingData_type) , intent(in) :: this
    integer                 , intent(in) :: tstep
    !temporary variables
    real(r8) :: flux
    integer  :: index

    if (this%forcing_type == steady_state) then
       index = 1
    else
       index = tstep
    end if

    flux = this%qflx_infl(index)

  end function infiltration

end module ForcingDataType
