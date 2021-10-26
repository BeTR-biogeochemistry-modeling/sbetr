module BeTR_GridMod

  use bshr_kind_mod    , only : r8 => shr_kind_r8
  use babortutils      , only : endrun
  use bshr_log_mod     , only : errMsg => shr_log_errMsg
  use betr_constants   , only : betr_filename_length
  use betr_constants   , only : betr_string_length, betr_string_length_long
  use betr_constants   , only : betr_namelist_buffer_size
  use betr_constants   , only : stdout

  use betr_varcon, only : bspval


  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  character(len=*), parameter          :: uniform_str = 'uniform'
  character(len=*), parameter          :: clm_str = 'clm'
  character(len=*), parameter          :: dataset_str = 'dataset'

  integer, parameter :: uniform_grid = 1
  integer, parameter :: clm_grid = 2
  integer, parameter :: dataset_grid = 3

  integer, parameter :: num_columns = 1

  type, public :: betr_grid_type
     character(len=betr_filename_length) :: grid_data_filename
     character(len=betr_string_length)   :: grid_data_format ! file format: netcdf, namelist, csv, etc.
     character(len=betr_string_length)   :: grid_type_str ! uniform, clm
     integer :: grid_type

     integer           :: nlevgrnd
     real(r8)          :: delta_z
     real(r8), pointer :: zsoi(:)  => null() !soil depth, node center 1 : nlevsoi
     real(r8), pointer :: zisoi(:) => null() !soil depth, interface,  0 : nlevsoi
     real(r8), pointer :: dzsoi(:) => null() !soil layer thickness

     real(r8), pointer :: bsw(:) => null() ! clap-hornberg parameter
     real(r8), pointer :: watsat(:) => null() ! saturated volumetric water content
     real(r8), pointer :: sucsat(:)=> null()
     real(r8), pointer :: pctsand(:)=> null()
     real(r8), pointer :: pctclay(:)=> null()
     real(r8), pointer :: cellorg(:)=> null()
     real(r8), pointer :: msurf_OM(:) => null()
     real(r8), pointer :: KM_OM(:)  => null()
     real(r8), pointer :: bd(:)     => null()
     real(r8)          :: totzsoi
     integer           :: lithological_class
   contains
     procedure, public  :: Init
     procedure, public  :: ReadNamelist
     procedure, public  :: ReadNetCDFData
     procedure, public  :: UpdateGridConst
     procedure, private :: InitAllocate
     procedure, private :: uniform_vertical_grid
     procedure, private :: clm_exponential_vertical_grid
     procedure, private :: set_interface_depths
     procedure, private :: set_msurf_sorption
  end type betr_grid_type


contains

  ! ---------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer)

    implicit none

    class(betr_grid_type),                    intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in)    :: namelist_buffer

    this%nlevgrnd=15
    call this%ReadNameList(namelist_buffer)
    call this%InitAllocate()
    select case (trim(this%grid_type_str))
       case (clm_str)
          this%grid_type = clm_grid
          call this%clm_exponential_vertical_grid()
       case (uniform_str)
          this%grid_type = uniform_grid
          call this%uniform_vertical_grid()
       case (dataset_str)
          this%grid_type = dataset_grid
          ! grid will be initialized by reading the dataset
       case default
          this%grid_type = clm_grid
          call this%clm_exponential_vertical_grid()
          write(*, *) 'WARNING: no grid data type specified, using clm.'
    end select

    ! select read routine based on data format.
    call this%ReadNetCDFData()

    write(*, *) 'dump grid: '
    write(*, *) 'dzsoi   =  ', this%dzsoi
    write(*, *) 'zsoi    =  ', this%zsoi
    write(*, *) 'zisoi   =  ', this%zisoi
    write(*, *) 'bsw     =  ', this%bsw
    write(*, *) 'watsat  =  ', this%watsat
    write(*, *) 'sucsat  =  ', this%sucsat
    write(*, *) 'pctsand =  ', this%pctsand
    write(*, *) 'pctclay =  ', this%pctclay
    write(*, *) 'cellorg =  ', this%cellorg
  end subroutine Init

  ! ---------------------------------------------------------------------------
  subroutine InitAllocate(this)

    implicit none

    class(betr_grid_type), intent(inout) :: this

    allocate(this%zsoi(1:this%nlevgrnd))
    this%zsoi(:) = bspval

    allocate(this%zisoi(0:this%nlevgrnd))
    this%zisoi(:) = bspval

    allocate(this%dzsoi(1:this%nlevgrnd))
    this%dzsoi(:) = bspval

    allocate(this%bsw(1:this%nlevgrnd))
    this%bsw(:) = bspval

    allocate(this%watsat(1:this%nlevgrnd))
    this%watsat(:) = bspval

    allocate(this%bd(1:this%nlevgrnd))
    this%bd(:) = bspval

    allocate(this%sucsat(1:this%nlevgrnd))
    this%sucsat(:) = bspval

    allocate(this%pctsand(1:this%nlevgrnd))
    this%pctsand(:) = bspval

    allocate(this%pctclay(1:this%nlevgrnd))
    this%pctclay(:) = bspval

    allocate(this%cellorg(1:this%nlevgrnd))
    this%cellorg(:) = bspval

    allocate(this%msurf_OM(1:this%nlevgrnd))
    this%msurf_OM(:) = bspval

    allocate(this%KM_OM(1:this%nlevgrnd))
    this%KM_OM(:) = bspval

  end subroutine InitAllocate

  ! ---------------------------------------------------------------------------
  subroutine ReadNamelist(this, namelist_buffer)

    use betr_constants, only : betr_namelist_buffer_size
    use betr_constants, only : betr_filename_length
    use betr_constants, only : betr_string_length_long
    use elm_varpar,     only : nlevtrc_soil, nlevsoi, nlevgrnd
    implicit none

    class(betr_grid_type),                    intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in)    :: namelist_buffer

    integer                                :: nml_error
    character(len=*), parameter            :: subname = 'ReadNameList'
    character(len=betr_string_length)      :: grid_data_format, grid_type_str
    character(len=betr_filename_length)    :: grid_data_filename
    character(len=betr_string_length_long) :: ioerror_msg
    real(r8)                               :: delta_z

    !-----------------------------------------------------------------------

    namelist / betr_grid /                                    &
         grid_type_str, grid_data_filename, grid_data_format, &
         nlevgrnd, delta_z, nlevtrc_soil

    grid_data_format   = ''
    grid_data_filename = ''
    grid_type_str      = clm_str
    nlevgrnd           = 15 ! default to clm grid
    nlevtrc_soil       = 10
    delta_z            = bspval

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=betr_grid, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          write(*, *) ioerror_msg
          call endrun(msg="ERROR reading betr_grid namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr grid :'
       write(stdout, *)
       write(stdout, *) ' betr_grid namelist settings :'
       write(stdout, *)
       write(stdout, betr_grid)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    this%grid_type_str = trim(grid_type_str)
    this%grid_data_format = trim(grid_data_format)
    this%grid_data_filename = trim(grid_data_filename)
    this%nlevgrnd = nlevgrnd
    this%delta_z = delta_z
    nlevsoi            = nlevtrc_soil

  end subroutine ReadNamelist

  !------------------------------------------------------------------------
  subroutine ReadNetCDFData(this)
    !DESCRIPTION
    !read netcdf data
    use bncdio_pio, only : file_desc_t
    use bncdio_pio, only : ncd_nowrite
    use bncdio_pio, only : ncd_pio_openfile
    use bncdio_pio, only : ncd_getvar, Var_desc_t
    use bncdio_pio, only : ncd_pio_closefile, check_var
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variables
    character(len=250)    :: ncf_in_filename
    type(file_desc_t)     :: ncf_in
    real(r8), allocatable :: data(:,:)
    integer               :: j
    type(Var_desc_t)  :: vardesc
    logical :: readvar

    ncf_in_filename = trim(this%grid_data_filename)
    call ncd_pio_openfile(ncf_in, ncf_in_filename, mode=ncd_nowrite)

    allocate(data(num_columns, this%nlevgrnd))

    call ncd_getvar(ncf_in, 'WATSAT', data)
    do j = 1, this%nlevgrnd
       this%watsat(j) = data(num_columns, j)
       this%bd(j)     = (1._r8 - this%watsat(j))*2.7e3_r8  !kg /m3
    enddo

    call ncd_getvar(ncf_in, 'BSW', data)
    do j = 1, this%nlevgrnd
       this%bsw(j) = data(num_columns, j)
    enddo

    call ncd_getvar(ncf_in, 'SUCSAT', data)
    do j = 1, this%nlevgrnd
       this%sucsat(j) = data(num_columns, j)
    enddo

    call ncd_getvar(ncf_in, 'PCTSAND', data)
    do j = 1, this%nlevgrnd
      this%pctsand(j) = data(num_columns, j)
    enddo

    call ncd_getvar(ncf_in, 'PCTCLAY', data)
    do j = 1, this%nlevgrnd
      this%pctclay(j) = data(num_columns, j)
    enddo

    call ncd_getvar(ncf_in, 'CELLORG', data)
    do j = 1, this%nlevgrnd
       this%cellorg(j) = data(num_columns, j)
    enddo

    if (this%grid_type == dataset_grid) then
       call ncd_getvar(ncf_in, 'DZSOI', data)
       this%totzsoi=0._r8
       do j = 1, this%nlevgrnd
          this%dzsoi(j) = data(num_columns, j)
          this%totzsoi = this%totzsoi + this%dzsoi(j)
       end do

       call ncd_getvar(ncf_in, 'ZSOI', data)
       do j = 1, this%nlevgrnd
          this%zsoi(j) = data(num_columns, j)
       end do

       call this%set_interface_depths()

    end if
    call check_var(ncf_in, 'lithological_class', vardesc, readvar)
    if(readvar)then
      call ncd_getvar(ncf_in,'lithological_class',this%lithological_class)
    else
      this%lithological_class=3   !this is a randomly assigned value
    endif
    call ncd_pio_closefile(ncf_in)

    deallocate(data)

    call this%set_msurf_sorption()
  end subroutine ReadNetCDFData

  ! ---------------------------------------------------------------------------
  subroutine uniform_vertical_grid(this)
    !DESCRIPTION
    !set uniform thickness grid
    use bshr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variable
    integer :: j

    if (abs(this%delta_z-bspval)<1.e-10_r8) then
       call endrun(msg="ERROR reading betr_grid namelist must specify delta_z. "//errmsg(mod_filename, __LINE__))
    end if

    ! thickness b/n two interfaces
    this%dzsoi(:) = this%delta_z

    ! node depths
    this%totzsoi = 0._r8
    do j = 1, this%nlevgrnd
       this%zsoi(j) = (real(j, r8) - 0.5_r8) * this%dzsoi(j)
       this%totzsoi = this%totzsoi + this%dzsoi(j)
    enddo

    call this%set_interface_depths()

  end subroutine uniform_vertical_grid

  ! ---------------------------------------------------------------------------
  subroutine clm_exponential_vertical_grid(this)
    !
    ! DESCRIPTION
    ! initialize the clm exporential vertical grid for computation
    !
    use bshr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variables
    real(r8)                             :: scalez = 0.025_r8 ! Soil layer thickness discretization (m)
    integer                              :: j

    ! node depths
    do j = 1, this%nlevgrnd
       this%zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)
    enddo

    ! thickness b/n two interfaces
    this%dzsoi(1) = 0.5_r8*(this%zsoi(1) + this%zsoi(2))

    do j = 2,this%nlevgrnd-1
       this%dzsoi(j) = 0.5_r8*(this%zsoi(j+1) - this%zsoi(j-1))

    enddo
    this%dzsoi(this%nlevgrnd) = this%zsoi(this%nlevgrnd) - this%zsoi(this%nlevgrnd-1)
    this%totzsoi = 0._r8
    do j = 1, this%nlevgrnd
      this%totzsoi = this%totzsoi + this%dzsoi(j)
    enddo
    call this%set_interface_depths()

  end subroutine clm_exponential_vertical_grid

  ! ---------------------------------------------------------------------------
  subroutine set_interface_depths(this)
    !DESCRIPTION
    !set node depth
    !USES
    use bshr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variables
    integer :: j

    this%zisoi(0) = 0._r8
    do j = 1, this%nlevgrnd-1
       this%zisoi(j) = 0.5_r8*(this%zsoi(j) + this%zsoi(j+1))
    enddo
    this%zisoi(this%nlevgrnd) = this%zsoi(this%nlevgrnd) + 0.5_r8*this%dzsoi(this%nlevgrnd)
    !print*,'zisoi',this%nlevgrnd,size(this%zisoi)
  end subroutine set_interface_depths

  ! ---------------------------------------------------------------------------
  subroutine UpdateGridConst(this, bounds, lbj, ubj, numf, filter, soilstate_vars, cnstate_vars)
  !
  !DESCRIPTION
  !update grid constant variables
  use SoilStateType, only : soilstate_type
  use decompMod         , only : bounds_type
  use CNStateType       , only : cnstate_type
  implicit none
  class(betr_grid_type), intent(inout) :: this
  type(bounds_type)        , intent(in)    :: bounds
  integer                  , intent(in)    :: numf
  integer                  , intent(in)    :: filter(:)
  integer                  , intent(in)    :: lbj, ubj
  type(soilstate_type)     , intent(inout) :: soilstate_vars
  type(cnstate_type)       , intent(inout) :: cnstate_vars
  integer :: c,fc, j

  do j = 1, ubj
    do fc = 1, numf
      c = filter(fc)
      soilstate_vars%watsat_col(c,j)   = this%watsat(j)
      soilstate_vars%bsw_col(c,j)      = this%bsw(j)
      soilstate_vars%sucsat_col(c,j)   = this%sucsat(j)
      soilstate_vars%cellsand_col(c,j) = this%pctsand(j)
      soilstate_vars%cellclay_col(c,j) = this%pctclay(j)
      soilstate_vars%cellorg_col(c,j)  = this%cellorg(j)
      soilstate_vars%bd_col(c,j)       = this%bd(j)
      cnstate_vars%pdep_prof_col(c,j)  = this%dzsoi(j)/this%totzsoi
    enddo
  enddo
  do fc = 1, numf
    c = filter(fc)
    cnstate_vars%lithoclass_col(c) = this%lithological_class
  enddo
  end subroutine UpdateGridConst
  ! ---------------------------------------------------------------------------
  subroutine set_msurf_sorption(this)
  use tracer_varcon, only : catomw
  implicit none
  class(betr_grid_type), intent(inout) :: this
  integer :: j
  do j = 1, this%nlevgrnd
    this%msurf_OM(j) = 1300._r8/catomw
    this%KM_OM(j) = 10._r8/catomw
  enddo
  end subroutine set_msurf_sorption

end module BeTR_GridMod
