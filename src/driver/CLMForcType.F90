module CLMForcType
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_const_mod       , only : rhoice => SHR_CONST_RHOICE
  use shr_const_mod       , only : rhoh2o => SHR_CONST_RHOFW

  implicit none
  private

  character(len=*), private, parameter :: module_filename = '__FILE__'
  
  type, public :: clmforc_type
     character(len=256) :: grid_filename
     character(len=256) :: forcing_filename
    real(r8), pointer :: t_soi(:,:)
    real(r8), pointer :: h2osoi_liqvol(:,:)
    real(r8), pointer :: h2osoi_liq(:,:)
    real(r8), pointer :: h2osoi_ice(:,:)
    real(r8), pointer :: qflx_infl(:)       !surface infiltration, mm/s
    real(r8), pointer :: qflx_tran_dep(:,:) !transpiration at depth, mm/s
    real(r8), pointer :: pbot(:)            !amtospheric pressure, Pa
    real(r8), pointer :: dzsoi(:)           !node thickness
    real(r8), pointer :: zsoi(:)            !node depth of each numerical node
    real(r8), pointer :: bsw(:)             !clap-hornberg parameter
    real(r8), pointer :: watsat(:)          !saturated volumetric water content
    real(r8), pointer :: h2osoi_icevol(:,:)
    integer           :: nlev
    integer           :: ntsnap
  contains
     procedure, private  :: Init
     procedure, private :: InitAllocate
     procedure, public  :: LoadForcingData
     procedure, private :: ReadNameList
  end type clmforc_type

  type(clmforc_type), public :: clmforc_vars

contains
  subroutine Init(this, diml, dimt)

  class(clmforc_type) :: this

  integer, intent(in) :: diml, dimt

  call this%InitAllocate(diml, dimt)

  end subroutine Init


  !------------------------------------------------------------------------
  subroutine InitAllocate(this, diml, dimt)

  class(clmforc_type) :: this
  !at this moment the variable size is fixed

  integer, intent(in) :: diml, dimt

  allocate(this%t_soi(dimt,diml))
  allocate(this%h2osoi_liqvol(dimt,diml))
  allocate(this%h2osoi_icevol(dimt,diml))
  allocate(this%h2osoi_liq(dimt,diml))
  allocate(this%h2osoi_ice(dimt,diml))
  allocate(this%qflx_infl(dimt))
  allocate(this%qflx_tran_dep(dimt,diml))
  allocate(this%zsoi(diml))
  allocate(this%dzsoi(diml))
  allocate(this%pbot(dimt))
  allocate(this%watsat(diml))
  allocate(this%bsw(diml))
  this%nlev=diml
  this%ntsnap=dimt
  end subroutine InitAllocate


  !------------------------------------------------------------------------
  subroutine LoadForcingData(this)
  use ncdio_pio
  class(clmforc_type) :: this
  character(len=250) :: ncf_in_filename_grid
  character(len=250) :: ncf_in_filename_forc
  type(file_desc_t)  :: ncf_in_grid
  type(file_desc_t)  :: ncf_in_forc
  real(r8), allocatable :: data_2d(:,:,:,:)
  real(r8), allocatable :: data_1d(:,:,:)
  real(r8), allocatable :: data1(:,:)
  integer :: dimlenl, dimlent
  integer :: j1, j2

  call this%ReadNameList('example.nl')
  ncf_in_filename_grid=trim(this%grid_filename)
  call ncd_pio_openfile(ncf_in_grid, ncf_in_filename_grid, mode=ncd_nowrite)

  ncf_in_filename_forc=trim(this%forcing_filename)

  call ncd_pio_openfile(ncf_in_forc, ncf_in_filename_forc, mode=ncd_nowrite)
  dimlenl = get_dim_len(ncf_in_forc,'levgrnd')
  dimlent = get_dim_len(ncf_in_forc,'time')

  !allocate memory
  call this%Init( dimlenl, dimlent)

  !load grid
  allocate(data1(1,dimlenl))
  call ncd_getvar(ncf_in_grid,'DZSOI',data1)
  do j2 = 1, dimlenl
    this%dzsoi(j2) = data1(1,j2)
  enddo

  call ncd_getvar(ncf_in_grid,'ZSOI',data1)
  do j2 = 1, dimlenl
    this%zsoi(j2) = data1(1,j2)
  enddo

  call ncd_getvar(ncf_in_grid,'WATSAT',data1)
  do j2 = 1, dimlenl
    this%watsat(j2) = data1(1,j2)
  enddo

  call ncd_getvar(ncf_in_grid,'BSW',data1)
  do j2 = 1, dimlenl
    this%bsw(j2) = data1(1,j2)
  enddo

  !read data
  allocate(data_2d(1,1,dimlenl,1:dimlent))
  allocate(data_1d(1,1,dimlent))
  call ncd_getvar(ncf_in_forc, 'TSOI', data_2d)

  !now assign the data

  do j2 = 1, dimlenl
  do j1 = 1, dimlent
    this%t_soi(j1,j2) = data_2d(1,1,j2,j1)
  enddo
  enddo


  !read h2osoi
  call ncd_getvar(ncf_in_forc, 'H2OSOI', data_2d)

  !now assign the data
  do j2 = 1, dimlenl
  do j1 = 1, dimlent
    this%h2osoi_liqvol(j1,j2) = data_2d(1,1,j2,j1)
  enddo
  enddo

  !read ice
  call ncd_getvar(ncf_in_forc, 'SOILICE', data_2d)

  !now assign the data
  do j2 = 1, dimlenl
  do j1 = 1, dimlent

    this%h2osoi_icevol(j1,j2) = data_2d(1,1,j2,j1)/rhoice/this%dzsoi(j2)
    this%h2osoi_liqvol(j1,j2) = this%h2osoi_liqvol(j1,j2) - this%h2osoi_icevol(j1,j2)
    this%h2osoi_ice(j1,j2) = data_2d(1,1,j2,j1)
    this%h2osoi_liq(j1,j2) = this%h2osoi_liqvol(j1,j2)*this%dzsoi(j2)*rhoh2o
  enddo
  enddo


  call ncd_getvar(ncf_in_forc,'QINFL',data_1d)

  do j1 =1, dimlent
    this%qflx_infl(j1) = data_1d(1,1,j1)
  enddo

  call ncd_getvar(ncf_in_forc,'PBOT',data_1d)

  do j1 =1, dimlent
    this%pbot(j1) = data_1d(1,1,j1)
  enddo

  call ncd_getvar(ncf_in_forc, 'QVEGT_DEP', data_2d)

  !now assign the data
  do j2 = 1, dimlenl
  do j1 = 1, dimlent
    this%qflx_tran_dep(j1,j2) = data_2d(1,1,j2,j1)
  enddo
  enddo

  call ncd_pio_closefile(ncf_in_forc)
  call ncd_pio_closefile(ncf_in_grid)

  deallocate(data_2d)
  deallocate(data_1d)
  deallocate(data1)


end subroutine LoadForcingData

! ----------------------------------------------------------------------

  subroutine ReadNameList(this, filename)
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
    
    implicit none
    ! !ARGUMENTS:
    class(clmforc_type) :: this
    character(len=*), intent(in) :: filename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                    ! error code
    integer :: unitn                   ! unit for namelist file
    character(len=32) :: subname = 'ReadNameList'
    character(len=256) :: grid_filename, forcing_filename

    !-----------------------------------------------------------------------

    namelist / forcing_inparm / grid_filename, forcing_filename

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in forcing_inparm  namelist'
       call opnfil (filename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'forcing_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, forcing_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading forcing_inparm namelist"//errmsg(module_filename, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if
    ! Broadcast namelist variables read in
    this%grid_filename = grid_filename
    this%forcing_filename = forcing_filename
    call shr_mpi_bcast(grid_filename, mpicom)
    call shr_mpi_bcast(forcing_filename, mpicom)

  end subroutine ReadNameList



end module CLMForcType
