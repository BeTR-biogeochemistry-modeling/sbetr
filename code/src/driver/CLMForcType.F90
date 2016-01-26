module CLMForcType
  use shr_kind_mod        , only : r8 => shr_kind_r8
implicit none


  type, public :: clmforc_type
    real(r8), pointer :: t_soi(:,:)
    real(r8), pointer :: h2o_soi(:,:)
    real(r8), pointer :: qflx_infl(:)
    real(r8), pointer :: qflx_tran_dep(:,:)
    real(r8), pointer :: dzsoi(:)     !node thickness
    real(r8), pointer :: zsoi(:)      !node depth of each numerical node
    integer           :: nlev
    integer           :: ntsnap
  contains
     procedure, private  :: Init
     procedure, private :: InitAllocate
     procedure, public  :: loadForc
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
  allocate(this%h2o_soi(dimt,diml))
  allocate(this%qflx_infl(dimt))
  allocate(this%qflx_tran_dep(dimt,diml))
  allocate(this%zsoi(diml))
  allocate(this%dzsoi(diml))

  this%nlev=diml
  this%ntsnap=dimt
  end subroutine InitAllocate


  !------------------------------------------------------------------------
  subroutine loadForc(this)
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


  ncf_in_filename_grid='/Users/jinyuntang/work/data_collection/clm_output/hourly_data_example/sierra_grid.nc'
  call ncd_pio_openfile(ncf_in_grid, ncf_in_filename_grid, mode=ncd_nowrite)


  ncf_in_filename_forc='/Users/jinyuntang/work/data_collection/clm_output/hourly_data_example/sierra_clmdef_hhour.clm2.h0.0021-01-01-00000.nc'

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


  call ncd_getvar(ncf_in_forc, 'H2OSOI', data_2d)

  !now assign the data
  do j2 = 1, dimlenl
  do j1 = 1, dimlent
    this%h2o_soi(j1,j2) = data_2d(1,1,j2,j1)
  enddo
  enddo

  call ncd_getvar(ncf_in_forc,'QINFL',data_1d)

  do j1 =1, dimlent
    this%qflx_infl(j1) = data_1d(1,1,j1)
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


  end subroutine loadForc

end module CLMForcType
