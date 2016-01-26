module CLMForcType
  use shr_kind_mod        , only : r8 => shr_kind_r8
implicit none


  type, public :: clmforc_type
    real(r8), pointer :: t_soi(:,:)
    real(r8), pointer :: h2o_soi(:,:)
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


  end subroutine InitAllocate


  !------------------------------------------------------------------------
  subroutine loadForc(this)
  use ncdio_pio
  class(clmforc_type) :: this
  character(len=250) :: ncf_in_filename
  type(file_desc_t)  :: ncf_in
  real(r8), allocatable :: data_2d(:,:,:,:)
  integer :: dimlenl, dimlent
  integer :: j1, j2

  ncf_in_filename='/Users/jinyuntang/work/data_collection/clm_output/hourly_data_example/sierra_clmdef_hhour.clm2.h0.0021-01-01-00000.nc'

  call ncd_pio_openfile(ncf_in, ncf_in_filename, mode=ncd_nowrite)
  dimlenl = get_dim_len(ncf_in,'levgrnd')
  dimlent = get_dim_len(ncf_in,'time')

  !allocate memory
  call this%Init( dimlenl, dimlent)

  !read data
  allocate(data_2d(1,1,dimlenl,1:dimlent))
  call ncd_getvar(ncf_in, 'TSOI', data_2d)

  !now assign the data
  do j2 = 1, dimlenl
  do j1 = 1, dimlent
    this%t_soi(j1,j2) = data_2d(1,1,j2,j1)
  enddo
  enddo


  call ncd_getvar(ncf_in, 'H2OSOI', data_2d)

  !now assign the data
  do j2 = 1, dimlenl
  do j1 = 1, dimlent
    this%h2o_soi(j1,j2) = data_2d(1,1,j2,j1)
  enddo
  enddo

  call ncd_pio_closefile(ncf_in)




  deallocate(data_2d)
  end subroutine loadForc

end module CLMForcType
