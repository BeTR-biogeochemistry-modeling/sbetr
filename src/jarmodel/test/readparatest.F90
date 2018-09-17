program readparatest


  use ecacnpParaType, only : create_jarpars_ecacnp, ecacnp_para_type
  use cdomParaType, only : cdomPara_type, create_jarpars_cdom
  use ncdio_pio   , only : file_desc_t, ncd_io
  use BetrStatusType   , only : betr_status_type
  use ncdio_pio   , only : ncd_pio_closefile, ncd_pio_openfile, ncd_nowrite
implicit none

  class(ecacnp_para_type), pointer :: centpara
  class(cdomPara_type), pointer :: cdompara
  type(betr_status_type) :: bstatus
  type(file_desc_t)  :: ncid  ! pio netCDF file id
  character(len=*), parameter :: fname1='/Users/jinyuntang/work/github/ACME-Climate/sbetr/tools/jarmodel.ecacnp_pars.03072018.nc'
  character(len=*), parameter :: fname2='/Users/jinyuntang/work/github/ACME-Climate/sbetr/tools/jarmodel.cdomcnp_pars.03092018.nc'

  allocate(centpara, source=create_jarpars_ecacnp())
  allocate(cdompara, source=create_jarpars_cdom())

  call centpara%Init(bstatus=bstatus)
  call cdompara%Init(bstatus=bstatus)

  call ncd_pio_openfile(ncid, fname1, ncd_nowrite)

  call centpara%readPars(ncid, bstatus)
  if(bstatus%check_status())print*, bstatus%print_msg()
  call ncd_pio_closefile(ncid)

  call centpara%printPars()
  print*,'-----------------------------------------'
  call ncd_pio_openfile(ncid, fname2, ncd_nowrite)

  call cdompara%readPars(ncid, bstatus)
  if(bstatus%check_status())print*, bstatus%print_msg()
  call ncd_pio_closefile(ncid)

  call cdompara%printPars()


end program readparatest
