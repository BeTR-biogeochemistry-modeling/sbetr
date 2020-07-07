
program main


  use histMod, only : histf_type
  use histMod, only : hist_var_str_len, hist_unit_str_len, hist_freq_str_len
  use anaforc, only : anal_tempf
  use shr_kind_mod, only : r8 => shr_kind_r8
  use  BeTR_TimeMod, only : betr_time_type
  use betr_varcon, only : var_flux_type, var_state_type
implicit none
  character(len=*), parameter :: modfile=__FILE__

  character(len=hist_var_str_len) :: varl(3)
  character(len=hist_unit_str_len) :: unitl(3)
  character(len=hist_freq_str_len) :: freql(3)
  integer :: vartype(3)
  real(r8) :: yval(1:1,1:3)
  integer :: ncols
  type(histf_type) :: hist
  type(betr_time_type) :: timer
  real(r8) :: dt, tsec
  print*,'netcdf io test in '//modfile

  varl(1)='tmp1'; unitl(1)='K'; freql(1)='day'; vartype(1)=var_state_type
  varl(2)='tmp2'; unitl(2)='K'; freql(2)='month'; vartype(2)=var_state_type
  varl(3)='tmp3'; unitl(3)='K'; freql(3)='hour'; vartype(3)=var_state_type
  ncols=1
  call timer%Init(namelist_buffer='junk data')
  dt=timer%get_step_size()
  call hist%init(ncols,varl, unitl, vartype, freql, 'iotest')

  do
    tsec=timer%get_cur_time()+dt*0.5_r8
    call anal_tempf(yval(1,1), tsec)
    yval(1,2)=yval(1,1)
    yval(1,3)=yval(1,1)
    call timer%update_time_stamp()
    call hist%hist_wrap(yval, timer)
    if(timer%its_time_to_exit())exit
  enddo


end program main
