program main

  use  BeTR_TimeMod, only : betr_time_type
implicit none


  type(betr_time_type) :: timer

  integer :: counter

  call timer%Init(namelist_buffer='junk data')


  counter = 0
  do
    counter = counter + 1
    call timer%update_time_stamp
!    call timer%print_cur_time

!    if(mod(counter, 48)==0)exit
!    if(mod(counter, 48*7)==0)exit
!    if(mod(counter, 48*31)==0)exit
!    if(mod(counter, 48*365)==0)exit
    if(mod(counter, 48*365*2)==0)exit
  enddo
  print*,'its_a_new_day',timer%its_a_new_day()
  print*,'its_a_new_week',timer%its_a_new_week()
  print*,'its_a_new_month',timer%its_a_new_month()
  print*,'its_a_new_year',timer%its_a_new_year()

end program main
