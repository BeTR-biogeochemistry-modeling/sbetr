module BeTR_TimeMod
!
! DESCRIPTION
! the module contains subroutine to march the time

  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use bshr_kind_mod  , only : i8 => shr_kind_i8
  use betr_constants , only : stdout, betr_string_length_long

  implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  type, public:: betr_time_type
     ! NOTE(bja, 201603) all real variables have units of seconds!
     real(r8) :: delta_time
     real(r8) :: stop_time
     real(r8) :: time
     real(r8) :: time0
     real(r8) :: timef
     real(r8) :: toy
     real(r8) :: restart_dtime
     integer  :: tstep
     integer  :: nelapstep
     integer  :: dow, dom, doy
     integer  :: moy, cyears, cdays
     real(r8) :: tod
     integer  :: hist_freq   !negative number, steps, positive number, 1: day, 30:mon, 365:year
     integer  :: stop_opt
   contains
     procedure, public :: Init
     procedure, public :: its_time_to_write_restart
     procedure, public :: its_time_to_exit
     procedure, public :: update_time_stamp
     procedure, public :: set_nstep
     procedure, public :: get_nstep
     procedure, public :: set_time_offset
     procedure, public :: get_days_per_year
     procedure, public :: get_step_size
     procedure, public :: get_cur_time
     procedure, public :: get_cur_timef
     procedure, private:: proc_nextstep
     procedure, public :: proc_initstep
     procedure, public :: print_cur_time
     procedure, public :: its_time_to_histflush
     procedure, private:: ReadNamelist
     procedure, public :: setClock
     procedure, public :: its_a_new_hour
     procedure, public :: its_a_new_day
     procedure, public :: its_a_new_week
     procedure, public :: its_a_new_month
     procedure, public :: its_a_new_year
     procedure, public :: get_ymdhs
     procedure, public :: get_cur_year
     procedure, public :: get_cur_day
     procedure, public :: is_first_step
     procedure, public :: print_model_time_stamp
  end type betr_time_type

  integer, parameter, private :: daz(12)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
contains

  subroutine setClock(this, dtime, nelapstep)

  implicit none
  class(betr_time_type), intent(inout) :: this
  real(r8), intent(in) :: dtime
  integer, intent(in) :: nelapstep

  this%delta_time = dtime
  this%nelapstep = nelapstep
  end subroutine setClock
  !-------------------------------------------------------------------------------
  function get_cur_time(this)result(ans)

  implicit none
  class(betr_time_type), intent(in) :: this
  real(r8) :: ans

  ans = this%time
  end function get_cur_time
  !-------------------------------------------------------------------------------
  function get_cur_timef(this)result(ans)

  implicit none
  class(betr_time_type), intent(in) :: this
  real(r8) :: ans

  ans = this%time+this%time0
  end function get_cur_timef
  !-------------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, masterproc)

    use betr_constants, only : betr_namelist_buffer_size

    implicit none

    class(betr_time_type), intent(inout) :: this
    character(len=*), optional, intent(in) :: namelist_buffer
    logical, optional, intent(in) :: masterproc
    this%tstep = 1
    this%time0 = 0._r8
    this%time  = 0._r8
    this%tod   = 0._r8
    this%toy   = 0._r8
    this%cyears = 0
    this%cdays  = 0
    this%dow    = 0
    this%dom    = 0
    this%doy    = 0
    this%moy    = 1
    this%hist_freq=-1
    this%nelapstep=0
    if(present(namelist_buffer))then
      if(present(masterproc))then
        call this%ReadNamelist(namelist_buffer, masterproc)
      else
        call this%ReadNamelist(namelist_buffer)
      endif
    endif
  end subroutine Init

  ! ----------------------------------------------------------------------

  subroutine ReadNamelist(this, namelist_buffer, masterproc)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use babortutils    , only : endrun
    use bshr_log_mod   , only : errMsg => shr_log_errMsg
    use betr_constants , only : stdout, betr_string_length_long, betr_namelist_buffer_size

    implicit none
    ! !ARGUMENTS:
    class(betr_time_type), intent(inout) :: this
    character(len=*), intent(in) :: namelist_buffer
    logical, optional, intent(in) :: masterproc
    !
    ! !LOCAL VARIABLES:
    integer :: nml_error
    character(len=*), parameter :: subname = 'betr_time%ReadNamelist'
    character(len=betr_string_length_long) :: ioerror_msg
    real(r8) :: delta_time         !model time step
    real(r8) :: stop_time          !when to stop
    integer  :: stop_n
    integer  :: hist_freq
    character(len=8)  :: stop_option        !1: step, 2:day, 3:year
    real(r8) :: restart_dtime      !when to write restart file
    logical :: masterproc_loc
    !-----------------------------------------------------------------------

    namelist / betr_time / delta_time, stop_n, stop_option, restart_dtime, hist_freq

    ! FIXME(bja, 201603) Only reading time variables in seconds!
    ! Should enable other values with unit coversions.

    ! FIXME(bja, 201603) assign some defaults, should eventually remove
    ! when all input files are updated.
    masterproc_loc=.true.
    if(present(masterproc))masterproc_loc=masterproc
    delta_time = 1800._r8                !half hourly time step
    stop_n=2       !by default 2 cycle
    stop_option='nyears'  !by default years
    this%stop_opt=1
    hist_freq=-1   !write every time step
    restart_dtime = -1._r8

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( index(trim(namelist_buffer),'betr_time')/=0 )then
       ioerror_msg=''
       read(namelist_buffer, nml=betr_time, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading sbetr_driver namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (masterproc_loc) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr time :'
       write(stdout, *)
       write(stdout, *) ' betr_time namelist settings :'
       write(stdout, *)
       write(stdout, betr_time)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif
    this%hist_freq = hist_freq
    this%delta_time = delta_time
    this%stop_time = delta_time*stop_n
    select case (trim(stop_option))
    case ('ndays')
      !day
      this%stop_time= stop_n * 86400._r8
      this%stop_opt=2
    case ('nyears')
      !year
      this%stop_time= stop_n * 86400._r8 * 365._r8
      this%stop_opt=3
    case ('nsteps')
      this%stop_time = stop_n * this%delta_time
      this%stop_opt=1
    case default
      call endrun(msg="ERROR setting up stop_option "//errmsg(mod_filename, __LINE__))
    end select
    if(restart_dtime < 0._r8)then
      this%restart_dtime  = this%stop_time
    else
      this%restart_dtime = restart_dtime
    endif
  end subroutine ReadNamelist

  !-------------------------------------------------------------------------------
    function its_time_to_write_restart(this)result(ans)
     !
     ! DESCRIPTION
     ! decide if to write restart file
     !

     implicit none
     class(betr_time_type), intent(in) :: this
     logical :: ans

     character(len=80) :: subname = 'its_time_to_write_restart'

     ans = (mod(this%time,this%restart_dtime) == 0)
     end function its_time_to_write_restart

  !-------------------------------------------------------------------------------
  function its_time_to_exit(this) result(ans)
    !
    ! DESCRIPTION
    ! decide if to exit the loop
    !

    implicit none
    class(betr_time_type), intent(in) :: this
    logical :: ans

    character(len=80) :: subname = 'its_time_to_exit'

    ans = (this%time >= this%stop_time)

  end function its_time_to_exit

  !-------------------------------------------------------------------------------
  subroutine update_time_stamp(this)
    !
    ! DESCRIPTION
    !

    implicit none

    class(betr_time_type), intent(inout) :: this

    character(len=80) :: subname='update_time_stamp'

    real(r8), parameter :: secpyear= 86400._r8*365._r8

    this%time = this%time + this%delta_time
    this%toy  = this%toy + this%delta_time

    this%tstep = this%tstep + 1
    !
    ! reset the clock every year, and assuming the time step
    ! size is always
    if(mod(this%toy, secpyear) == 0) then
       this%tstep = 1
    end if

    !update time of the day
    this%tod=this%tod+this%delta_time

    !update varaibles when it is to start a new day
    if(this%its_a_new_day())then
      this%tod=0._r8
      this%dom=this%dom+1
      this%dow = mod(this%dow + 1, 7)
      this%doy = this%doy + 1
      this%cdays= this%cdays + 1
    endif

    if(this%its_a_new_month())then
      this%moy=this%moy+1
      if(this%moy==13)then
        this%moy=1
        this%doy=mod(this%doy,365)
      endif
      this%dom=0
    endif
    if(this%its_a_new_year())this%cyears=this%cyears+1
    call this%proc_nextstep()
  end subroutine update_time_stamp

  !-------------------------------------------------------------------------------
  real(r8) function get_step_size(this) result(rc)

    ! Return the step size in seconds.

    implicit none

    class(betr_time_type), intent(in) :: this

    rc = this%delta_time
    return

  end function get_step_size

  !-------------------------------------------------------------------------------
  integer function get_nstep(this)

    ! Return the timestep number.
    implicit none

    class(betr_time_type), intent(in) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'

    get_nstep = this%nelapstep

  end function get_nstep


  !-------------------------------------------------------------------------------
  subroutine set_nstep(this, nstep)

    ! Return the timestep number.
    implicit none
    class(betr_time_type), intent(inout) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'
    integer, intent(in) :: nstep

    this%nelapstep = nstep

    if(this%its_a_new_year())then
      this%tstep = 1
    endif

  end subroutine set_nstep


  !-------------------------------------------------------------------------------
  subroutine set_time_offset(this, nstep, continue_run)

    ! Return the timestep number.
    implicit none
    class(betr_time_type), intent(inout) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'
    integer, intent(in) :: nstep
    logical, intent(in) :: continue_run

    this%time0  = nstep*this%delta_time
    if(continue_run)then
      this%nelapstep = nstep
    else
      this%nelapstep = 0
    endif
    if(this%its_a_new_year())then
      this%tstep = 1
    endif
  end subroutine set_time_offset

  !-------------------------------------------------------------------------------
  subroutine proc_nextstep(this)

    implicit none

    class(betr_time_type), intent(inout) :: this

    this%nelapstep = this%nelapstep + 1
  end subroutine proc_nextstep

  !-------------------------------------------------------------------------------
  subroutine proc_initstep(this)

    implicit none
    class(betr_time_type), intent(inout) :: this

    this%nelapstep = 0
  end subroutine proc_initstep

  !-------------------------------------------------------------------------------
  integer function get_days_per_year(this, offset)

    implicit none

    class(betr_time_type), intent(in) :: this
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.

    ! Positive for future times, negative
    ! for previous times.

    ! remove unused dummy arg compiler warning
    if (offset > 0) continue
    if (this%tstep > 0) continue

    !hardwire at the moment
    get_days_per_year = 365
  end function get_days_per_year

  !-------------------------------------------------------------------------------

  subroutine print_cur_time(this)
  implicit none
  class(betr_time_type), intent(in) :: this

  print*,'time=',this%time
  print*,'tod =',this%tod
  print*,'dow =',this%dow
  print*,'dom =',this%dom
  print*,'doy =',this%doy
  print*,'moy =',this%moy
  print*,'cyears=',this%cyears

  end subroutine print_cur_time


  !-------------------------------------------------------------------------------
  subroutine get_ymdhs(this, ymdhs)

  implicit none
  class(betr_time_type), intent(in) :: this
  character(len=*), intent(out) :: ymdhs

  write(ymdhs,'(I4.4,I2.2,I2.2,I2.2,I4.4)')this%cyears,this%moy,this%dom,&
     int(this%tod/3600.0),int(mod(this%tod,3600._r8))
  end subroutine get_ymdhs
  !-------------------------------------------------------------------------------
  function its_a_new_hour(this)result(yesno)
  implicit none
  class(betr_time_type), intent(in) :: this
  logical :: yesno


  yesno= abs(mod(this%tod, 3600._r8))<1.e-3_r8

  end function its_a_new_hour

  !-------------------------------------------------------------------------------

  function its_a_new_day(this)result(yesno)
  implicit none
  class(betr_time_type), intent(in) :: this
  logical :: yesno


  yesno= abs(mod(this%tod, 86400._r8))<1.e-3_r8

  end function its_a_new_day

  !-------------------------------------------------------------------------------

  function its_a_new_week(this)result(yesno)
  implicit none
  class(betr_time_type), intent(in) :: this
  logical :: yesno


  yesno= (this%dow == 0 .and. this%tod<1.e-3_r8)

  end function its_a_new_week

  !-------------------------------------------------------------------------------

  function its_a_new_month(this)result(yesno)
  implicit none
  class(betr_time_type), intent(in) :: this
  logical :: yesno


  yesno = ((this%dom == daz(this%moy) .or. this%dom==0) .and. this%tod < 1.e-3_r8)

  end function its_a_new_month

  !-------------------------------------------------------------------------------
  function its_a_new_year(this)result(yesno)
  implicit none
  class(betr_time_type), intent(in) :: this
  logical :: yesno

  yesno = (this%moy == 1 .and. this%dom == 0 .and. this%tod < 1.e-3_r8)

  end function its_a_new_year
  !-------------------------------------------------------------------------------

  function its_time_to_histflush(this)result(yesno)

  implicit none
  class(betr_time_type), intent(in) :: this
  logical :: yesno

  if(this%hist_freq<0)then
    !by time steps
    yesno = (mod(this%nelapstep,this%hist_freq)==0)
  elseif(this%hist_freq==0)then
    !no history file output until the last time step
    yesno=this%its_time_to_exit()
  elseif(this%hist_freq==1)then
    !by day
    yesno=this%its_a_new_day()
  elseif(this%hist_freq==30)then
    !by month
    yesno=this%its_a_new_month()
  elseif(this%hist_freq==365)then
    !by year
    yesno=this%its_a_new_year()
  elseif(this%hist_freq==9999)then
    yesno=.false.
  endif
  end function its_time_to_histflush
  !-------------------------------------------------------------------------------
  function get_cur_year(this)result(ans)
  implicit none
  class(betr_time_type), intent(in) :: this
  integer :: ans
  ans = this%cyears
  end function get_cur_year
  !-------------------------------------------------------------------------------
  function get_cur_day(this)result(ans)
  implicit none
  class(betr_time_type), intent(in) :: this
  integer :: ans
  ans = this%cdays
  end function get_cur_day
  !-------------------------------------------------------------------------------
  function is_first_step(this)result(ans)
  implicit none
  class(betr_time_type), intent(in) :: this
  logical :: ans
  ans = (this%nelapstep==0)
  end function is_first_step

  !-------------------------------------------------------------------------------
  subroutine print_model_time_stamp(this, iulog)

  implicit none
  class(betr_time_type), intent(in) :: this
  integer, intent(in) :: iulog

  select case(this%stop_opt)
  case (1)
      write(iulog,*)'step', this%get_nstep()
  case (2)
    if (this%its_a_new_day()) then
      write(iulog,*)'day', this%get_cur_day()
    end if
  case (3)
    if (this%its_a_new_year()) then
      write(iulog,*)'year', this%get_cur_year()
    end if
  end select

  end subroutine print_model_time_stamp
  !-------------------------------------------------------------------------------
end module BeTR_TimeMod
