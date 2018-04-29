module histMod

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8, i4 => shr_kind_i4
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use ncdio_pio   , only : file_desc_t
implicit none
 private

  character(len=*), parameter :: mod_filename=__FILE__
  integer, parameter :: nclocks=5
  integer, parameter :: clock_hour =1
  integer, parameter :: clock_day  =2
  integer, parameter :: clock_week =3
  integer, parameter :: clock_month=4
  integer, parameter :: clock_year =5

  integer, public, parameter :: hist_var_str_len=40
  integer, public, parameter :: hist_unit_str_len=16
  integer, public, parameter :: hist_freq_str_len=6
  integer, private, parameter:: hist_fname_str_len=256
  type, public :: histf_type
    character(len=40), pointer :: varnames(:) => null()
    character(len=6) , pointer :: hrfreq(:)   => null() !'hour','day','week','month','year'
    character(len=12), pointer :: units(:)    => null()
    real(r8)         , pointer :: counter(:)  => null()
    character(len=hist_fname_str_len), pointer :: ncfname(:)     => null()
    real(r8)         , pointer :: yvals(:)    => null()
    integer :: nvars
    integer :: nh_vars
    integer :: nd_vars
    integer :: nw_vars
    integer :: nm_vars
    integer :: ny_vars
    integer, pointer :: nh_varid(:) => null()
    integer, pointer :: nd_varid(:) => null()
    integer, pointer :: nw_varid(:) => null()
    integer, pointer :: nm_varid(:) => null()
    integer, pointer :: ny_varid(:) => null()
    integer, pointer :: record(:)   => null()
    real(r8) :: dtime
  contains
    procedure, public  :: init
    procedure, public  :: hist_wrap
    procedure, public  :: histrst
    procedure, private :: initAlloc
    procedure, private :: hist_create
    procedure, private :: hist_write
    procedure, private :: hist_add_var
    procedure, private :: proc_counter
    procedure, private :: reset_counter
    procedure, private :: proc_record

  end type histf_type

contains

  subroutine histrst(this, gfname, rwflag, yymmddhhss)

  !DESCRIPTION
  !writes file for restart
  use ncdio_pio, only : ncd_pio_openfile, ncd_pio_createfile
  use ncdio_pio, only : ncd_pio_closefile, ncd_defdim, ncd_defvar
  use ncdio_pio, only : ncd_putvar, ncd_enddef, ncd_unlimited
  use ncdio_pio, only : ncd_nowrite, ncd_float, ncd_getvar
  use shr_const_mod, only : spval => SHR_CONST_SPVAL
  implicit none
  class(histf_type), intent(inout):: this
  character(len=*), intent(in) :: gfname
  character(len=*), intent(in) :: rwflag
  character(len=*), intent(in) :: yymmddhhss
  character(len=128) :: fname
  type(file_desc_t) :: ncf
  real(r8), pointer :: ptr1d(:)
  integer :: recordDimID

  write(fname,'(A)')trim(gfname)//'.hr.'//trim(yymmddhhss)//'.nc'
  if(trim(rwflag)=='read')then
    call ncd_pio_openfile(ncf, fname, ncd_nowrite)
    ptr1d => this%yvals
    call ncd_getvar(ncf,'vars',ptr1d)
    ptr1d => this%counter
    call ncd_getvar(ncf,'counters',ptr1d)
    call ncd_pio_closefile(ncf)
  else
    call ncd_pio_createfile(ncf, trim(fname))
    !the temporal dimension is infinite
    call ncd_defdim(ncf,'clocks',5,recordDimID)
    call ncd_defdim(ncf,'numvar',ncd_unlimited,recordDimID)

    call ncd_defvar(ncf, 'counters', ncd_float,              &
            dim1name='clocks',long_name='counters',          &
            units='none', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'vars', ncd_float,              &
            dim1name='numvar',long_name='vars',          &
            units='none', missing_value=spval, fill_value=spval)

    call ncd_enddef(ncf)
    call ncd_putvar(ncf,'counters',this%counter)
    call ncd_putvar(ncf,'vars',this%yvals)
    call ncd_pio_closefile(ncf)
  endif

  end subroutine histrst

!--------------------------------------------------------
  subroutine init(this, varlist, unitlist, hrfreq, gfname, dtime)
  use ncdio_pio, only : ncd_enddef, ncd_pio_closefile
  implicit none
  class(histf_type), intent(inout):: this
  character(len=hist_var_str_len), intent(in) :: varlist(:)
  character(len=hist_unit_str_len), intent(in) :: unitlist(:)
  character(len=hist_freq_str_len), optional, intent(in) :: hrfreq(:)
  character(len=*), optional,  intent(in) :: gfname    !generic filename
  real(r8), optional, intent(in) :: dtime
  character(len=256) :: loc_gfname
  integer :: n, nh, nd, nw, nm, ny
  logical, target :: yes_hour, yes_day, yes_mon, yes_week, yes_year
  type(file_desc_t) :: ncid(nclocks)
  integer :: clock_id
  logical, pointer :: yes_flag

  SHR_ASSERT_ALL((size(varlist)   == size(unitlist)), errMsg(mod_filename,__LINE__))

  if(present(hrfreq))then
    SHR_ASSERT_ALL((size(varlist)   == size(hrfreq)), errMsg(mod_filename,__LINE__))
  endif
  if(present(dtime))then
    this%dtime = dtime
  else
    this%dtime = 1._r8
  endif
  do n = 1, nclocks
    ncid(n)%fh=-1
  enddo
  yes_hour=.false.; yes_day=.false.; yes_mon=.false.; yes_week=.false.; yes_year=.false.

  if(present(gfname))then
    write(loc_gfname,'(A)')trim(gfname)
  else
    loc_gfname='sompack_hist'
  endif
  this%nvars=size(varlist)

  this%nh_vars = 0
  this%nd_vars = 0
  this%nw_vars = 0
  this%nm_vars = 0
  this%ny_vars = 0
  if(present(hrfreq))then
    do n = 1, this%nvars
      select case (trim(hrfreq(n)))
      case ('hour')
        this%nh_vars=this%nh_vars+1
      case ('day')
        this%nd_vars=this%nd_vars+1
      case ('week')
        this%nw_vars=this%nw_vars+1
      case ('year')
        this%ny_vars=this%ny_vars+1
      case default
        this%nm_vars=this%nm_vars+1
      end select
    enddo
  else
    this%nm_vars=this%nvars
  endif

  call this%initAlloc()

  this%varnames(1:this%nvars) = varlist(1:this%nvars)
  this%units(1:this%nvars) = unitlist(1:this%nvars)

  nh=0
  nd=0
  nw=0
  nm=0
  ny=0
  if(present(hrfreq))then
    this%hrfreq(1:this%nvars) = hrfreq(1:this%nvars)
    do n = 1, this%nvars
      select case (trim(this%hrfreq(n)))
      case ('hour')
        clock_id=clock_hour
        yes_flag => yes_hour
        nh=nh+1
        this%nh_varid(nh) = n
      case ('day')
        clock_id=clock_day
        yes_flag => yes_day
        nd=nd+1
        this%nd_varid(nd) = n
      case ('week')
        clock_id=clock_week
        yes_flag => yes_week
        nw=nw+1
        this%nw_varid(nw) = n
      case ('year')
        clock_id=clock_year
        yes_flag => yes_year
        ny=ny+1
        this%ny_varid(ny) = n
      case default
        clock_id=clock_month
        yes_flag => yes_mon
        nm=nm+1
        this%nm_varid(nm) = n
      end select
      call this%hist_create(loc_gfname, yes_flag,trim(this%hrfreq(n)),ncid(clock_id))
      call this%hist_add_var(ncid(clock_id),this%varnames(n),this%units(n), trim(this%hrfreq(n)))
    enddo

  else
    !default output is by month
    do n = 1, this%nvars
      call this%hist_create(loc_gfname, yes_mon, 'month',ncid(clock_month))
      call this%hist_add_var(ncid(clock_month),this%varnames(n), this%units(n), 'month')
      nm=nm+1
      this%nm_varid(nm) = n
    enddo
  endif

  do n = 1, nclocks
    if(ncid(n)%fh>0)then
      call ncd_enddef(ncid(n))
      call ncd_pio_closefile(ncid(n))
    endif
  enddo
  end subroutine init

!--------------------------------------------------------

  subroutine initAlloc(this)
  implicit none
  class(histf_type), intent(inout):: this
  integer :: n

  allocate(this%varnames(this%nvars));this%varnames(:)=''
  allocate(this%hrfreq(this%nvars)); this%hrfreq(:)=''
  allocate(this%counter(nclocks)); this%counter(:)=0._r8
  allocate(this%units(this%nvars))
  allocate(this%yvals(this%nvars)); this%yvals(:) = 0._r8
  allocate(this%ncfname(nclocks)); this%ncfname(:)=''
  allocate(this%record(nclocks));  this%record(:)=0
  if(this%nh_vars>0)allocate(this%nh_varid(this%nh_vars))
  if(this%nd_vars>0)allocate(this%nd_varid(this%nd_vars))
  if(this%nw_vars>0)allocate(this%nw_varid(this%nw_vars))
  if(this%nm_vars>0)allocate(this%nm_varid(this%nm_vars))
  if(this%ny_vars>0)allocate(this%ny_varid(this%ny_vars))

  end subroutine initAlloc

!--------------------------------------------------------
  subroutine hist_wrap(this, yval, timer)
  use  BeTR_TimeMod, only : betr_time_type
  implicit none
  class(histf_type), intent(inout):: this
  real(r8), dimension(:), intent(in) :: yval
  type(betr_time_type), intent(in) :: timer

  integer :: clockid, id

  SHR_ASSERT_ALL((size(yval)   == this%nvars), errMsg(mod_filename,__LINE__))

  ! use daxpy - compute y := alpha * x + y
  ! SUBROUTINE DAXPY(N, ALPHA, X, INCX, Y, INCY)
  call daxpy(this%nvars, 1._r8, yval, 1, this%yvals, 1)

  call this%proc_counter()

  if(timer%its_a_new_hour())then
    call this%hist_write(clock_hour, this%nh_vars, this%nh_varid)
  endif

  if(timer%its_a_new_day())then
    call this%hist_write(clock_day, this%nd_vars, this%nd_varid)
  endif

  if(timer%its_a_new_week())then
    call this%hist_write(clock_week, this%nw_vars, this%nw_varid)
  endif

  if(timer%its_a_new_month())then
    call this%hist_write(clock_month, this%nm_vars, this%nm_varid)
  endif

  if(timer%its_a_new_year())then
    call this%hist_write(clock_year, this%ny_vars, this%ny_varid)
  endif

  end subroutine hist_wrap

!--------------------------------------------------------

  subroutine hist_create(this, gname, yesno, freq, ncid)
  !
  !DESCRIPTION
  ! create history file for writting
  use ncdio_pio, only : ncd_pio_createfile
  use betr_varcon, only : spval  => bspval
  use ncdio_pio, only : ncd_defvar
  use ncdio_pio, only : ncd_defdim, ncd_unlimited, ncd_float
  implicit none
  class(histf_type) , intent(inout):: this
  character(len=*)    , intent(in) :: gname
  logical          , intent(inout) :: yesno
  character(len=*) , intent(in) :: freq
  type(file_desc_t), intent(out):: ncid

  integer :: id
  integer :: recordDimID

  if(yesno)return
  yesno=.true.
  select case(freq)
    case ('hour')
      id = clock_hour
    case ('day')
      id = clock_day
    case ('week')
      id = clock_week
    case ('year')
      id = clock_year
    case default
      id = clock_month
  end select
  write(this%ncfname(id),'(A)')trim(gname)//'.hist.'//trim(freq)//'.nc'
  call ncd_pio_createfile(ncid, this%ncfname(id))

  !the temporal dimension is infinite
  call ncd_defdim(ncid,trim(freq),ncd_unlimited,recordDimID)

  end subroutine hist_create
!--------------------------------------------------------
  subroutine hist_add_var(this, ncid, varname, units, freq)

  use ncdio_pio,  only : ncd_defvar,  ncd_float, file_desc_t
  use betr_varcon, only : spval => bspval
  implicit none
  class(histf_type) , intent(inout):: this
  character(len=*), intent(in) :: varname
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: freq
  type(file_desc_t), intent(inout):: ncid

  call ncd_defvar(ncid, varname, ncd_float,              &
        dim1name=trim(freq),long_name=varname,               &
        units=units, missing_value=spval, fill_value=spval)

  end subroutine hist_add_var

!--------------------------------------------------------
  subroutine proc_counter(this)
  implicit none
  class(histf_type) , intent(inout):: this

  integer :: n

  do n = 1, nclocks
    this%counter(n)=this%counter(n)+1._r8
  enddo
  end subroutine proc_counter

!--------------------------------------------------------
  subroutine reset_counter(this, clockid)
  implicit none
  class(histf_type) , intent(inout):: this
  integer , intent(in) :: clockid

  this%counter(clockid)=0._r8
  end subroutine reset_counter

!--------------------------------------------------------
  subroutine proc_record(this, clockid)
  implicit none
  class(histf_type) , intent(inout):: this
  integer , intent(in) :: clockid

  this%record(clockid)=this%record(clockid)+1
  end subroutine proc_record

!--------------------------------------------------------
  subroutine hist_write(this, clockid, nvars, varid)
  use ncdio_pio, only : file_desc_t, ncd_pio_closefile, ncd_putvar
  use ncdio_pio, only : ncd_pio_openfile_for_write
  implicit none
  class(histf_type) , intent(inout):: this
  integer, intent(in) :: clockid
  integer, intent(in) :: nvars
  integer, intent(in) :: varid(1:nvars)

  integer :: n, id
  type(file_desc_t) :: ncid

  if(nvars==0)return
  call this%proc_record(clockid)
  call ncd_pio_openfile_for_write(ncid, this%ncfname(clockid))
  do n =1, nvars
    id = varid(n)
    this%yvals(id) = this%yvals(id)/this%counter(clockid)/this%dtime
    call ncd_putvar(ncid,this%varnames(id),this%record(clockid),this%yvals(id))
    this%yvals(id) = 0._r8
  enddo
  call this%reset_counter(clockid)
  call ncd_pio_closefile(ncid)
  end subroutine hist_write
end module histMod
