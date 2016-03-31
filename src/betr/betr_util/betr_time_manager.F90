module betr_time_manager

   use bshr_kind_mod, only: r8 => shr_kind_r8


   implicit none
   private

   ! Public methods

   public ::&
        set_nstep               , &
        get_nstep               , & !
        get_days_per_year       , & !
        get_step_size              !return step size in seconds

   public :: proc_nextstep, &
        proc_initstep

   integer, save :: nelapstep
  contains

  integer function get_step_size()result(rc)

    ! Return the step size in seconds.
    implicit none


    rc = 1800
    return


  end function get_step_size

  !=========================================================================================
  integer function get_nstep()

    ! Return the timestep number.

    character(len=*), parameter :: sub = 'betr::get_nstep'


    get_nstep = nelapstep


  end function get_nstep


  !=========================================================================================
  subroutine set_nstep(nstep)

    ! Return the timestep number.
  implicit none
    character(len=*), parameter :: sub = 'betr::get_nstep'
    integer, intent(in) :: nstep

    nelapstep = nstep


  end subroutine set_nstep



  subroutine proc_nextstep()

  implicit none

  nelapstep = nelapstep + 1
  end subroutine proc_nextstep

  subroutine proc_initstep()

  implicit none

  nelapstep = 0
  end subroutine proc_initstep



  integer function get_days_per_year( offset )

  implicit none
  ! Arguments
  integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
  ! Positive for future times, negative
  ! for previous times.

  ! remove unused dummy arg compiler warning
  if (offset > 0) continue

  !hardwire at the moment
  get_days_per_year = 365
  end function get_days_per_year
end module betr_time_manager
