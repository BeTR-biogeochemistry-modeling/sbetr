module shistFileMod
!
! module contains subroutines to define field and dimensions

implicit none
  private
  save
  public :: hist_file_create
  public :: hist_def_fld2d
  public :: hist_def_fld1d

contains

  subroutine hist_file_create(ncid, nlevgrnd, ncol)
  !
  ! DESCRIPTIONS
  ! 
  ! Create a netcdf file
  !
  use netcdf
  use clm_varcon, only : spval
  use ncdio_pio, only : check_ret, ncd_defvar, file_desc_t, ncd_defdim
  implicit none
  class(file_desc_t), intent(inout) :: ncid       
  integer,            intent(in) :: nlevgrnd  ! number of vertical layers
  integer,            intent(in) :: ncol      ! number of columns
  !returning variable
  
  !local variables
  character(len=255) :: subname='hist_file_creat'
  integer :: recordDimID

! define the dimensions
  !the temporal dimension is infinite
  call ncd_defdim(ncid,'time',nf90_unlimited,recordDimID)
    
  !number of vertical layers
  call ncd_defdim(ncid, 'levgrnd', nlevgrnd, recordDimID)

  !number of columns
  call ncd_defdim(ncid, 'ncol', ncol, recordDimID)

  !define the time dimension
  call ncd_defvar(ncid, 'time',nf90_float,dim1name='time',  &
      long_name='time', units = 'none',   &
      missing_value=spval, fill_value=spval)

  end subroutine hist_file_create
  
!--------------------------------------------------------------------------------  
  subroutine hist_def_fld2d(ncid, varname, nf90_type, dim1name, dim2name, long_name, units)
  !
  ! DESCRIPTION
  ! define 1d field
  !
  use ncdio_pio,  only : ncd_defvar,file_desc_t
  use clm_varcon, only : spval
  implicit none
  class(file_desc_t), intent(inout) :: ncid                    !file id
  character(len=*), intent(in) :: varname        !variable name
  integer,          intent(in) :: nf90_type      !data type
  character(len=*), intent(in) :: dim1name       !name of the 1st dim
  character(len=*), intent(in) :: dim2name       !name of the 2nd dim
  character(len=*), intent(in) :: long_name      !long name of the variable
  character(len=*), intent(in) :: units          !variable units
  
  call ncd_defvar(ncid, varname, nf90_type,                  &
        dim1name=dim1name,dim2name=dim2name,                 &
        dim3name="time",long_name=long_name,                 &
        units=units, missing_value=spval, fill_value=spval)
  
  end subroutine hist_def_fld2d

!--------------------------------------------------------------------------------
  subroutine hist_def_fld1d(ncid, varname, nf90_type, dim1name, long_name, units)
  !
  ! DESCRIPTION
  ! define 1d field
  !
  use netcdf
  use ncdio_pio, only : ncd_defvar, file_desc_t
  use clm_varcon, only : spval 
  implicit none
  class(file_desc_t), intent(inout) :: ncid                   !file id
  character(len=*), intent(in) :: varname       !variable name
  integer,          intent(in) :: nf90_type     !data type
  character(len=*), intent(in) :: dim1name      !name of the 1st dim
  character(len=*), intent(in) :: long_name     !long name of the variable
  character(len=*), intent(in) :: units         !variable units
  
  call ncd_defvar(ncid, varname, nf90_type,                 & 
        dim1name=dim1name,                                  &
        dim2name="time",long_name=long_name,                &
        units=units, missing_value=spval, fill_value=spval)
  
  end subroutine hist_def_fld1d
  

end module shistFileMod
