module ForcDataType


  use JarBgcForcType, only : JarBGC_forc_type
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use ncdio_pio   , only : file_desc_t, ncd_io
  use ncdio_pio    , only : ncd_nowrite
  use BetrStatusType   , only : betr_status_type
  use ncdio_pio   , only : ncd_pio_closefile, ncd_pio_openfile, ncd_nowrite
implicit none
  private
  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: init_forc
  public :: load_forc

  interface load_forc
    module procedure load_forc_const, load_forc_transient
  end interface load_forc


  type, private :: forc_data_type
   real(r8), pointer :: pbot(:)            !PBOT
   real(r8), pointer :: h2osoi_liq(:)      !SOILLIQ
   real(r8), pointer :: h2osoi_vol(:)      !H2OSOI
   real(r8), pointer :: air_vol(:)         !watsat - h2osoi_vol
   real(r8), pointer :: air_temp(:)        !TBOT
   real(r8), pointer :: temp(:)            !TSOI
   real(r8), pointer :: finudated(:)       !not provided, set to zero
   real(r8) :: dzsoi                       !DZSOI
   real(r8) :: zsoi                        !ZSOI
   real(r8) :: pct_sand
   real(r8) :: pct_clay
   real(r8) :: cellorg
   real(r8) :: pH
   integer  :: data_len
  contains
    procedure, public :: init
    procedure, private :: ReadForcingData
  end type forc_data_type

  type(forc_data_type) :: forc_data


contains

  subroutine init(this, data_len)

  implicit none
  class(forc_data_type), intent(inout) :: this
  integer, intent(in) :: data_len

  allocate(this%pbot(data_len))      ; this%pbot(:) = 1.01325e5_r8
  allocate(this%h2osoi_liq(data_len)); this%h2osoi_liq(:) = 3._r8
  allocate(this%h2osoi_vol(data_len)); this%h2osoi_vol(:) = 0.3_r8
  allocate(this%air_vol(data_len))   ; this%air_vol(:)  = 0.1_r8
  allocate(this%air_temp(data_len))  ; this%air_temp(:) = 298._r8
  allocate(this%finudated(data_len)) ; this%finudated(:) = 0._r8
  allocate(this%temp(data_len))      ; this%temp(:) = 283._r8
  this%dzsoi=0.1_r8
  this%zsoi=0.05_r8
  this%pct_sand=30._r8
  this%pct_clay=15._r8
  this%cellorg=50._r8
  this%pH = 6._r8
  end subroutine init
  !--------------------------------------------------------------------
  subroutine ReadForcingData(this, forcing_filename)
  use ncdio_pio    , only : ncd_getvar
  implicit none
  class(forc_data_type), intent(inout) :: this
  character(len=*), intent(in) :: forcing_filename

  type(file_desc_t) :: ncf_in_forc

  call ncd_pio_openfile(ncf_in_forc, forcing_filename, mode=ncd_nowrite)

  call ncd_getvar(ncf_in_forc, 'h2osoi_liq', this%h2osoi_liq)
  call ncd_getvar(ncf_in_forc, 'h2osoi_vol', this%h2osoi_vol)
  call ncd_getvar(ncf_in_forc, 'air_vol', this%air_vol)
  call ncd_getvar(ncf_in_forc, 'tair', this%air_temp)
  call ncd_getvar(ncf_in_forc, 'tsoi', this%temp)
  call ncd_getvar(ncf_in_forc, 'patm', this%pbot)
  call ncd_getvar(ncf_in_forc, 'dzsoi', this%dzsoi)
  call ncd_getvar(ncf_in_forc, 'zsoi', this%zsoi)
  call ncd_getvar(ncf_in_forc, 'pct_sand', this%pct_sand)
  call ncd_getvar(ncf_in_forc, 'pct_clay', this%pct_clay)
  call ncd_getvar(ncf_in_forc, 'pH', this%pH)
  call ncd_getvar(ncf_in_forc, 'cellorg', this%cellorg)

  call ncd_pio_closefile(ncf_in_forc)

  end subroutine ReadForcingData
  !--------------------------------------------------------------------
  subroutine init_forc(namelist_buffer)

  use ncdio_pio    , only : get_dim_len
  use betr_constants , only : betr_filename_length, betr_string_length_long
  use babortutils    , only : endrun
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  implicit none
  character(len=*), intent(in) :: namelist_buffer


  ! !LOCAL VARIABLES:
  integer                                :: nml_error
  character(len=betr_string_length_long) :: ioerror_msg
  character(len=betr_filename_length)    :: forcing_format, forcing_type_name, forcing_filename

  namelist / forcing_inparm / forcing_type_name, forcing_filename, forcing_format

  integer :: data_len
  type(file_desc_t) :: ncf_in_forc
  data_len=1

  if ( .true. )then
     ioerror_msg=''
     read(namelist_buffer, nml=forcing_inparm, iostat=nml_error, iomsg=ioerror_msg)
     if (nml_error /= 0) then
        call endrun(msg="ERROR reading forcing_inparm namelist "//errmsg(mod_filename, __LINE__))
     end if
  end if

  call ncd_pio_openfile(ncf_in_forc, forcing_filename, mode=ncd_nowrite)
  data_len=get_dim_len(ncf_in_forc,'time')

  call ncd_pio_closefile(ncf_in_forc)
  !determine data length

  call forc_data%init(data_len)
  forc_data%data_len = data_len

  call forc_data%ReadForcingData(forcing_filename)
  end subroutine init_forc

  !--------------------------------------------------------------------
  subroutine load_forc_const(soil_forc, atm_forc)
  !DESCRIPTION
  !load constant soil data

  use SoilForcType, only : soil_forc_type
  use AtmForcType , only : atm_forc_type
  implicit none
  class(soil_forc_type), intent(inout) :: soil_forc
  class(atm_forc_type) , intent(inout) :: atm_forc

  soil_forc%pct_clay = forc_data%pct_clay
  soil_forc%pct_sand = forc_data%pct_sand
  soil_forc%cellorg  = forc_data%cellorg
  soil_forc%dzsoi    = forc_data%dzsoi
  soil_forc%depz     = forc_data%zsoi
  soil_forc%pH       = forc_data%pH
  call soil_forc%init()
  call atm_forc%init()
  end subroutine load_forc_const
  !--------------------------------------------------------------------

  subroutine load_forc_transient(om_forc, nut_forc, atm_forc, soil_forc, nstep)

  !DESCRIPTION
  !load transient forcing data
  use SoilForcType, only : soil_forc_type
  use AtmForcType , only : atm_forc_type
  use OMForcType  , only : om_forc_type
  use NutForcType , only : nut_forc_type
  implicit none
  type(om_forc_type), intent(inout) :: om_forc
  type(nut_forc_type), intent(inout):: nut_forc
  type(atm_forc_type), intent(inout):: atm_forc
  type(soil_forc_type),intent(inout):: soil_forc
  integer, intent(in) :: nstep

  integer :: rec

  rec=mod(nstep,forc_data%data_len)
  if(rec==0)rec=forc_data%data_len

  atm_forc%patm_pascal = forc_data%pbot(rec)
  atm_forc%air_temp    = forc_data%air_temp(rec)
  soil_forc%h2osoi_liq = forc_data%h2osoi_liq(rec)
  soil_forc%h2osoi_vol = forc_data%h2osoi_vol(rec)
  soil_forc%air_vol    = max(forc_data%air_vol(rec),0._r8)
  soil_forc%temp       = forc_data%temp(rec)
  soil_forc%finundated = forc_data%finudated(rec)
  call soil_forc%update()

  om_forc%cflx_input_litr_met = 1.e-8_r8
  om_forc%cflx_input_litr_cel = 1.e-8_r8
  om_forc%cflx_input_litr_lig = 1.e-8_r8
  om_forc%cflx_input_litr_cwd = 1.e-8_r8
  om_forc%cflx_input_litr_fwd = 0._r8
  om_forc%cflx_input_litr_lwd = 0._r8
  nut_forc%sflx_minn_input_nh4= 1.e-10_r8
  nut_forc%sflx_minn_input_no3= 0._r8
  nut_forc%sflx_minp_input_po4= 1.e-12_r8

!  print*,'patm=',atm_forc%patm_pascal
!  print*,'air_temp=',atm_forc%air_temp
!  print*,'h2osoi_liq',soil_forc%h2osoi_vol
!  print*,'air_vol=',soil_forc%air_vol
!  print*,'tsoi=',soil_forc%temp
!  print*,'inudated=',soil_forc%finundated

  end subroutine load_forc_transient

end module ForcDataType
