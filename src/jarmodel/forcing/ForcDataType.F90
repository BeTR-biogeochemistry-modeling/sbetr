module ForcDataType


  use JarBgcForcType, only : JarBGC_forc_type
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use ncdio_pio   , only : file_desc_t, ncd_io
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
   real(r8), pointer :: pbot(:)
   real(r8), pointer :: h2osoi_liq(:)
   real(r8), pointer :: h2osoi_vol(:)
   real(r8), pointer :: air_vol(:)
   real(r8), pointer :: air_temp(:)
   real(r8), pointer :: temp(:)
   real(r8), pointer :: finudated(:)
   real(r8) :: dzsoi
   real(r8) :: zsoi
   real(r8) :: pct_sand
   real(r8) :: pct_clay
   real(r8) :: cellorg
   real(r8) :: pH
   integer  :: data_len
  contains
    procedure, public :: init
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
  allocate(this%finudated(data_len)); this%finudated(:) = 0._r8
  allocate(this%temp(data_len))      ; this%temp(:) = 283._r8
  this%dzsoi=0.1_r8
  this%zsoi=0.05_r8
  this%pct_sand=30._r8
  this%pct_clay=15._r8
  this%cellorg=50._r8

  end subroutine init

  !--------------------------------------------------------------------
  subroutine init_forc(nc_forc_file)

  implicit none
  character(len=*), intent(in) :: nc_forc_file

  integer :: data_len
  data_len=1
  !determine data length
  call forc_data%init(data_len)
  forc_data%data_len = data_len

  end subroutine init_forc
  !--------------------------------------------------------------------
  subroutine load_forc_const(soil_forc)
  !DESCRIPTION
  !load constant soil data

  use SoilForcType, only : soil_forc_type
  implicit none
  class(soil_forc_type), intent(inout) :: soil_forc


  soil_forc%pct_clay = forc_data%pct_clay
  soil_forc%pct_sand = forc_data%pct_sand
  soil_forc%cellorg  = forc_data%cellorg
  soil_forc%dzsoi    = forc_data%dzsoi
  soil_forc%depz     = forc_data%zsoi
  soil_forc%pH       = forc_data%pH
  call soil_forc%init()
  end subroutine load_forc_const
  !--------------------------------------------------------------------

  subroutine load_forc_transient(om_forc, atm_forc, soil_forc, nstep)

  !DESCRIPTION
  !load transient forcing data
  use SoilForcType, only : soil_forc_type
  use AtmForcType , only : atm_forc_type
  use OMForcType  , only : om_forc_type

  implicit none
  type(om_forc_type), intent(inout) :: om_forc
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
  soil_forc%air_vol    = forc_data%air_vol(rec)
  soil_forc%temp       = forc_data%temp(rec)
  soil_forc%finundated = forc_data%finudated(rec)


  end subroutine load_forc_transient



end module ForcDataType
