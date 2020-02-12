module SoilForcType

  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: soil_forc_type
  real(r8) :: temp                   !temperature
  real(r8) :: h2osoi_vol             !volumetric soil moisture, [m3/m3]
  real(r8) :: h2osoi_liq             !liquid water content, [kg/m2]
  real(r8) :: h2osoi_liqvol          !volumetric liquid water content, [m3/m3]
  real(r8) :: air_vol                !volumetric air, [m3/m3]
  real(r8) :: finundated             !innudated fraction of the soil, [none]
  real(r8) :: soilpsi                !soilwater pontential  [MPa], positive
  real(r8) :: depz                   !depth of the soil sampled, [m]
  real(r8) :: dzsoi                  !soil thickness, [m]
  real(r8) :: sucsat                 ! minimum soil suction [mm]
  real(r8) :: bsw                    !
  real(r8) :: bd                     !bulk density, [kg /m3]
  real(r8) :: pct_sand               ! percent of sand in soil
  real(r8) :: pct_clay               ! percent of clay in soil
  real(r8) :: watsat                 ! saturated water content [m3/m3]
  real(r8) :: watfc                  ! field capacity [m3/m3]
  real(r8) :: cellorg                ! soil organic matter content, [kg/m3]
  real(r8) :: pH                     ! soil pH
  real(r8) :: hksat                  ! saturated hydraulic conductivity of mineral soil, [mm/s]
  real(r8) :: hk                     ! hydraulic conductivity [mm/s]
  real(r8) :: Diff_Darcy             ! Darcy diffusivity [m2/s]
  real(r8) :: tauaqu                 ! aqueous tortuosity
  real(r8) :: taugas                 ! gaseous tortuosity
  contains
    procedure, public  :: init
    procedure, private :: get_tortuosity
    procedure, public  :: update
  end type soil_forc_type


contains


  subroutine init(this)

  implicit none
  class(soil_forc_type), intent(inout) :: this

  real(r8) :: om_frac
  real(r8), parameter :: organic_max = 160._r8

  om_frac= (this%cellorg/organic_max)**2._r8

  call set_soil_hydro_property(this%pct_sand, this%pct_clay, om_frac, this%depz, this%bd, &
       this%watsat,this%bsw, this%sucsat, this%hksat)

  end subroutine init
  !-----------------------------------------------------------------------

  subroutine get_tortuosity(this)

  implicit none
  class(soil_forc_type), intent(inout) :: this

  real(r8) :: ice_vol
  real(r8) :: eff_por

  ice_vol = max(this%h2osoi_vol-this%h2osoi_liqvol,0._r8)
  eff_por = this%watsat-ice_vol

  this%taugas = (this%air_vol/eff_por)**(3._r8/this%bsw)*this%air_vol
  this%tauaqu = (min(this%h2osoi_liqvol/eff_por,1._r8))**(this%bsw/3._r8-1._r8)*this%h2osoi_liqvol
  end subroutine get_tortuosity

  !-----------------------------------------------------------------------

  subroutine update(this)

  use betr_varcon , only : grav => bgrav
  use betr_varcon , only : hfus => bhfus,tfrz => btfrz
  use betr_varcon , only : denh2o  => bdenh2o
  implicit none
  class(soil_forc_type), intent(inout) :: this

  real(r8), parameter :: e_ice=6._r8
  real(r8) :: icefrac
  real(r8) :: s, imped
  real(r8) :: dsmpds

  this%h2osoi_liqvol = this%h2osoi_liq/(denh2o*this%dzsoi)

  s= max(this%h2osoi_liqvol/this%watsat,1.e-3_r8)
  icefrac = (this%h2osoi_vol-this%h2osoi_liqvol)/this%watsat
  imped = 10._r8**(-e_ice*icefrac)

  call soil_hk(this%hksat, imped, s, this%bsw, this%hk)

  if(this%temp<tfrz)then
    this%soilpsi = -hfus*(tfrz-this%temp)/(grav*this%temp) * 1000._r8   ![mm]
    dsmpds = 0._r8
  else
    call soil_suction(this%sucsat, s, this%bsw, this%soilpsi, dsmpds)
  endif
  this%soilpsi = max(this%soilpsi * grav*1.e-6_r8, -15._r8)
  !Darcy diffusivity in m2/s
  this%Diff_Darcy = this%hk*dsmpds/this%watsat * 1.e-6_r8

  call this%get_tortuosity()
  end subroutine update
  !-----------------------------------------------------------------------
  subroutine set_soil_hydro_property(sand, clay, om_frac, zsoi, bd, watsat, bsw, sucsat, hksat)

  use FuncPedotransferMod, only : pedotransf, init_pedof, get_ipedof
  implicit none
  real(r8), intent(in) :: sand
  real(r8), intent(in) :: clay
  real(r8), intent(in) :: zsoi
  real(r8), intent(in) :: om_frac
  real(r8), intent(out):: bd
  real(r8), intent(out):: watsat
  real(r8), intent(out):: bsw
  real(r8), intent(out):: sucsat
  real(r8), intent(out):: hksat

  integer :: ipedof
  real(r8):: om_watsat
  real(r8):: xksat
  real(r8):: om_b, om_sucsat, om_hksat
  real(r8), parameter :: zsapric   = 0.5_r8       ! depth (m) that organic matter takes on characteristics of sapric peat
  real(r8), parameter :: pcalpha   = 0.5_r8       ! percolation threshold
  real(r8), parameter :: pcbeta    = 0.139_r8     ! percolation exponent
  real(r8) :: perc_norm, perc_frac
  real(r8) :: uncon_frac, uncon_hksat
  call init_pedof
  ipedof = get_ipedof(0)

  call pedotransf(ipedof, sand, clay, watsat, bsw, sucsat, xksat)

  om_watsat  = max(0.93_r8 - 0.1_r8   *(zsoi/zsapric), 0.83_r8)
  om_b       = min(2.7_r8  + 9.3_r8   *(zsoi/zsapric), 12.0_r8)
  om_sucsat  = min(10.3_r8 - 0.2_r8   *(zsoi/zsapric), 10.1_r8)
  om_hksat   = max(0.28_r8 - 0.2799_r8*(zsoi/zsapric), 0.0001_r8)

  bd  = (1._r8 - watsat)*2.7e3_r8
  watsat    = (1._r8 - om_frac) * watsat + om_watsat*om_frac

  bsw       = (1._r8-om_frac) * (2.91_r8 + 0.159_r8*clay) + om_frac*om_b
  sucsat    = (1._r8-om_frac) * sucsat + om_sucsat*om_frac


  ! perc_frac is zero unless perf_frac greater than percolation threshold
  if (om_frac > pcalpha) then
    perc_norm=(1._r8 - pcalpha)**(-pcbeta)
    perc_frac=perc_norm*(om_frac - pcalpha)**pcbeta
  else
    perc_frac=0._r8
  endif

  ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
  uncon_frac=(1._r8-om_frac)+(1._r8-perc_frac)*om_frac

  ! uncon_hksat is series addition of mineral/organic conductivites
  if (om_frac < 1._r8) then
    uncon_hksat=uncon_frac/((1._r8-om_frac)/xksat &
         +((1._r8-perc_frac)*om_frac)/om_hksat)
  else
    uncon_hksat = 0._r8
  end if
  hksat  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat

  end subroutine set_soil_hydro_property

  !-----------------------------------------------------------------------
  subroutine soil_hk(hksat, imped, s, bsw, hk, dhkds)
    !
    ! !DESCRIPTION:
    ! Compute hydraulic conductivity
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hksat    !saturated hydraulic conductivity [mm/s]
    real(r8), intent(in) :: imped    !ice impedance
    real(r8), intent(in) :: s        !reletive saturation, [0, 1]
    real(r8), intent(in) :: bsw      !shape parameter
    real(r8), intent(out):: hk       !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out):: dhkds    !d[hk]/ds   [mm/s]
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'soil_hk'
    !-----------------------------------------------------------------------

    !compute hydraulic conductivity
    hk=imped*hksat*s**(2._r8*bsw+3._r8)

    !compute the derivative
    if(present(dhkds))then
       dhkds=(2._r8*bsw+3._r8)*hk/s
    endif

  end subroutine soil_hk

  !-----------------------------------------------------------------------
  subroutine soil_suction(smpsat, s, bsw, smp, dsmpds)
    !
    ! !DESCRIPTION:
    ! Compute soil suction potential
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: smpsat   !minimum soil suction, positive [mm]
    real(r8), intent(in)            :: s        !reletive saturation, [0, 1]
    real(r8), intent(in)            :: bsw      !shape parameter
    real(r8), intent(out)           :: smp      !soil suction, negative, [mm]
    real(r8), optional, intent(out) :: dsmpds   !d[smp]/ds, [mm]
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'soil_suction'
    !-----------------------------------------------------------------------

    !compute soil suction potential, negative
    smp = -smpsat*s**(-bsw)

    !compute derivative
    if(present(dsmpds))then
       dsmpds=-bsw*smp/s
    endif

  end subroutine soil_suction



end module SoilForcType
