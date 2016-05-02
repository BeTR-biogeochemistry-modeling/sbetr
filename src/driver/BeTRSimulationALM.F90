module BeTRSimulationALM
  !
  ! !DESCRIPTION:
  !  API for using BeTR in ALM
  !
  ! !USES:
  !
#include "shr_assert.h"
  use abortutils        , only : endrun
  use clm_varctl        , only : iulog
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_kind_mod      , only : r8 => shr_kind_r8

  use BeTRSimulation    , only : betr_simulation_type
  use decompMod         , only : bounds_type

  use BeTRSimulation    , only : betr_simulation_type
  use BeTR_TimeMod      , only : betr_time_type
  use EcophysConType    , only : ecophyscon_type
  use tracer_varcon     , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use BeTR_PatchType    , only : betr_pft
  use BeTR_ColumnType   , only : betr_col
  use BeTR_LandunitType , only : betr_lun
  use betr_decompMod    , only : betr_bounds_type
  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_alm_type
     type(ecophyscon_type) :: ecophyscon
   contains
     procedure :: Init                              => ALMInit
     procedure, public :: StepWithoutDrainage       => ALMStepWithoutDrainage
     procedure, public :: StepWithDrainage          => ALMStepWithDrainage
     !unique subroutines
     procedure, public :: DiagnoseDtracerFreezeThaw => ALMDiagnoseDtracerFreezeThaw
     procedure, public :: CalcDewSubFlux            => ALMCalcDewSubFlux
     procedure, public :: SoilFluxStateRecv         => ALMBetrSoilFluxStateRecv
     procedure, public :: CalcSmpL                  => ALMCalcSmpL
     procedure, public :: PlantSoilBGCSend          => ALMBetrPlantSoilBGCSend
     procedure, public :: PlantSoilBGCRecv          => ALMBetrPlantSoilBGCRecv
  end type betr_simulation_alm_type

  public :: create_betr_simulation_alm

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_alm()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_simulation_alm_type), pointer :: create_betr_simulation_alm
    class(betr_simulation_alm_type), pointer :: simulation

    allocate(simulation)
    create_betr_simulation_alm => simulation

  end function create_betr_simulation_alm

!-------------------------------------------------------------------------------

  subroutine ALMInit(this, base_filename, namelist_buffer, &
       bounds, waterstate, cnstate)
    !DESCRIPTION
    !Initialize BeTR for ALM
    !
    !USES
    !data types from alm
    use PatchType       , only : pft
    use ColumnType      , only : col
    use LandunitType    , only : lun
    use pftvarcon       , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use WaterStateType  , only : waterstate_type
    use CNStateType     , only : cnstate_type
    use landunit_varcon , only : istcrop, istice, istsoil
    use clm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil

    !betr types
    use BeTRSimulation      , only : BeTRSimulationInit
    use betr_constants      , only : betr_filename_length
    use betr_constants      , only : betr_namelist_buffer_size
    use BeTR_PatchType      , only : betr_pft
    use BeTR_ColumnType     , only : betr_col
    use BeTR_LandunitType   , only : betr_lun
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use BeTR_landvarconType , only : betr_landvarcon
    use BeTR_decompMod      , only : betr_bounds_type

    implicit none
    class(betr_simulation_alm_type)          , intent(inout) :: this
    character(len=betr_filename_length)      , intent(in)    :: base_filename
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    type(bounds_type)                        , intent(in)    :: bounds
    type(waterstate_type)                    , intent(inout) :: waterstate
    type(cnstate_type)                       , intent(inout) :: cnstate

    !temporary variables
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj

    !grid size
    betr_nlevsoi      = nlevsoi
    betr_nlevsno      = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    !grid horizontal bounds
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

    betr_pft%wtcol                     => pft%wtcol
    betr_pft%column                    => pft%column
    betr_pft%itype                     => pft%itype
    betr_pft%landunit                  => pft%landunit

    betr_col%landunit                  => col%landunit
    betr_col%gridcell                  => col%gridcell
    betr_col%snl                       => col%snl
    betr_col%dz                        => col%dz
    betr_col%zi                        => col%zi
    betr_col%z                         => col%z

    betr_lun%itype                     => lun%itype
    betr_lun%ifspecial                 => lun%ifspecial

    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice


    ! now call the base simulation init to continue initialization
    call this%BeTRInit(base_filename, namelist_buffer, &
         betr_bounds, waterstate, cnstate)

  end subroutine ALMInit
!-------------------------------------------------------------------------------
  subroutine ALMStepWithoutDrainage(this, betr_time, bounds,  col ,   &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)
   !DESCRIPTION
   !march one time step without doing drainage
   !
   !USES
    use SoilStateType     , only : soilstate_type
    use WaterStateType    , only : Waterstate_Type
    use TemperatureType   , only : temperature_type
    use ChemStateType     , only : chemstate_type
    use WaterfluxType     , only : waterflux_type
    use ColumnType        , only : column_type
    use atm2lndType       , only : atm2lnd_type
    use SoilHydrologyType , only : soilhydrology_type
    use CNStateType       , only : cnstate_type
    use CNCarbonFluxType  , only : carbonflux_type
    use CanopyStateType   , only : canopystate_type
    use PatchType         , only : pft
    use LandunitType      , only : lun
    use pftvarcon         , only : crop
    use clm_varpar        , only : nlevsno, nlevsoi, nlevtrc_soil
    use tracer_varcon     , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    implicit none
    ! !ARGUMENTS :
    class(betr_simulation_alm_type) , intent(inout) :: this
    class(betr_time_type)           , intent(in)    :: betr_time
    type(bounds_type)               , intent(in)    :: bounds ! bounds
    type(column_type)               , intent(in)    :: col ! column type
    type(Waterstate_Type)           , intent(in)    :: waterstate_vars ! water state variables
    type(soilstate_type)            , intent(in)    :: soilstate_vars ! column physics variable
    type(temperature_type)          , intent(in)    :: temperature_vars ! energy state variable
    type(chemstate_type)            , intent(in)    :: chemstate_vars
    type(atm2lnd_type)              , intent(in)    :: atm2lnd_vars
    type(soilhydrology_type)        , intent(in)    :: soilhydrology_vars
    type(cnstate_type)              , intent(inout) :: cnstate_vars
    type(canopystate_type)          , intent(in)    :: canopystate_vars
    type(carbonflux_type)           , intent(in)    :: carbonflux_vars
    type(waterflux_type)            , intent(inout) :: waterflux_vars

    !TEMPORARY VARIABLES
    type(betr_bounds_type)     :: betr_bounds

    !pass necessary data for correct subroutine call
    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    betr_pft%wtcol     => pft%wtcol
    betr_pft%column    => pft%column
    betr_pft%itype     => pft%itype
    betr_pft%landunit  => pft%landunit
    betr_pft%crop      => crop

    betr_col%landunit  => col%landunit
    betr_col%gridcell  => col%gridcell
    betr_col%snl       => col%snl
    betr_col%dz        => col%dz
    betr_col%zi        => col%zi
    betr_col%z         => col%z

    betr_lun%itype     => lun%itype
    betr_lun%ifspecial => lun%ifspecial

    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

    call this%SetBiophysForcing(betr_bounds, &
      cnstate_vars = cnstate_vars,           &
      carbonflux_vars=carbonflux_vars,       &
      waterstate_vars=waterstate_vars,       &
      waterflux_vars=waterflux_vars,         &
      temperature_vars=temperature_vars,     &
      soilhydrology_vars=soilhydrology_vars, &
      atm2lnd_vars=atm2lnd_vars,             &
      canopystate_vars=canopystate_vars,     &
      chemstate_vars=chemstate_vars,         &
      soilstate_vars=soilstate_vars)

    call this%betr%step_without_drainage(betr_time, betr_bounds,               &
         this%num_soilc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc, this%biogeo_flux, this%biogeo_state)

  end subroutine ALMStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine ALMStepWithDrainage(this, bounds,  col)
   !DESCRIPTION
   ! march one step with drainage
   !
   !USES
    use ColumnType     , only : column_type
    use WaterFluxType  , only : waterflux_type
    use LandunitType   , only : lun
    use clm_varpar     , only : nlevsno, nlevsoi, nlevtrc_soil

    use betr_decompMod , only : betr_bounds_type
    use tracer_varcon  , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_alm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    type(column_type)               , intent(in)    :: col ! column type

    !temporary variables
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0

    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    betr_col%landunit  => col%landunit
    betr_col%gridcell  => col%gridcell
    betr_col%snl       => col%snl
    betr_col%dz        => col%dz
    betr_col%zi        => col%zi
    betr_col%z         => col%z

    betr_lun%itype     => lun%itype
    betr_lun%ifspecial => lun%ifspecial

    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

    call this%betr%step_with_drainage(betr_bounds,      &
         this%num_soilc, this%filter_soilc, this%jtops, &
         this%biogeo_flux)

  end subroutine ALMStepWithDrainage

  !------------------------------------------------------------------------
  subroutine ALMBetrPlantSoilBGCSend(this, bounds, num_soilc,  filter_soilc)
  !this should be expanded
  implicit none
  class(betr_simulation_alm_type) :: this
  type(bounds_type) , intent(in)  :: bounds
  integer           , intent(in)  :: num_soilc
  integer           , intent(in)  :: filter_soilc(:)

  ! remove compiler warnings
  if (this%num_soilc > 0)     continue
  if (bounds%begc > 0)        continue
  if (num_soilc > 0)          continue
  if (size(filter_soilc) > 0) continue

  !pull in fluxes of external inputs and plant demand
  !call this%betr%plant_soilbgc%lsm_betr_plant_soilbgc_send(bounds, num_soilc, &
  !               filter_soilc, this%biogeo_state, this%biogeo_flux, ecophyscon_vars)

  !pull in all state variables and update tracers
  end subroutine ALMBetrPlantSoilBGCSend

  !------------------------------------------------------------------------
  subroutine ALMBetrPlantSoilBGCRecv(this, bounds, num_soilc,  filter_soilc)
  !this should be expanded
  implicit none
  class(betr_simulation_alm_type) :: this
  type(bounds_type) , intent(in)  :: bounds
  integer           , intent(in)  :: num_soilc
  integer           , intent(in)  :: filter_soilc(:)

  ! remove compiler warnings
  if (this%num_soilc > 0)     continue
  if (bounds%begc > 0)        continue
  if (num_soilc > 0)          continue
  if (size(filter_soilc) > 0) continue

  !pull in fluxes of external inputs and plant demand
  !call this%betr%plant_soilbgc%lsm_betr_plant_soilbgc_recv(bounds, num_soilc, &
  !               filter_soilc, this%biogeo_flux)
  !reset the biogeo_flux
  !pull in all state variables and update tracers
  end subroutine ALMBetrPlantSoilBGCRecv
  !------------------------------------------------------------------------

  subroutine ALMDiagnoseDtracerFreezeThaw(this, bounds, num_nolakec, filter_nolakec, col, lun, &
    waterstate_vars)
  !
  ! DESCRIPTION
  ! aqueous tracer partition based on freeze-thaw
  !
  ! USES
  use ColumnType            , only : column_type
  use LandunitType          , only : landunit_type
  use WaterStateType        , only : waterstate_type
  implicit none
  !
  ! Arguments
  class(betr_simulation_alm_type)    :: this
  type(bounds_type)     , intent(in) :: bounds
  integer               , intent(in) :: num_nolakec                        ! number of column non-lake points in column filter
  integer               , intent(in) :: filter_nolakec(:)                  ! column filter for non-lake points
  type(landunit_type)   , intent(in) :: lun
  type(waterstate_type) , intent(in) :: waterstate_vars
  type(column_type)     , intent(in) :: col                                ! column type

  !temporary variables
  type(betr_bounds_type)     :: betr_bounds

  betr_col%landunit  => col%landunit
  betr_col%gridcell  => col%gridcell
  betr_col%snl       => col%snl
  betr_col%dz        => col%dz
  betr_col%zi        => col%zi
  betr_col%z         => col%z

  betr_lun%itype     => lun%itype
  betr_lun%ifspecial => lun%ifspecial

  betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
  betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
  betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
  betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
  betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

  call this%SetBiophysForcing(betr_bounds, &
      waterstate_vars=waterstate_vars)

  call this%betr%diagnose_dtracer_freeze_thaw(betr_bounds, num_nolakec, filter_nolakec,  &
    this%biophys_forc)

  end subroutine ALMDiagnoseDtracerFreezeThaw

  !------------------------------------------------------------------------
  subroutine ALMCalcDewSubFlux(this, betr_time, &
       bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars)
   !DESCRIPTION
    ! Calculate tracer flux from dew or/and sublimation
    !External interface called by ALM
    use ColumnType      , only : col
    use LandunitType    , only : lun
    use WaterfluxType   , only : waterflux_type
    use WaterstateType  , only : waterstate_type
    use clm_varcon      , only : denh2o,spval
    use landunit_varcon , only : istsoil, istcrop
    use betr_decompMod  , only : betr_bounds_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_alm_type) , intent(inout) :: this
    class(betr_time_type)           , intent(in)    :: betr_time
    type(betr_bounds_type)          , intent(in)    :: bounds
    integer                         , intent(in)    :: num_hydrologyc ! number of column soil points in column filter_soilc
    integer                         , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
    type(waterstate_type)           , intent(in)    :: waterstate_vars
    type(waterflux_type)            , intent(in)    :: waterflux_vars

    call this%SetBiophysForcing(bounds, &
      waterstate_vars=waterstate_vars,  &
      waterflux_vars=waterflux_vars)

    call this%betr%calc_dew_sub_flux(betr_time,           &
         bounds, num_hydrologyc, filter_soilc_hydrologyc, &
        this%biophys_forc, this%betr%tracers, this%betr%tracerfluxes, this%betr%tracerstates)

  end subroutine ALMCalcDewSubFlux

  !------------------------------------------------------------------------
  subroutine ALMBetrSoilFluxStateRecv(this, bounds, num_soilc, filter_soilc)
  !this should be expanded and called after tracer update with drainage
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_alm_type) :: this
  type(bounds_type) , intent(in)  :: bounds
  integer           , intent(in)  :: num_soilc
  integer           , intent(in)  :: filter_soilc(:)

  ! remove compiler warnings
  if (this%num_soilc > 0)     continue
  if (bounds%begc > 0)        continue
  if (num_soilc > 0)          continue
  if (size(filter_soilc) > 0) continue

!x  call this%betr%bgc_reaction%betr_lsm_flux_state_sendback(bounds, &
!       num_soilc, filter_soilc, this%biogeo_state, this%biogeo_flux)

  end subroutine ALMBetrSoilFluxStateRecv

  !------------------------------------------------------------------------
  subroutine ALMCalcSmpL(this, bounds, lbj, ubj, numf, filter, t_soisno, &
     soilstate_vars, waterstate_vars, soil_water_retention_curve)
  !DESCRIPTION
  ! calculate soil suction potential
  !
  !USES
  use SoilStateType              , only : soilstate_type
  use WaterStateType             , only : waterstate_type
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
  use clm_varcon                 , only : grav,hfus,tfrz
  implicit none
  !ARGUMENTS
  class(betr_simulation_alm_type)                        :: this
  type(bounds_type)                      , intent(in)    :: bounds  ! bounds
  integer                                , intent(in)    :: lbj, ubj                                          ! lower and upper bounds, make sure they are > 0
  integer                                , intent(in)    :: numf                                              ! number of columns in column filter
  integer                                , intent(in)    :: filter(:)                                         ! column filter
  real(r8)                               , intent(in)    :: t_soisno(bounds%begc: , lbj: )                    ! soil temperature
  type(soilstate_type)                   , intent(in)    :: soilstate_vars
  type(waterstate_type)                  , intent(inout) :: waterstate_vars
  class(soil_water_retention_curve_type) , intent(in)    :: soil_water_retention_curve

  !local variables
  real(r8) :: s_node
  integer  :: fc, c, j

  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(mod_filename,__LINE__))

  ! remove compiler warnings
  if (this%num_soilc > 0) continue

  associate(                                                     & !
    h2osoi_vol        =>    waterstate_vars%h2osoi_vol_col     , & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
    smp_l             =>    waterstate_vars%smp_l_col          , & ! Output: [real(r8) (:,:) ]  soil suction (mm)
    bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
    watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
    sucsat            =>    soilstate_vars%sucsat_col            & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
  )

  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(j>=1)then
        if(t_soisno(c,j)<tfrz)then
          smp_l(c,j)= hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
        else
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_l(c,j))
        endif

      endif
    enddo
  enddo
  end associate
  end subroutine ALMCalcSmpL

end module BeTRSimulationALM
