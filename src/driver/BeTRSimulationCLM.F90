module BeTRSimulationCLM
  !
  ! !DESCRIPTION:
  !  factory to load the specific betr simulator
  !
  ! !USES:
  !
#include "shr_assert.h"

  use shr_kind_mod        , only : r8 => shr_kind_r8

  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_log_mod, only : errMsg => shr_log_errMsg

  use BeTRSimulation, only : betr_simulation_type

  use decompMod, only : bounds_type
  use EcophysConType, only : ecophyscon_type
  use BeTR_PatchType, only : betr_pft
  use BeTR_ColumnType, only : betr_col
  use BeTR_LandunitType,only : betr_lun
  use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_clm_type

     type(ecophyscon_type) :: ecophyscon

     ! NOTE(bja, 201603) CLM stubb types go here!

   contains
     procedure, public :: Init => CLMInit
     procedure, public :: StepWithoutDrainage => CLMStepWithoutDrainage
     procedure, public :: StepWithDrainage => CLMStepWithDrainage
     procedure, public :: ConsistencyCheck => betr_clm_h2oiso_consistency_check
     procedure, public :: diagnose_dtracer_freeze_thaw_clm
     procedure, public :: calc_dew_sub_flux_clm
     procedure, public :: clm_betr_flux_state_receive
     procedure, public :: calc_smp_l_clm
  end type betr_simulation_clm_type

  public :: create_betr_simulation_clm

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_clm()

    implicit none

    class(betr_simulation_clm_type), pointer :: create_betr_simulation_clm
    class(betr_simulation_clm_type), pointer :: simulation

    allocate(simulation)

    create_betr_simulation_clm => simulation

  end function create_betr_simulation_clm


  !-------------------------------------------------------------------------------

  subroutine CLMInit(this, base_filename, namelist_buffer, bounds, waterstate, cnstate)

    use BeTRSimulation, only : BeTRSimulationInit
    use betr_constants, only : betr_namelist_buffer_size
    use betr_constants, only : betr_filename_length

    use BeTR_PatchType, only : betr_pft
    use BeTR_ColumnType, only : betr_col
    use BeTR_LandunitType,only : betr_lun
    use BeTR_pftvarconType, only : betr_pftvarcon
    use PatchType, only : pft
    use ColumnType, only : col
    use LandunitType,only : lun
    use pftvarcon, only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass

    use WaterStateType, only : waterstate_type
    use BeTR_WaterStateType, only : betr_waterstate_type

    use CNStateType, only : cnstate_type
    use BeTR_CNStateType, only : betr_cnstate_type

    use landunit_varcon
    use BeTR_landvarconType, only : betr_landvarcon
    use BeTR_decompMod              , only : betr_bounds_type
    use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    character(len=betr_filename_length), intent(in) :: base_filename
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    type(bounds_type), intent(in) :: bounds
    type(waterstate_type), intent(inout) :: waterstate
    type(cnstate_type), intent(inout) :: cnstate

    type(betr_bounds_type)     :: betr_bounds
    type(betr_waterstate_type) :: betr_waterstate
    type(betr_cnstate_type) :: betr_cnstate
    integer :: lbj, ubj

    betr_nlevsoi = nlevsoi
    betr_nlevsno = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    betr_pft%wtcol                        => pft%wtcol
    betr_pft%column                       => pft%column
    betr_pft%itype                        => pft%itype
    betr_pft%landunit                     => pft%landunit

    betr_col%landunit                     => col%landunit
    betr_col%gridcell                     => col%gridcell
    betr_col%snl                          => col%snl
    betr_col%dz                           => col%dz
    betr_col%zi                           => col%zi
    betr_col%z                            => col%z

    betr_lun%itype                        => lun%itype
    betr_lun%ifspecial                    => lun%ifspecial

    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice

    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    betr_waterstate%h2osoi_liq_col => waterstate%h2osoi_liq_col
    betr_waterstate%h2osoi_ice_col => waterstate%h2osoi_ice_col

    betr_cnstate%isoilorder  => cnstate%isoilorder

    ! allocate the reaction types that may only be known to this
    ! simulation type.

    ! now call the base simulation init to continue initialization
    call this%BeTRInit(base_filename, namelist_buffer, &
         betr_bounds, betr_waterstate, betr_cnstate)

  end subroutine CLMInit


  !---------------------------------------------------------------------------------
  subroutine CLMStepWithoutDrainage(this, bounds, col, &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)

    use SoilStateType, only : soilstate_type
    use WaterStateType, only : Waterstate_Type
    use TemperatureType, only : temperature_type
    use ChemStateType, only : chemstate_type
    use WaterfluxType, only : waterflux_type
    use ColumnType, only : column_type
    use BGCReactionsMod, only : bgc_reaction_type
    use atm2lndType, only : atm2lnd_type
    use SoilHydrologyType, only : soilhydrology_type
    use BeTR_CarbonFluxType, only : betr_carbonflux_type
    use CNStateType, only : cnstate_type
    use CNCarbonFluxType, only : carbonflux_type
    use CanopyStateType, only : canopystate_type


    use BeTR_WaterstateType   , only : betr_waterstate_type
    use BeTR_WaterfluxType    , only : betr_waterflux_type
    use BeTR_TemperatureType, only : betr_temperature_type
    use BeTR_SoilHydrologyType, only : betr_soilhydrology_type
    use BeTR_atm2lndType, only : betr_atm2lnd_type
    use betr_decompMod    , only : betr_bounds_type
    use BeTR_CanopyStateType, only : betr_canopystate_type
    use BeTR_ChemStateType, only : betr_chemstate_type
    use BeTR_SoilStateType, only : betr_soilstate_type
    use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use PatchType, only : pft
    use LandunitType,only : lun
    use pftvarcon, only : crop

    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this

    ! !ARGUMENTS :
    type(bounds_type), intent(in) :: bounds ! bounds
    type(column_type), intent(in) :: col ! column type
    type(Waterstate_Type), intent(in) :: waterstate_vars ! water state variables
    type(soilstate_type), intent(in) :: soilstate_vars ! column physics variable
    type(temperature_type), intent(in) :: temperature_vars ! energy state variable
    type(chemstate_type), intent(in) :: chemstate_vars
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    type(cnstate_type), intent(inout) :: cnstate_vars
    type(canopystate_type), intent(in) :: canopystate_vars
    type(carbonflux_type), intent(in) :: carbonflux_vars
    type(waterflux_type), intent(inout) :: waterflux_vars

    !temporary variables
    !temporary variables
    type(betr_waterflux_type)  :: betr_waterflux_vars
    type(betr_waterstate_type)  :: betr_waterstate_vars
    type(betr_temperature_type):: betr_temperature_vars
    type(betr_soilhydrology_type) :: betr_soilhydrology_vars
    type(betr_atm2lnd_type) :: betr_atm2lnd_vars
    type(betr_canopystate_type)  :: betr_canopystate_vars
    type(betr_chemstate_type) :: betr_chemstate_vars
    type(betr_soilstate_type) :: betr_soilstate_vars ! column physics variable
    type(betr_bounds_type)     :: betr_bounds


    integer  :: lbj, ubj ! lower and upper bounds, make sure they are > 0

    !pass necessary data for correct subroutine call
    betr_nlevsoi = nlevsoi
    betr_nlevsno = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    this%betr%cnstates%isoilorder => cnstate_vars%isoilorder

    this%betr%carbonfluxes%annsum_npp_patch => carbonflux_vars%annsum_npp_patch
    this%betr%carbonfluxes%agnpp_patch => carbonflux_vars%agnpp_patch
    this%betr%carbonfluxes%bgnpp_patch => carbonflux_vars%bgnpp_patch

    betr_pft%wtcol => pft%wtcol
    betr_pft%column => pft%column
    betr_pft%itype => pft%itype
    betr_pft%landunit => pft%landunit
    betr_pft%crop => crop

    betr_col%landunit => col%landunit
    betr_col%gridcell => col%gridcell
    betr_col%snl => col%snl
    betr_col%dz => col%dz
    betr_col%zi => col%zi
    betr_col%z => col%z

    betr_lun%itype => lun%itype
    betr_lun%ifspecial => lun%ifspecial

    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    !assign waterstate

    betr_waterstate_vars%h2osoi_liq_col    => waterstate_vars%h2osoi_liq_col
    betr_waterstate_vars%h2osoi_ice_col    => waterstate_vars%h2osoi_ice_col

    betr_waterstate_vars%h2osoi_liq_old    => waterstate_vars%h2osoi_liq_old
    betr_waterstate_vars%h2osoi_ice_old    => waterstate_vars%h2osoi_ice_old
    betr_waterstate_vars%h2osoi_liqvol_col => waterstate_vars%h2osoi_liqvol_col
    betr_waterstate_vars%h2osoi_icevol_col => waterstate_vars%h2osoi_icevol_col
    betr_waterstate_vars%h2osoi_vol_col    => waterstate_vars%h2osoi_vol_col
    betr_waterstate_vars%air_vol_col       => waterstate_vars%air_vol_col
    betr_waterstate_vars%finundated_col    => waterstate_vars%finundated_col
    betr_waterstate_vars%rho_vap           => waterstate_vars%rho_vap
    betr_waterstate_vars%rhvap_soi         => waterstate_vars%rhvap_soi
    betr_waterstate_vars%smp_l_col         => waterstate_vars%smp_l_col
    betr_waterstate_vars%frac_h2osfc_col   => waterstate_vars%frac_h2osfc_col

    betr_waterflux_vars%qflx_adv_col       => waterflux_vars%qflx_adv_col
    betr_waterflux_vars%qflx_infl_col      => waterflux_vars%qflx_infl_col
    betr_waterflux_vars%qflx_surf_col      => waterflux_vars%qflx_surf_col
    betr_waterflux_vars%qflx_rootsoi       => waterflux_vars%qflx_rootsoi
    betr_waterflux_vars%qflx_gross_evap_soil_col  => waterflux_vars%qflx_gross_evap_soil_col
    betr_waterflux_vars%qflx_gross_infl_soil_col  => waterflux_vars%qflx_gross_infl_soil_col
    betr_waterflux_vars%qflx_rootsoi_col        => waterflux_vars%qflx_rootsoi_col
    betr_waterflux_vars%qflx_drain_vr_col       => waterflux_vars%qflx_drain_vr_col
    betr_waterflux_vars%qflx_totdrain_col       => waterflux_vars%qflx_totdrain_col
    betr_waterflux_vars%qflx_dew_grnd_col       => waterflux_vars%qflx_dew_grnd_col
    betr_waterflux_vars%qflx_dew_snow_col       => waterflux_vars%qflx_dew_snow_col
    betr_waterflux_vars%qflx_sub_snow_vol_col   => waterflux_vars%qflx_sub_snow_vol_col
    betr_waterflux_vars%qflx_sub_snow_col       => waterflux_vars%qflx_sub_snow_col
    betr_waterflux_vars%qflx_h2osfc2topsoi_col  => waterflux_vars%qflx_h2osfc2topsoi_col
    betr_waterflux_vars%qflx_snow2topsoi_col    => waterflux_vars%qflx_snow2topsoi_col
    betr_waterflux_vars%qflx_tran_veg_patch     => waterflux_vars%qflx_tran_veg_patch

    betr_temperature_vars%t_soisno_col        => temperature_vars%t_soisno_col
    betr_temperature_vars%t_soi_10cm          => temperature_vars%t_soi_10cm
    betr_temperature_vars%t_veg_patch         => temperature_vars%t_veg_patch

    betr_soilhydrology_vars%fracice_col    => soilhydrology_vars%fracice_col
    betr_soilhydrology_vars%zwts_col       => soilhydrology_vars%zwts_col
    betr_soilhydrology_vars%qflx_bot_col   => soilhydrology_vars%qflx_bot_col

    betr_atm2lnd_vars%forc_pbot_downscaled_col => atm2lnd_vars%forc_pbot_downscaled_col

    betr_canopystate_vars%altmax_col      => canopystate_vars%altmax_col
    betr_canopystate_vars%altmax_lastyear_col   => canopystate_vars%altmax_lastyear_col
    betr_canopystate_vars%lbl_rsc_h2o_patch     => canopystate_vars%lbl_rsc_h2o_patch
    betr_canopystate_vars%elai_patch      => canopystate_vars%elai_patch

    betr_chemstate_vars%soil_pH => chemstate_vars%soil_pH

    betr_soilstate_vars%bsw_col  => soilstate_vars%bsw_col
    betr_soilstate_vars%watsat_col => soilstate_vars%watsat_col
    betr_soilstate_vars%eff_porosity_col => soilstate_vars%eff_porosity_col
    betr_soilstate_vars%soilpsi_col  => soilstate_vars%soilpsi_col
    betr_soilstate_vars%cellorg_col  => soilstate_vars%cellorg_col
    betr_soilstate_vars%cellclay_col  => soilstate_vars%cellclay_col
    betr_soilstate_vars%cellsand_col  => soilstate_vars%cellsand_col
    betr_soilstate_vars%bd_col   => soilstate_vars%bd_col
    betr_soilstate_vars%watfc_col  => soilstate_vars%watfc_col
    betr_soilstate_vars%sucsat_col => soilstate_vars%sucsat_col
    betr_soilstate_vars%rootfr_patch => soilstate_vars%rootfr_patch

    call this%betr%step_without_drainage(betr_bounds,   &
         this%num_soilc, this%filter_soilc, this%num_soilp, this%filter_soilp,  &
         betr_atm2lnd_vars, betr_soilhydrology_vars, betr_soilstate_vars, &
         betr_waterstate_vars, betr_waterflux_vars, betr_temperature_vars,  &
         betr_chemstate_vars, this%betr%cnstates, betr_canopystate_vars, &
         this%betr%carbonfluxes)


  end subroutine CLMStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine CLMStepWithDrainage(this, bounds, &
       waterflux_vars, col)

    use ColumnType, only : column_type
    use MathfuncMod, only : safe_div
    use WaterFluxType, only : waterflux_type

    use betr_decompMod, only : betr_bounds_type
    use LandunitType,only : lun
    use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use BeTR_WaterfluxType    , only : betr_waterflux_type
    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil

    implicit none

    ! !ARGUMENTS:
    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(waterflux_type), intent(in) :: waterflux_vars
    type(column_type), intent(in) :: col ! column type

    type(betr_waterflux_type)  :: betr_waterflux_vars
    type(betr_bounds_type)     :: betr_bounds
    integer  :: lbj, ubj ! lower and upper bounds, make sure they are > 0

    betr_nlevsoi = nlevsoi
    betr_nlevsno = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    betr_col%landunit => col%landunit
    betr_col%gridcell => col%gridcell
    betr_col%snl => col%snl
    betr_col%dz => col%dz
    betr_col%zi => col%zi
    betr_col%z => col%z

    betr_lun%itype => lun%itype
    betr_lun%ifspecial => lun%ifspecial

    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

    betr_waterflux_vars%qflx_adv_col       => waterflux_vars%qflx_adv_col
    betr_waterflux_vars%qflx_infl_col      => waterflux_vars%qflx_infl_col
    betr_waterflux_vars%qflx_surf_col      => waterflux_vars%qflx_surf_col
    betr_waterflux_vars%qflx_rootsoi       => waterflux_vars%qflx_rootsoi
    betr_waterflux_vars%qflx_gross_evap_soil_col  => waterflux_vars%qflx_gross_evap_soil_col
    betr_waterflux_vars%qflx_gross_infl_soil_col  => waterflux_vars%qflx_gross_infl_soil_col
    betr_waterflux_vars%qflx_rootsoi_col        => waterflux_vars%qflx_rootsoi_col
    betr_waterflux_vars%qflx_drain_vr_col       => waterflux_vars%qflx_drain_vr_col
    betr_waterflux_vars%qflx_totdrain_col       => waterflux_vars%qflx_totdrain_col
    betr_waterflux_vars%qflx_dew_grnd_col       => waterflux_vars%qflx_dew_grnd_col
    betr_waterflux_vars%qflx_dew_snow_col       => waterflux_vars%qflx_dew_snow_col
    betr_waterflux_vars%qflx_sub_snow_vol_col   => waterflux_vars%qflx_sub_snow_vol_col
    betr_waterflux_vars%qflx_sub_snow_col       => waterflux_vars%qflx_sub_snow_col
    betr_waterflux_vars%qflx_h2osfc2topsoi_col  => waterflux_vars%qflx_h2osfc2topsoi_col
    betr_waterflux_vars%qflx_snow2topsoi_col    => waterflux_vars%qflx_snow2topsoi_col
    betr_waterflux_vars%qflx_tran_veg_patch     => waterflux_vars%qflx_tran_veg_patch


    call this%betr%step_with_drainage(betr_bounds,   &
         this%num_soilc, this%filter_soilc, this%jtops, &
         betr_waterflux_vars)

  end subroutine CLMStepWithDrainage

  !------------------------------------------------------------------------

  subroutine diagnose_dtracer_freeze_thaw_clm(this, bounds, num_nolakec, filter_nolakec, col, lun, &
    waterstate_vars)
    !
    ! DESCRIPTION
    ! aqueous tracer partition based on freeze-thaw
    !
    ! USES
    !
    use ColumnType, only : column_type
    use LandunitType, only : landunit_type
    use BeTR_WaterStateType, only : betr_waterstate_type
    use landunit_varcon, only : istsoil
    use betr_decompMod    , only : betr_bounds_type
    use WaterStateType, only : waterstate_type
    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(betr_bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_nolakec ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:) ! column filter for non-lake points
    type(column_type), intent(in) :: col ! column type
    type(landunit_type), intent(in) :: lun
    type(waterstate_type), intent(in) :: waterstate_vars


    !temporary variables
    type(betr_waterstate_type)  :: betr_waterstate_vars
    type(betr_bounds_type)     :: betr_bounds

    betr_col%landunit => col%landunit
    betr_col%gridcell => col%gridcell
    betr_col%snl => col%snl
    betr_col%dz => col%dz
    betr_col%zi => col%zi
    betr_col%z => col%z

    betr_lun%itype => lun%itype
    betr_lun%ifspecial => lun%ifspecial

    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

    !assign waterstate
    betr_waterstate_vars%h2osoi_liq_col    => waterstate_vars%h2osoi_liq_col
    betr_waterstate_vars%h2osoi_ice_col    => waterstate_vars%h2osoi_ice_col

    betr_waterstate_vars%h2osoi_liq_old    => waterstate_vars%h2osoi_liq_old
    betr_waterstate_vars%h2osoi_ice_old    => waterstate_vars%h2osoi_ice_old
    betr_waterstate_vars%h2osoi_liqvol_col => waterstate_vars%h2osoi_liqvol_col
    betr_waterstate_vars%h2osoi_icevol_col => waterstate_vars%h2osoi_icevol_col
    betr_waterstate_vars%h2osoi_vol_col    => waterstate_vars%h2osoi_vol_col
    betr_waterstate_vars%air_vol_col       => waterstate_vars%air_vol_col
    betr_waterstate_vars%finundated_col    => waterstate_vars%finundated_col
    betr_waterstate_vars%rho_vap           => waterstate_vars%rho_vap
    betr_waterstate_vars%rhvap_soi         => waterstate_vars%rhvap_soi
    betr_waterstate_vars%smp_l_col         => waterstate_vars%smp_l_col
    betr_waterstate_vars%frac_h2osfc_col   => waterstate_vars%frac_h2osfc_col

    call this%betr%diagnose_dtracer_freeze_thaw(bounds, num_nolakec, filter_nolakec,  &
      betr_waterstate_vars)

  end subroutine diagnose_dtracer_freeze_thaw_clm

  !------------------------------------------------------------------------
  subroutine calc_dew_sub_flux_clm(this, bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars)

    ! External interface called by CLM

    use ColumnType, only : col
    use LandunitType, only : lun
    use BeTR_WaterfluxType, only : betr_waterflux_type
    use BeTR_WaterstateType, only : betr_waterstate_type
    use clm_varcon, only : denh2o,spval
    use landunit_varcon, only : istsoil, istcrop
    use betr_decompMod    , only : betr_bounds_type

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(betr_bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc ! number of column soil points in column filter_soilc
    integer, intent(in) :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
    type(betr_waterstate_type), intent(in) :: waterstate_vars
    type(betr_waterflux_type), intent(in) :: waterflux_vars

  call this%betr%calc_dew_sub_flux(bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars, this%betr%tracers, this%betr%tracerfluxes, this%betr%tracerstates)

  end subroutine calc_dew_sub_flux_clm

  !------------------------------------------------------------------------
  subroutine clm_betr_flux_state_receive(this, bounds, num_soilc, filter_soilc)

    use betr_decompMod    , only : betr_bounds_type

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(betr_bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)

    call this%betr%bgc_reaction%lsm_betr_flux_state_receive(bounds, &
         num_soilc, filter_soilc, &
         this%betr%tracerstates, this%betr%tracerfluxes,  this%betr%tracers)

  end subroutine clm_betr_flux_state_receive

  !------------------------------------------------------------------------
  subroutine calc_smp_l_clm(this, bounds, lbj, ubj, &
       numf, filter, t_soisno, soilstate_vars, waterstate_vars, soil_water_retention_curve)

    use SoilStateType, only : soilstate_type
    use WaterStateType, only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use betr_decompMod    , only : betr_bounds_type

    use clm_varcon, only : grav,hfus,tfrz

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(betr_bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer, intent(in) :: numf                                              ! number of columns in column filter
    integer, intent(in) :: filter(:)                                         ! column filter
    real(r8), intent(in) :: t_soisno(bounds%begc:, lbj: )                    ! soil temperature
    type(soilstate_type), intent(in) :: soilstate_vars
    type(waterstate_type), intent(inout) :: waterstate_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

    !local variables
    real(r8) :: s_node
    integer :: fc, c, j

    SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(mod_filename,__LINE__))


    associate(                                                     & !
         h2osoi_vol =>    waterstate_vars%h2osoi_vol_col, & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
         smp_l =>    waterstate_vars%smp_l_col, & ! Output: [real(r8) (:,:) ]  soil suction (mm)
         bsw =>    soilstate_vars%bsw_col, & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         watsat =>    soilstate_vars%watsat_col, & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         sucsat =>    soilstate_vars%sucsat_col            & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         )


      do j = lbj, ubj
         do fc = 1, numf
            c = filter(fc)
            if(j>=1)then

               if(t_soisno(c,j)<tfrz)then
                  smp_l(c,j)= hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
               else
                  s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
                  ! FIXME(bja, 201603) this depends on the SWRC
                  ! implemnted by the driving LSM. Doesn't agree with
                  ! stub version...
                  ! the following call is CLM specific, jyt, Mar 13, 2016
                  !Xcall soil_water_retention_curve%soil_suction(c, j, s_node, soilstate_vars, smp_l(c,j))
               endif

            endif
         enddo
      enddo
    end associate
  end subroutine calc_smp_l_clm

  !------------------------------------------------------------------------
  subroutine betr_clm_readParams(this, ncid)

    use ncdio_pio, only : file_desc_t

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(file_desc_t), intent(inout) :: ncid  ! pio netCDF file id
    call this%betr%bgc_reaction%readParams(ncid, this%betr%tracers)
  end subroutine betr_clm_readParams

  !------------------------------------------------------------------------
  subroutine betr_clm_h2oiso_consistency_check(this, &
       bounds, ubj, num_soilc, filter_soilc, waterstate_vars)
    !
    ! check the overall water mass consistency between betr and clm
    !use clm_instMod
    use MathfuncMod, only : dot_sum
    use clm_varcon, only : denh2o

    use WaterStateType, only : waterstate_type

    use clm_time_manager, only : get_nstep

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: ubj
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)
    type(waterstate_type), intent(in) :: waterstate_vars
    real(r8), allocatable :: eyev(:)
    integer :: fc, c
    real(r8):: totwater, err

    associate( &
         h2osoi_ice => waterstate_vars%h2osoi_ice_col, & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq => waterstate_vars%h2osoi_liq_col, & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)
         end_tracer_molarmass => this%betr%tracerstates%end_tracer_molarmass_col, &
         id_trc_o18_h2o => this%betr%tracers%id_trc_o18_h2o           &
         )

      allocate(eyev(1:ubj))
      eyev=1._r8
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         totwater=dot_sum(h2osoi_ice(c,1:ubj),eyev) + dot_sum(h2osoi_liq(c,1:ubj),eyev)
         err = totwater-end_tracer_molarmass(c,id_trc_o18_h2o)
         print*,get_nstep(),'diff',c, totwater, end_tracer_molarmass(c,id_trc_o18_h2o),err, err/totwater
      enddo
      deallocate(eyev)

    end associate
  end subroutine betr_clm_h2oiso_consistency_check
end module BeTRSimulationCLM
