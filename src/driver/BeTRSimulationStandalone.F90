module BeTRSimulationStandalone
  !
  ! !DESCRIPTION:
  !  BeTR standalone simulation class.
  !
  !  BeTR simulation class are API definitions, mapping data
  !  structures from a specific LSM, e.g. CLM, ALM, into BeTR data
  !  structures. The standalone class use BeTR data structures
  !  natively, so is mapping BeTR to BeTR. This class shouldn't be
  !  doing much.
  !
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use betr_decompMod    , only : betr_bounds_type
  use decompMod, only : bounds_type
  use BeTRSimulation, only : betr_simulation_type
  use BeTR_CNStateType, only : betr_cnstate_type
  use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use EcophysConType, only : ecophyscon_type

  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_standalone_type
     ! NOTE(bja, 201603) LSM specific types here!
     type(ecophyscon_type) :: ecophyscon

     ! NOTE(bja, 201603) most (all?) BeTR types go into the base
     ! class.

   contains
     procedure, public :: Init => StandaloneInit
     procedure, public :: StepWithoutDrainage => StandaloneStepWithoutDrainage
     procedure, public :: StepWithDrainage => StandaloneStepWithDrainage
  end type betr_simulation_standalone_type

  public :: create_betr_simulation_standalone

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_standalone()

    implicit none

    class(betr_simulation_standalone_type), pointer :: create_betr_simulation_standalone
    class(betr_simulation_standalone_type), pointer :: simulation

    allocate(simulation)

    create_betr_simulation_standalone => simulation

  end function create_betr_simulation_standalone

  !-------------------------------------------------------------------------------

  subroutine StandaloneInit(this, reaction_method, bounds, waterstate)

    use BeTRSimulation, only : BeTRSimulationInit
    use ReactionsFactory, only : create_bgc_reaction_type, &
         create_plant_soilbgc_type

    use BeTR_PatchType, only : betr_pft
    use BeTR_ColumnType, only : betr_col
    use BeTR_LandunitType,only : betr_lun
    use BeTR_pftvarconType, only : betr_pftvarcon
    use PatchType, only : pft
    use ColumnType, only : col
    use LandunitType,only : lun
    use pftvarcon, only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use clm_instMod, only : cnstate_vars
    use WaterStateType, only : waterstate_type
    use BeTR_WaterStateType, only : betr_waterstate_type
    use landunit_varcon
    use BeTR_landvarconType, only : betr_landvarcon
    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil

    implicit none

    class(betr_simulation_standalone_type), intent(inout) :: this
    character(len=*), intent(in) :: reaction_method
    type(bounds_type)    , intent(in) :: bounds

    type(waterstate_type), intent(inout) :: waterstate

    type(betr_bounds_type)     :: betr_bounds
    integer  :: lbj, ubj


    betr_nlevsoi = nlevsoi
    betr_nlevsno = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    !pass necessary data for correct subroutine call
    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg

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

    ! allocate the reaction types that may only be known to this
    ! simulation type.
    allocate(this%betr%bgc_reaction, source=create_bgc_reaction_type(reaction_method))
    allocate(this%betr%plant_soilbgc, source=create_plant_soilbgc_type(reaction_method))

    ! now call the base simulation init to continue initialization
    call BeTRSimulationInit(this, reaction_method, bounds, waterstate)

    !pass necessary data
    this%betr%cnstates%isoilorder  => cnstate_vars%isoilorder

    call this%betr%bgc_reaction%init_betr_lsm_bgc_coupler(betr_bounds, this%betr%plant_soilbgc, &
         this%betr%tracers, this%betr%tracerstates, this%betr%cnstates, &
         this%ecophyscon)

  end subroutine StandaloneInit


  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithoutDrainage(this, bounds,  &
       num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)

    use SoilStateType, only : soilstate_type
    use WaterStateType, only : Waterstate_Type
    use TemperatureType, only : temperature_type
    use ChemStateType, only : chemstate_type
    use WaterfluxType, only : waterflux_type
    use ColumnType, only : column_type
    use atm2lndType, only : atm2lnd_type
    use SoilHydrologyType, only : soilhydrology_type
    use CNStateType, only : cnstate_type
    use CNCarbonFluxType, only : carbonflux_type
    use CanopyStateType, only : canopystate_type

    use BGCReactionsMod, only : bgc_reaction_type
    use BeTR_CarbonFluxType, only : betr_carbonflux_type
    use BeTR_PatchType, only : betr_pft
    use BeTR_ColumnType, only : betr_col
    use BeTR_LandunitType,only : betr_lun
    use BeTR_WaterstateType   , only : betr_waterstate_type
    use BeTR_WaterfluxType    , only : betr_waterflux_type
    use PatchType, only : pft
    use LandunitType,only : lun
    use pftvarcon, only : crop

    implicit none

    class(betr_simulation_standalone_type), intent(inout) :: this

    ! !ARGUMENTS :
    type(bounds_type), intent(in) :: bounds ! bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: num_soilp
    integer, intent(in) :: filter_soilp(:) ! pft filter

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
    type(betr_carbonflux_type) :: betr_carbonflux_vars
    type(betr_waterflux_type)  :: betr_waterflux_vars
    type(betr_waterstate_type)  :: betr_waterstate_vars
    type(betr_bounds_type)     :: betr_bounds
    integer  :: lbj, ubj ! lower and upper bounds, make sure they are > 0

    !pass necessary data for correct subroutine call
    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    this%betr%cnstates%isoilorder          => cnstate_vars%isoilorder

    betr_carbonflux_vars%annsum_npp_patch => carbonflux_vars%annsum_npp_patch
    betr_carbonflux_vars%agnpp_patch      => carbonflux_vars%agnpp_patch
    betr_carbonflux_vars%bgnpp_patch      => carbonflux_vars%bgnpp_patch

    betr_pft%wtcol                        => pft%wtcol
    betr_pft%column                       => pft%column
    betr_pft%itype                        => pft%itype
    betr_pft%landunit                     => pft%landunit
    betr_pft%crop                         => crop

    betr_col%landunit                     => col%landunit
    betr_col%gridcell                     => col%gridcell
    betr_col%snl                          => col%snl
    betr_col%dz                           => col%dz
    betr_col%zi                           => col%zi
    betr_col%z                            => col%z

    betr_lun%itype                        => lun%itype
    betr_lun%ifspecial                    => lun%ifspecial

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

    call this%betr%step_without_drainage(betr_bounds, lbj, ubj, &
         num_soilc, filter_soilc, num_soilp, filter_soilp,  &
         atm2lnd_vars, soilhydrology_vars, soilstate_vars, &
         betr_waterstate_vars, temperature_vars, betr_waterflux_vars, &
         chemstate_vars, this%betr%cnstates, canopystate_vars, &
         betr_carbonflux_vars, this%betr%tracers, this%betr%bgc_reaction, &
         this%betr_aerecond_vars, this%betr%tracerboundaryconds, this%betr%tracercoeffs, &
         this%betr%tracerstates, this%betr%tracerfluxes, this%betr%plant_soilbgc)

  end subroutine StandaloneStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithDrainage(this, bounds,  &
       num_soilc, filter_soilc, jtops, waterflux_vars, col)

    use ColumnType    , only : column_type
    use MathfuncMod   , only : safe_div
    use WaterFluxType , only : waterflux_type
    use BeTR_ColumnType, only : betr_col
    use BeTR_LandunitType,only : betr_lun
    use LandunitType,only : lun
    use BeTR_WaterfluxType    , only : betr_waterflux_type
    implicit none

    ! !ARGUMENTS:
    class(betr_simulation_standalone_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: jtops(bounds%begc: )
    type(waterflux_type)    , intent(in) :: waterflux_vars
    type(column_type), intent(in) :: col ! column type

    type(betr_waterflux_type)  :: betr_waterflux_vars
    type(betr_bounds_type)     :: betr_bounds
    integer  :: lbj, ubj ! lower and upper bounds, make sure they are > 0

    !pass necessary data for correct subroutine call
    !set lbj and ubj
    betr_bounds%lbj  = 1          ; betr_bounds%ubj  = betr_nlevsoi
    betr_bounds%begp = bounds%begp; betr_bounds%endp = bounds%endp
    betr_bounds%begc = bounds%begc; betr_bounds%endc = bounds%endc
    betr_bounds%begl = bounds%begl; betr_bounds%endl = bounds%endl
    betr_bounds%begg = bounds%begg; betr_bounds%endg = bounds%endg
    lbj = betr_bounds%lbj; ubj = betr_bounds%ubj

    betr_col%landunit                     => col%landunit
    betr_col%gridcell                     => col%gridcell
    betr_col%snl                          => col%snl
    betr_col%dz                           => col%dz
    betr_col%zi                           => col%zi
    betr_col%z                            => col%z

    betr_lun%itype                        => lun%itype
    betr_lun%ifspecial                    => lun%ifspecial

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

    call this%betr%step_with_drainage(betr_bounds, lbj, ubj, &
         num_soilc, filter_soilc, &
         jtops, betr_waterflux_vars, this%betr%tracers, this%betr%tracercoeffs, &
         this%betr%tracerstates,  this%betr%tracerfluxes)

  end subroutine StandaloneStepWithDrainage


end module BeTRSimulationStandalone
