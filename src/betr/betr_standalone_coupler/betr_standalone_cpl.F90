module betr_standalone_cpl

  use decompMod, only : bounds_type
  use betr_instMod
  use betr_initializeMod
 public :: betr_initialize_standalone
 public :: run_betr_one_step_without_drainage_standalone
 public :: run_betr_one_step_with_drainage_standalone
contains

  !-------------------------------------------------------------------------------
  subroutine run_betr_one_step_without_drainage_standalone(bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
     atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, temperature_vars, waterflux_vars, chemstate_vars, &
     cnstate_vars, canopystate_vars, carbonflux_vars)

     !
     ! !DESCRIPTION:
     ! run betr code one time step forward, without drainage calculation

     ! !USES:
     use BetrBGCMod, only : run_betr_one_step_without_drainage
     use SoilStateType, only : soilstate_type
     use WaterStateType, only : Waterstate_Type
     use TemperatureType, only : temperature_type
     use ChemStateType, only : chemstate_type
     use WaterfluxType, only : waterflux_type
     use ColumnType, only : column_type
     use BGCReactionsMod, only : bgc_reaction_type
     use atm2lndType, only : atm2lnd_type
     use SoilHydrologyType, only : soilhydrology_type
     use BeTR_CNStateType, only : betr_cnstate_type
     use BeTR_CarbonFluxType, only : betr_carbonflux_type
     use CNStateType, only : cnstate_type
     use CNCarbonFluxType, only : carbonflux_type
     use CanopyStateType, only : canopystate_type
     use BeTR_PatchType, only : betr_pft
     use BeTR_ColumnType,only : betr_col
     use PatchType, only : pft
     use pftvarcon, only : crop

     !
     ! !ARGUMENTS :
     type(bounds_type)                , intent(in)    :: bounds                     ! bounds
     integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
     integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
     integer                          , intent(in)    :: num_soilp
     integer                          , intent(in)    :: filter_soilp(:)            ! pft filter
     integer                          , intent(in)    :: lbj, ubj                   ! lower and upper bounds, make sure they are > 0

     type(column_type)                , intent(in)    :: col                        ! column type
     type(Waterstate_Type)            , intent(in)    :: waterstate_vars            ! water state variables
     type(soilstate_type)             , intent(in)    :: soilstate_vars             ! column physics variable
     type(temperature_type)           , intent(in)    :: temperature_vars           ! energy state variable
     type(chemstate_type)             , intent(in)    :: chemstate_vars
     type(atm2lnd_type)               , intent(in)    :: atm2lnd_vars
     type(soilhydrology_type)         , intent(in)    :: soilhydrology_vars
     type(cnstate_type)               , intent(inout) :: cnstate_vars
     type(canopystate_type)           , intent(in)    :: canopystate_vars
     type(carbonflux_type)            , intent(in)    :: carbonflux_vars
     type(waterflux_type)             , intent(inout) :: waterflux_vars



     !temporary variables
     type(betr_cnstate_type)        :: betr_cnstate_vars
     type(betr_carbonflux_type)     :: betr_carbonflux_vars

     !pass necessary data for correct subroutine call

     betr_cnstate_vars%isoilorder          => cnstate_vars%isoilorder

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

     call run_betr_one_step_without_drainage(bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp,    &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, temperature_vars, waterflux_vars, chemstate_vars, &
       betr_cnstate_vars, canopystate_vars, betr_carbonflux_vars, betrtracer_vars, bgc_reaction, betr_aerecond_vars,   &
       tracerboundarycond_vars, tracercoeff_vars, tracerstate_vars, tracerflux_vars, plant_soilbgc)

  end subroutine run_betr_one_step_without_drainage_standalone

  !-------------------------------------------------------------------------------
  subroutine betr_initialize_standalone(bounds, lbj, ubj)

  use clm_instMod
  use betr_initializeMod, only : betr_initialize
  use BeTR_CNStateType, only : betr_cnstate_type
  use EcophysConType, only : ecophyscon
  use BeTR_PatchType, only : betr_pft
  use PatchType, only : pft
  use pftvarcon, only : nc3_arctic_grass, crop, nc3_nonarctic_grass, nc4_grass, noveg
  use BeTR_pftvarconType, only : betr_pftvarcon
  implicit none
  type(bounds_type)    , intent(in) :: bounds
  integer              , intent(in) :: lbj, ubj

  !temporary variables
  type(betr_cnstate_type)   :: betr_cnstate_vars

  betr_pft%wtcol                        => pft%wtcol
  betr_pft%column                       => pft%column
  betr_pft%itype                        => pft%itype
  betr_pft%landunit                     => pft%landunit

  call betr_initialize(bounds, lbj, ubj, waterstate_vars)

  !pass necessary data
  betr_cnstate_vars%isoilorder          => cnstate_vars%isoilorder
  betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
  betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
  betr_pftvarcon%nc4_grass           = nc4_grass
  betr_pftvarcon%noveg               = noveg

  call bgc_reaction%init_betr_lsm_bgc_coupler(bounds, plant_soilbgc, &
       betrtracer_vars, tracerstate_vars, betr_cnstate_vars, ecophyscon)

  end subroutine betr_initialize_standalone


  !-------------------------------------------------------------------------------
  subroutine run_betr_one_step_with_drainage_standalone(bounds, lbj, ubj, num_soilc, filter_soilc, &
       jtops, waterflux_vars, col)

  ! !USES:

  use BetrBGCMod, only : run_betr_one_step_with_drainage
  use ColumnType, only : column_type
  use MathfuncMod, only : safe_div
  use WaterFluxType, only : waterflux_type
  use BeTR_ColumnType,only : betr_col
  implicit none
  ! !ARGUMENTS:
  type(bounds_type),        intent(in)    :: bounds
  integer,                  intent(in)    :: lbj, ubj
  integer,                  intent(in)    :: num_soilc                          ! number of columns in column filter_soilc
  integer,                  intent(in)    :: filter_soilc(:)                    ! column filter_soilc
  integer,                  intent(in)    :: jtops(bounds%begc: )
  type(waterflux_type)    , intent(in)    :: waterflux_vars
  type(column_type),        intent(in)    :: col                                ! column type


  betr_col%landunit                     => col%landunit
  betr_col%gridcell                     => col%gridcell
  betr_col%snl                          => col%snl
  betr_col%dz                           => col%dz
  betr_col%zi                           => col%zi
  betr_col%z                            => col%z

  call run_betr_one_step_with_drainage(bounds, lbj, ubj, num_soilc, filter_soilc, &
       jtops, waterflux_vars, betrtracer_vars, tracercoeff_vars, tracerstate_vars,  tracerflux_vars)

  end subroutine run_betr_one_step_with_drainage_standalone
end module betr_standalone_cpl
