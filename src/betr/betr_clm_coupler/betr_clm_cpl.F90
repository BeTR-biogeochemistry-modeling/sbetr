module betr_clm_cpl
#include "shr_assert.h"
  use decompMod          , only : bounds_type
  use betr_instMod
  use clm_time_manager   , only : get_nstep
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use clm_varpar         , only : nlevtrc_soil
implicit none
 public :: betr_initialize_clm
 public :: run_betr_one_step_without_drainage_clm
 public :: run_betr_one_step_with_drainage_clm
 public :: diagnose_dtracer_freeze_thaw_clm
 public :: calc_dew_sub_flux_clm
 public :: betr_clm_flux_statevar_feedback
 public :: calc_smp_l_clm
contains

  !-------------------------------------------------------------------------------
  subroutine run_betr_one_step_without_drainage_clm(bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
     atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, temperature_vars, waterflux_vars, chemstate_vars, &
     cnstate_vars, canopystate_vars, carbonflux_vars)

     !
     ! !DESCRIPTION:
     ! run betr code one time step forward, without drainage calculation

     ! !USES:
     use BetrBGCMod                   , only          : run_betr_one_step_without_drainage
     use SoilStateType                , only          : soilstate_type
     use WaterStateType               , only          : Waterstate_Type
     use TemperatureType              , only          : temperature_type
     use ChemStateType                , only          : chemstate_type
     use WaterfluxType                , only          : waterflux_type
     use ColumnType                   , only          : column_type
     use atm2lndType                  , only          : atm2lnd_type
     use SoilHydrologyType            , only          : soilhydrology_type
     use BeTR_CNStateType             , only          : betr_cnstate_type
     use BeTR_CarbonFluxType          , only          : betr_carbonflux_type
     use CNVegStateType               , only          : cnstate_type => cnveg_state_type
     use CNVegCarbonFluxType          , only          : carbonflux_type  => cnveg_carbonflux_type
     use CanopyStateType              , only          : canopystate_type
     use BeTR_PatchType               , only          : betr_pft
     use PatchType                    , only          : pft => patch
     use pftconMod                    , only          : pftcon


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
     betr_pft%crop                         => pftcon%crop

     call run_betr_one_step_without_drainage(bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp, col ,   &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, temperature_vars, waterflux_vars, chemstate_vars, &
       betr_cnstate_vars, canopystate_vars, betr_carbonflux_vars, betrtracer_vars, bgc_reaction, betr_aerecond_vars,   &
       tracerboundarycond_vars, tracercoeff_vars, tracerstate_vars, tracerflux_vars, plant_soilbgc)

  end subroutine run_betr_one_step_without_drainage_clm

  !-------------------------------------------------------------------------------
  subroutine betr_initialize_clm(bounds, lbj, ubj)

  use clm_instMod
  use betr_initializeMod, only : betr_initialize
  use pftconMod         , only : nc3_arctic_grass, nc3_nonarctic_grass, nc4_grass, noveg
  use BeTR_CNStateType  , only : betr_cnstate_type
  use EcophysConType    , only : ecophyscon
  use BeTR_PatchType    , only : betr_pft
  use BeTR_pftvarconType, only : betr_pftvarcon
  use PatchType         , only : pft => patch
  implicit none
  type(bounds_type)    , intent(in) :: bounds
  integer              , intent(in) :: lbj, ubj

  !temporary variables
  type(betr_cnstate_type)   :: betr_cnstate_vars

  betr_pft%wtcol                        => pft%wtcol
  betr_pft%column                       => pft%column
  betr_pft%itype                        => pft%itype
  betr_pft%landunit                     => pft%landunit

  call betr_initialize(bounds, lbj, ubj, waterstate_inst)

  !pass necessary data
  betr_cnstate_vars%isoilorder          => cnveg_state_inst%isoilorder
  betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
  betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
  betr_pftvarcon%nc4_grass           = nc4_grass
  betr_pftvarcon%noveg               = noveg

  call bgc_reaction%init_betr_lsm_bgc_coupler(bounds, plant_soilbgc, &
       betrtracer_vars, tracerstate_vars, betr_cnstate_vars, ecophyscon)

  end subroutine betr_initialize_clm


  !-------------------------------------------------------------------------------
  subroutine run_betr_one_step_with_drainage_clm(bounds, lbj, ubj, num_soilc, filter_soilc, &
       jtops, waterflux_vars, col)

  ! !USES:

  use BetrBGCMod            , only : run_betr_one_step_with_drainage
  use ColumnType            , only : column_type
  use MathfuncMod           , only : safe_div
  use WaterFluxType         , only : waterflux_type
  implicit none
  ! !ARGUMENTS:
  type(bounds_type),        intent(in)    :: bounds
  integer,                  intent(in)    :: lbj, ubj
  integer,                  intent(in)    :: num_soilc                          ! number of columns in column filter_soilc
  integer,                  intent(in)    :: filter_soilc(:)                    ! column filter_soilc
  integer,                  intent(in)    :: jtops(bounds%begc: )
  type(waterflux_type)    , intent(in)    :: waterflux_vars
  type(column_type),        intent(in)    :: col                                ! column type

  call run_betr_one_step_with_drainage(bounds, lbj, ubj, num_soilc, filter_soilc, &
       jtops, waterflux_vars, col, betrtracer_vars, tracercoeff_vars, tracerstate_vars,  tracerflux_vars)

  end subroutine run_betr_one_step_with_drainage_clm

  !------------------------------------------------------------------------

  subroutine diagnose_dtracer_freeze_thaw_clm(bounds, num_nolakec, filter_nolakec, col, lun, &
    waterstate_vars)
  !
  ! DESCRIPTION
  ! aqueous tracer partition based on freeze-thaw
  !
  ! USES
  !
  use ColumnType            , only : column_type
  use LandunitType          , only : landunit_type
  use WaterStateType        , only : waterstate_type
  use landunit_varcon       , only : istsoil
  use TracerParamsMod       , only : diagnose_dtracer_freeze_thaw
  !
  ! Arguments
  implicit none
  type(bounds_type)      , intent(in)    :: bounds
  integer                , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
  integer                , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points
  type(landunit_type)    , intent(in)    :: lun
  type(waterstate_type)  , intent(in)    :: waterstate_vars
  type(column_type)      , intent(in)    :: col                                ! column type


  call diagnose_dtracer_freeze_thaw(bounds, num_nolakec, filter_nolakec, lun, &
    waterstate_vars, betrtracer_vars, tracerstate_vars)

  end subroutine diagnose_dtracer_freeze_thaw_clm


  !------------------------------------------------------------------------
  subroutine calc_dew_sub_flux_clm(bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars)

  use ColumnType            , only : col
  use LandunitType          , only : lun
  use WaterfluxType         , only : waterflux_type
  use BetrBGCMod            , only : calc_dew_sub_flux, check_mass_err
  use WaterstateType        , only : waterstate_type
  use clm_varcon            , only : denh2o,spval
  use landunit_varcon       , only : istsoil, istcrop

  ! !ARGUMENTS:
  type(bounds_type)         , intent(in)    :: bounds
  integer                   , intent(in)    :: num_hydrologyc             ! number of column soil points in column filter_soilc
  integer                   , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
  type(waterstate_type)     , intent(in)    :: waterstate_vars
  type(waterflux_type)      , intent(in)    :: waterflux_vars


  call calc_dew_sub_flux(bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars, betrtracer_vars, tracerflux_vars, tracerstate_vars)


  end subroutine calc_dew_sub_flux_clm

  !------------------------------------------------------------------------
  subroutine betr_clm_flux_statevar_feedback(bounds, num_soilc, filter_soilc)

  implicit none
  ! !ARGUMENTS:
  type(bounds_type)         , intent(in)   :: bounds
  integer                   , intent(in)   :: num_soilc
  integer                   , intent(in)   :: filter_soilc(:)

  call bgc_reaction%betr_lsm_flux_statevar_feedback(bounds, &
       num_soilc, filter_soilc,                    &
       tracerstate_vars, tracerflux_vars,  betrtracer_vars)

  end subroutine betr_clm_flux_statevar_feedback



  !------------------------------------------------------------------------
  subroutine calc_smp_l_clm(bounds, lbj, ubj, numf, filter, t_soisno, soilstate_vars, waterstate_vars, soil_water_retention_curve)

  use SoilStateType              , only : soilstate_type
  use WaterStateType             , only : waterstate_type
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
  use clm_varcon                 , only : grav,hfus,tfrz

  implicit none
  type(bounds_type)         , intent(in)    :: bounds  ! bounds
  integer                   , intent(in)    :: lbj, ubj                                          ! lower and upper bounds, make sure they are > 0
  integer                   , intent(in)    :: numf                                              ! number of columns in column filter
  integer                   , intent(in)    :: filter(:)                                         ! column filter
  real(r8)                  , intent(in)    :: t_soisno(bounds%begc: , lbj: )                    ! soil temperature
  type(soilstate_type)      , intent(in)    :: soilstate_vars
  type(waterstate_type)     , intent(inout) :: waterstate_vars
  class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

  !local variables
  real(r8) :: s_node
  integer  :: fc, c, j

  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(__FILE__,__LINE__))


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
          call soil_water_retention_curve%soil_suction(c,j, s_node, soilstate_vars, smp_l(c,j))
        endif

      endif
    enddo
  enddo
  end associate
  end subroutine calc_smp_l_clm


  !------------------------------------------------------------------------
  subroutine h2oiso_consistency_check()
  implicit none

   
  end subroutine h2oiso_consistency_check
end module betr_clm_cpl
