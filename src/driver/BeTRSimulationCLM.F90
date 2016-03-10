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

     procedure, public :: diagnose_dtracer_freeze_thaw_clm
     procedure, public :: calc_dew_sub_flux_clm
     procedure, public :: clm_betr_flux_state_receive
     procedure, public :: calc_smp_l_clm
     procedure, public :: betr_clm_h2oiso_consistency_check
  end type betr_simulation_clm_type

  public :: create_betr_simulation_clm
  
contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_clm()

    implicit none

    class(betr_simulation_clm_type), pointer :: create_betr_simulation_clm
    class(betr_simulation_clm_type), pointer :: simulation

    allocate(simulation)

    call simulation%Init()

    create_betr_simulation_clm => simulation

  end function create_betr_simulation_clm

  
  !-------------------------------------------------------------------------------

  subroutine CLMInit(this, reaction_method, bounds, lbj, ubj, waterstate)

    use BeTRSimulation, only : BeTRSimulationInit
    use ReactionsFactoryStandalone, only : create_standalone_bgc_reaction_type, &
         create_standalone_plant_soilbgc_type

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

    use landunit_varcon
    use BeTR_landvarconType, only : betr_landvarcon
    use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    character(len=*), intent(in) :: reaction_method
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: lbj, ubj
    type(waterstate_type), intent(inout) :: waterstate

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

    ! allocate the reaction types that may only be known to this
    ! simulation type.
    allocate(this%bgc_reaction, source=create_standalone_bgc_reaction_type(reaction_method))
    allocate(this%plant_soilbgc, source=create_standalone_plant_soilbgc_type(reaction_method))

    ! now call the base simulation init to continue initialization
    call BeTRSimulationInit(this, reaction_method, bounds, lbj, ubj, waterstate)

    !pass necessary data
    this%betr_cnstate_vars%isoilorder  => cnstate_vars%isoilorder

    call this%bgc_reaction%init_betr_lsm_bgc_coupler(bounds, this%plant_soilbgc, &
         this%betrtracer_vars, this%tracerstate_vars, this%betr_cnstate_vars, &
         this%ecophyscon)

  end subroutine CLMInit


  !---------------------------------------------------------------------------------
  subroutine CLMStepWithoutDrainage(this, bounds, lbj, ubj, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, col,   &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)

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
    use BeTR_CarbonFluxType, only : betr_carbonflux_type
    use CNStateType, only : cnstate_type
    use CNCarbonFluxType, only : carbonflux_type
    use CanopyStateType, only : canopystate_type
    use BeTR_PatchType, only : betr_pft
    use BeTR_ColumnType, only : betr_col
    use BeTR_LandunitType,only : betr_lun

    use PatchType, only : pft
    use LandunitType,only : lun
    use pftvarcon, only : crop

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this

    ! !ARGUMENTS :
    type(bounds_type), intent(in) :: bounds ! bounds
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: num_soilp
    integer, intent(in) :: filter_soilp(:) ! pft filter
    integer, intent(in) :: lbj, ubj ! lower and upper bounds, make sure they are > 0

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

    !pass necessary data for correct subroutine call

    this%betr_cnstate_vars%isoilorder => cnstate_vars%isoilorder

    betr_carbonflux_vars%annsum_npp_patch => carbonflux_vars%annsum_npp_patch
    betr_carbonflux_vars%agnpp_patch => carbonflux_vars%agnpp_patch
    betr_carbonflux_vars%bgnpp_patch => carbonflux_vars%bgnpp_patch

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

    call run_betr_one_step_without_drainage(bounds, lbj, ubj, &
         num_soilc, filter_soilc, num_soilp, filter_soilp,  &
         atm2lnd_vars, soilhydrology_vars, soilstate_vars, &
         waterstate_vars, temperature_vars, waterflux_vars, &
         chemstate_vars, this%betr_cnstate_vars, canopystate_vars, &
         betr_carbonflux_vars, this%betrtracer_vars, this%bgc_reaction, &
         this%betr_aerecond_vars, this%tracerboundarycond_vars, this%tracercoeff_vars, &
         this%tracerstate_vars, this%tracerflux_vars, this%plant_soilbgc)

  end subroutine CLMStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine CLMStepWithDrainage(this, bounds, lbj, ubj, &
       num_soilc, filter_soilc, jtops, waterflux_vars, col)

    use BetrBGCMod, only : run_betr_one_step_with_drainage
    use ColumnType, only : column_type
    use MathfuncMod, only : safe_div
    use WaterFluxType, only : waterflux_type
    use BeTR_ColumnType, only : betr_col
    use BeTR_LandunitType,only : betr_lun
    use LandunitType,only : lun
    implicit none

    ! !ARGUMENTS:
    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: lbj, ubj
    integer, intent(in) :: num_soilc ! number of columns in column filter_soilc
    integer, intent(in) :: filter_soilc(:) ! column filter_soilc
    integer, intent(in) :: jtops(bounds%begc: )
    type(waterflux_type), intent(in) :: waterflux_vars
    type(column_type), intent(in) :: col ! column type


    betr_col%landunit => col%landunit
    betr_col%gridcell => col%gridcell
    betr_col%snl => col%snl
    betr_col%dz => col%dz
    betr_col%zi => col%zi
    betr_col%z => col%z

    betr_lun%itype => lun%itype
    betr_lun%ifspecial => lun%ifspecial

    call run_betr_one_step_with_drainage(bounds, lbj, ubj, &
         num_soilc, filter_soilc, &
         jtops, waterflux_vars, this%betrtracer_vars, this%tracercoeff_vars, &
         this%tracerstate_vars,  this%tracerflux_vars)

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
    use WaterStateType, only : waterstate_type
    use landunit_varcon, only : istsoil
    use TracerParamsMod, only : diagnose_dtracer_freeze_thaw

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_nolakec ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:) ! column filter for non-lake points
    type(landunit_type), intent(in) :: lun
    type(waterstate_type), intent(in) :: waterstate_vars
    type(column_type), intent(in) :: col ! column type


    call diagnose_dtracer_freeze_thaw(bounds, num_nolakec, filter_nolakec, lun, &
         waterstate_vars, betrtracer_vars, tracerstate_vars)

  end subroutine diagnose_dtracer_freeze_thaw_clm

  !------------------------------------------------------------------------
  subroutine calc_dew_sub_flux_clm(this, bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars)

    use ColumnType, only : col
    use LandunitType, only : lun
    use WaterfluxType, only : waterflux_type
    use BetrBGCMod, only : calc_dew_sub_flux, check_mass_err
    use WaterstateType, only : waterstate_type
    use clm_varcon, only : denh2o,spval
    use landunit_varcon, only : istsoil, istcrop

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc ! number of column soil points in column filter_soilc
  integer, intent(in) :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
  type(waterstate_type), intent(in) :: waterstate_vars
  type(waterflux_type), intent(in) :: waterflux_vars

  call calc_dew_sub_flux(bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars, betrtracer_vars, tracerflux_vars, tracerstate_vars)

  end subroutine calc_dew_sub_flux_clm

  !------------------------------------------------------------------------
  subroutine clm_betr_flux_state_receive(this, bounds, num_soilc, filter_soilc)

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)

    call this%bgc_reaction%lsm_betr_flux_state_receive(bounds, &
         num_soilc, filter_soilc, &
         tracerstate_vars, tracerflux_vars,  betrtracer_vars)

  end subroutine clm_betr_flux_state_receive

  !------------------------------------------------------------------------
  subroutine calc_smp_l_clm(this, bounds, lbj, ubj, &
       numf, filter, t_soisno, soilstate_vars, waterstate_vars, soil_water_retention_curve)

    use SoilStateType, only : soilstate_type
    use WaterStateType, only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use clm_varcon, only : grav,hfus,tfrz

    implicit none

    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  ! bounds
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
                  call soil_water_retention_curve%soil_suction(c,j, s_node, soilstate_vars, smp_l(c,j))
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
    
    call this%bgc_reaction%readParams(ncid, betrtracer_vars)
  end subroutine betr_clm_readParams
  
  !------------------------------------------------------------------------
  subroutine betr_clm_h2oiso_consistency_check(this, &
       bounds, ubj, num_soilc, filter_soilc)
    !
    ! check the overall water mass consistency between betr and clm
    use clm_instMod
    use MathfuncMod, only : dot_sum
    use clm_varcon, only : denh2o
    
    implicit none
    
    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: ubj
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)
    real(r8), allocatable :: eyev(:)
    integer :: fc, c
    real(r8):: totwater, err
    
    associate( &
         h2osoi_ice =>    waterstate_inst%h2osoi_ice_col, & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq =>    waterstate_inst%h2osoi_liq_col, & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)
         end_tracer_molarmass => tracerstate_vars%end_tracer_molarmass_col, &
         id_trc_o18_h2o => betrtracer_vars%id_trc_o18_h2o           &
         )
      
      allocate(eyev(1:ubj))
      eyev=1._r8
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         totwater=dot_sum(h2osoi_ice(c,1:ubj),eyev) + dot_sum(h2osoi_liq(c,1:ubj),eyev)
         err = totwater-end_tracer_molarmass(c,id_trc_o18_h2o)
         print*,get_nstep(),'diff',c, totwater, err, err/totwater
      enddo
      deallocate(eyev)
      
    end associate
  end subroutine betr_clm_h2oiso_consistency_check
end module BeTRSimulationCLM
