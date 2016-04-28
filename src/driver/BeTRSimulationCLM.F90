module BeTRSimulationCLM
  !
  ! !DESCRIPTION:
  !  CLM-BeTR interface
  !
  ! !USES:
  !
#include "shr_assert.h"

  use shr_kind_mod        , only : r8 => shr_kind_r8

  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_log_mod, only : errMsg => shr_log_errMsg

  use decompMod, only : bounds_type
  use EcophysConType, only : ecophyscon_type
  use BeTRSimulation, only : betr_simulation_type
  use betr_decompMod    , only : betr_bounds_type
  use BeTR_TimeMod, only : betr_time_type
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
     !clm unique subroutines
     procedure, public :: ConsistencyCheck => clm_h2oiso_consistency_check
     procedure, public :: DiagnoseDtracerFreezeThaw => CLMDiagnoseDtracerFreezeThaw
     procedure, public :: CalcDewSubFlux => CLMCalcDewSubFlux
     procedure, public :: BetrFluxStateReceive      => CLMBetrFluxStateReceive
     procedure, public :: CalcSmpL                  => CLMCalcSmpL
  end type betr_simulation_clm_type

  public :: create_betr_simulation_clm

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_clm()
    !DESCRIPTION
    ! constructor
    implicit none

    class(betr_simulation_clm_type), pointer :: create_betr_simulation_clm
    class(betr_simulation_clm_type), pointer :: simulation

    allocate(simulation)
    create_betr_simulation_clm => simulation

  end function create_betr_simulation_clm

  !-------------------------------------------------------------------------------

  subroutine CLMInit(this, base_filename, namelist_buffer, bounds, waterstate, cnstate)
    !DESCRIPTION
    !initialize interface
    !
    !USES
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
    use CNStateType, only : cnstate_type
    use landunit_varcon, only : istcrop, istice, istsoil
    use BeTR_landvarconType, only : betr_landvarcon
    use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type), intent(inout) :: this
    character(len=betr_filename_length), intent(in) :: base_filename
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    type(bounds_type), intent(in) :: bounds
    type(waterstate_type), intent(inout) :: waterstate
    type(cnstate_type), intent(inout) :: cnstate
    !TEMPORARY VARIABLES
    type(betr_bounds_type)     :: betr_bounds
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

    ! allocate the reaction types that may only be known to this
    ! simulation type.
    ! now call the base simulation init to continue initialization
    call this%BeTRInit(base_filename, namelist_buffer, &
         betr_bounds, waterstate, cnstate)

  end subroutine CLMInit


  !---------------------------------------------------------------------------------
  subroutine CLMStepWithoutDrainage(this, betr_time, bounds, col, &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, &
       temperature_vars, waterflux_vars, chemstate_vars, &
       cnstate_vars, canopystate_vars, carbonflux_vars)
   !DESCRIPTION
   !march one step without drainage
   !
   !USES
    use SoilStateType, only : soilstate_type
    use WaterStateType, only : Waterstate_Type
    use TemperatureType, only : temperature_type
    use ChemStateType, only : chemstate_type
    use WaterfluxType, only : waterflux_type
    use ColumnType, only : column_type
    use BGCReactionsMod, only : bgc_reaction_type
    use atm2lndType, only : atm2lnd_type
    use SoilHydrologyType, only : soilhydrology_type
    use CNStateType, only : cnstate_type
    use CNCarbonFluxType, only : carbonflux_type
    use CanopyStateType, only : canopystate_type
    use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use PatchType, only : pft
    use LandunitType,only : lun
    use pftvarcon, only : crop
    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    ! !ARGUMENTS :
    class(betr_simulation_clm_type), intent(inout) :: this
    class(betr_time_type), intent(in) :: betr_time
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
    type(betr_bounds_type)     :: betr_bounds
    integer  :: lbj, ubj ! lower and upper bounds, make sure they are > 0

    !pass necessary data for correct subroutine call
    betr_nlevsoi = nlevsoi
    betr_nlevsno = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

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

    call this%SetBiophysForcing(betr_bounds, &
      cnstate_vars = cnstate_vars, &
      carbonflux_vars=carbonflux_vars, &
      waterstate_vars=waterstate_vars, &
      waterflux_vars=waterflux_vars, &
      temperature_vars=temperature_vars,&
      soilhydrology_vars=soilhydrology_vars,&
      atm2lnd_vars=atm2lnd_vars,&
      canopystate_vars=canopystate_vars, &
      chemstate_vars=chemstate_vars, &
      soilstate_vars=soilstate_vars)

    call this%betr%step_without_drainage(betr_time, betr_bounds,   &
         this%num_soilc, this%filter_soilc, this%num_soilp, this%filter_soilp,  &
         this%biophys_forc, this%biogeo_flux, this%biogeo_state)

  end subroutine CLMStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine CLMStepWithDrainage(this, bounds, col)
   !DESCRIPTION
   !march one step with drainage
  !
  !USES
    use ColumnType, only : column_type
    use MathfuncMod, only : safe_div
    use WaterFluxType, only : waterflux_type
    use betr_decompMod, only : betr_bounds_type
    use LandunitType,only : lun
    use tracer_varcon, only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use clm_varpar, only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(column_type), intent(in) :: col ! column type

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

    call this%betr%step_with_drainage(betr_bounds,   &
         this%num_soilc, this%filter_soilc, this%jtops, &
         this%biogeo_flux)

  end subroutine CLMStepWithDrainage

  !------------------------------------------------------------------------

  subroutine CLMDiagnoseDtracerFreezeThaw(this, bounds, num_nolakec, filter_nolakec, col, lun, &
    waterstate_vars)
    !
    ! DESCRIPTION
    ! aqueous tracer partition based on freeze-thaw
    !
    ! USES
    use ColumnType, only : column_type
    use LandunitType, only : landunit_type
    use landunit_varcon, only : istsoil
    use betr_decompMod    , only : betr_bounds_type
    use WaterStateType, only : waterstate_type
    implicit none
    !!ARGUMENTS
    class(betr_simulation_clm_type), intent(inout) :: this
    type(betr_bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_nolakec ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:) ! column filter for non-lake points
    type(column_type), intent(in) :: col ! column type
    type(landunit_type), intent(in) :: lun
    type(waterstate_type), intent(in) :: waterstate_vars

    !temporary variables
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

    call this%SetBiophysForcing(betr_bounds, &
      waterstate_vars=waterstate_vars)

    call this%betr%diagnose_dtracer_freeze_thaw(bounds, num_nolakec, filter_nolakec,  &
      this%biophys_forc)

  end subroutine CLMDiagnoseDtracerFreezeThaw

  !------------------------------------------------------------------------
  subroutine CLMCalcDewSubFlux(this, betr_time, &
       bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars)
    !DESCRIPTION
    ! External interface called by CLM
    !
    !USES
    use ColumnType, only : col
    use LandunitType, only : lun
    use WaterfluxType, only : waterflux_type
    use WaterstateType, only : waterstate_type
    use clm_varcon, only : denh2o,spval
    use landunit_varcon, only : istsoil, istcrop
    use betr_decompMod    , only : betr_bounds_type
    use BeTR_TimeMod, only : betr_time_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type), intent(inout) :: this
    class(betr_time_type), intent(in) :: betr_time
    type(betr_bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc ! number of column soil points in column filter_soilc
    integer, intent(in) :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
    type(waterstate_type), intent(in) :: waterstate_vars
    type(waterflux_type), intent(in) :: waterflux_vars

    call this%SetBiophysForcing(bounds, &
      waterstate_vars=waterstate_vars, &
      waterflux_vars=waterflux_vars)

   call this%betr%calc_dew_sub_flux(betr_time, bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       this%biophys_forc, this%betr%tracers, this%betr%tracerfluxes, this%betr%tracerstates)

  end subroutine CLMCalcDewSubFlux

  !------------------------------------------------------------------------
  subroutine CLMBetrFluxStateReceive(this, bounds, num_soilc, filter_soilc)
   !DESCRIPTION
   !this is to expaneded
   !
   !USES
    use betr_decompMod    , only : betr_bounds_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type), intent(inout) :: this
    type(betr_bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)

    call this%betr%bgc_reaction%lsm_betr_flux_state_receive(bounds, &
         num_soilc, filter_soilc, &
         this%betr%tracerstates, this%betr%tracerfluxes,  this%betr%tracers)

  end subroutine CLMBetrFluxStateReceive

  !------------------------------------------------------------------------
  subroutine CLMCalcSmpL(this, bounds, lbj, ubj, &
       numf, filter, t_soisno, soilstate_vars, waterstate_vars, soil_water_retention_curve)
   !DESCRIPTION
   ! calculate water suction potential
   !
   !USES
    use SoilStateType, only : soilstate_type
    use WaterStateType, only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use betr_decompMod    , only : betr_bounds_type
    use clm_varcon, only : grav,hfus,tfrz
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer, intent(in) :: numf                                              ! number of columns in column filter
    integer, intent(in) :: filter(:)                                         ! column filter
    real(r8), intent(in) :: t_soisno(bounds%begc:, lbj: )                    ! soil temperature
    type(soilstate_type), intent(in) :: soilstate_vars
    type(waterstate_type), intent(inout) :: waterstate_vars
    class(soil_water_retention_curve_type), intent(in), optional :: soil_water_retention_curve

    !temporary variables
    real(r8) :: s_node
    integer :: fc, c, j

    SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(mod_filename,__LINE__))

    ! humor the compiler about unused variables
    if (this%num_soilc > 0) continue
    if (present(soil_water_retention_curve)) continue

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
  end subroutine CLMCalcSmpL

!X!  !------------------------------------------------------------------------
!X!  subroutine betr_clm_readParams(this, ncid)
!X!
!X!    use ncdio_pio, only : file_desc_t
!X!
!X!    implicit none
!X!
!X!    class(betr_simulation_clm_type), intent(inout) :: this
!X!    type(file_desc_t), intent(inout) :: ncid  ! pio netCDF file id
!X!    call this%betr%bgc_reaction%readParams(ncid, this%betr%tracers)
!X!  end subroutine betr_clm_readParams

  !------------------------------------------------------------------------
  subroutine clm_h2oiso_consistency_check(this, &
       bounds, ubj, num_soilc, filter_soilc, waterstate_vars)
    !DESCRIPTION
    ! check the overall water mass consistency between betr and clm
    !
    !USES
    use MathfuncMod, only : dot_sum
    use clm_varcon, only : denh2o
    use WaterStateType, only : waterstate_type
    use clm_time_manager, only : get_nstep
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: ubj
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)
    type(waterstate_type), intent(in) :: waterstate_vars
    !TEMPORARY VARIABLES
    real(r8), allocatable :: eyev(:)
    integer :: fc, c
    real(r8):: totwater, err

    ! humor the compiler about unused variables
    if (bounds%begc > 0) continue

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
  end subroutine clm_h2oiso_consistency_check

end module BeTRSimulationCLM
