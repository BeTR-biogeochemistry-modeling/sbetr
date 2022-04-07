module BeTRSimulationStandalone
  !
  ! !DESCRIPTION:
  !  BeTR standalone simulation class.
  !
  !  BeTR simulation class are API definitions, mapping data
  !  structures from a specific LSM, e.g. CLM, ELM, into BeTR data
  !  structures. The standalone class use BeTR data structures
  !  natively, so is mapping BeTR to BeTR. This class shouldn't be
  !  doing much.
  !
#include "shr_assert.h"
  use abortutils          , only : endrun
  use elm_varctl          , only : iulog
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use betr_decompMod      , only : betr_bounds_type
  use decompMod           , only : bounds_type
  use BeTRSimulation      , only : betr_simulation_type
  use tracer_varcon       , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use EcophysConType      , only : ecophyscon_type
  use betr_varcon         , only : betr_maxpatch_pft
  use CNCarbonStateType  , only : carbonstate_type
  use CNNitrogenStateType, only : nitrogenstate_type
  use pfNitrogenStateType, only : pf_nitrogenstate_type
  use PhosphorusStateType, only : phosphorusstate_type
  use CNCarbonFluxType   , only : carbonflux_type
  use PfCarbonFluxType  , only : pf_carbonflux_type
  use CNNitrogenFluxType , only : nitrogenflux_type
  use PfNitrogenFluxType , only : pf_nitrogenflux_type
  use PhosphorusFluxType , only : phosphorusflux_type
  use PfPhosphorusFluxType , only : pf_phosphorusflux_type
  use WaterStateType  , only : waterstate_type
  use WaterfluxType     , only : waterflux_type
  use TemperatureType   , only : temperature_type
  use PfTemperatureType   , only : pf_temperature_type
  use PfWaterfluxType     , only : pf_waterflux_type
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_standalone_type
     ! NOTE(bja, 201603) LSM specific types here!

     ! NOTE(bja, 201603) most (all?) BeTR types go into the base
     ! class.

   contains
     procedure, public :: Init                => StandaloneInitOffline
     procedure, public :: StepWithoutDrainage => StandaloneStepWithoutDrainage
     procedure, public :: StepWithDrainage    => StandaloneStepWithDrainage
     procedure, public :: SetBiophysForcing   => StandaloneSetBiophysForcing
     procedure, public :: PlantSoilBGCSend    => StandalonePlantSoilBGCSend
     procedure, public :: PlantSoilBGCRecv    => StandalonePlantSoilBGCRecv
     procedure, public :: CalcSmpL            => StandaloneCalcSmpL
     procedure, private :: set_transient_kinetics_par
  end type betr_simulation_standalone_type

  public :: create_betr_simulation_standalone

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_standalone()
  ! DESCRIPTION
  ! constructor
    implicit none

    class(betr_simulation_standalone_type), pointer :: create_betr_simulation_standalone
    class(betr_simulation_standalone_type), pointer :: simulation

    allocate(simulation)
    create_betr_simulation_standalone => simulation

  end function create_betr_simulation_standalone



  !-------------------------------------------------------------------------------

  subroutine StandaloneInitOffline(this, bounds, lun, col, pft, waterstate, namelist_buffer, base_filename, case_id)

    !DESCRIPTION
    !initialize standalone betr
    !
    !USES
    use betr_constants      , only : betr_namelist_buffer_size
    use betr_constants      , only : betr_filename_length
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use PatchType           , only : patch_type
    use ColumnType          , only : column_type
    use LandunitType        , only : landunit_type
    use pftvarcon           , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use WaterStateType      , only : waterstate_type
    use CNStateType         , only : cnstate_type
    use landunit_varcon     , only : istcrop, istice, istsoil
    use BeTR_landvarconType , only : betr_landvarcon
    use elm_varpar          , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    !ARGUMENTS
    class(betr_simulation_standalone_type)   , intent(inout) :: this
    character(len=*)                         , intent(in)    :: base_filename
    character(len=*)                         , intent(in)    :: case_id
    character(len=*)                         , intent(in)    :: namelist_buffer
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout)    :: col
    type(patch_type)                         , intent(in)    :: pft
    type(waterstate_type)                    , intent(inout) :: waterstate

    betr_nlevsoi      = nlevsoi
    betr_nlevsno      = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    !pass necessary data for correct subroutine call

    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice

    ! allocate the reaction types that may only be known to this
    ! simulation type by deallocating and overriding methods created
    ! in betr%Init().

    ! now call the base simulation init to continue initialization
    call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer, base_filename, case_id)

    !pass necessary data

  end subroutine StandaloneInitOffline


  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithoutDrainage(this, bounds, col, pft)
    !DESCRIPTION
    !march one step without drainage
    !
    !USES
    use ColumnType        , only : column_type
    use PatchType         , only : patch_type
    use BeTR_TimeMod      , only : betr_time_type
    use LandunitType      , only : lun
    use pftvarcon         , only : crop
    implicit none
    !ARGUMENTS
    class(betr_simulation_standalone_type) , intent(inout) :: this
    type(bounds_type)                      , intent(in)    :: bounds ! bounds
    type(column_type)                      , intent(in)    :: col ! column type
    type(patch_type)                       , intent(in)    :: pft
    !temporary variables
    type(betr_bounds_type)     :: betr_bounds

    integer  :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer  :: c
    !pass necessary data for correct subroutine call
    !set lbj and ubj

    call this%BeTRSetBounds(betr_bounds)

    call this%bsimstatus%reset()

    call this%BeTRSetcps(bounds, col, pft)

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle

      call this%biophys_forc(c)%frac_normalize(this%betr_pft(c)%npfts, 1, betr_nlevtrc_soil)

      call this%betr(c)%step_without_drainage(this%betr_time, betr_bounds, this%betr_col(c), &
         this%betr_pft(c), this%num_surfc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc(c), this%biogeo_flux(c), this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine StandaloneStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithDrainage(this, bounds, col)
    !DESCRIPTION
    !march one step with drainage
    !
    !USES
    use ColumnType        , only : column_type
    use MathfuncMod       , only : safe_div
    use WaterFluxType     , only : waterflux_type
    use LandunitType      , only : lun
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_standalone_type) , intent(inout) :: this
    type(bounds_type)                      , intent(in)    :: bounds
    type(column_type)                      , intent(in)    :: col ! column type

    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer :: c, c_l, begc_l, endc_l

    !pass necessary data for correct subroutine call
    !set lbj and ubj

    call this%bsimstatus%reset()
    call this%BeTRSetBounds(betr_bounds)
    call this%BeTRSetcps(bounds, col)
    c_l = 1; begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle

      call this%betr(c)%step_with_drainage(betr_bounds, this%betr_col(c), &
         this%num_surfc, this%filter_soilc,        &
         this%jtops, this%biogeo_flux(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo

  end subroutine StandaloneStepWithDrainage

  !------------------------------------------------------------------------
  subroutine StandaloneSetBiophysForcing(this, bounds, col, pft,  carbonflux_vars,pf_carbonflux_vars,&
    waterstate_vars, waterflux_vars, pf_waterflux_vars, temperature_vars, pf_temperature_vars, &
    soilhydrology_vars, atm2lnd_vars, canopystate_vars, chemstate_vars, soilstate_vars, cnstate_vars)
  !DESCRIPTION
  !pass in biogeophysical variables for running betr
  !USES
  use ColumnType        , only : column_type
  use PatchType         , only : patch_type
  use SoilStateType     , only : soilstate_type
  use WaterStateType    , only : Waterstate_Type
  use TemperatureType   , only : temperature_type
  use ChemStateType     , only : chemstate_type
  use WaterfluxType     , only : waterflux_type
  use atm2lndType       , only : atm2lnd_type
  use SoilHydrologyType , only : soilhydrology_type
  use CNStateType       , only : cnstate_type
  use CNCarbonFluxType  , only : carbonflux_type
  use CanopyStateType   , only : canopystate_type
  use elm_varpar        , only : nlevsno, nlevsoi
  implicit none
  !ARGUMENTS
  class(betr_simulation_standalone_type) , intent(inout)        :: this
  type(bounds_type)                      , intent(in)           :: bounds
  type(column_type)                      , intent(in)    :: col ! column type
  type(patch_type)                       , intent(in) :: pft
  type(cnstate_type)                     , optional, intent(in) :: cnstate_vars
  type(carbonflux_type)       , optional, intent(in) :: carbonflux_vars
  type(pf_carbonflux_type)    , optional, intent(in) :: pf_carbonflux_vars
  type(Waterstate_Type)       , optional, intent(in) :: Waterstate_vars
  type(waterflux_type)        , optional, intent(in) :: waterflux_vars
  type(pf_waterflux_type)     , optional, intent(in) :: pf_waterflux_vars
  type(temperature_type)      , optional, intent(in) :: temperature_vars
  type(pf_temperature_type)   , optional, intent(in) :: pf_temperature_vars
  type(soilhydrology_type)    , optional, intent(in) :: soilhydrology_vars
  type(atm2lnd_type)          , optional, intent(in) :: atm2lnd_vars
  type(canopystate_type)      , optional, intent(in) :: canopystate_vars
  type(chemstate_type)        , optional, intent(in) :: chemstate_vars
  type(soilstate_type)        , optional, intent(in) :: soilstate_vars
  integer :: j, c, c_l


  call this%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, &
    carbonflux_vars, pf_carbonflux_vars, waterstate_vars,  waterflux_vars, pf_waterflux_vars, &
    temperature_vars, pf_temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars)

  !the following will be standalone bgc specific
  end subroutine StandaloneSetBiophysForcing

  !------------------------------------------------------------------------
  subroutine StandalonePlantSoilBGCSend(this, bounds, col, pft, num_surfc,  filter_soilc, cnstate_vars, &
    carbonflux_vars,  c13_cflx_vars, c14_cflx_vars, nitrogenflux_vars, phosphorusflux_vars, &
    PlantMicKinetics_vars)

  !read in biogeochemical fluxes from elm for soil bgc modeling
  !these are C, N and P fluxes from root input, surface litter input
  !atmospheric deposition, fire (negative), and fertilization
  !Because of possible harvest activity that is
  !related to dynamic land use, input profiles are computed in elm.
  !
  use ColumnType         , only : column_type
  use PatchType          , only : patch_type
  use CNCarbonFluxType   , only : carbonflux_type
  use CNNitrogenFluxType , only : nitrogenflux_type
  use PhosphorusFluxType , only : phosphorusflux_type
  use CNStateType        , only : cnstate_type
  use elm_varpar         , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  use mathfuncMod        , only : apvb,bisnan
  use tracer_varcon      , only : use_c13_betr, use_c14_betr
  implicit none
  class(betr_simulation_standalone_type) , intent(inout)        :: this
  type(bounds_type) , intent(in)  :: bounds
  type(column_type) , intent(in)  :: col ! column type
  type(patch_type)  , intent(in)  :: pft ! pft type
  integer           , intent(in)  :: num_surfc
  integer           , intent(in)  :: filter_soilc(:)
  type(cnstate_type), intent(in)  :: cnstate_vars
  type(carbonflux_type), intent(in):: carbonflux_vars
  type(carbonflux_type), intent(in):: c13_cflx_vars
  type(carbonflux_type), intent(in):: c14_cflx_vars
  type(nitrogenflux_type), intent(in):: nitrogenflux_vars
  type(phosphorusflux_type), intent(in):: phosphorusflux_vars
  type(PlantMicKinetics_type), intent(in) :: PlantMicKinetics_vars

  integer :: c_l, c, fc, j
  type(betr_bounds_type)   :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)

  !reset and prepare for retrieval
  do c = bounds%begc, bounds%endc
    if(.not. this%active_col(c))cycle
    call this%biogeo_state(c)%reset(value_column=0._r8, active_soibgc=this%do_soibgc())
    call this%biogeo_flux(c)%reset(value_column=0._r8, active_soibgc=this%do_soibgc())
  enddo
  if(.not. this%do_soibgc())return

  c_l = 1
  !set kinetic parameters
  call this%set_transient_kinetics_par(betr_bounds, col, pft, num_surfc, filter_soilc, PlantMicKinetics_vars)

  do fc = 1, num_surfc
    c = filter_soilc(fc)
    call this%biophys_forc(c)%reset(value_column=0._r8)
    this%biophys_forc(c)%isoilorder(c_l) = 1
    call this%biophys_forc(c)%c12flx%reset(value_column=0._r8)
    call this%biophys_forc(c)%n14flx%reset(value_column=0._r8)
    call this%biophys_forc(c)%p31flx%reset(value_column=0._r8)
    this%biophys_forc(c)%lithotype_col(c_l) = cnstate_vars%lithoclass_col(c)
    if(use_c13_betr)then
      call this%biophys_forc(c)%c13flx%reset(value_column=0._r8)
    endif

    if(use_c14_betr)then
      call this%biophys_forc(c)%c14flx%reset(value_column=0._r8)
    endif

    this%biophys_forc(c)%frac_loss_lit_to_fire_col(c_l) = 0._r8
    this%biophys_forc(c)%frac_loss_cwd_to_fire_col(c_l)=0._r8
  enddo

  do j = betr_bounds%lbj, betr_bounds%ubj
    do fc = 1, num_surfc
      c = filter_soilc(fc)
      this%biophys_forc(c)%c12flx%rt_vr_col(c_l,j) = carbonflux_vars%rr_vr(c,j)
      this%biophys_forc(c)%n14flx%nflx_minn_input_nh4_vr_col(c_l,j)=1.e-10_r8
      this%biophys_forc(c)%p31flx%pflx_minp_input_po4_vr_col(c_l,j)=1.e-11_r8
      this%biophys_forc(c)%pweath_prof_col(c_l,j) = cnstate_vars%pdep_prof_col(c,j)
      this%biophys_forc(c)%c12flx%cflx_input_litr_met_vr_col(c_l,j) = carbonflux_vars%cflx_input_litr_met_vr(c,j)
      this%biophys_forc(c)%c12flx%cflx_input_litr_cel_vr_col(c_l,j) = carbonflux_vars%cflx_input_litr_cel_vr(c,j)
      this%biophys_forc(c)%c12flx%cflx_input_litr_lig_vr_col(c_l,j) = carbonflux_vars%cflx_input_litr_lig_vr(c,j)
      this%biophys_forc(c)%c12flx%cflx_input_litr_cwd_vr_col(c_l,j) = carbonflux_vars%cflx_input_litr_cwd_vr(c,j)

      this%biophys_forc(c)%n14flx%nflx_input_litr_met_vr_col(c_l,j) = nitrogenflux_vars%nflx_input_litr_met_vr(c,j)
      this%biophys_forc(c)%n14flx%nflx_input_litr_cel_vr_col(c_l,j) = nitrogenflux_vars%nflx_input_litr_cel_vr(c,j)
      this%biophys_forc(c)%n14flx%nflx_input_litr_lig_vr_col(c_l,j) = nitrogenflux_vars%nflx_input_litr_lig_vr(c,j)
      this%biophys_forc(c)%n14flx%nflx_input_litr_cwd_vr_col(c_l,j) = nitrogenflux_vars%nflx_input_litr_cwd_vr(c,j)
      this%biophys_forc(c)%n14flx%nflx_minn_input_nh4_vr_col(c_l,j) = nitrogenflux_vars%nflx_minn_input_nh4_vr(c,j)
      this%biophys_forc(c)%n14flx%nflx_minn_input_no3_vr_col(c_l,j) = nitrogenflux_vars%nflx_minn_input_no3_vr(c,j)

      this%biophys_forc(c)%p31flx%pflx_input_litr_met_vr_col(c_l,j) = phosphorusflux_vars%pflx_input_litr_met_vr(c,j)
      this%biophys_forc(c)%p31flx%pflx_input_litr_cel_vr_col(c_l,j) = phosphorusflux_vars%pflx_input_litr_cel_vr(c,j)
      this%biophys_forc(c)%p31flx%pflx_input_litr_lig_vr_col(c_l,j) = phosphorusflux_vars%pflx_input_litr_lig_vr(c,j)
      this%biophys_forc(c)%p31flx%pflx_input_litr_cwd_vr_col(c_l,j) = phosphorusflux_vars%pflx_input_litr_cwd_vr(c,j)
      this%biophys_forc(c)%p31flx%pflx_minp_input_po4_vr_col(c_l,j) = phosphorusflux_vars%pflx_minp_input_po4_vr(c,j)

      this%biophys_forc(c)%biochem_pmin_vr(c_l,j) =0._r8
    enddo
  enddo

  end subroutine StandalonePlantSoilBGCSend

  !------------------------------------------------------------------------
  subroutine StandalonePlantSoilBGCRecv(this, bounds, col, pft, num_surfc,  filter_soilc,&
   c12state_vars, c12flux_vars, c13state_vars, c13flux_vars, c14state_vars, c14flux_vars, &
   n14state_vars, n14flux_vars, p31state_vars, p31flux_vars)

  !DESCRIPTION
  !retrieve the variables

  use ColumnType         , only : column_type
  use PatchType          , only : patch_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use elm_time_manager    , only : get_nstep
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use tracer_varcon       , only : use_c13_betr, use_c14_betr
  use pftvarcon           , only : noveg
  use MathfuncMod         , only : safe_div
  use tracer_varcon       , only : reaction_method
  implicit none
  class(betr_simulation_standalone_type) , intent(inout) :: this
  type(bounds_type) , intent(in)  :: bounds
  type(patch_type)            , intent(in) :: pft
  type(column_type)           , intent(in)    :: col ! column type
  integer           , intent(in)  :: num_surfc
  integer           , intent(in)  :: filter_soilc(:)
  type(carbonstate_type), intent(inout) :: c12state_vars
  type(carbonstate_type), intent(inout) :: c13state_vars
  type(carbonstate_type), intent(inout) :: c14state_vars
  type(nitrogenstate_type), intent(inout) :: n14state_vars
  type(phosphorusstate_type), intent(inout) :: p31state_vars
  type(carbonflux_type)  , intent(inout):: c12flux_vars    !return carbon fluxes through DON?
  type(carbonflux_type)  , intent(inout):: c13flux_vars    !return carbon fluxes through DON?
  type(carbonflux_type)  , intent(inout):: c14flux_vars    !return carbon fluxes through DON?
  type(nitrogenflux_type), intent(inout):: n14flux_vars
  type(phosphorusflux_type), intent(inout):: p31flux_vars

  integer :: c, fc, p, pi, c_l

  !TEMPORARY VARIABLES
  type(betr_bounds_type)     :: betr_bounds
  integer :: begc_l, endc_l

  !summarize the fluxes and state variables
  c_l = 1
  call this%BeTRSetBounds(betr_bounds)
  begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

  !retrieve and return
  do c = bounds%begc, bounds%endc
    if(.not. this%active_col(c))cycle

    call this%betr(c)%retrieve_biostates(betr_bounds,  1, betr_nlevsoi, &
       this%num_surfc, this%filter_soilc, this%jtops, this%biogeo_state(c),this%bstatus(c))

    if(this%bstatus(c)%check_status())then
      call this%bsimstatus%setcol(c)
      call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
      exit
    endif
    call this%betr(c)%retrieve_biofluxes(this%num_surfc, this%filter_soilc, &
      this%num_soilp, this%filter_soilp, this%biogeo_flux(c))

    call this%biogeo_state(c)%summary(betr_bounds, 1, betr_nlevtrc_soil,&
         this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil), &
         this%betr_col(c)%zi(begc_l:endc_l,1:betr_nlevtrc_soil),this%do_soibgc())

    call this%biogeo_flux(c)%summary(betr_bounds, 1, betr_nlevtrc_soil, &
         this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil),this%do_soibgc())

    call this%biophys_forc(c)%summary(betr_bounds, 1, betr_nlevtrc_soil, &
         this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil), this%do_soibgc())
  enddo

  if(.not. this%do_soibgc())return

  do fc = 1, num_surfc
    c = filter_soilc(fc)

    !recollect soil respirations, fire and hydraulic loss
    c12flux_vars%hr(c) = this%biogeo_flux(c)%c12flux_vars%hr_col(c_l)
    c12flux_vars%fire_decomp_closs(c) = this%biogeo_flux(c)%c12flux_vars%fire_decomp_closs_col(c_l)
    c12flux_vars%som_c_leached(c) = &
        this%biogeo_flux(c)%c12flux_vars%som_c_leached_col(c_l) + &
        this%biogeo_flux(c)%c12flux_vars%som_c_qdrain_col(c_l)
    c12flux_vars%som_c_runoff(c) = this%biogeo_flux(c)%c12flux_vars%som_c_runoff_col(c_l)

    !the following is for consistency with the ELM definitation, which computes
    !som_c_leached_col as a numerical roundoff
    c12flux_vars%som_c_leached(c)=-c12flux_vars%som_c_leached(c)
    c12flux_vars%co2_soi_flx(c) = this%biogeo_flux(c)%c12flux_vars%co2_soi_flx_col(c_l)
    if(use_c13_betr)then
      c13flux_vars%hr(c) = this%biogeo_flux(c)%c13flux_vars%hr_col(c_l)
      c13flux_vars%fire_decomp_closs(c) = this%biogeo_flux(c)%c13flux_vars%fire_decomp_closs_col(c_l)
    endif

    if(use_c14_betr)then
      c14flux_vars%hr(c) = this%biogeo_flux(c)%c14flux_vars%hr_col(c_l)
      c14flux_vars%fire_decomp_closs(c) = this%biogeo_flux(c)%c14flux_vars%fire_decomp_closs_col(c_l)
    endif

    !recollect  nitrifications, nitrifier-N2O loss, denitrifications
    n14flux_vars%f_nit(c) = this%biogeo_flux(c)%n14flux_vars%f_nit_col(c_l)
    n14flux_vars%f_denit(c)= this%biogeo_flux(c)%n14flux_vars%f_denit_col(c_l)
    n14flux_vars%denit(c)= n14flux_vars%f_denit(c)
    n14flux_vars%f_n2o_nit(c)=this%biogeo_flux(c)%n14flux_vars%f_n2o_nit_col(c_l)

    !hydraulic loss
    n14flux_vars%smin_no3_leached(c)= &
        this%biogeo_flux(c)%n14flux_vars%smin_no3_leached_col(c_l) + &
        this%biogeo_flux(c)%n14flux_vars%smin_no3_qdrain_col(c_l)
    n14flux_vars%smin_nh4_leached(c)= &
        this%biogeo_flux(c)%n14flux_vars%smin_nh4_leached_col(c_l) + &
        this%biogeo_flux(c)%n14flux_vars%smin_nh4_qdrain_col(c_l)


    n14flux_vars%som_n_leached(c) = &
        this%biogeo_flux(c)%n14flux_vars%som_n_leached_col(c_l) + &
        this%biogeo_flux(c)%n14flux_vars%som_n_qdrain_col(c_l)
    n14flux_vars%nh3_soi_flx(c) = this%biogeo_flux(c)%n14flux_vars%nh3_soi_flx_col(c_l)
    n14flux_vars%som_n_runoff(c) = this%biogeo_flux(c)%n14flux_vars%som_n_runoff_col(c_l)
    !the following is for consistency with the ELM definitation, which computes
    !som_n_leached_col as a numerical roundoff
    n14flux_vars%som_n_leached(c) = - n14flux_vars%som_n_leached(c)

    n14flux_vars%smin_no3_runoff(c)=this%biogeo_flux(c)%n14flux_vars%smin_no3_runoff_col(c_l)
    n14flux_vars%smin_nh4_runoff(c)=this%biogeo_flux(c)%n14flux_vars%smin_nh4_runoff_col(c_l)
    !fire loss
    n14flux_vars%fire_decomp_nloss(c) = this%biogeo_flux(c)%n14flux_vars%fire_decomp_nloss_col(c_l)
    n14flux_vars%supplement_to_sminn(c) = this%biogeo_flux(c)%n14flux_vars%supplement_to_sminn_col(c_l)
    !no nh4 volatilization and runoff/leaching loss at this moment

    !recollect mineral phosphorus loss
    !Remark: now hydraulic mineral p loss lumps all three fluxes, Jinyun Tang
    p31flux_vars%sminp_runoff(c)=this%biogeo_flux(c)%p31flux_vars%sminp_runoff_col(c_l)
    p31flux_vars%sminp_leached(c) = &
       this%biogeo_flux(c)%p31flux_vars%sminp_leached_col(c_l) + &
       this%biogeo_flux(c)%p31flux_vars%sminp_qdrain_col(c_l)

    p31flux_vars%supplement_to_sminp(c) = this%biogeo_flux(c)%p31flux_vars%supplement_to_sminp_col(c_l)
    p31flux_vars%secondp_to_occlp(c) = this%biogeo_flux(c)%p31flux_vars%secondp_to_occlp_col(c_l)
    p31flux_vars%fire_decomp_ploss(c) = this%biogeo_flux(c)%p31flux_vars%fire_decomp_ploss_col(c_l)
    p31flux_vars%som_p_leached(c) = &
        this%biogeo_flux(c)%p31flux_vars%som_p_leached_col(c_l) + &
        this%biogeo_flux(c)%p31flux_vars%som_p_qdrain_col(c_l)
    p31flux_vars%som_p_runoff(c) = this%biogeo_flux(c)%p31flux_vars%som_p_runoff_col(c_l)
    !the following is for consistency with the ELM definitation, which computes
    !som_p_leached_col as a numerical roundoff
    p31flux_vars%som_p_leached(c) = -p31flux_vars%som_p_leached(c)

    c12state_vars%beg_totomc(c) = c12state_vars%totomc(c)

    !recollect soil organic carbon, soil organic nitrogen, and soil organic phosphorus
    c12state_vars%cwdc(c) = this%biogeo_state(c)%c12state_vars%cwdc_col(c_l)
    c12state_vars%totlitc(c) = this%biogeo_state(c)%c12state_vars%totlitc_col(c_l)
    c12state_vars%totsomc(c) = this%biogeo_state(c)%c12state_vars%totsomc_col(c_l)
    c12state_vars%totlitc_1m(c) = this%biogeo_state(c)%c12state_vars%totlitc_1m_col(c_l)
    c12state_vars%totsomc_1m(c) = this%biogeo_state(c)%c12state_vars%totsomc_1m_col(c_l)
    c12state_vars%totomc(c) = this%biogeo_state(c)%c12state_vars%totomc(c_l)
    c12state_vars%cmass_residual(c) = c12state_vars%totomc(c) - c12state_vars%beg_totomc(c) &
       - (this%biophys_forc(c)%c12flx%cflx_input_col(c_l)-c12flux_vars%hr(c)  &
       - c12flux_vars%fire_decomp_closs(c)+c12flux_vars%som_c_leached(c) &
       - c12flux_vars%som_c_runoff(c)) * this%betr_time%get_step_size()

    if(use_c13_betr)then
      c13state_vars%beg_totomc(c) = c13state_vars%totomc(c)
      c13state_vars%cwdc(c) = this%biogeo_state(c)%c13state_vars%cwdc_col(c_l)
      c13state_vars%totlitc(c) = this%biogeo_state(c)%c13state_vars%totlitc_col(c_l)
      c13state_vars%totsomc(c) = this%biogeo_state(c)%c13state_vars%totsomc_col(c_l)
      c13state_vars%totlitc_1m(c) = this%biogeo_state(c)%c13state_vars%totlitc_1m_col(c_l)
      c13state_vars%totsomc_1m(c) = this%biogeo_state(c)%c13state_vars%totsomc_1m_col(c_l)
      c13state_vars%totomc(c) = this%biogeo_state(c)%c13state_vars%totomc(c_l)
      c13state_vars%cmass_residual(c) = c13state_vars%totomc(c) - c13state_vars%beg_totomc(c) &
       - (this%biophys_forc(c)%c13flx%cflx_input_col(c_l)-c13flux_vars%hr(c)  &
       - c13flux_vars%fire_decomp_closs(c)+c13flux_vars%som_c_leached(c) &
       - c13flux_vars%som_c_runoff(c)) * this%betr_time%get_step_size()

    endif
    if(use_c14_betr)then
      c14state_vars%beg_totomc(c) = c14state_vars%totomc(c)
      c14state_vars%cwdc(c) = this%biogeo_state(c)%c14state_vars%cwdc_col(c_l)
      c14state_vars%totlitc(c) = this%biogeo_state(c)%c14state_vars%totlitc_col(c_l)
      c14state_vars%totsomc(c) = this%biogeo_state(c)%c14state_vars%totsomc_col(c_l)
      c14state_vars%totlitc_1m(c) = this%biogeo_state(c)%c14state_vars%totlitc_1m_col(c_l)
      c14state_vars%totsomc_1m(c) = this%biogeo_state(c)%c14state_vars%totsomc_1m_col(c_l)
      c14state_vars%totomc(c) = this%biogeo_state(c)%c14state_vars%totomc(c_l)
      c14state_vars%cmass_residual(c) = c14state_vars%totomc(c) - c14state_vars%beg_totomc(c) &
       - (this%biophys_forc(c)%c14flx%cflx_input_col(c_l)-c14flux_vars%hr(c)  &
       - c14flux_vars%fire_decomp_closs(c)+c14flux_vars%som_c_leached(c) &
       - c14flux_vars%som_c_runoff(c)) * this%betr_time%get_step_size()
    endif
    n14state_vars%beg_totsoin(c) = n14state_vars%totsoin(c)
    n14state_vars%cwdn(c) = this%biogeo_state(c)%n14state_vars%cwdn_col(c_l)
    n14state_vars%totlitn(c) = this%biogeo_state(c)%n14state_vars%totlitn_col(c_l)
    n14state_vars%totsomn(c) = this%biogeo_state(c)%n14state_vars%totsomn_col(c_l)
    n14state_vars%totlitn_1m(c) = this%biogeo_state(c)%n14state_vars%totlitn_1m_col(c_l)
    n14state_vars%totsomn_1m(c) = this%biogeo_state(c)%n14state_vars%totsomn_1m_col(c_l)
    n14state_vars%totsoin(c) = this%biogeo_state(c)%n14state_vars%totsoin(c_l)
    n14state_vars%nmass_residual(c) = n14state_vars%totsoin(c) - n14state_vars%beg_totsoin(c) &
       - (this%biophys_forc(c)%n14flx%nflx_input_col(c_l) &
       + this%biophys_forc(c)%n14flx%nflx_minninput_col(c_l)  &
       - n14flux_vars%fire_decomp_nloss(c)+n14flux_vars%som_n_leached(c) &
       - n14flux_vars%som_n_runoff(c)-n14flux_vars%nh3_soi_flx(c) &
       - n14flux_vars%smin_no3_leached(c) - n14flux_vars%smin_nh4_leached(c) &
       - n14flux_vars%smin_no3_runoff(c) - n14flux_vars%smin_nh4_runoff(c) &
       - n14flux_vars%f_n2o_nit(c)-n14flux_vars%denit(c)) * this%betr_time%get_step_size()


    p31state_vars%beg_totsoip(c) = p31state_vars%totsoip(c)
    p31state_vars%cwdp(c) = this%biogeo_state(c)%p31state_vars%cwdp_col(c_l)
    p31state_vars%totlitp(c) = this%biogeo_state(c)%p31state_vars%totlitp_col(c_l)
    p31state_vars%totsomp(c) = this%biogeo_state(c)%p31state_vars%totsomp_col(c_l)
    p31state_vars%totlitp_1m(c) = this%biogeo_state(c)%p31state_vars%totlitp_1m_col(c_l)
    p31state_vars%totsomp_1m(c) = this%biogeo_state(c)%p31state_vars%totsomp_1m_col(c_l)

    !recollect inorganic nitrogen (smin_nh4, smin_no3), and inorganic phosphorus (disolvable and protected)
    n14state_vars%sminn(c) = this%biogeo_state(c)%n14state_vars%sminn_col(c_l)
    n14state_vars%smin_nh4(c)=this%biogeo_state(c)%n14state_vars%sminn_nh4_col(c_l)
    n14state_vars%smin_no3(c)=this%biogeo_state(c)%n14state_vars%sminn_no3_col(c_l)

    p31state_vars%sminp(c) = this%biogeo_state(c)%p31state_vars%sminp_col(c_l)
    p31state_vars%occlp(c) = this%biogeo_state(c)%p31state_vars%occlp_col(c_l)
    !print*,'smin_nh4',n14state_vars%smin_nh4_col(c)

    p31state_vars%totsoip(c) = this%biogeo_state(c)%p31state_vars%totsoip(c_l)
    p31state_vars%pmass_residual(c) = p31state_vars%totsoip(c) - p31state_vars%beg_totsoip(c) &
       - (this%biophys_forc(c)%p31flx%pflx_input_col(c_l) &
       + this%biophys_forc(c)%p31flx%pminp_input_col(c_l)  &
       - p31flux_vars%fire_decomp_ploss(c)+p31flux_vars%som_p_leached(c) &
       -p31flux_vars%sminp_leached(c) - p31flux_vars%sminp_runoff(c)) &
       * this%betr_time%get_step_size()

    c12state_vars%som1c(c) = this%biogeo_state(c)%c12state_vars%som1c_col(c_l)
    c12state_vars%som2c(c) = this%biogeo_state(c)%c12state_vars%som2c_col(c_l)
    c12state_vars%som3c(c) = this%biogeo_state(c)%c12state_vars%som3c_col(c_l)
    c12state_vars%domc(c)  = this%biogeo_state(c)%c12state_vars%domc_col(c_l)

    !print*,'som1c',c12state_vars%som1c_col(c)
    if(use_c13_betr)then
      c13state_vars%som1c(c) = this%biogeo_state(c)%c13state_vars%som1c_col(c_l)
      c13state_vars%som2c(c) = this%biogeo_state(c)%c13state_vars%som2c_col(c_l)
      c13state_vars%som3c(c) = this%biogeo_state(c)%c13state_vars%som3c_col(c_l)
      c13state_vars%domc(c)  = this%biogeo_state(c)%c13state_vars%domc_col(c_l)
    endif
    if(use_c14_betr)then
      c14state_vars%som1c(c) = this%biogeo_state(c)%c14state_vars%som1c_col(c_l)
      c14state_vars%som2c(c) = this%biogeo_state(c)%c14state_vars%som2c_col(c_l)
      c14state_vars%som3c(c) = this%biogeo_state(c)%c14state_vars%som3c_col(c_l)
      c14state_vars%domc(c)  = this%biogeo_state(c)%c14state_vars%domc_col(c_l)
    endif

    if(index(trim(reaction_method),'ecacnp')/=0 .or. &
      index(trim(reaction_method),'cdom')/=0 .or. &
      index(trim(reaction_method),'keca')/=0 .or. &
      index(trim(reaction_method),'ch4soil')/=0)then
      n14state_vars%som1n(c) = this%biogeo_state(c)%n14state_vars%som1n_col(c_l)
      n14state_vars%som2n(c) = this%biogeo_state(c)%n14state_vars%som2n_col(c_l)
      n14state_vars%som3n(c) = this%biogeo_state(c)%n14state_vars%som3n_col(c_l)
      n14state_vars%domn(c) = this%biogeo_state(c)%n14state_vars%domn_col(c_l)

      p31state_vars%som1p(c) = this%biogeo_state(c)%p31state_vars%som1p_col(c_l)
      p31state_vars%som2p(c) = this%biogeo_state(c)%p31state_vars%som2p_col(c_l)
      p31state_vars%som3p(c) = this%biogeo_state(c)%p31state_vars%som3p_col(c_l)
      p31state_vars%domp(c) = this%biogeo_state(c)%p31state_vars%domp_col(c_l)
    endif
  enddo
  end subroutine StandalonePlantSoilBGCRecv

  !------------------------------------------------------------------------
  subroutine set_transient_kinetics_par(this, betr_bounds, col, pft, num_surfc, filter_soilc, PlantMicKinetics_vars)
  !DESCRIPTION
  !set kinetic parameters for column c
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  use tracer_varcon      , only : reaction_method,natomw,patomw
  use pftvarcon             , only : noveg
  use PatchType          , only : patch_type
  use ColumnType        , only : column_type
  implicit none
  class(betr_simulation_standalone_type) , intent(inout) :: this
  type(betr_bounds_type), intent(in) :: betr_bounds
  type(column_type)     , intent(in)    :: col ! column type
  type(patch_type)      , intent(in) :: pft
  integer, intent(in) :: num_surfc
  integer, intent(in) :: filter_soilc(:)
  type(PlantMicKinetics_type), intent(in) :: PlantMicKinetics_vars

  integer :: j, fc, c, p, pi, pp, c_l

  associate(      &
    plant_nh4_vmax_vr_patch => PlantMicKinetics_vars%plant_nh4_vmax_vr_patch, &
    plant_no3_vmax_vr_patch => PlantMicKinetics_vars%plant_no3_vmax_vr_patch, &
    plant_p_vmax_vr_patch   => PlantMicKinetics_vars%plant_p_vmax_vr_patch, &
    plant_nh4_km_vr_patch   => PlantMicKinetics_vars%plant_nh4_km_vr_patch, &
    plant_no3_km_vr_patch   => PlantMicKinetics_vars%plant_no3_km_vr_patch, &
    plant_p_km_vr_patch     => PlantMicKinetics_vars%plant_p_km_vr_patch , &
    plant_eff_ncompet_b_vr_patch => PlantMicKinetics_vars%plant_eff_ncompet_b_vr_patch , &
    plant_eff_pcompet_b_vr_patch => PlantMicKinetics_vars%plant_eff_pcompet_b_vr_patch , &
    minsurf_nh4_compet_vr_col => PlantMicKinetics_vars%minsurf_nh4_compet_vr_col, &
    minsurf_p_compet_vr_col => PlantMicKinetics_vars%minsurf_p_compet_vr_col , &
    plant_eff_frootc_vr_patch => PlantMicKinetics_vars%plant_eff_frootc_vr_patch &
  )
  c_l = 1
  do fc = 1, num_surfc
    c = filter_soilc(fc)
    pp = 0

    do pi = 1, betr_maxpatch_pft
      if (pi <= col%npfts(c)) then
        p = col%pfti(c) + pi - 1

        if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
          pp = pp + 1
          do j =1, betr_bounds%ubj
            this%betr(c)%plantNutkinetics%plant_nh4_vmax_vr_patch(pp,j) = 0._r8
            this%betr(c)%plantNutkinetics%plant_no3_vmax_vr_patch(pp,j) = 0._r8
            this%betr(c)%plantNutkinetics%plant_p_vmax_vr_patch(pp,j) = 0._r8
            this%betr(c)%plantNutkinetics%plant_nh4_km_vr_patch(pp,j) = 0._r8
            this%betr(c)%plantNutkinetics%plant_no3_km_vr_patch(pp,j) = 0._r8
            this%betr(c)%plantNutkinetics%plant_p_km_vr_patch(pp,j) = 0._r8
            this%betr(c)%plantNutkinetics%plant_eff_ncompet_b_vr_patch(pp,j)= 0._r8
            this%betr(c)%plantNutkinetics%plant_eff_pcompet_b_vr_patch(pp,j)= 0._r8
            this%betr(c)%plantNutkinetics%plant_eff_frootc_vr_patch(pp,j) = 0._r8
          enddo
        endif
      endif
    enddo
    this%betr(c)%nactpft = pp
    do j = 1, betr_bounds%ubj
      this%betr(c)%plantNutkinetics%minsurf_p_compet_vr_col(c_l,j) = 1._r8
      this%betr(c)%plantNutkinetics%minsurf_nh4_compet_vr_col(c_l,j) = 1._r8
      this%betr(c)%plantNutkinetics%minsurf_dom_compet_vr_col(c_l,j) = PlantMicKinetics_vars%minsurf_dom_compet_vr_col(c,j)
    enddo
  enddo

  !the following parameters are specific to ECACNP, and I assume they are
  !grid specific as they currently used in elm-cnp.
  do j =1, betr_bounds%ubj
    do fc = 1, num_surfc
      c = filter_soilc(fc)
      this%betr(c)%plantNutkinetics%km_minsurf_p_vr_col(c_l,j) = 1._r8
      this%betr(c)%plantNutkinetics%km_minsurf_nh4_vr_col(c_l,j)=1._r8
      this%betr(c)%plantNutkinetics%km_minsurf_dom_vr_col(c_l,j)=PlantMicKinetics_vars%km_minsurf_dom_vr_col(c,j)
    enddo
  enddo

  end associate
  end subroutine set_transient_kinetics_par

  !------------------------------------------------------------------------
  subroutine StandaloneCalcSmpL(this, bounds, lbj, ubj, numf, filter, t_soisno, &
     soilstate_vars, waterstate_vars, soil_water_retention_curve)
  !DESCRIPTION
  ! calculate soil suction potential
  !
  !USES
  use SoilStateType              , only : soilstate_type
  use WaterStateType             , only : waterstate_type
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
  use elm_varcon                 , only : grav,hfus,tfrz
  implicit none
  !ARGUMENTS
  class(betr_simulation_standalone_type) , intent(inout) :: this
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
  integer  :: fc, c, j, c_l
  real(r8) :: dsmpds_top
  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(mod_filename,__LINE__))

  ! remove compiler warnings
  if (this%num_surfc > 0) continue

  associate(                                                     & !
    h2osoi_vol        =>    waterstate_vars%h2osoi_vol         , & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
    smp_l             =>    waterstate_vars%smp_l              , & ! Output: [real(r8) (:,:) ]  soil suction (mm)
    bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
    watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
    sucsat            =>    soilstate_vars%sucsat_col            & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
  )

  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(.not. this%active_col(c))cycle
      if(j==1)then
        if(t_soisno(c,j)<tfrz)then
          smp_l(c,j)= -hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
        else
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_l(c,j), dsmpds_top)
!  the following will be implemented later.
!          call soil_water_retention_curve%soil_hk(hksat, 1._r8, s_node, bsw(c,j), hk)
        endif

!        this%biophys_forc(c)%Dw_hk(c_l) = hk *
      else
        if(t_soisno(c,j)<tfrz)then
          smp_l(c,j)= -hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
        else
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_l(c,j))
        endif
      endif
    enddo
  enddo
  end associate
  end subroutine StandaloneCalcSmpL

end module BeTRSimulationStandalone
