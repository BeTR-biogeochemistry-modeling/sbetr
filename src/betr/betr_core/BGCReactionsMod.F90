module BGCReactionsMod
  !
  ! !DESCRIPTION:
  ! template for doing bgc reaction in betr
  !
  ! !USES:
  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type

  implicit none

  private

  public ::  bgc_reaction_type

  type, abstract :: bgc_reaction_type
     private
     ! dummy var to remove compiler warnings
     logical, public :: dummy_compiler_warning
   contains
     !initialize betr bgc
     procedure(Init_betrbgc_interface)                    , deferred :: Init_betrbgc

     !doing bgc reaction
     procedure(calc_bgc_reaction_interface)               , deferred :: calc_bgc_reaction

     !set boundary condition for related tracer transport
     procedure(set_boundary_conditions_interface)         , deferred :: set_boundary_conditions

     procedure(init_boundary_condition_type_interface)    , deferred :: init_boundary_condition_type

     !do equilibrium tracer chemistry
     procedure(do_tracer_equilibration_interface )        , deferred :: do_tracer_equilibration

     !do cold initialization of different tracers
     procedure(initCold_interface)                        , deferred :: initCold

     !read in implementation specific parameters
     procedure(readParams_interface)                      , deferred :: readParams

     !send back soil state flux variables to other parts of lsm
     procedure(lsm_betr_flux_state_receive_interface)     , deferred :: lsm_betr_flux_state_receive

  end type bgc_reaction_type

  abstract interface
     !----------------------------------------------------------------------
     subroutine Init_betrbgc_interface(this, bounds, lbj, ubj, betrtracer_vars, bstatus)
       !
       ! !DESCRIPTION:
       ! template for init_betrbgc
       !
       ! !USES:
       use BeTRTracerType  , only : BeTRtracer_type
       use BeTR_decompMod  , only : betr_bounds_type
       use BetrStatusType  , only : betr_status_type
       !
       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type) , intent(inout)    :: this
       type(betr_bounds_type)   , intent(in)    :: bounds
       integer                  , intent(in)    :: lbj, ubj
       type(BeTRtracer_type )   , intent(inout) :: betrtracer_vars
       type(betr_status_type)   , intent(out)   :: bstatus

     end subroutine Init_betrbgc_interface
     !----------------------------------------------------------------------
     subroutine calc_bgc_reaction_interface(this, bounds, col, lbj, ubj, num_soilc, filter_soilc, &
          num_soilp,filter_soilp, jtops, dtime, betrtracer_vars, tracercoeff_vars,  biophysforc,    &
          tracerstate_vars, tracerflux_vars,  tracerboundarycond_vars, plant_soilbgc, betr_status)
       !
       ! !DESCRIPTION:
       ! template for calc_bgc_reaction
       !
       ! !USES:
       use bshr_kind_mod            , only : r8 => shr_kind_r8
       use TracerBoundaryCondType   , only : tracerboundarycond_type
       use tracerfluxType           , only : tracerflux_type
       use tracerstatetype          , only : tracerstate_type
       use tracercoeffType          , only : tracercoeff_type
       use BeTR_decompMod           , only : betr_bounds_type
       use BeTRTracerType           , only : BeTRTracer_Type
       use PlantSoilBGCMod          , only : plant_soilbgc_type
       use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
       use BetrStatusType           , only : betr_status_type
       use betr_columnType          , only : betr_column_type
       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)         , intent(inout)   :: this
       type(betr_bounds_type)           , intent(in)    :: bounds                      ! bounds
       type(betr_column_type)           , intent(in)    :: col
       integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
       integer                          , intent(in)    :: filter_soilc(:)             ! column filter
       integer                          , intent(in)    :: num_soilp
       integer                          , intent(in)    :: filter_soilp(:)
       integer                          , intent(in)    :: jtops( : )                  ! top index of each column
       integer                          , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
       real(r8)                         , intent(in)    :: dtime                       ! model time step
       type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
       type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
       type(tracercoeff_type)           , intent(in)    :: tracercoeff_vars
       type(tracerstate_type)           , intent(inout) :: tracerstate_vars
       type(tracerflux_type)            , intent(inout) :: tracerflux_vars
       type(tracerboundarycond_type)    , intent(inout) :: tracerboundarycond_vars !
       class(plant_soilbgc_type)        , intent(inout) :: plant_soilbgc
       type(betr_status_type)           , intent(out)   :: betr_status

     end subroutine calc_bgc_reaction_interface
     !----------------------------------------------------------------------

     subroutine set_boundary_conditions_interface(this, bounds, num_soilc, filter_soilc, dz_top, &
          betrtracer_vars, biophysforc, biogeo_flux, tracerboundarycond_vars, betr_status)

       ! !DESCRIPTION:
       ! template for set_boundary_conditions
       !
       ! !USES:
       use TracerBoundaryCondType   , only : tracerboundarycond_type
       use BeTR_decompMod           , only : betr_bounds_type
       use BeTRTracerType           , only : BeTRTracer_Type
       use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
       use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
       use bshr_kind_mod            , only : r8 => shr_kind_r8
       use BetrStatusType           , only : betr_status_type
       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)         , intent(inout)    :: this                       !
       type(betr_bounds_type)           , intent(in)    :: bounds                     !
       integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter
       integer                          , intent(in)    :: filter_soilc(:)            ! column filter
       type(betrtracer_type)            , intent(in)    :: betrtracer_vars            !
       real(r8)                         , intent(in)    :: dz_top( : )                !
       type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
       type(betr_biogeo_flux_type)      , intent(in)    :: biogeo_flux
       type(tracerboundarycond_type)    , intent(inout) :: tracerboundarycond_vars !
       type(betr_status_type)           , intent(out)   :: betr_status

     end subroutine set_boundary_conditions_interface

     !----------------------------------------------------------------------

     subroutine init_boundary_condition_type_interface(this, bounds, &
          betrtracer_vars, tracerboundarycond_vars )
       !
       ! !DESCRIPTION:
       ! template for init_boundary_condition
       !
       ! !USES:
       use BeTRTracerType        , only : betrtracer_type
       use TracerBoundaryCondType, only : tracerboundarycond_type
       use BeTR_decompMod        , only : betr_bounds_type

       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)      , intent(inout) :: this
       type(BeTRtracer_type )        , intent(in) :: betrtracer_vars
       type(betr_bounds_type)        , intent(in) :: bounds
       type(tracerboundarycond_type) , intent(in) :: tracerboundarycond_vars

     end subroutine init_boundary_condition_type_interface

     !-------------------------------------------------------------------------------
     subroutine do_tracer_equilibration_interface(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
          betrtracer_vars, tracercoeff_vars, tracerstate_vars, betr_status)
       !
       ! !DESCRIPTION:
       ! template for do_tracer_equilibration
       ! !USES:
       !
       use tracerstatetype       , only : tracerstate_type
       use tracercoeffType       , only : tracercoeff_type
       use BeTRTracerType        , only : BeTRTracer_Type
       use BeTR_decompMod        , only : betr_bounds_type
       use BetrStatusType        , only : betr_status_type
       ! !ARGUMENTS:
       import :: bgc_reaction_type

       class(bgc_reaction_type)   , intent(inout)    :: this
       type(betr_bounds_type)     , intent(in)    :: bounds
       integer                    , intent(in)    :: lbj, ubj
       integer                    , intent(in)    :: jtops( : )        ! top label of each column
       integer                    , intent(in)    :: num_soilc
       integer                    , intent(in)    :: filter_soilc(:)
       type(betrtracer_type)      , intent(in)    :: betrtracer_vars
       type(tracercoeff_type)     , intent(in)    :: tracercoeff_vars
       type(tracerstate_type)     , intent(inout) :: tracerstate_vars
       type(betr_status_type)     , intent(out)   :: betr_status

     end subroutine do_tracer_equilibration_interface

     !-------------------------------------------------------------------------------
     subroutine InitCold_interface(this, bounds, col, betrtracer_vars, biophysforc, tracerstate_vars)
       !
       ! !DESCRIPTION:
       ! template for InitCold
       ! !USES:
       !
       use BeTRTracerType           , only : BeTRTracer_Type
       use tracerstatetype          , only : tracerstate_type
       use BeTR_decompMod           , only : betr_bounds_type
       use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
       use betr_columnType          , only : betr_column_type
       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)         , intent(inout) :: this
       type(betr_bounds_type)           , intent(in)    :: bounds
       type(betr_column_type)           , intent(in)    :: col
       type(BeTRTracer_Type)            , intent(in)    :: betrtracer_vars
       type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
       type(tracerstate_type)           , intent(inout) :: tracerstate_vars


     end subroutine InitCold_interface

     !-------------------------------------------------------------------------------
     subroutine readParams_interface(this, name_list_buffer, betrtracer_vars)
       !
       ! !DESCRIPTION:
       ! template for readParams
       ! !USES:
       use BeTRTracerType           , only : BeTRTracer_Type

       ! !ARGUMENTS:
       import :: bgc_reaction_type

       class(bgc_reaction_type)          , intent(inout)    :: this
       character(len=*)                  , intent(in)  :: name_list_buffer
       type(BeTRTracer_Type)             , intent(inout) :: betrtracer_vars

     end subroutine readParams_interface

     !-------------------------------------------------------------------------------
     subroutine lsm_betr_flux_state_receive_interface(this, bounds,num_soilc, filter_soilc, &
          tracerstate_vars, tracerflux_vars,  betrtracer_vars)

       ! !DESCRIPTION:
       ! template for lsm_betr_flux_state_receive
       ! !USES:
       use BeTR_decompMod           , only : betr_bounds_type
       use tracerfluxType           , only : tracerflux_type
       use tracerstatetype          , only : tracerstate_type
       use BeTRTracerType           , only : BeTRTracer_Type
       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)   , intent(inout) :: this                  !
       type(betr_bounds_type)     , intent(in) :: bounds                ! bounds
       integer                    , intent(in) :: num_soilc
       integer                    , intent(in) :: filter_soilc(:)
       type(betrtracer_type)      , intent(in) :: betrtracer_vars       ! betr configuration information
       type(tracerstate_type)     , intent(in) :: tracerstate_vars      !
       type(tracerflux_type)      , intent(in) :: tracerflux_vars       !

     end subroutine lsm_betr_flux_state_receive_interface

  end interface
end module BGCReactionsMod
