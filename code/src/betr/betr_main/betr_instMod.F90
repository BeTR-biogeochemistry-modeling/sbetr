module betr_instMod
  !
  ! !DESCRIPTION:
  !  subroutines to initialize betr based on namelist

  ! !USES:
  use BeTRTracerType            , only : BeTRtracer_type
  use TracerCoeffType           , only : TracerCoeff_type
  use TracerFluxType            , only : TracerFlux_type
  use TracerStateType           , only : TracerState_type
  use tracerboundarycondType    , only : tracerboundarycond_type
  use BeTR_aerocondType         , only : betr_aerecond_type
  implicit none
  save
  private   ! By default everything is public


  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  type(BeTRtracer_type)                , public :: betrtracer_vars
  type(TracerCoeff_type)               , public :: tracercoeff_vars
  type(TracerFlux_type)                , public :: tracerflux_vars
  type(TracerState_type)               , public :: tracerState_vars
  type(tracerboundarycond_type)        , public :: tracerboundarycond_vars
  type(betr_aerecond_type)              ,public :: betr_aerecond_vars

end module betr_instMod
