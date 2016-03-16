module H2OIsotopePlantSoilBGCType


  use PlantSoilBGCMod , only : plant_soilbgc_type
  implicit none
  save
  private


  public :: plant_soilbgc_h2oiso_run_type

  type, extends(plant_soilbgc_type) :: &
    plant_soilbgc_h2oiso_run_type
  private
    contains
    procedure :: Init_plant_soilbgc
    procedure :: plant_soilbgc_summary
    procedure :: integrate_vr_flux_to_2D
    procedure :: lsm_betr_plant_soilbgc_send
  end type plant_soilbgc_h2oiso_run_type


  interface plant_soilbgc_h2oiso_run_type
    module procedure constructor
  end interface plant_soilbgc_h2oiso_run_type

  contains

  !-------------------------------------------------------------------------------
  type(plant_soilbgc_h2oiso_run_type) function constructor()
  !
  ! !DESCRIPTION:
  ! create an object of type plant_soilbgc_h2oiso_run_type.
  ! Right now it is purposely empty

  end function constructor

  !-------------------------------------------------------------------------------
  subroutine Init_plant_soilbgc(this, bounds, lbj, ubj)

  !
  ! !DESCRIPTION:
  ! template for init_betrbgc
  !
  ! !USES:
  use BeTR_decompMod             , only : betr_bounds_type

  ! !ARGUMENTS:
  class(plant_soilbgc_h2oiso_run_type) , intent(in) :: this
  type(betr_bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj

  end subroutine Init_plant_soilbgc


  !----------------------------------------------------------------------
  subroutine plant_soilbgc_summary(this,bounds, lbj, ubj, numf, &
       filter, dz, betrtracer_vars, tracerflux_vars)

  ! !USES:
  use BeTRTracerType        , only : BeTRtracer_type
  use tracerfluxType        , only : tracerflux_type
  use BeTR_decompMod             , only : betr_bounds_type
  use bshr_kind_mod          , only : r8 => shr_kind_r8

  ! !ARGUMENTS:

  class(plant_soilbgc_h2oiso_run_type) , intent(in) :: this
  type(betr_bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)
  real(r8)                  , intent(in) :: dz(bounds%begc:bounds%endc,1:ubj)
  type(BeTRtracer_type )    , intent(in) :: betrtracer_vars
  type(tracerflux_type)     , intent(in) :: tracerflux_vars


  end subroutine plant_soilbgc_summary


  !----------------------------------------------------------------------

  subroutine integrate_vr_flux_to_2D(this, bounds, numf, filter)

  use BeTR_decompMod             , only : betr_bounds_type
  ! !ARGUMENTS:

  class(plant_soilbgc_h2oiso_run_type) , intent(in) :: this
  type(betr_bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)


  end subroutine integrate_vr_flux_to_2D

  !----------------------------------------------------------------------

  subroutine lsm_betr_plant_soilbgc_send(this, bounds, numf, filter)

  use BeTR_decompMod             , only : betr_bounds_type
  ! !ARGUMENTS:

  class(plant_soilbgc_h2oiso_run_type) , intent(in) :: this
  type(betr_bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)

  end subroutine lsm_betr_plant_soilbgc_send
end module H2OIsotopePlantSoilBGCType
