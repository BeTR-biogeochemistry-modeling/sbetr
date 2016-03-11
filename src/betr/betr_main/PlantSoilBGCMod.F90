module PlantSoilBGCMod

!
! !DESCRIPTION:
! template for doing plant and soil bgc coupling
!
! !USES:

type, abstract :: plant_soilbgc_type
   private
 contains

   !initialize Init_plant_soilbgc
   procedure(Init_plant_soilbgc_interface)                    , deferred :: Init_plant_soilbgc

   procedure(lsm_betr_plant_soilbgc_send_interface)           , deferred :: lsm_betr_plant_soilbgc_send

   procedure(plant_soilbgc_summary_interface)                 , deferred :: plant_soilbgc_summary

   procedure(integrate_vr_flux_to_2D_interface)               , deferred :: integrate_vr_flux_to_2D
end type plant_soilbgc_type


  abstract interface
  !----------------------------------------------------------------------
  subroutine Init_plant_soilbgc_interface(this, bounds, lbj, ubj)

  !
  ! !DESCRIPTION:
  ! template for init_betrbgc
  !
  ! !USES:
  use BeTR_decompMod         , only : bounds_type  => betr_bounds_type

  ! !ARGUMENTS:
  import :: plant_soilbgc_type
  class(plant_soilbgc_type) , intent(in) :: this
  type(bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj

  end subroutine Init_plant_soilbgc_interface


  !----------------------------------------------------------------------
  subroutine plant_soilbgc_summary_interface(this,bounds, lbj, ubj, numf, &
       filter, dz, betrtracer_vars, tracerflux_vars)

  ! !USES:
  use BeTRTracerType        , only : BeTRtracer_type
  use tracerfluxType        , only : tracerflux_type
  use BeTR_decompMod        , only : bounds_type  => betr_bounds_type
  use shr_kind_mod          , only : r8 => shr_kind_r8

  ! !ARGUMENTS:
  import :: plant_soilbgc_type

  class(plant_soilbgc_type) , intent(in) :: this
  type(bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)
  real(r8)                  , intent(in) :: dz(bounds%begc:bounds%endc,1:ubj)
  type(BeTRtracer_type )    , intent(in) :: betrtracer_vars
  type(tracerflux_type)     , intent(in) :: tracerflux_vars


  end subroutine plant_soilbgc_summary_interface


  !----------------------------------------------------------------------

  subroutine integrate_vr_flux_to_2D_interface(this, bounds, numf, filter)
  !
  ! !DESCRIPTIONS
  ! integrate 3d fluxes into 2d fluxes
  use BeTR_decompMod         , only : bounds_type  => betr_bounds_type
  ! !ARGUMENTS:
  import :: plant_soilbgc_type

  class(plant_soilbgc_type) , intent(in) :: this
  type(bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)


  end subroutine integrate_vr_flux_to_2D_interface
  !----------------------------------------------------------------------

  subroutine lsm_betr_plant_soilbgc_send_interface(this, bounds, numf, filter)
  !
  ! send lsm variables into betr
  ! the interface will be further revised for plant soilbgc coupling
  !
  use BeTR_decompMod         , only : bounds_type  => betr_bounds_type
  ! !ARGUMENTS:
  import :: plant_soilbgc_type

  class(plant_soilbgc_type) , intent(in) :: this
  type(bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)

  end subroutine lsm_betr_plant_soilbgc_send_interface
  end interface
end module PlantSoilBGCMod
