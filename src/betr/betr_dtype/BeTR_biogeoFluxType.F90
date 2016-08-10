module BeTR_biogeoFluxType
  !DESCRIPTION
  !module for flux data exchange between lsm and betr
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  type betr_biogeo_flux_type
    real(r8), pointer :: qflx_adv_col             (:,:) => null() !advection velocity from one layer to another, (0:nlevgrnd), positive downward
    real(r8), pointer :: qflx_gross_evap_soil_col (:)   => null() ! col gross infiltration from soil, this satisfies the relationship qflx_infl_col = qflx_gross_infl_soil_col-qflx_gross_evap_soil_col
    real(r8), pointer :: qflx_gross_infl_soil_col (:)   => null() ! col gross infiltration, before considering the evaporation, mm/s
    real(r8), pointer :: qflx_infl_col            (:)   => null()  !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_drain_vr_col        (:,:) => null() ! col liquid water losted as drainage (m /time step)
    real(r8), pointer :: qflx_totdrain_col        (:)   => null() ! col total liquid water drainage  (m/time step), updated in betr

    real(r8), pointer :: hr_vr_col(:,:)       => null()  !vertically resolved heterotrophic respiration, g C/m2/s
    real(r8), pointer :: f_denit_vr_col(:,:)  => null()  !vertically resolved denitrification, g N /m2/s
    real(r8), pointer :: f_nit_vr_col(:,:)    => null()  !vertically resolved nitrification, g N/m2/s
    real(r8), pointer :: f_n2o_nit_vr(:,:)    => null()  !vertically resolved n2o production from nitrification, gN/m2/s
    contains
      procedure, public  :: Init
      procedure, private :: InitAllocate
  end type betr_biogeo_flux_type

  public :: create_betr_biogeoFlux

contains

  function create_betr_biogeoFlux()
  ! DESCRIPTION
  ! constructor
  implicit none
  class(betr_biogeo_flux_type), pointer :: create_betr_biogeoFlux
  class(betr_biogeo_flux_type), pointer :: biogeflux

  allocate(biogeflux)
  create_betr_biogeoFlux => biogeflux

  end function create_betr_biogeoFlux

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_biogeo_flux_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_biogeo_flux_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  Allocate(this%qflx_adv_col             (begc:endc,lbj-1:ubj) ) !advection velocity from one layer to another, (0:nlevgrnd), positive downward
  allocate(this%qflx_gross_evap_soil_col (begc:endc)  ) ! col gross infiltration from soil, this satisfies the relationship qflx_infl_col = qflx_gross_infl_soil_col-qflx_gross_evap_soil_col
  allocate(this%qflx_gross_infl_soil_col (begc:endc) ) ! col gross infiltration, before considering the evaporation, mm/s
  allocate(this%qflx_infl_col            (begc:endc)  )  !infiltration (mm H2O /s)
  allocate(this%qflx_drain_vr_col        (begc:endc,lbj:ubj) ) ! col liquid water losted as drainage (m /time step)
  allocate(this%qflx_totdrain_col        (begc:endc)   ) ! col total liquid water drainage  (m/time step), updated in betr

  allocate(this%hr_vr_col(begc:endc,lbj:ubj))
  allocate(this%f_denit_vr_col(begc:endc,lbj:ubj))
  allocate(this%f_nit_vr_col(begc:endc,lbj:ubj))
  allocate(this%f_n2o_nit_vr(begc:endc,lbj:ubj))

  end subroutine InitAllocate



end module BeTR_biogeoFluxType
