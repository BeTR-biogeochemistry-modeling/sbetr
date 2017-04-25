module BeTR_biogeoFluxType
  !DESCRIPTION
  !module for flux data exchange between lsm and betr
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
  use tracer_varcon, only : use_c13_betr, use_c14_betr
  use BeTR_carbonfluxRecvType, only : betr_carbonflux_recv_type
  use BeTR_nitrogenfluxRecvType, only : betr_nitrogenflux_recv_type
  use BeTR_phosphorusfluxRecvType, only : betr_phosphorusflux_recv_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type betr_biogeo_flux_type
    real(r8), pointer :: qflx_adv_col             (:,:) => null() !advection velocity from one layer to another, (0:nlevgrnd), positive downward
    real(r8), pointer :: qflx_gross_evap_soil_col (:)   => null() ! col gross infiltration from soil, this satisfies the relationship qflx_infl_col = qflx_gross_infl_soil_col-qflx_gross_evap_soil_col
    real(r8), pointer :: qflx_gross_infl_soil_col (:)   => null() ! col gross infiltration, before considering the evaporation, mm/s
    real(r8), pointer :: qflx_infl_col            (:)   => null()  !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_drain_vr_col        (:,:) => null() ! col liquid water losted as drainage (m /time step)
    real(r8), pointer :: qflx_totdrain_col        (:)   => null() ! col total liquid water drainage  (m/time step), updated in betr


    type(betr_carbonflux_recv_type) :: c12flux_vars
    type(betr_carbonflux_recv_type) :: c13flux_vars
    type(betr_carbonflux_recv_type) :: c14flux_vars
    type(betr_nitrogenflux_recv_type) :: n14flux_vars
    type(betr_phosphorusflux_recv_type) :: p31flux_vars

    contains
      procedure, public  :: Init
      procedure, private :: InitAllocate
      procedure, public  :: reset
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

  call this%c12flux_vars%Init(bounds)
  if(use_c13_betr)then
    call this%c13flux_vars%Init(bounds)
  endif
  if(use_c14_betr)then
    call this%c14flux_vars%Init(bounds)
  endif
  call this%n14flux_vars%Init(bounds)
  call this%p31flux_vars%Init(bounds)

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

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_biogeo_flux_type)  :: this
  real(r8), intent(in) :: value_column

  call this%c12flux_vars%reset(value_column)
  if(use_c13_betr)then
    call this%c13flux_vars%reset(value_column)
  endif
  if(use_c14_betr)then
    call this%c14flux_vars%reset(value_column)
  endif
  call this%n14flux_vars%reset(value_column)
  call this%p31flux_vars%reset(value_column)

  end subroutine reset

end module BeTR_biogeoFluxType
