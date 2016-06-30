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

    real(r8), pointer :: cflx_input_litr_met_vr_col(:,:) => null() ! metabolic litter input
    real(r8), pointer :: cflx_input_litr_cel_vr_col(:,:) => null() ! cellulose litter input
    real(r8), pointer :: cflx_input_litr_lig_vr_col(:,:) => null() ! lignin litter input
    real(r8), pointer :: cflx_input_litr_cwd_vr_col(:,:) => null() ! coarse woody debries input
    !The only loss is through fire and no som is lost through burning
    real(r8), pointer :: cflx_output_litr_met_vr_col(:,:) => null() ! metabolic litter input
    real(r8), pointer :: cflx_output_litr_cel_vr_col(:,:) => null() ! cellulose litter input
    real(r8), pointer :: cflx_output_litr_lig_vr_col(:,:) => null() ! lignin litter input
    real(r8), pointer :: cflx_output_litr_cwd_vr_col(:,:) => null() ! coarse woody debries input

    real(r8), pointer :: nflx_input_litr_met_vr_col(:,:) => null() ! metabolic litter input
    real(r8), pointer :: nflx_input_litr_cel_vr_col(:,:) => null() ! cellulose litter input
    real(r8), pointer :: nflx_input_litr_lig_vr_col(:,:) => null() ! lignin litter input
    real(r8), pointer :: nflx_input_litr_cwd_vr_col(:,:) => null() ! coarse woody debries input
    !The only loss is through fire and no som is lost through burning
    real(r8), pointer :: nflx_output_litr_met_vr_col(:,:) => null() ! metabolic litter input
    real(r8), pointer :: nflx_output_litr_cel_vr_col(:,:) => null() ! cellulose litter input
    real(r8), pointer :: nflx_output_litr_lig_vr_col(:,:) => null() ! lignin litter input
    real(r8), pointer :: nflx_output_litr_cwd_vr_col(:,:) => null() ! coarse woody debries input

    real(r8), pointer :: pflx_input_litr_met_vr_col(:,:) => null() ! metabolic litter input
    real(r8), pointer :: pflx_input_litr_cel_vr_col(:,:) => null() ! cellulose litter input
    real(r8), pointer :: pflx_input_litr_lig_vr_col(:,:) => null() ! lignin litter input
    real(r8), pointer :: pflx_input_litr_cwd_vr_col(:,:) => null() ! coarse woody debries input
    !The only loss is through fire and no som is lost through burning
    real(r8), pointer :: pflx_output_litr_met_vr_col(:,:) => null() ! metabolic litter input
    real(r8), pointer :: pflx_output_litr_cel_vr_col(:,:) => null() ! cellulose litter input
    real(r8), pointer :: pflx_output_litr_lig_vr_col(:,:) => null() ! lignin litter input
    real(r8), pointer :: pflx_output_litr_cwd_vr_col(:,:) => null() ! coarse woody debries input

    real(r8), pointer :: sflx_minn_input_nh4_vr_col(:,:)  => null() !mineral nh4 input through deposition & fertilization
    real(r8), pointer :: sflx_minn_input_no3_vr_col(:,:)  => null() !mineral no3 input through deposition & fertilization
    real(r8), pointer :: sflx_minp_input_po4_vr_col(:,:)  => null() !mineral phosphorus input through deposition & fertilization
    real(r8), pointer :: sflx_minp_weathering_po4_vr_col(:,:)  => null() !mineral phosphorus input through weathering


    real(r8), pointer :: sflx_minn_nh4_fix_vr_col(:,:) => null()    !nitrogen fixation
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

  allocate(this%cflx_input_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%cflx_input_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%cflx_input_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%cflx_input_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%cflx_output_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%cflx_output_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%cflx_output_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%cflx_output_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input

  allocate(this%nflx_input_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%nflx_input_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%nflx_input_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%nflx_input_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%nflx_output_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%nflx_output_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%nflx_output_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%nflx_output_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input

  allocate(this%pflx_input_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%pflx_input_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%pflx_input_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%pflx_input_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%pflx_output_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%pflx_output_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%pflx_output_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%pflx_output_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input

  allocate(this%sflx_minn_input_nh4_vr_col(begc:endc,lbj:ubj)) !mineral nh4 input through deposition & fertilization
  allocate(this%sflx_minn_input_no3_vr_col(begc:endc,lbj:ubj)) !mineral no3 input through deposition & fertilization
  allocate(this%sflx_minp_input_po4_vr_col(begc:endc,lbj:ubj)) !mineral phosphorus input through weathering, deposition & fertilization
  allocate(this%sflx_minn_nh4_fix_vr_col(begc:endc,lbj:ubj))   !nh4 from fixation
  allocate(this%sflx_minp_weathering_po4_vr_col(begc:endc,lbj:ubj)) !p from weathering
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_biogeo_flux_type)       :: this
  real(r8), intent(in) :: value_column


  this%cflx_input_litr_met_vr_col(:,:) = value_column
  this%cflx_input_litr_cel_vr_col(:,:) = value_column
  this%cflx_input_litr_lig_vr_col(:,:)= value_column
  this%cflx_input_litr_cwd_vr_col(:,:)= value_column
  this%cflx_output_litr_met_vr_col(:,:)= value_column
  this%cflx_output_litr_cel_vr_col(:,:)= value_column
  this%cflx_output_litr_lig_vr_col(:,:)= value_column
  this%cflx_output_litr_cwd_vr_col(:,:)= value_column

  this%nflx_input_litr_met_vr_col(:,:)= value_column
  this%nflx_input_litr_cel_vr_col(:,:)= value_column
  this%nflx_input_litr_lig_vr_col(:,:)= value_column
  this%nflx_input_litr_cwd_vr_col(:,:)= value_column
  this%nflx_output_litr_met_vr_col(:,:)= value_column
  this%nflx_output_litr_cel_vr_col(:,:)= value_column
  this%nflx_output_litr_lig_vr_col(:,:)= value_column
  this%nflx_output_litr_cwd_vr_col(:,:)= value_column

  this%pflx_input_litr_met_vr_col(:,:)= value_column
  this%pflx_input_litr_cel_vr_col(:,:)= value_column
  this%pflx_input_litr_lig_vr_col(:,:)= value_column
  this%pflx_input_litr_cwd_vr_col(:,:)= value_column
  this%pflx_output_litr_met_vr_col(:,:)= value_column
  this%pflx_output_litr_cel_vr_col(:,:)= value_column
  this%pflx_output_litr_lig_vr_col(:,:)= value_column
  this%pflx_output_litr_cwd_vr_col(:,:)= value_column

  this%sflx_minn_input_nh4_vr_col(:,:)= value_column
  this%sflx_minn_input_no3_vr_col(:,:)= value_column
  this%sflx_minp_input_po4_vr_col(:,:)= value_column
  this%sflx_minn_nh4_fix_vr_col(:,:)= value_column
  this%sflx_minp_weathering_po4_vr_col(:,:) = value_column
  end subroutine reset

end module BeTR_biogeoFluxType
