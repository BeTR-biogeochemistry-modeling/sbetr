module BeTR_biogeophysInputType
  !DESCRIPTION
  !module for input forcing data from lsm
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
  use BeTR_carbonfluxType, only : betr_carbonflux_type
  use BeTR_nitrogenfluxType, only : betr_nitrogenflux_type
  use BeTR_phosphorusfluxType, only : betr_phosphorusflux_type

implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type betr_biogeophys_input_type
    ! cnstate_vars
    integer, pointer :: isoilorder(:)           => null()  ! soil order
    real(r8), pointer:: frac_loss_lit_to_fire_col(:) => null() !fraction of litter cnp loss through fire
    real(r8), pointer:: frac_loss_cwd_to_fire_col(:) => null() !fraction of cwd cnp loss through fire
    integer, pointer :: lithoclass(:)           => null()
    !carbon flux
    real(r8), pointer :: annsum_npp_patch(:)    => null()  !annual npp
    real(r8), pointer :: agnpp_patch(:)         => null()
    real(r8), pointer :: bgnpp_patch(:)         => null()
    !waterstate
    real(r8), pointer :: h2osoi_liq_col(:,:)    => null()    !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_col(:,:)    => null()    !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_liqvol_col(:,:) => null()   !volumetric liquid water content
    real(r8), pointer :: h2osoi_icevol_col(:,:) => null()   !volumetric ice water content
    real(r8), pointer :: h2osoi_vol_col(:,:)    => null()    !volumetric water content, total
    real(r8), pointer :: air_vol_col(:,:)       => null()    !volume possessed by air
    real(r8), pointer :: finundated_col(:)      => null()    ! fraction of column that is inundated, this is for bgc caclulation in betr
    real(r8), pointer :: rho_vap(:,:)           => null()    !concentration of bulk water vapor, assume in equilibrium with the liquid water
    real(r8), pointer :: rhvap_soi(:,:)         => null()    !soil relative humidity
    real(r8), pointer :: smp_l_col     (:,:)    => null()    ! col liquid phase soil matric potential, mm
    real(r8), pointer :: frac_h2osfc_col (:)    => null()    ! col fractional area with surface water greater than zero
    !waterflux
    real(r8), pointer :: qflx_surf_col            (:)      => null()  !surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_rootsoi_col         (:,:)    => null() ! col root and soil water exchange [mm H2O/s] [+ into root]
    real(r8), pointer :: qflx_rootsoi_frac_patch  (:,:)    => null()
    real(r8), pointer :: qflx_dew_grnd_col        (:)      => null() ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
    real(r8), pointer :: qflx_dew_snow_col        (:)      => null() ! col surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_vol_col    (:)      => null()
    real(r8), pointer :: qflx_sub_snow_col        (:)      => null() ! col sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_h2osfc2topsoi_col   (:)      => null() ! col liquid water coming from surface standing water top soil (mm H2O/s)
    real(r8), pointer :: qflx_snow2topsoi_col     (:)      => null() ! col liquid water coming from residual snow to topsoil (mm H2O/s)
    real(r8), pointer :: qflx_tran_veg_patch      (:)      => null()
    real(r8), pointer :: qflx_runoff_col          (:)      => null() !col total runoff
    !temperature
    real(r8), pointer :: t_soisno_col(:,:)                 => null()      !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: t_soi_10cm(:)                     => null()  !soil temperature in top 10cm of soil (Kelvin)
    real(r8), pointer :: t_veg_patch              (:)      => null()  ! patch vegetation temperature (Kelvin)
    !soilhydrology
    real(r8), pointer :: fracice_col       (:,:)           => null() ! col fractional impermeability (-)
    real(r8), pointer :: qflx_bot_col      (:)             => null() ! bottom of soil col flux, (mm/s)
    !atm2lnd
    real(r8), pointer :: forc_pbot_downscaled_col      (:) => null() ! downscaled atm pressure (Pa)
    real(r8), pointer :: forc_t_downscaled_col         (:) => null() ! downscaled atm temperature (Kelvin)
    real(r8), pointer :: n2_ppmv_col                   (:) => null() ! n2 concentration in ppmv
    real(r8), pointer :: o2_ppmv_col                   (:) => null() ! o2 concentration in ppmv
    real(r8), pointer :: ar_ppmv_col                   (:) => null() ! ar concentration in ppmv
    real(r8), pointer :: co2_ppmv_col                  (:) => null() ! co2 concentration in ppmv
    real(r8), pointer :: ch4_ppmv_col                  (:) => null() ! ch4 concentration in ppmv
    real(r8), pointer :: n2o_ppmv_col                  (:) => null() ! n2o concentration in ppmv
    real(r8), pointer :: no_ppmv_col                   (:) => null() ! no concentration in ppmv
    real(r8), pointer :: nh3_ppmv_col                  (:) => null() ! nh3 concentration in ppmv
    !canopystate
    real(r8) , pointer :: altmax_col               (:) => null() ! col maximum annual depth of thaw
    real(r8) , pointer :: altmax_lastyear_col      (:) => null() ! col prior year maximum annual depth of thaw
    real(r8),  pointer :: lbl_rsc_h2o_patch        (:) => null() ! laminar boundary layer resistance for water over dry leaf (s/m)
    real(r8) , pointer :: elai_patch               (:) => null() ! patch canopy one-sided leaf area index with burying by snow
    !chemstate
    real(r8), pointer :: soil_pH(:,:)    ! soil pH (-nlevsno+1:nlevgrnd)
    !soilstate
    real(r8), pointer :: bsw_col(:,:)                  => null() !Clapp and Hornberger "b" (nlevgrnd)
    real(r8), pointer :: watsat_col(:,:)               => null() !volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: eff_porosity_col(:,:)         => null() !effective porosity = porosity - vol_ice (nlevgrnd)
    real(r8), pointer :: soilpsi_col          (:,:)    => null() ! col soil water potential in each soil layer (MPa) (CN)
    real(r8), pointer :: cellorg_col          (:,:)    => null() ! col organic matter for gridcell containing column (1:nlevsoi)
    real(r8), pointer :: cellclay_col         (:,:)    => null() ! clay value for gridcell containing column (1:nlevsoi)
    real(r8), pointer :: cellsand_col         (:,:)    => null() ! sand value for gridcell containing column (1:nlevsoi)
    real(r8), pointer :: bd_col               (:,:)    => null() ! col bulk density of dry soil material [kg/m^3] (CN)
    real(r8), pointer :: watfc_col            (:,:)    => null() ! col volumetric soil water at field capacity (nlevsoi)
    real(r8), pointer :: sucsat_col           (:,:)    => null() ! col minimum soil suction (mm) (nlevgrnd)
    real(r8), pointer :: rootfr_patch         (:,:)    => null() ! patch fraction of roots in each soil layer (nlevgrnd)
    real(r8), pointer :: rr_patch(:,:)  => null()
    real(r8), pointer :: froot_prof_patch(:,:) => null()
    real(r8), pointer :: frootc_patch(:) => null()
    real(r8), pointer :: cn_scalar_patch(:) => null()
    real(r8), pointer :: cp_scalar_patch(:) => null()
    real(r8), pointer :: dic_prod_vr_col(:,:) => null()
    real(r8), pointer :: doc_prod_vr_col(:,:) => null()
    real(r8), pointer :: biochem_pmin_vr(:,:) => null()
    real(r8), pointer :: scalaravg_col(:) => null()
    real(r8), pointer :: dom_scalar_col(:) => null()
    real(r8), pointer :: solutionp_vr_col(:,:) => null()
    real(r8), pointer :: labilep_vr_col(:,:) => null()
    real(r8), pointer :: secondp_vr_col(:,:) => null()
    real(r8), pointer :: occlp_vr_col(:,:) => null()
    integer, pointer  :: lithotype_col(:) => null()
    !carbon fluxes
    type(betr_carbonflux_type) :: c12flx
    type(betr_carbonflux_type) :: c13flx
    type(betr_carbonflux_type) :: c14flx
    type(betr_nitrogenflux_type):: n14flx
    type(betr_phosphorusflux_type) :: p31flx
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, private :: InitCold
    procedure, public  :: summary
    procedure, public  :: reset
    procedure, public  :: frac_normalize
  end type betr_biogeophys_input_type

  public :: create_betr_biogeophys_input
contains

  function create_betr_biogeophys_input()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_biogeophys_input_type), pointer :: create_betr_biogeophys_input
    class(betr_biogeophys_input_type), pointer :: biogeophys_input

    allocate(biogeophys_input)
    create_betr_biogeophys_input => biogeophys_input

  end function create_betr_biogeophys_input
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
  use tracer_varcon, only : use_c13_betr, use_c14_betr
  implicit none
  class(betr_biogeophys_input_type)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  call this%c12flx%Init(bounds)
  call this%n14flx%Init(bounds)
  call this%p31flx%Init(bounds)
  if(use_c13_betr)call this%c13flx%Init(bounds)
  if(use_c14_betr)call this%c14flx%Init(bounds)

  call this%InitCold()
  end subroutine Init
  !------------------------------------------------------------------------


  subroutine InitCold(this)

  implicit none
  class(betr_biogeophys_input_type)  :: this

  this%n2_ppmv_col(:) = 7.8e5_r8
  this%o2_ppmv_col(:) = 2.1e5_r8
  this%ar_ppmv_col(:) = 0.9e5_r8
  this%co2_ppmv_col(:)= 400._r8
  this%ch4_ppmv_col(:)= 1.8_r8
  this%n2o_ppmv_col(:)= 250.e-3_r8
  this%no_ppmv_col(:) = 0._r8
  this%nh3_ppmv_col(:)= 0._r8
  this%forc_pbot_downscaled_col(:) = 1.013e5_r8
  this%forc_t_downscaled_col(:) = 288._r8
  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  implicit none
  class(betr_biogeophys_input_type)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  ! cnstate_vars
  allocate(this%isoilorder(begc:endc))  ! soil order
  allocate(this%lithoclass(begc:endc))  ! lithology class
  allocate(this%frootc_patch(begp:endp))
  allocate(this%cn_scalar_patch(begp:endp))
  allocate(this%cp_scalar_patch(begp:endp))
  !carbon flux
  allocate (this%annsum_npp_patch(  begp:endp))  !annual npp
  allocate (this%agnpp_patch(       begp:endp))
  allocate (this%bgnpp_patch(       begp:endp))
  allocate(this%rr_patch (begp:endp, lbj:ubj))
  allocate(this%froot_prof_patch(begp:endp, lbj:ubj))
  allocate(this%frac_loss_lit_to_fire_col(begc:endc)) !fraction of litter cnp loss through fire
  allocate(this%frac_loss_cwd_to_fire_col(begc:endc)) !fraction of cwd cnp loss through fire
  !waterstate
  allocate (this%frac_h2osfc_col (  begc:endc         ) ) ! col fractional area with surface water greater than zero
  allocate (this%finundated_col(    begc:endc         ) ) ! fraction of column that is inundated, this is for bgc caclulation in betr
  allocate (this%h2osoi_liq_col(    begc:endc,lbj:ubj ) ) !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
  allocate (this%h2osoi_ice_col(    begc:endc,lbj:ubj ) ) !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
  allocate (this%h2osoi_liqvol_col( begc:endc,lbj:ubj ) ) !volumetric liquid water content
  allocate (this%h2osoi_icevol_col( begc:endc,lbj:ubj ) ) !volumetric ice water content
  allocate (this%h2osoi_vol_col(    begc:endc,lbj:ubj ) ) !volumetric water content, total
  allocate (this%air_vol_col(       begc:endc,lbj:ubj ) ) !volume possessed by air
  allocate (this%rho_vap(           begc:endc,lbj:ubj ) ) !concentration of bulk water vapor, assume in equilibrium with the liquid water
  allocate (this%rhvap_soi(         begc:endc,lbj:ubj ) ) !soil relative humidity
  allocate (this%smp_l_col     (    begc:endc,lbj:ubj ) ) ! col liquid phase soil matric potential, mm

  !waterflux
  allocate(this%qflx_surf_col            (begc:endc         ) )  !surface runoff (mm H2O /s)
  allocate(this%qflx_dew_grnd_col        (begc:endc         ) ) ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
  allocate(this%qflx_dew_snow_col        (begc:endc         ) ) ! col surface dew added to snow pack (mm H2O /s) [+]
  allocate(this%qflx_sub_snow_vol_col    (begc:endc         ) )
  allocate(this%qflx_sub_snow_col        (begc:endc         ) ) ! col sublimation rate from snow pack (mm H2O /s) [+]
  allocate(this%qflx_h2osfc2topsoi_col   (begc:endc         ) ) ! col liquid water coming from surface standing water top soil (mm H2O/s)
  allocate(this%qflx_snow2topsoi_col     (begc:endc         ) ) ! col liquid water coming from residual snow to topsoil (mm H2O/s)
  allocate(this%qflx_tran_veg_patch      (begp:endp         ) )
  allocate(this%qflx_rootsoi_col         (begc:endc,lbj:ubj ) ) ! col root and soil water exchange [mm H2O/s] [+ into root]
  allocate(this%qflx_rootsoi_frac_patch  (begp:endp,lbj:ubj ) ) ! col root and soil water exchange [mm H2O/s] [+ into root]
  allocate(this%qflx_runoff_col          (begc:endc))
  !temperature
  allocate(this%t_soi_10cm               (begc:endc)) !soil temperature in top 10cm of soil (Kelvin)
  allocate(this%t_veg_patch              (begp:endp)) ! patch vegetation temperature (Kelvin)
  allocate(this%t_soisno_col             (begc:endc,lbj:ubj)) !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

  !soilhydrology
  allocate(this%qflx_bot_col      (begc:endc))        ! bottom of soil col flux, (mm/s)
  allocate(this%fracice_col       (begc:endc,lbj:ubj)) ! col fractional impermeability (-)

  !atm2lnd
  allocate(this%forc_pbot_downscaled_col      ( begc:endc)  ) ! downscaled atm pressure (Pa)
  allocate(this%forc_t_downscaled_col         ( begc:endc)  ) ! downscaled atm temperature (Kelvin)
  allocate(this%n2_ppmv_col(begc:endc))
  allocate(this%o2_ppmv_col(begc:endc))
  allocate(this%ar_ppmv_col(begc:endc))
  allocate(this%co2_ppmv_col(begc:endc))
  allocate(this%ch4_ppmv_col(begc:endc))
  allocate(this%n2o_ppmv_col(begc:endc))
  allocate(this%no_ppmv_col(begc:endc))
  allocate(this%nh3_ppmv_col(begc:endc))

  !canopystate
  allocate(this%altmax_col               (      begc:endc )  ) ! col maximum annual depth of thaw
  allocate(this%altmax_lastyear_col      (      begc:endc )  ) ! col prior year maximum annual depth of thaw
  allocate(this%lbl_rsc_h2o_patch        (      begp:endp )  ) ! laminar boundary layer resistance for water over dry leaf (s/m)
  allocate(this%elai_patch               (      begp:endp ) )  ! patch canopy one-sided leaf area index with burying by snow

  !chemstate
  allocate(this%soil_pH(begc:endc,lbj:ubj))    ! soil pH (-nlevsno+1:nlevgrnd)

  !soilstate
  allocate(this%bsw_col(begc:endc,lbj:ubj)         )       !Clapp and Hornberger "b" (nlevgrnd)
  allocate(this%watsat_col(begc:endc,lbj:ubj)       )      !volumetric soil water at saturation (porosity) (nlevgrnd)
  allocate(this%eff_porosity_col(begc:endc,lbj:ubj)   )    !effective porosity = porosity - vol_ice (nlevgrnd)
  allocate(this%soilpsi_col          (begc:endc,lbj:ubj) ) ! col soil water potential in each soil layer (MPa) (CN)
  allocate(this%cellorg_col          (begc:endc,lbj:ubj) ) ! col organic matter for gridcell containing column (1:nlevsoi)
  allocate(this%cellclay_col         (begc:endc,lbj:ubj) ) ! clay value for gridcell containing column (1:nlevsoi)
  allocate(this%cellsand_col         (begc:endc,lbj:ubj) ) ! sand value for gridcell containing column (1:nlevsoi)
  allocate(this%bd_col               (begc:endc,lbj:ubj) ) ! col bulk density of dry soil material [kg/m^3] (CN)
  allocate(this%watfc_col            (begc:endc,lbj:ubj) ) ! col volumetric soil water at field capacity (nlevsoi)
  allocate(this%sucsat_col           (begc:endc,lbj:ubj) ) ! col minimum soil suction (mm) (nlevgrnd)
  allocate(this%rootfr_patch         (begp:endp,lbj:ubj) ) ! patch fraction of roots in each soil layer (nlevgrnd)
  allocate(this%lithotype_col        (begc:endc))
  allocate(this%solutionp_vr_col(begc:endc,lbj:ubj))
  allocate(this%labilep_vr_col(begc:endc,lbj:ubj))
  allocate(this%secondp_vr_col(begc:endc,lbj:ubj))
  allocate(this%occlp_vr_col(begc:endc,lbj:ubj))

  allocate(this%dic_prod_vr_col(begc:endc,lbj:ubj))
  allocate(this%doc_prod_vr_col(begc:endc,lbj:ubj))
  allocate(this%biochem_pmin_vr(begc:endc,lbj:ubj))
  allocate(this%scalaravg_col(begc:endc))
  allocate(this%dom_scalar_col(begc:endc)); this%dom_scalar_col(:) = 1._r8
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  use tracer_varcon, only : use_c13_betr, use_c14_betr
  implicit none
  class(betr_biogeophys_input_type)  :: this
  real(r8), intent(in) :: value_column


  call this%c12flx%reset(value_column)
  call this%n14flx%reset(value_column)
  call this%p31flx%reset(value_column)
  if(use_c13_betr)call this%c13flx%reset(value_column)
  if(use_c14_betr)call this%c14flx%reset(value_column)

  end subroutine reset

  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz)
  use tracer_varcon, only : use_c13_betr, use_c14_betr
  implicit none
  class(betr_biogeophys_input_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)

  !integrate
  call this%c12flx%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj))

  if(use_c13_betr)then
     call this%c13flx%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj))
  endif

  if(use_c14_betr)then
     call this%c14flx%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj))
  endif

  call this%n14flx%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj))

  call this%p31flx%summary(bounds, lbj, ubj, dz)
  end subroutine summary

  !------------------------------------------------------------------------
  subroutine frac_normalize(this, npfts, lbj, ubj)
  !
  ! DESCRIPTION
  ! normalize some fractional data
  implicit none
  class(betr_biogeophys_input_type),intent(inout)  :: this
  integer , intent(in) :: npfts
  integer , intent(in) :: lbj, ubj

  integer :: j, p
  real(r8) :: nrmscal

  do j = lbj, ubj
    nrmscal= sum(this%qflx_rootsoi_frac_patch(1:npfts,j))
    if(nrmscal/=0._r8)then
      do p = 1, npfts
        this%qflx_rootsoi_frac_patch(p,j) = this%qflx_rootsoi_frac_patch(p,j) / nrmscal
      end do
    endif
  enddo

  end subroutine frac_normalize
end module BeTR_biogeophysInputType
