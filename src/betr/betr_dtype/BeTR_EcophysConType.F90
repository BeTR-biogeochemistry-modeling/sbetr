module BeTR_EcophysConType

  ! !USES:
  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use betr_decompMod  , only : betr_bounds_type
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  !
  implicit none
  private
  !
  !
  ! !PUBLIC TYPES:
  type, public :: betr_ecophyscon_type
     integer , pointer :: noveg         (:)        => null() ! value for not vegetated
     integer , pointer :: tree          (:)        => null() ! tree or not?
     real(r8), pointer :: smpso         (:)        => null() ! soil water potential at full stomatal opening (mm)
     real(r8), pointer :: smpsc         (:)        => null() ! soil water potential at full stomatal closure (mm)
     real(r8), pointer :: fnitr         (:)        => null() ! foliage nitrogen limitation factor (-)
     real(r8), pointer :: foln          (:)        => null() ! foliage nitrogen (%)
     real(r8), pointer :: dleaf         (:)        => null() ! characteristic leaf dimension (m)
     real(r8), pointer :: c3psn         (:)        => null() ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), pointer :: xl            (:)        => null() ! leaf/stem orientation index
     real(r8), pointer :: z0mr          (:)        => null() ! ratio of momentum roughness length to canopy top height (-)
     real(r8), pointer :: displar       (:)        => null() ! ratio of displacement height to canopy top height (-)
     real(r8), pointer :: roota_par     (:)        => null() ! rooting distribution parameter [1/m]
     real(r8), pointer :: rootb_par     (:)        => null() ! rooting distribution parameter [1/m]
     real(r8), pointer :: rootprof_beta (:)        => null() ! rooting distribution parameter for C and N inputs [unitless]
     real(r8), pointer :: dwood         (:)        => null() ! wood density (gC/m3)
     real(r8), pointer :: slatop        (:)        => null() ! specific leaf area at top of canopy, projected area basis [m^2/gC]
     real(r8), pointer :: dsladlai      (:)        => null() ! dSLA/dLAI, projected area basis [m^2/gC]
     real(r8), pointer :: leafcn        (:)        => null() ! leaf C:N (gC/gN)
     real(r8), pointer :: flnr          (:)        => null() ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
     real(r8), pointer :: woody         (:)        => null() ! binary flag for woody lifeform (1=woody, 0=not woody)
     real(r8), pointer :: lflitcn       (:)        => null() ! leaf litter C:N (gC/gN)
     real(r8), pointer :: frootcn       (:)        => null() ! fine root C:N (gC/gN)
     real(r8), pointer :: livewdcn      (:)        => null() ! live wood (phloem and ray parenchyma) C:N (gC/gN)
     real(r8), pointer :: deadwdcn      (:)        => null() ! dead wood (xylem and heartwood) C:N (gC/gN)
     real(r8), pointer :: graincn       (:)        => null() ! grain C:N (gC/gN) for prognostic crop model
     real(r8), pointer :: froot_leaf    (:)        => null() ! allocation parameter: new fine root C per new leaf C (gC/gC)
     real(r8), pointer :: stem_leaf     (:)        => null() ! allocation parameter: new stem c per new leaf C (gC/gC)
     real(r8), pointer :: croot_stem    (:)        => null() ! allocation parameter: new coarse root C per new stem C (gC/gC)
     real(r8), pointer :: flivewd       (:)        => null() ! allocation parameter: fraction of new wood that is live
                                                             ! (phloem and ray parenchyma) (no units)
     real(r8), pointer :: fcur          (:)        => null() ! allocation parameter: fraction of allocation that goes
                                                             ! to currently displayed growth, remainder to storage
     real(r8), pointer :: lf_flab       (:)        => null() ! leaf litter labile fraction
     real(r8), pointer :: lf_fcel       (:)        => null() ! leaf litter cellulose fraction
     real(r8), pointer :: lf_flig       (:)        => null() ! leaf litter lignin fraction
     real(r8), pointer :: fr_flab       (:)        => null() ! fine root litter labile fraction
     real(r8), pointer :: fr_fcel       (:)        => null() ! fine root litter cellulose fraction
     real(r8), pointer :: fr_flig       (:)        => null() ! fine root litter lignin fraction
     real(r8), pointer :: leaf_long     (:)        => null() ! leaf longevity (yrs)
     real(r8), pointer :: evergreen     (:)        => null() ! binary flag for evergreen leaf habit (0 or 1)
     real(r8), pointer :: stress_decid  (:)        => null() ! binary flag for stress-deciduous leaf habit (0 or 1)
     real(r8), pointer :: season_decid  (:)        => null() ! binary flag for seasonal-deciduous leaf habit (0 or 1)
     real(r8), pointer :: cc_leaf       (:)        => null() ! combustion completeness factor for leaf (0 to 1)
     real(r8), pointer :: cc_lstem      (:)        => null() ! combustion completeness factor for live stem (0 to 1)
     real(r8), pointer :: cc_dstem      (:)        => null() ! combustion completeness factor for dead stem (0 to 1)
     real(r8), pointer :: cc_other      (:)        => null() ! combustion completeness factor for other plant tissues (0 to 1)
     real(r8), pointer :: fm_leaf       (:)        => null() ! fire-related mortality factor for leaf (0 to 1)
     real(r8), pointer :: fm_lstem      (:)        => null() ! fire-related mortality factor for live stem (0 to 1)
     real(r8), pointer :: fm_dstem      (:)        => null() ! fire-related mortality factor for dead stem (0 to 1)
     real(r8), pointer :: fm_other      (:)        => null() ! fire-related mortality factor for other plant tissues (0 to 1)
     real(r8), pointer :: fm_root       (:)        => null() ! fire-related mortality factor for fine roots (0 to 1)
     real(r8), pointer :: fm_lroot      (:)        => null() ! fire-related mortality factor for live roots (0 to 1)
     real(r8), pointer :: fm_droot      (:)        => null() ! fire-related mortality factor for dead roots (0 to 1)
     real(r8), pointer :: fertnitro     (:)        => null() ! fertilizer applied (crop)
     real(r8), pointer :: fleafcn       (:)        => null() ! C:N during grain fill; leaf (crop)
     real(r8), pointer :: ffrootcn      (:)        => null() ! C:N during grain fill; froot (crop)
     real(r8), pointer :: fstemcn       (:)        => null() ! C:N during grain fill; stem (crop)
     real(r8), pointer :: presharv      (:)        => null() ! porportion of residue harvested (crop)
     real(r8), pointer :: convfact      (:)        => null() ! converstion factor to bu/acre (crop)
     real(r8), pointer :: fyield        (:)        => null() ! fraction of grain that is actually harvested (crop)

     real(r8), pointer :: leafcp        (:)        => null() ! leaf C:P (gC/gP)
     real(r8), pointer :: lflitcp       (:)        => null() ! leaf litter C:P (gC/gP)
     real(r8), pointer :: frootcp       (:)        => null() ! fine root C:P (gC/gP)
     real(r8), pointer :: livewdcp      (:)        => null() ! live wood (phloem and ray parenchyma) C:P (gC/gP)
     real(r8), pointer :: deadwdcp      (:)        => null() ! dead wood (xylem and heartwood) C:P (gC/gP)
     real(r8), pointer :: graincp       (:)        => null() ! grain C:P (gC/gP) for prognostic crop model

     ! pft dependent parameters for phosphorus for nutrient competition
     real(r8), pointer :: vmax_plant_nh4(:)        => null() ! vmax for plant nh4 uptake
     real(r8), pointer :: vmax_plant_no3(:)        => null() ! vmax for plant no3 uptake
     real(r8), pointer :: vmax_plant_p(:)          => null() ! vmax for plant p uptake
     real(r8), pointer :: vmax_minsurf_p_vr(:,:)   => null() ! vmax for p adsorption
     real(r8), pointer :: km_plant_nh4(:)          => null() ! km for plant nh4 uptake
     real(r8), pointer :: km_plant_no3(:)          => null() ! km for plant no3 uptake
     real(r8), pointer :: km_plant_p(:)            => null() ! km for plant p uptake
     real(r8), pointer :: km_minsurf_p_vr(:,:)     => null() ! km for p adsorption
     real(r8)          :: km_decomp_nh4                      ! km for microbial decomposer nh4 uptake
     real(r8)          :: km_decomp_no3                      ! km for microbial decomposer no3 uptake
     real(r8)          :: km_decomp_p                        ! km for microbial decomposer p uptake
     real(r8)          :: km_nit                             ! km for nitrifier nh4 uptake
     real(r8)          :: km_den                             ! km for denitrifier no3 uptake
     real(r8), pointer :: decompmicc_patch_vr(:,:) => null() ! microbial decomposer biomass gc/m3
     real(r8)          :: vmax_nfix                          ! vmax of symbiotic n2 fixation
     real(r8)          :: km_nfix                            ! km of symbiotic n2 fixation
     real(r8), pointer :: vmax_ptase_vr(:)         => null() ! vmax of biochemical p production
     real(r8)          :: km_ptase                           ! km of biochemical p production
     real(r8)          :: lamda_ptase                        ! critical value that incur biochemical production
     real(r8), pointer :: i_vc(:)                  => null() ! intercept of photosynthesis vcmax ~ leaf n content regression model
     real(r8), pointer :: s_vc(:)                  => null() ! slope of photosynthesis vcmax ~ leaf n content regression model
   contains
     procedure, public  :: Init
     procedure, private :: InitAllocate
  end type betr_ecophyscon_type


contains

  subroutine Init(this, bounds, numpft, nsoilorder)
  implicit none
  class(betr_Ecophyscon_type)             :: this
  integer,                intent(in) :: numpft
  type(betr_bounds_type), intent(in) :: bounds
  integer,                intent(in) :: nsoilorder

  call this%InitAllocate(numpft, bounds%ubj, nsoilorder)

  end subroutine Init
  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, numpft, ubj, nsoilorder)
    !
    ! !USES:

    implicit none
    class(betr_Ecophyscon_type) :: this
    integer, intent(in) :: numpft
    integer, intent(in) :: ubj
    integer, intent(in) :: nsoilorder
    !
    ! !LOCAL VARIABLES:
    integer :: m, ib, j
    !------------------------------------------------------------------------

    allocate(this%noveg         (0:numpft)) ; this%noveg        (:)   =huge(1)
    allocate(this%tree          (0:numpft)) ; this%tree         (:)   =huge(1)
    allocate(this%smpso         (0:numpft)) ; this%smpso        (:)   =nan
    allocate(this%smpsc         (0:numpft)) ; this%smpsc        (:)   =nan
    allocate(this%fnitr         (0:numpft)) ; this%fnitr        (:)   =nan
    allocate(this%foln          (0:numpft)) ; this%foln         (:)   =nan
    allocate(this%dleaf         (0:numpft)) ; this%dleaf        (:)   =nan
    allocate(this%c3psn         (0:numpft)) ; this%c3psn        (:)   =nan
    allocate(this%xl            (0:numpft)) ; this%xl           (:)   =nan
    allocate(this%z0mr          (0:numpft)) ; this%z0mr         (:)   =nan
    allocate(this%displar       (0:numpft)) ; this%displar      (:)   =nan
    allocate(this%roota_par     (0:numpft)) ; this%roota_par    (:)   =nan
    allocate(this%rootb_par     (0:numpft)) ; this%rootb_par    (:)   =nan
    allocate(this%slatop        (0:numpft)) ; this%slatop       (:)   =nan
    allocate(this%dsladlai      (0:numpft)) ; this%dsladlai     (:)   =nan
    allocate(this%leafcn        (0:numpft)) ; this%leafcn       (:)   =nan
    allocate(this%flnr          (0:numpft)) ; this%flnr         (:)   =nan
    allocate(this%woody         (0:numpft)) ; this%woody        (:)   =nan
    allocate(this%lflitcn       (0:numpft)) ; this%lflitcn      (:)   =nan
    allocate(this%frootcn       (0:numpft)) ; this%frootcn      (:)   =nan
    allocate(this%livewdcn      (0:numpft)) ; this%livewdcn     (:)   =nan
    allocate(this%deadwdcn      (0:numpft)) ; this%deadwdcn     (:)   =nan
    allocate(this%graincn       (0:numpft)) ; this%graincn      (:)   =nan
    allocate(this%froot_leaf    (0:numpft)) ; this%froot_leaf   (:)   =nan
    allocate(this%stem_leaf     (0:numpft)) ; this%stem_leaf    (:)   =nan
    allocate(this%croot_stem    (0:numpft)) ; this%croot_stem   (:)   =nan
    allocate(this%flivewd       (0:numpft)) ; this%flivewd      (:)   =nan
    allocate(this%fcur          (0:numpft)) ; this%fcur         (:)   =nan
    allocate(this%lf_flab       (0:numpft)) ; this%lf_flab      (:)   =nan
    allocate(this%lf_fcel       (0:numpft)) ; this%lf_fcel      (:)   =nan
    allocate(this%lf_flig       (0:numpft)) ; this%lf_flig      (:)   =nan
    allocate(this%fr_flab       (0:numpft)) ; this%fr_flab      (:)   =nan
    allocate(this%fr_fcel       (0:numpft)) ; this%fr_fcel      (:)   =nan
    allocate(this%fr_flig       (0:numpft)) ; this%fr_flig      (:)   =nan
    allocate(this%leaf_long     (0:numpft)) ; this%leaf_long    (:)   =nan
    allocate(this%evergreen     (0:numpft)) ; this%evergreen    (:)   =nan
    allocate(this%stress_decid  (0:numpft)) ; this%stress_decid (:)   =nan
    allocate(this%season_decid  (0:numpft)) ; this%season_decid (:)   =nan
    allocate(this%dwood         (0:numpft)) ; this%dwood        (:)   =nan
    allocate(this%rootprof_beta (0:numpft)) ; this%rootprof_beta(:)   =nan
    allocate(this%fertnitro     (0:numpft)) ; this%fertnitro    (:)   =nan
    allocate(this%fleafcn       (0:numpft)) ; this%fleafcn      (:)   =nan
    allocate(this%ffrootcn      (0:numpft)) ; this%ffrootcn     (:)   =nan
    allocate(this%fstemcn       (0:numpft)) ; this%fstemcn      (:)   =nan
    allocate(this%presharv      (0:numpft)) ; this%presharv     (:)   =nan
    allocate(this%convfact      (0:numpft)) ; this%convfact     (:)   =nan
    allocate(this%fyield        (0:numpft)) ; this%fyield       (:)   =nan


    allocate(this%leafcp        (0:numpft))               ; this%leafcp       (:)   =nan
    allocate(this%lflitcp       (0:numpft))               ; this%lflitcp      (:)   =nan
    allocate(this%frootcp       (0:numpft))               ; this%frootcp      (:)   =nan
    allocate(this%livewdcp      (0:numpft))               ; this%livewdcp     (:)   =nan
    allocate(this%deadwdcp      (0:numpft))               ; this%deadwdcp     (:)   =nan
    allocate(this%graincp       (0:numpft))               ; this%graincp      (:)   =nan

    allocate( this%vmax_plant_nh4(0:numpft))              ; this%vmax_plant_nh4(:)        =nan
    allocate( this%vmax_plant_no3(0:numpft))              ; this%vmax_plant_no3(:)        =nan
    allocate( this%vmax_plant_p(0:numpft))                ; this%vmax_plant_p(:)          =nan
    allocate( this%vmax_minsurf_p_vr(0:nsoilorder,1:ubj)) ; this%vmax_minsurf_p_vr(:,:)   =nan
    allocate( this%km_plant_nh4(0:numpft))                ; this%km_plant_nh4(:)          =nan
    allocate( this%km_plant_no3(0:numpft))                ; this%km_plant_no3(:)          =nan
    allocate( this%km_plant_p(0:numpft))                  ; this%km_plant_p(:)            =nan
    allocate( this%km_minsurf_p_vr(0:nsoilorder,1:ubj))   ; this%km_minsurf_p_vr(:,:)     =nan
    allocate( this%decompmicc_patch_vr(0:numpft,1:ubj))   ; this%decompmicc_patch_vr(:,:) =nan
    allocate( this%vmax_ptase_vr(1:ubj))                  ; this%vmax_ptase_vr(:)         =nan
    allocate( this%i_vc(0:numpft))                        ; this%i_vc(:)                  =nan
    allocate( this%s_vc(0:numpft))                        ; this%s_vc(:)                  =nan


  end subroutine InitAllocate

end module BeTR_EcophysConType
