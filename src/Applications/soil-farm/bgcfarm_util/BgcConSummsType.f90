module BgcConSummsType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
implicit none

 private

 type, public :: BgcConSumms_type

  ! Microbial parameters
    real(r8) :: gmax_mic
    real(r8) :: yld_mic
    real(r8) :: yld_enz
    real(r8) :: yld_res
    real(r8) :: fenz2poly
    real(r8) :: minsite
    real(r8) :: mic_transp
    real(r8) :: decay_mic0
    real(r8) :: decay_enz
    real(r8) :: pmax_enz
    real(r8) :: ea_vmax_mic
    real(r8) :: ea_vmax_enz
    real(r8) :: ea_kaff_mono_mic
    real(r8) :: ea_kaff_enz_poly
    real(r8) :: ea_mr_mic
    real(r8) :: ea_kappa_mic
    real(r8) :: ea_kaff_mono_msurf
    real(r8) :: ea_kaff_enz_msurf
    real(r8) :: ref_vmax_enz
    real(r8) :: ref_kaff_enz_poly
    real(r8) :: ref_kaff_enz_msurf
    real(r8) :: ref_kaff_mono_msurf
    real(r8) :: ref_mr_mic
    real(r8) :: ref_kappa_mic
    real(r8) :: ref_kaff_mono_mic
    real(r8) :: ref_vmax_mic

  !decomposition
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding

  real(r8) :: rf_l1s1_bgc !CO2 from metabolic C decomp.
  real(r8) :: rf_l2s1_bgc
  real(r8) :: rf_l3s2_bgc
  real(r8) :: rf_s2s1_bgc
  real(r8) :: rf_s3s1_bgc

  real(r8) :: cwd_fcel_bgc
  real(r8) :: cwd_flig_bgc
  real(r8) :: lwd_fcel_bgc
  real(r8) :: lwd_flig_bgc
  real(r8) :: fwd_fcel_bgc
  real(r8) :: fwd_flig_bgc

  !nitrification-denitrification
  real(r8) :: nitrif_n2o_loss_frac
  real(r8) :: organic_max
  real(r8) :: rij_kro_a
  real(r8) :: rij_kro_alpha
  real(r8) :: rij_kro_beta
  real(r8) :: rij_kro_gamma
  real(r8) :: rij_kro_delta
  real(r8) :: surface_tension_water

  real(r8) :: init_cn_met
  real(r8) :: init_cn_cel
  real(r8) :: init_cn_lig
  real(r8) :: init_cn_cwd
  real(r8) :: init_cn_lwd
  real(r8) :: init_cn_fwd
  real(r8) :: init_cn_mic
  real(r8) :: init_cn_res
  real(r8) :: init_cn_enz
  real(r8) :: init_cn_mono
  real(r8) :: init_cn_poly

  real(r8) :: init_cp_met
  real(r8) :: init_cp_cel
  real(r8) :: init_cp_lig
  real(r8) :: init_cp_cwd
  real(r8) :: init_cp_lwd
  real(r8) :: init_cp_fwd
  real(r8) :: init_cp_mic
  real(r8) :: init_cp_res
  real(r8) :: init_cp_enz
  real(r8) :: init_cp_mono
  real(r8) :: init_cp_poly

  real(r8) :: init_cc13_met
  real(r8) :: init_cc13_cel
  real(r8) :: init_cc13_lig
  real(r8) :: init_cc13_cwd
  real(r8) :: init_cc13_lwd
  real(r8) :: init_cc13_fwd
  real(r8) :: init_cc13_mic
  real(r8) :: init_cc13_res
  real(r8) :: init_cc13_enz
  real(r8) :: init_cc13_mono
  real(r8) :: init_cc13_poly

  real(r8) :: init_cc14_met
  real(r8) :: init_cc14_cel
  real(r8) :: init_cc14_lig
  real(r8) :: init_cc14_cwd
  real(r8) :: init_cc14_lwd
  real(r8) :: init_cc14_fwd
  real(r8) :: init_cc14_mic
  real(r8) :: init_cc14_res
  real(r8) :: init_cc14_enz
  real(r8) :: init_cc14_mono
  real(r8) :: init_cc14_poly
  real(r8) :: c14decay_const
  real(r8) :: c14decay_som_const
  real(r8) :: c14decay_dom_const
  real(r8) :: c14decay_bm_const

  logical :: use_c13
  logical :: use_c14
  logical :: nop_limit                                              !switch for P limitation
  logical :: non_limit
  !ECA nutrient competition
  real(r8), pointer :: vmax_minp_soluble_to_secondary(:)  => null() !maximum conversion rate of soluble P into secondary P

  !inorganic phosphorus cycling
  real(r8), pointer :: frac_p_sec_to_sol(:)    => null()   !fraction of released secondary phosphorus that goes into soluble form
  real(r8), pointer :: minp_secondary_decay(:) => null()   !decay rate of secondary phosphorus
  !real(r8), pointer :: spinup_factor(:)
 contains
   procedure, public  :: Init
   procedure, private :: InitAllocate
   procedure, private :: set_defpar_default
   !procedure, public  :: apply_spinup_factor
   !procedure, private :: ReadNamelist
 end type BgcConSumms_type

 type(BgcConSumms_type), public :: bgc_con_summs
contains
  !--------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  !use betr_ctrl     , only : betr_spinup_state
  implicit none
  class(BgcConSumms_type), intent(inout) :: this
  character(len=betr_namelist_buffer_size_ext) , intent(in)    :: namelist_buffer
  type(betr_status_type)                   , intent(out) :: bstatus

  call bstatus%reset()

  call this%InitAllocate()
  call this%set_defpar_default()

  !update parameter from namelist - add back later
  !call this%ReadNamelist(namelist_buffer, bstatus)
  !if(betr_spinup_state/=0)then
  !  call this%apply_spinup_factor()
  !endif
  end subroutine Init
  !--------------------------------------------------------------------
  subroutine InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(BgcConSumms_type), intent(inout) :: this

  allocate(this%minp_secondary_decay(betr_max_soilorder))
  allocate(this%vmax_minp_soluble_to_secondary(betr_max_soilorder))
  allocate(this%frac_p_sec_to_sol(betr_max_soilorder))
  !allocate(this%spinup_factor(9))
  !the following will be actually calculated from CNP bgc
  end subroutine InitAllocate
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)
  implicit none
  class(BgcConSumms_type), intent(inout) :: this
  real(r8) :: half_life

  half_life = 5568._r8 ! yr
  half_life = half_life * 86400._r8 * 365._r8
  this%c14decay_const = - log(0.5_r8) / half_life
  this%c14decay_som_const  =this%c14decay_const
  this%c14decay_dom_const  =this%c14decay_const
  this%c14decay_Bm_const  =this%c14decay_const
  ! Parameters
    this%gmax_mic  = 0.1025             ! Maximum microbial growth rate (1/day)
    this%yld_mic = 0.8                  ! Growth efficiency of microbes (g mic/g res)
    this%yld_enz = 0.8                  ! Growth efficiency of enzymes (g enz/g res)
    this%yld_res = 0.5                  ! Assimilation efficiency from monomer uptake (g res/g mono)
    this%fenz2poly = 0.2                ! Proportion of degraded exoenzyme into polymers (g poly/g enz)
    this%minsite = 1000                 ! Abundance of mineral surface              (g C surface/m3)
    this%mic_transp = 0.05              ! Scaling factor between transporter and microbial structural biomass
    this%decay_mic0 = 0.01314           ! Reference microbial death rate (1/day)
    this%decay_enz = 0.0061              ! Enzyme turnover tate (1/day)
    this%pmax_enz = 0.0019              ! Maximum enzyme production rate (1/day)

  ! Set up parameters for activation energy of different processes
    this%ea_vmax_mic            = 45000     ! Ea for maximum rate of monomer uptake (K)
    this%ea_vmax_enz            = 45000     ! Ea for maximum rate of polymer degradation (K)
    this%ea_kaff_mono_mic       = 1804.086  ! Ea for monomer-microbe affinity (K)
    this%ea_kaff_enz_poly       = 1804.086  ! Ea for polymer-enzyme affinity (K)
    this%ea_mr_mic              = 60000     ! Ea for maintenance (K)
    this%ea_kappa_mic           = 60000     ! Ea for reserve export (K)
    this%ea_kaff_mono_msurf     = 10000     ! Ea for monomer-mineral affinity (K)
    this%ea_kaff_enz_msurf      = 10000     ! Ea for enzyme-mineral affinity (K)

  ! Enzymes
    this%ref_vmax_enz           = 2.4133    ! Maximum rate of polymer degradation (1/day)
    this%ref_kaff_enz_poly      = 200       ! Affinity parameter for enzymatic polymer degradation (g enzymes/m3)
    this%ref_kaff_enz_msurf     = 50        ! Affinity parameter for surface adsorption of enzymes (g enzymes/m3)
  
  ! Monomer adsorption
    this%ref_kaff_mono_msurf   = 25        ! Affinity parameter or mineral surface adsorption of monomers (g monomers/m3)

  ! Microbes
    this%ref_mr_mic             = 0.0231      ! Microbial maintenance rate (1/day)
    this%ref_kappa_mic          = 0.0537      ! Reserve turnover rate (1/day)
    this%ref_kaff_mono_mic      = 1         ! Affinity parameter for microbial monomer uptake (g monomers/m3)
    this%ref_vmax_mic           = 10.9343    ! Maximum rate of monomer assimilation (1/day) 

  !decomposition
  this%Q10                   = 2._r8
  this%froz_q10              = 10._r8
  this%decomp_depth_efolding = 1._r8

  !following is based on Table 15.4 in CLM4.5 tech note
  this%rf_l1s1_bgc           = 0.55_r8
  this%rf_l2s1_bgc           = 0.5_r8
  this%rf_l3s2_bgc           = 0.5_r8
  this%cwd_fcel_bgc          = 0.76_r8
  this%cwd_flig_bgc          = 0.24_r8

  this%init_cn_met  = 90._r8  !mass based
  this%init_cn_cel  = 90._r8  !mass based
  this%init_cn_lig  = 90._r8  !mass based
  this%init_cn_cwd  = 90._r8  !mass based
  this%init_cn_mic = 8._r8   !mass based
  this%init_cn_res = 8._r8   !mass based
  this%init_cn_enz = 3._r8   !mass based
  this%init_cn_mono = 11._r8  !mass based
  this%init_cn_poly = 11._r8  !mass based !*** check all these ratios!

  this%init_cp_met  = 1600._r8
  this%init_cp_cel  = 2000._r8
  this%init_cp_lig  = 2500._r8
  this%init_cp_cwd  = 4500._r8
  this%init_cp_mic = 110._r8 !mass based
  this%init_cp_res = 110._r8 !mass based
  this%init_cp_enz = 110._r8 !mass based
  this%init_cp_mono = 320._r8 !mass based
  this%init_cp_poly = 114._r8 !mass based

  !nitrification-denitrification
  this%nitrif_n2o_loss_frac  = 1.e-4_r8   !Arah and Vinten, 1995
  this%organic_max           = 160._r8    !organic matter content (kg/m3) where soil is assumed to act like peat
  this%rij_kro_a             = 1.5e-10_r8 ! Arah and Vinten, 1995
  this%rij_kro_alpha         = 1.26_r8    ! Arah and Vinten, 1995
  this%rij_kro_beta          = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_gamma         = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_delta         = 0.85_r8    ! Arah and Vinten, 1995
  this%surface_tension_water = 73.e-3_r8  ! (J/m^2), Arah and Vinten, 1995

  !ECA nutrient competition
  this%vmax_minp_soluble_to_secondary(:) = 1.e-5_r8  !1/s
  !inorganic phosphorus cycling
  !Note: (1._r8-frac_p_sec_to_sol)*minp_secondary_decay = occlusion rate
  this%frac_p_sec_to_sol(:)              = 0.9_r8    !fraction of released secondary phosphorus that goes into soluble form
  this%minp_secondary_decay(:)           = 1.e-5_r8  !decay rate of secondary phosphorus, 1/s

  this%use_c13 = .false.
  this%use_c14 = .false.
  this%nop_limit=.false.
  this%non_limit=.false.
  this%init_cc13_met = 0._r8
  this%init_cc13_cel = 0._r8
  this%init_cc13_lig = 0._r8
  this%init_cc13_cwd = 0._r8
  this%init_cc13_mic = 0._r8
  this%init_cc13_res = 0._r8
  this%init_cc13_enz = 0._r8
  this%init_cc13_mono= 0._r8
  this%init_cc13_poly= 0._r8

  this%init_cc14_met = 0._r8
  this%init_cc14_cel = 0._r8
  this%init_cc14_lig = 0._r8
  this%init_cc14_cwd = 0._r8
  this%init_cc14_mic = 0._r8
  this%init_cc14_res = 0._r8
  this%init_cc14_enz = 0._r8
  this%init_cc14_mono= 0._r8
  this%init_cc14_poly= 0._r8
  end subroutine set_defpar_default
!--------------------------------------------------------------------

end module BgcConSummsType
