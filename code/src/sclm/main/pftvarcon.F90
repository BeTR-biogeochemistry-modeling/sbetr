module pftvarcon
  use shr_kind_mod, only : r8 => shr_kind_r8
 use clm_varpar  , only : mxpft
implicit none


  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  real(r8), allocatable :: crop(:)        !crop pft: 0. = not crop, 1. = crop pft
  integer :: noveg                  !value for not vegetated

  real(r8), allocatable :: i_vc(:)             ! intercept of photosynthesis vcmax ~ leaf N content regression model
  real(r8), allocatable :: s_vc(:)             ! slope of photosynthesis vcmax ~ leaf N content regression model
  real(r8), allocatable :: leafcn_obs(:)       !leaf C:N [gC/gN]
  real(r8), allocatable :: frootcn_obs(:)      !fine root C:N (gC/gN)
  real(r8), allocatable :: livewdcn_obs(:)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8), allocatable :: deadwdcn_obs(:)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8), allocatable :: leafcp_obs(:)       !leaf C:P [gC/gP]
  real(r8), allocatable :: frootcp_obs(:)      !fine root C:P (gC/gP)
  real(r8), allocatable :: livewdcp_obs(:)     !live wood (phloem and ray parenchyma) C:P (gC/gP)
  real(r8), allocatable :: deadwdcp_obs(:)     !dead wood (xylem and heartwood) C:P (gC/gP)
  real(r8), allocatable :: VMAX_PTASE_vr(:)    ! VMAX of biochemical P production
  real(r8)              :: KM_PTASE            ! KM of biochemical P production
  real(r8)              :: lamda_ptase         ! critical value that incur biochemical production
  real(r8), allocatable :: VMAX_PLANT_NH4(:)   ! VMAX for plant NH4 uptake
  real(r8), allocatable :: VMAX_PLANT_NO3(:)   ! VMAX for plant NO3 uptake
  real(r8), allocatable :: VMAX_PLANT_P(:)     ! VMAX for plant P uptake
  real(r8), allocatable :: VMAX_MINSURF_P_vr(:,:)! VMAX for P adsorption -> move to soilorder_varcon
  real(r8), allocatable :: KM_PLANT_NH4(:)     ! KM for plant NH4 uptake
  real(r8), allocatable :: KM_PLANT_NO3(:)     ! KM for plant NO3 uptake
  real(r8), allocatable :: KM_PLANT_P(:)       ! KM for plant P uptake
  real(r8), allocatable :: KM_MINSURF_P_vr(:,:)! KM for P adsorption -> move to soilorder_varcon
  real(r8)              :: KM_DECOMP_NH4       ! KM for microbial decomposer NH4 uptake
  real(r8)              :: KM_DECOMP_NO3       ! KM for microbial decomposer NO3 uptake
  real(r8)              :: KM_DECOMP_P         ! KM for microbial decomposer P uptake
  real(r8)              :: KM_NIT              ! KM for nitrifier NH4 uptake
  real(r8)              :: KM_DEN              ! KM for denitrifier NO3 uptake
  real(r8), allocatable :: decompmicc_patch_vr(:,:) ! microbial decomposer biomass gC/m3
  real(r8)              :: VMAX_NFIX           ! VMAX of symbiotic N2 fixation
  real(r8)              :: KM_NFIX             ! KM of symbiotic N2 fixation
contains



!-----------------------------------------------------------------------
subroutine pftconrd

! !ARGUMENTS:
implicit none



    allocate( crop          (0:mxpft) )
end subroutine pftconrd
end module pftvarcon
