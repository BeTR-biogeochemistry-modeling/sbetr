module BgcCentCnpCompetType
!
! code to do ECA based competition
  ! !USES:
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use BiogeoConType , only : BiogeoCon_type
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: Compet_ECA_type
    real(r8), pointer :: vmax_minn_plant(:)     => null()   !number of maximum pft
    real(r8), pointer :: vmax_minp_plant(:)     => null()   !number of maximum pft
    real(r8), pointer :: vmax_minn_mic(:)      => null()
    real(r8), pointer :: vmax_minp_mic(:)      => null()
    real(r8), pointer :: vmax_den(:)           => null()
    real(r8), pointer :: vmax_nit(:) => null()
    real(r8), pointer :: kaff_minn_no3_plant(:) => null()   !number of maximum pft
    real(r8), pointer :: kaff_minn_nh4_plant(:)  => null()  !number of maximum pft
    real(r8), pointer :: kaff_minp_plant(:)    => null()    !number of maximum pft
    real(r8), pointer :: kaff_minn_nh4_mic(:) => null()     !number of maximum soil order
    real(r8), pointer :: kaff_minn_no3_mic(:)  => null()    !number of maximum soil order
    real(r8), pointer :: kaff_minp_mic(:)    => null()    !number of maximum pft
    real(r8), pointer :: kaff_minn_nh4_nit(:) => null()
    real(r8), pointer :: kaff_minn_no3_den(:) => null()
    real(r8), pointer :: kaff_minn_nh4_msurf(:) => null() !number of soil orders
    real(r8), pointer :: kaff_minp_msurf(:)   => null()   !number of soil orders
    real(r8), pointer :: vmax_minp_secondary_to_occlude(:) => null()
    real(r8), pointer :: vmax_minp_soluble_to_secondary(:) => null()

  contains
    procedure, public :: Init
    procedure, private:: InitAllocate
    procedure, private:: InitPar
    procedure, public :: run_compet_phosphorus
    procedure, public :: run_compet_nitrogen
  end type Compet_ECA_type

contains
  !-------------------------------------------------------------------------------
  subroutine Init(this, biogeo_con)

  implicit none
  class(Compet_ECA_type), intent(inout) :: this
  type(BiogeoCon_type),intent(in) :: biogeo_con

  call this%InitAllocate()

  call this%InitPar(biogeo_con)
  end subroutine Init
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(Compet_ECA_type), intent(inout) :: this

  allocate(this%vmax_minn_plant(betr_maxpatch_pft))
  allocate(this%vmax_minp_plant(betr_maxpatch_pft))
  allocate(this%kaff_minn_no3_plant(betr_maxpatch_pft))
  allocate(this%kaff_minn_nh4_plant(betr_maxpatch_pft))
  allocate(this%kaff_minp_plant(betr_maxpatch_pft))

  allocate(this%kaff_minn_nh4_msurf(betr_max_soilorder))
  allocate(this%kaff_minp_msurf(betr_max_soilorder))
  allocate(this%kaff_minn_nh4_mic(betr_max_soilorder))
  allocate(this%kaff_minn_no3_mic(betr_max_soilorder))
  allocate(this%kaff_minp_mic(betr_max_soilorder))
  allocate(this%kaff_minn_nh4_nit(betr_max_soilorder))
  allocate(this%kaff_minn_no3_den(betr_max_soilorder))
  allocate(this%vmax_minp_secondary_to_occlude(betr_max_soilorder))
  allocate(this%vmax_minp_soluble_to_secondary(betr_max_soilorder))
  allocate(this%vmax_nit(betr_max_soilorder))
  allocate(this%vmax_den(betr_max_soilorder))
  allocate(this%vmax_minp_mic(betr_max_soilorder))
  allocate(this%vmax_minn_mic(betr_max_soilorder))
  end subroutine InitAllocate
  !-------------------------------------------------------------------------------
  subroutine InitPar(this, biogeo_con)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(Compet_ECA_type) , intent(inout) :: this
  type(BiogeoCon_type)   , intent(in) :: biogeo_con

  integer :: j

  do j = 1, betr_maxpatch_pft
    this%vmax_minn_plant(j) = biogeo_con%vmax_minn_plant(j)
    this%vmax_minp_plant(j) = biogeo_con%vmax_minp_plant(j)
    this%kaff_minn_nh4_plant(j) = biogeo_con%kaff_minn_nh4_plant(j)
    this%kaff_minn_no3_plant(j) = biogeo_con%kaff_minn_no3_plant(j)
    this%kaff_minp_plant(j) = biogeo_con%kaff_minp_plant(j)
  enddo

  do j = 1, betr_max_soilorder
    this%kaff_minn_nh4_msurf(j) = biogeo_con%kaff_minn_nh4_msurf(j)
    this%kaff_minp_msurf(j) = biogeo_con%kaff_minp_msurf(j)
    this%kaff_minn_nh4_mic(j) = biogeo_con%kaff_minn_nh4_mic(j)
    this%kaff_minn_no3_mic(j) = biogeo_con%kaff_minn_no3_mic(j)
    this%kaff_minn_nh4_nit(j) = biogeo_con%kaff_minn_nh4_nit(j)
    this%kaff_minn_no3_den(j) = biogeo_con%kaff_minn_no3_den(j)
    this%kaff_minp_mic(j) = biogeo_con%kaff_minp_mic(j)
    this%vmax_nit(j) = biogeo_con%vmax_nit(j)
    this%vmax_den(j) = biogeo_con%vmax_den(j)
    this%vmax_minn_mic(j) = biogeo_con%vmax_minn_mic(j)
    this%vmax_minp_mic(j) = biogeo_con%vmax_minp_mic(j)
    this%vmax_minp_secondary_to_occlude(j) = biogeo_con%vmax_minp_secondary_to_occlude(j)
    this%vmax_minp_soluble_to_secondary(j) = biogeo_con%vmax_minp_soluble_to_secondary(j)
  enddo
  end subroutine InitPar
  !-------------------------------------------------------------------------------

  subroutine run_compet_nitrogen(this, smin_nh4, smin_no3, mic_pot_nn_flx, pot_f_nit, &
   pot_f_denit, plant_ntypes, plant_froot, plant_vtype, msurf_nh4, &
   soilorder, ECA_factor_nit, ECA_factor_den, ECA_factor_nh4_mic, &
    ECA_factor_no3_mic, ECA_flx_nh4_plants,ECA_flx_no3_plants, ECA_factor_msurf_nh4)

  use KineticsMod    , only : ecacomplex_cell_norm
  use BetrStatusType , only : betr_status_type
  implicit none
  class(Compet_ECA_type), intent(inout) :: this
  real(r8), intent(in) :: smin_nh4
  real(r8), intent(in) :: smin_no3
  real(r8), intent(in) :: mic_pot_nn_flx
  real(r8), intent(in) :: pot_f_nit
  real(r8), intent(in) :: pot_f_denit
  integer , intent(in) :: plant_ntypes
  real(r8), intent(in) :: plant_froot(plant_ntypes)
  integer,  intent(in) :: plant_vtype(plant_ntypes)
  real(r8), intent(in) :: msurf_nh4
  integer,  intent(in) :: soilorder
  real(r8), intent(out):: ECA_factor_nit
  real(r8), intent(out):: ECA_factor_den
  real(r8), intent(out):: ECA_factor_nh4_mic
  real(r8), intent(out):: ECA_factor_no3_mic
  real(r8), intent(out):: ECA_flx_nh4_plants(plant_ntypes)
  real(r8), intent(out):: ECA_flx_no3_plants(plant_ntypes)
  real(r8), intent(out):: ECA_factor_msurf_nh4
  !local variables
  real(r8), pointer :: kaff(:,:)
  real(r8), pointer :: substrate(:)
  real(r8), pointer :: entity(:)
  real(r8), pointer :: se_complex(:,:)
  real(r8) :: b_nit, b_den, b_mic
  integer :: tot_entity
  integer :: jj
  type(betr_status_type) :: bstatus

  !nit + denit + decomp_mic + msurf + plant
  tot_entity = 4 + plant_ntypes
  allocate(kaff(2,tot_entity))
  allocate(se_complex(2,tot_entity))
  allocate(substrate(2))
  allocate(entity(tot_entity))
  !form the affinity matrix
  !nh4
  kaff(1,:) = (/this%kaff_minn_nh4_nit(soilorder), 0._r8, &
      this%kaff_minn_nh4_mic(soilorder), this%kaff_minn_nh4_msurf(soilorder), &
      this%kaff_minn_nh4_plant(plant_vtype(1:plant_ntypes))/)
  !no3
  kaff(2,:) = (/0._r8, this%kaff_minn_no3_den(soilorder),&
      this%kaff_minn_no3_mic(soilorder), 0._r8, &
      this%kaff_minn_no3_plant(plant_vtype(1:plant_ntypes))/)
  !form the substrate vector
  substrate(:)=(/smin_nh4,smin_no3/)
  !form the competitor vector
  b_nit = pot_f_nit/this%vmax_nit(soilorder)
  b_den = pot_f_denit/this%vmax_den(soilorder)
  b_mic = mic_pot_nn_flx/this%vmax_minn_mic(soilorder)
  entity(:)=(/b_nit,b_den,b_mic,msurf_nh4,plant_froot(1:plant_ntypes)/)
  !do ECA calculation
  call ecacomplex_cell_norm(kaff,substrate,entity, se_complex, bstatus)
  ECA_factor_nit = se_complex(1,1)
  ECA_factor_den = se_complex(2,2)
  ECA_factor_nh4_mic = se_complex(1,3)
  ECA_factor_no3_mic = se_complex(2,3)
  ECA_factor_msurf_nh4 = se_complex(1,4)
  do jj = 1, plant_ntypes
    ECA_flx_nh4_plants(jj) = this%vmax_minn_plant(plant_vtype(jj)) * &
      plant_froot(plant_vtype(jj)) * se_complex(1,4+jj)

    ECA_flx_no3_plants(jj) = this%vmax_minn_plant(plant_vtype(jj)) * &
      plant_froot(plant_vtype(jj)) * se_complex(2,4+jj)
  enddo
  deallocate(kaff)
  deallocate(substrate)
  deallocate(entity)
  deallocate(se_complex)
  end subroutine run_compet_nitrogen
  !-------------------------------------------------------------------------------

  subroutine run_compet_phosphorus(this, sminp_soluble, mic_pot_np_flx, plant_ntypes, plant_froot, &
     plant_vtype, msurf_minp, soilorder, ECA_factor_phosphorus_mic, &
     ECA_flx_phosphorus_plants,ECA_flx_msurf_minp)
  use KineticsMod    , only : ecacomplex_cell_norm
  use BetrStatusType , only : betr_status_type
  implicit none
  class(Compet_ECA_type), intent(inout) :: this
  real(r8), intent(in) :: sminp_soluble
  real(r8), intent(in) :: mic_pot_np_flx
  integer , intent(in) :: plant_ntypes
  real(r8), intent(in) :: plant_froot(plant_ntypes)
  integer,  intent(in) :: plant_vtype(plant_ntypes)
  real(r8), intent(in) :: msurf_minp
  integer,  intent(in) :: soilorder
  real(r8), intent(out):: ECA_factor_phosphorus_mic
  real(r8), intent(out):: ECA_flx_phosphorus_plants(plant_ntypes)
  real(r8), intent(out):: ECA_flx_msurf_minp

  !local variables
  real(r8), pointer :: kaff(:,:)
  real(r8), pointer :: substrate(:)
  real(r8), pointer :: entity(:)
  real(r8), pointer :: se_complex(:,:)
  real(r8) :: b_mic
  integer :: tot_entity
  integer :: jj
  type(betr_status_type) :: bstatus

  !decomp_mic + msurf + plant
  tot_entity = 2 + plant_ntypes
  allocate(kaff(1,tot_entity))
  allocate(se_complex(1,tot_entity))
  allocate(substrate(1))
  allocate(entity(tot_entity))

  kaff(1,:) = (/this%kaff_minp_mic(soilorder), this%kaff_minp_msurf(soilorder), &
      this%kaff_minp_plant(plant_vtype(1:plant_ntypes))/)

  b_mic = mic_pot_np_flx/this%vmax_minp_mic(soilorder)
  entity(:)=(/b_mic,msurf_minp,plant_froot(1:plant_ntypes)/)
  substrate(:)=(/sminp_soluble/)

  !do ECA calculation
  call ecacomplex_cell_norm(kaff,substrate,entity, se_complex, bstatus)
  ECA_factor_phosphorus_mic = se_complex(1,1)
  do jj = 1, plant_ntypes
    ECA_flx_phosphorus_plants(jj) = this%vmax_minp_plant(plant_vtype(jj)) * &
      plant_froot(plant_vtype(jj)) * se_complex(1,2+jj)
  enddo
  ECA_flx_msurf_minp = this%vmax_minp_soluble_to_secondary(soilorder)*se_complex(1,2)*msurf_minp

  deallocate(kaff)
  deallocate(substrate)
  deallocate(entity)
  deallocate(se_complex)
  end subroutine run_compet_phosphorus

end module BgcCentCnpCompetType
