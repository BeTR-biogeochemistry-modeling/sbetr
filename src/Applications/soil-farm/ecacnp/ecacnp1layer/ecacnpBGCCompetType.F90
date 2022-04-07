module ecacnpBGCCompetType
!
! code to do ECA based competition
  ! !USES:
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use BetrStatusType      , only : betr_status_type
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: Compet_ECACNP_type
    real(r8), pointer :: mumax_minn_nh4_plant(:)     => null()   !number of maximum pft
    real(r8), pointer :: mumax_minn_no3_plant(:)     => null()   !number of maximum pft
    real(r8), pointer :: mumax_minp_plant(:)     => null()   !number of maximum pft
    real(r8), pointer :: kaff_minn_no3_plant(:) => null()   !number of maximum pft
    real(r8), pointer :: kaff_minn_nh4_plant(:)  => null()  !number of maximum pft
    real(r8), pointer :: kaff_minp_plant(:)    => null()    !number of maximum pft
    real(r8), pointer :: plant_froot_nn(:) => null()
    real(r8), pointer :: plant_froot_np(:) => null()
    real(r8), pointer :: plant_eff_frootc_patch(:) => null()
    real(r8) :: compet_bn_mic   !decomposer n competition transporter
    real(r8) :: compet_bp_mic   !decomposer p competition transporter
    real(r8) :: compet_bn_den
    real(r8) :: compet_bn_nit
    real(r8) :: kaff_minn_nh4_mic
    real(r8) :: kaff_minn_no3_mic
    real(r8) :: kaff_minp_mic
    real(r8) :: kaff_minn_nh4_nit
    real(r8) :: kaff_minn_no3_den
    real(r8) :: kaff_minn_nh4_msurf
    real(r8) :: kaff_minp_msurf
    real(r8) :: topt
    logical  :: debug
  contains
    procedure, public :: Init
    procedure, private:: InitAllocate
    procedure, public :: run_compet_phosphorus
    procedure, public :: run_compet_nitrogen
    procedure, public :: UpdateParas
  end type Compet_ECACNP_type

contains
  !-------------------------------------------------------------------------------
  subroutine Init(this, biogeo_con, bstatus)
  use BiogeoConType             , only : BiogeoCon_type
  use ecacnpParaType            , only : ecacnp_para_type
  implicit none
  class(Compet_ECACNP_type), intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)     , intent(out)   :: bstatus

  character(len=256) :: msg

  call this%InitAllocate()

  call bstatus%reset()
  select type(biogeo_con)
  type is(ecacnp_para_type)
    this%kaff_minn_nh4_mic = biogeo_con%km_decomp_nh4
    this%kaff_minn_no3_mic = biogeo_con%km_decomp_no3
    this%kaff_minp_mic     = biogeo_con%km_decomp_p
    this%kaff_minn_nh4_nit = biogeo_con%km_nit
    this%kaff_minn_no3_den = biogeo_con%km_den
  class default
    write(msg,'(A)')'Wrong parameter type passed in for Init in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine Init

!------------------------------------------
  subroutine UpdateParas(this,  biogeo_con)
  use BiogeoConType             , only : BiogeoCon_type
  use ecacnpParaType            , only : ecacnp_para_type
  implicit none
  class(Compet_ECACNP_type), intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con

  select type(biogeo_con)
  type is(ecacnp_para_type)
    this%kaff_minn_nh4_mic = biogeo_con%km_decomp_nh4
    this%kaff_minn_no3_mic = biogeo_con%km_decomp_no3
    this%kaff_minp_mic     = biogeo_con%km_decomp_p
    this%kaff_minn_nh4_nit = biogeo_con%km_nit
    this%kaff_minn_no3_den = biogeo_con%km_den
  end select
  end subroutine UpdateParas

  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(Compet_ECACNP_type), intent(inout) :: this

  allocate(this%mumax_minn_nh4_plant(betr_maxpatch_pft))
  allocate(this%mumax_minn_no3_plant(betr_maxpatch_pft))
  allocate(this%mumax_minp_plant(betr_maxpatch_pft))
  allocate(this%kaff_minn_no3_plant(betr_maxpatch_pft))
  allocate(this%kaff_minn_nh4_plant(betr_maxpatch_pft))
  allocate(this%kaff_minp_plant(betr_maxpatch_pft))
  allocate(this%plant_froot_nn(betr_maxpatch_pft));
  allocate(this%plant_froot_np(betr_maxpatch_pft));
  allocate(this%plant_eff_frootc_patch(betr_maxpatch_pft))
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------

  subroutine run_compet_nitrogen(this, non_limit, smin_nh4, smin_no3,  &
   plant_ntypes, msurf_nh4, ECA_factor_nit, ECA_factor_den, ECA_factor_nh4_mic, &
    ECA_factor_no3_mic, ECA_flx_nh4_plants,ECA_flx_no3_plants, ECA_factor_msurf_nh4)

  use KineticsMod    , only : ecacomplex_cell_norm
  use BetrStatusType , only : betr_status_type
  implicit none
  class(Compet_ECACNP_type), intent(inout) :: this
  logical , intent(in) :: non_limit
  real(r8), intent(in) :: smin_nh4
  real(r8), intent(in) :: smin_no3
  integer , intent(in) :: plant_ntypes
  real(r8), intent(in) :: msurf_nh4
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
  if(plant_ntypes>0)then
    kaff(1,:) = (/this%kaff_minn_nh4_nit, 0._r8, &
      this%kaff_minn_nh4_mic, this%kaff_minn_nh4_msurf, &
      this%kaff_minn_nh4_plant(1:plant_ntypes)/)
  !no3
    kaff(2,:) = (/0._r8, this%kaff_minn_no3_den,&
      this%kaff_minn_no3_mic, 0._r8, &
      this%kaff_minn_no3_plant(1:plant_ntypes)/)
  else
    kaff(1,:) = (/this%kaff_minn_nh4_nit, 0._r8, &
      this%kaff_minn_nh4_mic, this%kaff_minn_nh4_msurf/)

    kaff(2,:) = (/0._r8, this%kaff_minn_no3_den,&
      this%kaff_minn_no3_mic, 0._r8/)

  endif
  !form the substrate vector
  substrate(:)=(/smin_nh4,smin_no3/)

  !form the competitor vector
  b_nit = this%compet_bn_nit
  b_den = this%compet_bn_den
  b_mic = this%compet_bn_mic
  if(plant_ntypes>0)then
    entity(:)=(/b_nit,b_den,b_mic,msurf_nh4,this%plant_froot_nn(1:plant_ntypes)/)
  else
    entity(:)=(/b_nit,b_den,b_mic,msurf_nh4/)
  endif
  !do ECA calculation
  call ecacomplex_cell_norm(kaff,substrate,entity, se_complex, bstatus)
  if(non_limit)then
    ECA_factor_nit = 1._r8
    ECA_factor_den = 1._r8
    ECA_factor_nh4_mic = 1._r8
    ECA_factor_no3_mic = 1._r8
    do jj = 1, plant_ntypes
      ECA_flx_nh4_plants(jj) = this%mumax_minn_nh4_plant(jj) * &
        this%plant_froot_nn(jj)

      ECA_flx_no3_plants(jj) = this%mumax_minn_no3_plant(jj) * &
        this%plant_froot_nn(jj)
    enddo
  else
    ECA_factor_nit = se_complex(1,1)
    ECA_factor_den = se_complex(2,2)
    ECA_factor_nh4_mic = se_complex(1,3)
    ECA_factor_no3_mic = se_complex(2,3)
    do jj = 1, plant_ntypes
      ECA_flx_nh4_plants(jj) = this%mumax_minn_nh4_plant(jj) * &
        this%plant_froot_nn(jj) * se_complex(1,4+jj)

      ECA_flx_no3_plants(jj) = this%mumax_minn_no3_plant(jj) * &
        this%plant_froot_nn(jj) * se_complex(2,4+jj)
    enddo
  endif
  ECA_factor_msurf_nh4 = se_complex(1,4)

  deallocate(kaff)
  deallocate(substrate)
  deallocate(entity)
  deallocate(se_complex)
  end subroutine run_compet_nitrogen
  !-------------------------------------------------------------------------------

  subroutine run_compet_phosphorus(this, nop_lim,  sminp_soluble, plant_ntypes,&
     msurf_minp, ECA_factor_phosphorus_mic, ECA_factor_minp_msurf, ECA_flx_phosphorus_plants)

  !
  !DESCRIPTION
  ! do eca competition of phosphorus
  use KineticsMod    , only : ecacomplex_cell_norm
  use BetrStatusType , only : betr_status_type
  implicit none
  class(Compet_ECACNP_type), intent(inout) :: this
  real(r8), intent(in) :: sminp_soluble
  logical , intent(in) :: nop_lim               !logical indicator of P limitation
  integer , intent(in) :: plant_ntypes
  real(r8), intent(in) :: msurf_minp
  real(r8), intent(out):: ECA_factor_phosphorus_mic
  real(r8), intent(out):: ECA_flx_phosphorus_plants(plant_ntypes)
  real(r8), intent(out):: ECA_factor_minp_msurf

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

  if(plant_ntypes>0)then
    kaff(1,:) = (/this%kaff_minp_mic, this%kaff_minp_msurf, &
      this%kaff_minp_plant(1:plant_ntypes)/)
  else
    kaff(1,:) = (/this%kaff_minp_mic, this%kaff_minp_msurf/)
  endif
  b_mic = this%compet_bp_mic
  if(plant_ntypes>0)then
    entity(:)=(/b_mic,msurf_minp,this%plant_froot_np(1:plant_ntypes)/)
  else
    entity(:)=(/b_mic,msurf_minp/)
  endif
  substrate(:)=(/sminp_soluble/)

  !given P is under competitation by microbes, plants and mineral surfaces
  !the definitation of P limitation is quite ambiguous. Following what
  !is done in ELM-ECA-CNP and ELM-CNP, P limitation is only applied
  !to biological uptake. P-limitation to mineral sorption is not
  !considered.

  !do ECA calculation
  call ecacomplex_cell_norm(kaff,substrate,entity, se_complex, bstatus)
  ECA_factor_minp_msurf = se_complex(1,2)

  if(nop_lim)then
    !no P limitation is imposed on biological reactions
    ECA_factor_phosphorus_mic = 1._r8
    ECA_factor_minp_msurf = 0._r8
    do jj = 1, plant_ntypes
      ECA_flx_phosphorus_plants(jj) = this%mumax_minp_plant(jj) * &
        this%plant_froot_np(jj)
    enddo
  else
    ECA_factor_phosphorus_mic = se_complex(1,1)
    do jj = 1, plant_ntypes
      ECA_flx_phosphorus_plants(jj) = this%mumax_minp_plant(jj) * &
        this%plant_froot_np(jj) * se_complex(1,2+jj)
    enddo
  endif

  deallocate(kaff)
  deallocate(substrate)
  deallocate(entity)
  deallocate(se_complex)
  end subroutine run_compet_phosphorus

end module ecacnpBGCCompetType
