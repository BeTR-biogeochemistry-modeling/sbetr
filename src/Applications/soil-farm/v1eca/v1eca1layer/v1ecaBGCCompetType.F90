module v1ecaBGCCompetType
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

  type, public :: Compet_V1ECA_type
    real(r8), pointer :: vmax_minn_nh4_plant(:)     => null()   !number of maximum pft, already scaled with transporter and other factors
    real(r8), pointer :: vmax_minn_no3_plant(:)     => null()   !number of maximum pft, already scaled with transporter and other factors
    real(r8), pointer :: vmax_minp_plant(:)     => null()   !number of maximum pft
    real(r8), pointer :: kaff_minn_no3_plant(:) => null()   !number of maximum pft
    real(r8), pointer :: kaff_minn_nh4_plant(:)  => null()  !number of maximum pft
    real(r8), pointer :: kaff_minp_plant(:)    => null()    !number of maximum pft
    real(r8), pointer :: plant_froot_nn(:) => null()
    real(r8), pointer :: plant_froot_np(:) => null()
    real(r8), pointer :: plant_eff_frootc_patch(:) => null()
    real(r8) :: vmax_minsurf_p
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
    real(r8) :: dsolutionp_dt
    real(r8) :: dlabp_dt
    real(r8) :: bd
    real(r8) :: h2osoi_vol
    logical  :: debug
  contains
    procedure, public :: Init
    procedure, private:: InitAllocate
    procedure, public :: run_compet_phosphorus
    procedure, public :: run_compet_nitrogen
  end type Compet_V1ECA_type

contains
  !-------------------------------------------------------------------------------
  subroutine Init(this, biogeo_con, bstatus)
  use BiogeoConType             , only : BiogeoCon_type
  use v1ecaParaType            , only : v1eca_para_type
  implicit none
  class(Compet_V1ECA_type), intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)     , intent(out)   :: bstatus

  character(len=256) :: msg

  call this%InitAllocate()

  call bstatus%reset()
  select type(biogeo_con)
  type is(v1eca_para_type)
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
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(Compet_V1ECA_type), intent(inout) :: this

  allocate(this%vmax_minn_nh4_plant(betr_maxpatch_pft))
  allocate(this%vmax_minn_no3_plant(betr_maxpatch_pft))
  allocate(this%vmax_minp_plant(betr_maxpatch_pft))
  allocate(this%kaff_minn_no3_plant(betr_maxpatch_pft))
  allocate(this%kaff_minn_nh4_plant(betr_maxpatch_pft))
  allocate(this%kaff_minp_plant(betr_maxpatch_pft))
  allocate(this%plant_froot_nn(betr_maxpatch_pft));
  allocate(this%plant_froot_np(betr_maxpatch_pft));
  allocate(this%plant_eff_frootc_patch(betr_maxpatch_pft))
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------

  subroutine run_compet_nitrogen(this, non_limit, sol_smin_nh4, sol_smin_no3,  &
   plant_ntypes, ECA_factor_nit, ECA_factor_den, ECA_factor_nh4_mic, &
    ECA_factor_no3_mic, ECA_flx_nh4_plants,ECA_flx_no3_plants)

  use KineticsMod    , only : ecacomplex_cell_norm
  use BetrStatusType , only : betr_status_type
  implicit none
  class(Compet_V1ECA_type), intent(inout) :: this
  logical , intent(in) :: non_limit
  real(r8), intent(in) :: sol_smin_nh4
  real(r8), intent(in) :: sol_smin_no3
  integer , intent(in) :: plant_ntypes
  real(r8), intent(out):: ECA_factor_nit
  real(r8), intent(out):: ECA_factor_den
  real(r8), intent(out):: ECA_factor_nh4_mic
  real(r8), intent(out):: ECA_factor_no3_mic
  real(r8), intent(out):: ECA_flx_nh4_plants(plant_ntypes)
  real(r8), intent(out):: ECA_flx_no3_plants(plant_ntypes)
  !local variables
  real(r8) :: e_km_nh4, e_km_no3, e_km_p
  integer :: tot_entity
  integer :: jj
  type(betr_status_type) :: bstatus

  !nit + denit + decomp_mic + msurf + plant

  if(non_limit)then
    ECA_factor_nit = 1._r8
    ECA_factor_den = 1._r8
    ECA_factor_nh4_mic = 1._r8
    ECA_factor_no3_mic = 1._r8
    ECA_flx_nh4_plants(1:plant_ntypes)=1._r8
    ECA_flx_no3_plants(1:plant_ntypes)=1._r8
  else
    e_km_nh4 = 0._r8;e_km_no3=0._r8
    do jj = 1, plant_ntypes
      e_km_nh4=e_km_nh4+this%plant_froot_nn(jj)/this%kaff_minn_nh4_plant(jj)
      e_km_no3=e_km_no3+this%plant_froot_nn(jj)/this%kaff_minn_no3_plant(jj)
    enddo
    e_km_nh4=e_km_nh4+this%compet_bn_mic/this%kaff_minn_nh4_mic+this%compet_bn_nit/this%kaff_minn_nh4_nit
    e_km_no3=e_km_no3+this%compet_bn_mic/this%kaff_minn_no3_mic+this%compet_bn_den/this%kaff_minn_no3_den

    do jj = 1, plant_ntypes
      ECA_flx_nh4_plants(jj) = sol_smin_nh4/(this%kaff_minn_nh4_plant(jj)*(1._r8 + &
          sol_smin_nh4/this%kaff_minn_nh4_plant(jj)+e_km_nh4))
      ECA_flx_no3_plants(jj) = sol_smin_no3/(this%kaff_minn_no3_plant(jj)*(1._r8+ &
          sol_smin_no3/this%kaff_minn_no3_plant(jj)+e_km_no3))
    enddo
    ECA_factor_nh4_mic=sol_smin_nh4/(this%kaff_minn_nh4_mic*(1._r8+&
        sol_smin_nh4/this%kaff_minn_nh4_mic+e_km_nh4))
    ECA_factor_no3_mic=sol_smin_no3/(this%kaff_minn_no3_mic*(1._r8+&
        sol_smin_no3/this%kaff_minn_no3_mic+e_km_no3))
    ECA_factor_nit=sol_smin_nh4/(this%kaff_minn_nh4_nit*(1._r8+&
        sol_smin_nh4/this%kaff_minn_nh4_nit+e_km_nh4))
    ECA_factor_den=sol_smin_no3/(this%kaff_minn_no3_den*(1._r8+&
        sol_smin_no3/this%kaff_minn_no3_den+e_km_no3))
  endif
  do jj = 1, plant_ntypes
      ECA_flx_nh4_plants(jj) = this%vmax_minn_nh4_plant(jj) * &
        ECA_flx_nh4_plants(jj)

      ECA_flx_no3_plants(jj) = this%vmax_minn_no3_plant(jj) * &
        ECA_flx_no3_plants(jj)
  enddo

  end subroutine run_compet_nitrogen
  !-------------------------------------------------------------------------------

  subroutine run_compet_phosphorus(this, nop_lim, sol_sminp_soluble, plant_ntypes,&
     ECA_factor_phosphorus_mic, ECA_factor_minp_msurf, ECA_flx_phosphorus_plants)

  !
  !DESCRIPTION
  ! do eca competition of phosphorus
  use KineticsMod    , only : ecacomplex_cell_norm
  use BetrStatusType , only : betr_status_type
  implicit none
  class(Compet_V1ECA_type), intent(inout) :: this
  real(r8), intent(in) :: sol_sminp_soluble
  logical , intent(in) :: nop_lim               !logical indicator of P limitation
  integer , intent(in) :: plant_ntypes
  real(r8), intent(out):: ECA_factor_phosphorus_mic
  real(r8), intent(out):: ECA_flx_phosphorus_plants(plant_ntypes)
  real(r8), intent(out):: ECA_factor_minp_msurf

  !local variables

  integer :: jj
  real(r8):: e_km_p
  type(betr_status_type) :: bstatus

  !decomp_mic + msurf + plant

  !do ECA calculation

  if(nop_lim)then
    !no P limitation is imposed on biological reactions
    ECA_factor_phosphorus_mic = 1._r8
    ECA_factor_minp_msurf = 0._r8
    ECA_flx_phosphorus_plants(1:plant_ntypes) = 1._r8
  else
    e_km_p = 0._r8
    do jj = 1, plant_ntypes
      e_km_p = e_km_p + this%plant_froot_np(jj)/this%kaff_minp_plant(jj)
    enddo
    e_km_p=e_km_p + this%compet_bp_mic/this%kaff_minp_mic + &
       this%kaff_minp_msurf/this%kaff_minp_msurf

    do jj = 1, plant_ntypes
      ECA_flx_phosphorus_plants(jj) = sol_sminp_soluble/(this%kaff_minp_plant(jj) * &
        (1._r8+sol_sminp_soluble/this%kaff_minp_plant(jj)+e_km_p))
    enddo
    ECA_factor_phosphorus_mic = sol_sminp_soluble/(this%kaff_minp_mic * &
      (1._r8 + sol_sminp_soluble/this%kaff_minp_mic + e_km_p))
    ECA_factor_minp_msurf = sol_sminp_soluble/(this%kaff_minp_msurf* &
      (1._r8 + sol_sminp_soluble/this%kaff_minp_msurf)+e_km_p)
  endif

  do jj = 1, plant_ntypes
    ECA_flx_phosphorus_plants(jj) = this%vmax_minp_plant(jj) * ECA_flx_phosphorus_plants(jj)
  enddo

  end subroutine run_compet_phosphorus

end module v1ecaBGCCompetType
