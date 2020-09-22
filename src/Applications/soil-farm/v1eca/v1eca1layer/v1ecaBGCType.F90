module v1ecaBGCType
#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  ! subroutines for stoichiometric configuration of the century bgc
  ! !History, created by Jinyun Tang, Dec, 2014.
  ! the P adsorption and desorption are assumed to occur between labile and soluble
  ! phase. The model assumes these processes follow the Langmuir isotherm, and
  ! are not contributing to the mass change through the process of competition.
  ! Rather, these effects are accounted separately in the module PhosphorusDynamicsMod
  ! This single layer model should return total plant p uptake sminp_to_plant_vr
  ! and total p uptake to microbes, where the latter could be either positive or negative
  ! !USES:
  use bshr_assert_mod, only : shr_assert
  use bshr_assert_mod, only : shr_assert_all, shr_assert_all_ext
  use bshr_assert_mod, only : shr_assert_any
  use bshr_kind_mod             , only : r8 => shr_kind_r8
  use bshr_log_mod              , only : errMsg => shr_log_errMsg
  use betr_varcon               , only : spval => bspval
  use betr_ctrl                 , only : spinup_state => betr_spinup_state
  use gbetrType                 , only : gbetr_type
  use BiogeoConType             , only : BiogeoCon_type
  use BetrStatusType            , only : betr_status_type
  use BeTRJarModel              , only : jar_model_type
  use v1ecaBGCDecompType      , only : DecompCent_type
  use v1ecaBGCIndexType       , only : v1eca_bgc_index_type
  use v1ecaBGCNitDenType      , only : century_nitden_type
  use v1ecaBGCSOMType         , only : CentSom_type
  use v1ecaBGCCompetType      , only : Compet_V1ECA_type
  use v1ecaParaType           , only : v1eca_para_type
  use tracer_varcon           , only : natomw
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !Note:
  !Keeping v1eca_bgc_index as a private member is a workaround to call the ode solver
  !it increase the memory for each instance of v1eca_bgc_type, but enables
  !the ode function to be called by the ode solver

  type, extends(jar_model_type), public :: v1eca_bgc_type
    type(DecompCent_type),private        :: decompkf_eca
    type(century_nitden_type), private   :: nitden
    type(CentSom_type), private          :: censom
    type(Compet_V1ECA_type), public        :: competECA
    type(v1eca_bgc_index_type), private :: v1eca_bgc_index
    real(r8), pointer                    :: ystates0(:)
    real(r8), pointer                    :: ystates1(:)
    real(r8), pointer                    :: k_decay(:)
    real(r8), pointer                    :: cascade_matrix(:,:)
    real(r8), pointer                    :: cascade_matrixd(:,:)
    real(r8), pointer                    :: cascade_matrixp(:,:)
    real(r8), pointer                    :: cascade_matnh4(:)
    real(r8), pointer                    :: alpha_n(:)
    real(r8), pointer                    :: alpha_p(:)
    real(r8)                             :: pot_f_nit
    real(r8)                             :: pot_f_denit
    real(r8)                             :: rt_ar
    real(r8), pointer                    :: frac_p_sec_to_sol(:)
    real(r8), pointer                    :: minp_secondary_decay(:)
    real(r8), pointer                    :: mumax_minp_soluble_to_secondary(:)
    integer                              :: plant_ntypes
    logical                              :: nop_limit
    logical                              :: non_limit
    real(r8), pointer                    :: scal_f(:)
    real(r8), pointer                    :: conv_f(:)
    real(r8), pointer                    :: conc_f(:)
    integer                              :: soilorder
    real(r8), private                    :: c14decay_const
    real(r8), private                    :: c14decay_som_const
    real(r8), private                    :: c14decay_pom_const
    real(r8), private                    :: c14decay_Bm_const
    logical , private                    :: use_c13
    logical , private                    :: use_c14
    real(r8), private                    :: beg_c_mass, beg_c13_mass, beg_c14_mass
    real(r8), private                    :: beg_n_mass
    real(r8), private                    :: beg_p_mass
    real(r8), private                    :: c_inflx,n_inflx, p_inflx
    logical                              :: bgc_on
    logical , private                    :: batch_mode
  contains
    procedure, public  :: init          => init_v1eca
    procedure, public  :: initVL         => initVL_v1eca
    procedure, public  :: runbgc        => runbgc_v1eca
    procedure, public  :: UpdateParas   => UpdateParas_v1eca
    procedure, public  :: sumup_cnp_msflx => sumup_cnp_msflx_v1eca
    procedure, public  :: getvarllen    => getvarllen_v1eca
    procedure, public  :: getvarlist    => getvarlist_v1eca
    procedure, public  :: init_cold     => init_cold_v1eca
    procedure, private :: calc_cascade_matrix
    procedure, private :: init_states
    procedure, private :: InitAllocate
    procedure, private :: arenchyma_gas_transport
    procedure, private :: sumup_cnp_mass
    procedure, private :: bgc_integrate
    procedure, private :: c14decay
    procedure, private :: checksum_cascade
    procedure, private :: begin_massbal_check
    procedure, private :: end_massbal_check
    procedure, private :: sum_tot_store
    procedure, public  :: display_index
  end type v1eca_bgc_type

  public :: create_jarmodel_v1eca
contains

  function create_jarmodel_v1eca()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(v1eca_bgc_type), pointer :: create_jarmodel_v1eca
    class(v1eca_bgc_type), pointer :: bgc

    allocate(bgc)
    create_jarmodel_v1eca => bgc

  end function create_jarmodel_v1eca

  !-------------------------------------------------------------------------------
  function getvarllen_v1eca(this)result(ans)

  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  integer :: ans

  ans = this%v1eca_bgc_index%nstvars

  end function getvarllen_v1eca
  !-------------------------------------------------------------------------------
  subroutine getvarlist_v1eca(this, nstvars, varnames, varunits, vartypes)
  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  integer, intent(in) :: nstvars
  character(len=*), intent(out) :: varnames(1:nstvars)
  character(len=*), intent(out) :: varunits(1:nstvars)
  integer         , intent(out) :: vartypes(1:nstvars)
  integer :: n

  do n = 1, nstvars
    vartypes(n) = this%v1eca_bgc_index%vartypes(n)
    write(varnames(n),'(A)')trim(this%v1eca_bgc_index%varnames(n))
    write(varunits(n),'(A)')trim(this%v1eca_bgc_index%varunits(n))
  enddo
  end subroutine getvarlist_v1eca

  !-------------------------------------------------------------------------------
  subroutine begin_massbal_check(this)

  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  real(r8) :: c_mass0, n_mass0, p_mass0

  print*,'begin_massbal_check'
  call this%sumup_cnp_msflx(this%ystates0, c_mass0,  n_mass0, p_mass0)
  write(*, '(A,3(X,E20.10))')'cnp_mass beg0=',c_mass0,n_mass0,p_mass0
  print*,'co20',this%ystates0(this%v1eca_bgc_index%lid_co2_hr)
  call this%sumup_cnp_msflx(this%ystates1, this%beg_c_mass, &
    this%beg_n_mass, this%beg_p_mass)
  write(*, '(A,3(X,E20.10))')'cnp_mass beg1=',this%beg_c_mass,this%beg_n_mass,this%beg_p_mass
  print*,'co21',this%ystates1(this%v1eca_bgc_index%lid_co2_hr)
  end subroutine begin_massbal_check
  !-------------------------------------------------------------------------------

  subroutine end_massbal_check(this, header)

  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  character(len=*), intent(in) :: header
  real(r8) :: c_mass, n_mass, p_mass
  real(r8) :: c_flx,n_flx,p_flx
  type(betr_status_type) :: bstatus
  real(r8) :: dmass_c, dmass_n, dmass_p
  real(r8) :: emass_c, emass_n, emass_p
  real(r8) :: remass_c, remass_n, remass_p

  real(r8), parameter :: tiny_val=1.e-10_r8
  print*,trim(header)
  call this%sumup_cnp_msflx(this%ystates1, c_mass, n_mass, p_mass,c_flx,n_flx,p_flx)

  dmass_c=c_mass - this%beg_c_mass
  dmass_n=n_mass - this%beg_n_mass
  dmass_p=p_mass - this%beg_p_mass
  print*,'dmassc',dmass_c
  emass_c=dmass_c + c_flx - this%c_inflx
  emass_n=dmass_n + n_flx - this%n_inflx
  emass_p=dmass_p + p_flx - this%p_inflx

  remass_c=emass_c/max(abs(dmass_c),tiny_val)
  remass_n=emass_n/max(abs(dmass_n),tiny_val)
  remass_p=emass_p/max(abs(dmass_p),tiny_val)
  print*,'co2_hr',this%ystates1(this%v1eca_bgc_index%lid_co2_hr)
  write(*,'(A)')'------------------------------------------------------------------------------------------------------'
  write(*,'(A)')'type             beg_mass             end_mass            emass                inflx            outflx'
  write(*, '(A,5(X,E20.10))')'c_mass bal=',this%beg_c_mass, c_mass, emass_c, this%c_inflx, c_flx
  write(*, '(A,5(X,E20.10))')'n_mass bal=',this%beg_n_mass, n_mass, emass_n, this%n_inflx, n_flx
  write(*, '(A,5(X,E20.10))')'p_mass bal=',this%beg_p_mass, p_mass, emass_p, this%p_inflx, p_flx
  if(maxval((/abs(remass_c),abs(remass_n),abs(remass_p)/))>1.e-3_r8)stop
  end subroutine end_massbal_check

  !-------------------------------------------------------------------------------
  subroutine UpdateParas_v1eca(this,  biogeo_con, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  type(betr_status_type)     , intent(out)   :: bstatus
  integer :: sr
  character(len=256) :: msg

  call bstatus%reset()

  select type(biogeo_con)
  type is(v1eca_para_type)
    do sr = 1, betr_max_soilorder

      this%frac_p_sec_to_sol(sr) = biogeo_con%frac_p_sec_to_sol(sr)

      this%minp_secondary_decay(sr) = biogeo_con%minp_secondary_decay(sr)

      this%mumax_minp_soluble_to_secondary(sr) = biogeo_con%vmax_minp_soluble_to_secondary(sr)
    enddo

    call this%nitden%UpdateParas(biogeo_con)

    if(this%use_c14)then
      this%c14decay_const=biogeo_con%c14decay_const
      this%c14decay_som_const=biogeo_con%c14decay_som_const
      this%c14decay_pom_const=biogeo_con%c14decay_pom_const
      this%c14decay_Bm_const=biogeo_con%c14decay_Bm_const
    endif

    call this%censom%UpdateParas(this%v1eca_bgc_index, biogeo_con)

    call this%decompkf_eca%UpdateParas(biogeo_con)

  class default
    write(msg,'(A)')'Wrong parameter type passed in for UpdateParas in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine UpdateParas_v1eca
  !-------------------------------------------------------------------------------
  subroutine init_v1eca(this,  biogeo_con,  batch_mode, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  logical                   , intent(in) :: batch_mode
  type(betr_status_type)    , intent(out) :: bstatus

  character(len=256) :: msg
  write(this%jarname, '(A)')'v1eca'

  this%bgc_on=.true.
  this%batch_mode=batch_mode
  select type(biogeo_con)
  type is(v1eca_para_type)
    call bstatus%reset()
    call this%v1eca_bgc_index%Init(biogeo_con%use_c13, biogeo_con%use_c14, &
     biogeo_con%non_limit, biogeo_con%nop_limit, betr_maxpatch_pft, this%batch_mode)

    this%nop_limit=biogeo_con%nop_limit
    this%non_limit=biogeo_con%non_limit

    call this%InitAllocate(this%v1eca_bgc_index)

    call this%censom%Init(this%v1eca_bgc_index, biogeo_con, bstatus)

    if(bstatus%check_status())return

    call this%decompkf_eca%Init(biogeo_con)

    call this%nitden%Init(biogeo_con)

    call this%competECA%Init(biogeo_con, bstatus)
    if(bstatus%check_status())return
    this%use_c13 = biogeo_con%use_c13

    this%use_c14 = biogeo_con%use_c14

    call this%UpdateParas(biogeo_con, bstatus)
    if(bstatus%check_status())return
  class default
    call bstatus%reset()
    write(msg,'(A)')'Wrong parameter type passed in for init_v1eca in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine init_v1eca
  !-------------------------------------------------------------------------------
  subroutine initVL_v1eca(this,  biogeo_con,  batch_mode, v1eca_bgc_index, bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in) :: biogeo_con
  logical                   , intent(in) :: batch_mode
  type(v1eca_bgc_index_type), intent(in) :: v1eca_bgc_index
  type(betr_status_type)    , intent(out) :: bstatus

  character(len=256) :: msg
  write(this%jarname, '(A)')'v1eca'

  this%bgc_on=.true.
  this%batch_mode=batch_mode
  select type(biogeo_con)
  type is(v1eca_para_type)
    call bstatus%reset()
    call this%v1eca_bgc_index%hcopy(v1eca_bgc_index)

    this%nop_limit=biogeo_con%nop_limit
    this%non_limit=biogeo_con%non_limit

    call this%InitAllocate(this%v1eca_bgc_index)

    call this%censom%Init(this%v1eca_bgc_index, biogeo_con, bstatus)

    if(bstatus%check_status())return

    call this%decompkf_eca%Init(biogeo_con)

    call this%nitden%Init(biogeo_con)

    call this%competECA%Init(biogeo_con, bstatus)
    if(bstatus%check_status())return
    this%use_c13 = biogeo_con%use_c13

    this%use_c14 = biogeo_con%use_c14

    call this%UpdateParas(biogeo_con, bstatus)
    if(bstatus%check_status())return
  class default
    call bstatus%reset()
    write(msg,'(A)')'Wrong parameter type passed in for init_v1eca in ' &
      // errMsg(mod_filename,__LINE__)
    call bstatus%set_msg(msg,err=-1)
    return
  end select
  end subroutine initVL_v1eca
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, v1eca_bgc_index)
  use v1ecaBGCIndexType , only : v1eca_bgc_index_type
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(v1eca_bgc_type)   , intent(inout) :: this
  type(v1eca_bgc_index_type) , intent(in):: v1eca_bgc_index

  associate(                                   &
    nom_pools=> v1eca_bgc_index%nom_pools,    &
    nstvars => v1eca_bgc_index%nstvars  ,     &
    nreactions => v1eca_bgc_index%nreactions, &
    nprimvars => v1eca_bgc_index%nprimvars    &
  )

  allocate(this%ystates0(nstvars)); this%ystates0(:) = 0._r8
  allocate(this%ystates1(nstvars)); this%ystates1(:) = 0._r8
  allocate(this%k_decay(nom_pools)); this%k_decay(:) = 0._r8
  allocate(this%scal_f(nprimvars));  this%scal_f(:) = 0._r8
  allocate(this%conv_f(nprimvars));  this%conv_f(:) = 0._r8
  allocate(this%conc_f(nprimvars));  this%conc_f(:) = 0._r8

  allocate(this%cascade_matrix (1:nstvars  , 1:nreactions)); this%cascade_matrix(:,:) = 0._r8
  allocate(this%cascade_matrixd(1:nprimvars, 1:nreactions)); this%cascade_matrixd(:,:) = 0._r8
  allocate(this%cascade_matrixp(1:nprimvars, 1:nreactions)); this%cascade_matrixp(:,:) = 0._r8
  allocate(this%cascade_matnh4(nom_pools)); this%cascade_matnh4(:)=0._r8

  allocate(this%alpha_n(nom_pools)); this%alpha_n(:) = 0._r8
  allocate(this%alpha_p(nom_pools)); this%alpha_p(:) = 0._r8
  allocate(this%frac_p_sec_to_sol(betr_max_soilorder)); this%frac_p_sec_to_sol(:) = 0._r8
  allocate(this%minp_secondary_decay(betr_max_soilorder)); this%minp_secondary_decay(:) = 0._r8
  allocate(this%mumax_minp_soluble_to_secondary(betr_max_soilorder)); this%mumax_minp_soluble_to_secondary(:) = 0._r8
  end associate
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------


  subroutine checksum_cascade(this, v1eca_bgc_index)

  use v1ecaBGCIndexType       , only : v1eca_bgc_index_type

  implicit none
  ! !ARGUMENTS:
  class(v1eca_bgc_type)     , intent(in) :: this
  type(v1eca_bgc_index_type)   , intent(in) :: v1eca_bgc_index

  real(r8) :: resc, resn, resp
  integer  :: reac,jj

  associate(                                                   &
    cascade_matrix => this%cascade_matrix                    , &
    lit1      => v1eca_bgc_index%lit1                       , & !
    lit2      => v1eca_bgc_index%lit2                       , & !
    lit3      => v1eca_bgc_index%lit3                       , & !
    som1      => v1eca_bgc_index%som1                       , & !
    som2      => v1eca_bgc_index%som2                       , & !
    som3      => v1eca_bgc_index%som3                       , & !
    cwd       => v1eca_bgc_index%cwd                        , & !
    c_loc     => v1eca_bgc_index%c_loc                      , & !
    n_loc     => v1eca_bgc_index%n_loc                      , & !
    p_loc     => v1eca_bgc_index%p_loc                      , & !
    c13_loc   => v1eca_bgc_index%c13_loc                    , & !
    c14_loc   => v1eca_bgc_index%c14_loc                    , & !
    nelms     => v1eca_bgc_index%nelms                      , & !
    lid_o2    => v1eca_bgc_index%lid_o2                     , & !
    lid_co2   => v1eca_bgc_index%lid_co2                    , & !
    lid_nh4   => v1eca_bgc_index%lid_nh4                    , & !
    lid_c14_co2=> v1eca_bgc_index%lid_c14_co2               , & !
    lid_c13_co2=> v1eca_bgc_index%lid_c13_co2               , & !
    lid_co2_hr => v1eca_bgc_index%lid_co2_hr                , &
    lid_minn_nh4_immob=> v1eca_bgc_index%lid_minn_nh4_immob , &
    lid_minp_immob => v1eca_bgc_index%lid_minp_immob        , &
    lid_minp_soluble=> v1eca_bgc_index%lid_minp_soluble     , &
    lit1_dek_reac => v1eca_bgc_index%lit1_dek_reac          , &
    lit2_dek_reac => v1eca_bgc_index%lit2_dek_reac          , &
    lit3_dek_reac => v1eca_bgc_index%lit3_dek_reac          , &
    som1_dek_reac => v1eca_bgc_index%som1_dek_reac          , &
    som2_dek_reac => v1eca_bgc_index%som2_dek_reac          , &
    som3_dek_reac => v1eca_bgc_index%som3_dek_reac          , &
    cwd_dek_reac => v1eca_bgc_index%cwd_dek_reac            , &
    lid_n2   => v1eca_bgc_index%lid_n2                      , & !
    lid_n2o   => v1eca_bgc_index%lid_n2o                    , & !
    lid_no3   => v1eca_bgc_index%lid_no3                    , & !
    lid_nh4_nit_reac => v1eca_bgc_index%lid_nh4_nit_reac    , & !
    lid_no3_den_reac => v1eca_bgc_index%lid_no3_den_reac      & !
  )
  print*,'checksum'
  print*,'som1c ','som1n ','som1p'

  reac=lit1_dek_reac
  resc = cascade_matrix((lit1-1)*nelms + c_loc, reac) + cascade_matrix((som1-1)*nelms + c_loc, reac) + &
     cascade_matrix(lid_co2, reac)
  resn = cascade_matrix((lit1-1)*nelms + n_loc, reac) + cascade_matrix((som1-1)*nelms + n_loc, reac) + &
     cascade_matrix(lid_nh4, reac)
  resp = cascade_matrix((lit1-1)*nelms + p_loc, reac) + cascade_matrix((som1-1)*nelms + p_loc, reac) + &
     cascade_matrix(lid_minp_soluble, reac)
  write(*,'(A,3(X,E20.10))')'lit1 resc, resn, resp =',resc,resn, resp

  reac = lit2_dek_reac
  resc = cascade_matrix((lit2-1)*nelms + c_loc, reac) + cascade_matrix((som1-1)*nelms + c_loc, reac) + &
     cascade_matrix(lid_co2, reac)
  resn = cascade_matrix((lit2-1)*nelms + n_loc, reac) + cascade_matrix((som1-1)*nelms + n_loc, reac) + &
     cascade_matrix(lid_nh4, reac)
  resp = cascade_matrix((lit2-1)*nelms + p_loc, reac) + cascade_matrix((som1-1)*nelms + p_loc, reac) + &
     cascade_matrix(lid_minp_soluble, reac)
  write(*,'(A,3(X,E20.10))')'lit2 resc, resn, resp =',resc,resn, resp

  reac = lit3_dek_reac
  resc = cascade_matrix((lit3-1)*nelms + c_loc, reac) + cascade_matrix((som2-1)*nelms + c_loc, reac) + &
     cascade_matrix(lid_co2, reac)
  resn = cascade_matrix((lit3-1)*nelms + n_loc, reac) + cascade_matrix((som2-1)*nelms + n_loc, reac) + &
     cascade_matrix(lid_nh4, reac)
  resp = cascade_matrix((lit3-1)*nelms + p_loc, reac) + cascade_matrix((som2-1)*nelms + p_loc, reac) + &
     cascade_matrix(lid_minp_soluble, reac)
  write(*,'(A,3(X,E20.10))')'lit3 resc, resn, resp =',resc,resn, resp

  reac = som1_dek_reac
  resc = cascade_matrix((som1-1)*nelms + c_loc, reac) + cascade_matrix((som3-1)*nelms + c_loc, reac) + &
     cascade_matrix((som2-1)*nelms + c_loc, reac) + cascade_matrix(lid_co2, reac)
  resn = cascade_matrix((som1-1)*nelms + n_loc, reac) + cascade_matrix((som3-1)*nelms + n_loc, reac) + &
     cascade_matrix((som2-1)*nelms + n_loc, reac) + cascade_matrix(lid_nh4, reac)
  resp = cascade_matrix((som1-1)*nelms + p_loc, reac) + cascade_matrix((som3-1)*nelms + p_loc, reac) + &
     cascade_matrix((som2-1)*nelms + p_loc, reac) + cascade_matrix(lid_minp_soluble, reac)
  write(*,'(A,3(X,E20.10))')'som1 resc, resn, resp =',resc,resn, resp

  reac = som2_dek_reac
  resc = cascade_matrix((som2-1)*nelms + c_loc, reac) + cascade_matrix((som1-1)*nelms + c_loc, reac) + &
     cascade_matrix((som3-1)*nelms + c_loc, reac) + cascade_matrix(lid_co2, reac)
  resn = cascade_matrix((som2-1)*nelms + n_loc, reac) + cascade_matrix((som1-1)*nelms + n_loc, reac) + &
     cascade_matrix((som3-1)*nelms + n_loc, reac) + cascade_matrix(lid_nh4, reac)
  resp = cascade_matrix((som2-1)*nelms + p_loc, reac) + cascade_matrix((som1-1)*nelms + p_loc, reac) + &
     cascade_matrix((som3-1)*nelms + p_loc, reac) + cascade_matrix(lid_minp_soluble, reac)
  write(*,'(A,3(X,E20.10))')'som2 resc, resn, resp =',resc,resn, resp

  reac = som3_dek_reac
  resc = cascade_matrix((som3-1)*nelms + c_loc, reac) + cascade_matrix((som1-1)*nelms + c_loc, reac) + &
     cascade_matrix(lid_co2, reac)
  resn = cascade_matrix((som3-1)*nelms + n_loc, reac) + cascade_matrix((som1-1)*nelms + n_loc, reac) + &
     cascade_matrix(lid_nh4, reac)
  resp = cascade_matrix((som3-1)*nelms + p_loc, reac) + cascade_matrix((som1-1)*nelms + p_loc, reac) + &
     cascade_matrix(lid_minp_soluble, reac)
  write(*,'(A,3(X,E20.10))')'som3 resc, resn, resp =',resc,resn, resp

  reac = cwd_dek_reac
  resc = cascade_matrix((cwd-1)*nelms + c_loc, reac) + cascade_matrix((lit2-1)*nelms + c_loc, reac) + &
     cascade_matrix((lit3-1)*nelms + c_loc, reac) + cascade_matrix(lid_co2, reac)
  resn = cascade_matrix((cwd-1)*nelms + n_loc, reac) + cascade_matrix((lit2-1)*nelms + n_loc, reac) + &
     cascade_matrix((lit3-1)*nelms + n_loc, reac) + cascade_matrix(lid_nh4, reac)
  resp = cascade_matrix((cwd-1)*nelms + p_loc, reac) + cascade_matrix((lit2-1)*nelms + p_loc, reac) + &
     cascade_matrix((lit3-1)*nelms + p_loc, reac) + cascade_matrix(lid_minp_soluble, reac)
  write(*,'(A,3(X,E20.10))')'cwd resc, resn, resp =',resc,resn, resp

  reac = lid_nh4_nit_reac
  resn = cascade_matrix(lid_nh4, reac) + cascade_matrix(lid_no3, reac) + 2._r8 * cascade_matrix(lid_n2o, reac)
  write(*,'(A,(X,E20.10))')'nit, resn =',resn

  reac = lid_no3_den_reac
  resn = cascade_matrix(lid_no3, reac) + cascade_matrix(lid_n2o, reac) * 2._r8 + &
    cascade_matrix(lid_n2, reac) * 2._r8
  write(*,'(A,3(X,E20.10))')'den, resn =',resn
  end associate
  end subroutine checksum_cascade

  !-------------------------------------------------------------------------------
  subroutine runbgc_v1eca(this,  is_surflit, dtime, bgc_forc, nstates, ystates0, ystatesf, bstatus)

  !DESCRIPTION
  !do bgc model integration for one step
  use JarBgcForcType        , only : JarBGC_forc_type
  use MathfuncMod           , only : pd_decomp
  use BetrStatusType        , only : betr_status_type
  use MathfuncMod           , only : safe_div
  use tracer_varcon         , only : catomw, natomw, patomw
  implicit none
  class(v1eca_bgc_type)  , intent(inout) :: this
  logical                    , intent(in)    :: is_surflit
  real(r8)                   , intent(in)    :: dtime
  type(JarBGC_forc_type)     , intent(in)    :: bgc_forc
  integer                    , intent(in)    :: nstates
  real(r8)                   , intent(out)   :: ystates0(nstates)
  real(r8)                   , intent(out)   :: ystatesf(nstates)
  type(betr_status_type)     , intent(out)   :: bstatus

  !local variables
  real(r8) :: pot_om_decay_rates(this%v1eca_bgc_index%nom_pools)
  real(r8) :: pot_co2_hr
  real(r8) :: pot_f_nit_mol_per_sec
  real(r8) :: n2_n2o_ratio_denit
  real(r8) :: yf(this%v1eca_bgc_index%nstvars)
  real(r8) :: o2_decomp_depth
  real(r8) :: time
  real(r8) :: frc_c13, frc_c14
  real(r8) :: c_mass1, n_mass1, p_mass1
  real(r8) :: c_mass2, n_mass2, p_mass2
  real(r8) :: c_flx, n_flx, p_flx
  integer :: jj
  character(len=*),parameter :: subname = 'runbgc_v1eca'
  associate(                                            &
    pctsand   => bgc_forc%pct_sand             , &  !sand in %
    rt_ar     => bgc_forc%rt_ar                , &  !root autotrophic respiration
    rt_ar_c13 => bgc_forc%rt_ar_c13            , &  !root autotrophic respiration
    rt_ar_c14 => bgc_forc%rt_ar_c14            , &  !root autotrophic respiration
    lid_nh4   => this%v1eca_bgc_index%lid_nh4        , &  !position id of nh4
    lid_no3   => this%v1eca_bgc_index%lid_no3        , &  !
    lid_o2    => this%v1eca_bgc_index%lid_o2         , &  !
    nom_pools => this%v1eca_bgc_index%nom_pools      , &  !number om pools
    nom_tot_elms=> this%v1eca_bgc_index%nom_tot_elms , &
    nstvars   => this%v1eca_bgc_index%nstvars        , &
    nprimvars => this%v1eca_bgc_index%nprimvars      , &
    nreactions=> this%v1eca_bgc_index%nreactions     , &
    lid_plant_minn_nh4  => this%v1eca_bgc_index%lid_plant_minn_nh4       , &
    lid_plant_minn_no3  => this%v1eca_bgc_index%lid_plant_minn_no3       , &
    lid_n2o_nit => this%v1eca_bgc_index%lid_n2o_nit,&
    lid_no3_den => this%v1eca_bgc_index%lid_no3_den,&
    lid_minp_soluble=> this%v1eca_bgc_index%lid_minp_soluble, &
    cascade_matrix => this%cascade_matrix             , &
    cascade_matrixp=> this%cascade_matrixp            , &
    cascade_matrixd=> this%cascade_matrixd            , &
    ystates1 => this%ystates1                           &
  )
  this%v1eca_bgc_index%debug = bgc_forc%debug
  if(this%v1eca_bgc_index%debug)print*,'enter runbgc_v1eca'
  this%rt_ar = rt_ar
  if(this%use_c13)then
    frc_c13 = safe_div(rt_ar_c13,rt_ar);
  else
    frc_c13 = 1._r8
  endif
  if(this%use_c14)then
    frc_c14 = safe_div(rt_ar_c14,rt_ar)
  else
    frc_c14 = 1._r8
  endif
  call bstatus%reset()

  !initialize state variables
  call this%init_states(this%v1eca_bgc_index, bgc_forc)
  this%competECA%bd = bgc_forc%bd
  this%competECA%h2osoi_vol = bgc_forc%h2osoi_vol
  ystates0(:) = this%ystates0(:)

!  call this%begin_massbal_check()
  !initialize decomposition scaling factors
  call this%decompkf_eca%set_decompk_scalar(ystates1(lid_o2), bgc_forc)

  !initialize all entries to zero
  cascade_matrix(:,:) = 0._r8

  !calculate default stoichiometry entries
  call this%calc_cascade_matrix(this%v1eca_bgc_index, cascade_matrix, frc_c13, frc_c14)

  !run century decomposition, return decay rates, cascade matrix, potential hr
  call this%censom%run_decomp(is_surflit, this%v1eca_bgc_index, dtime, ystates1(1:nom_tot_elms),&
      this%decompkf_eca, bgc_forc%pct_sand, bgc_forc%pct_clay,this%alpha_n, this%alpha_p, &
      cascade_matrix, this%k_decay(1:nom_pools), pot_co2_hr, bstatus)
  if(bstatus%check_status())return

  call this%nitden%calc_pot_nitr(ystates1(lid_nh4), bgc_forc, this%decompkf_eca, pot_f_nit_mol_per_sec)

  !calculate potential o2 consumption
!  print*,'pot_co2_hr',pot_co2_hr
!  print*,'rt_ar', rt_ar
!  print*,'pot_f_nit', pot_f_nit_mol_per_sec * this%nitden%get_nit_o2_scef()
  o2_decomp_depth = pot_co2_hr + rt_ar + pot_f_nit_mol_per_sec * this%nitden%get_nit_o2_scef()

  !take a minimum > 0 to avoid singularity in calculating anaerobic fractions
  o2_decomp_depth = max(o2_decomp_depth,1.e-40_r8)

  !run nitrification-denitrification, returns cascade_matrix, decay rates
  call this%nitden%run_nitden(this%v1eca_bgc_index, bgc_forc, this%decompkf_eca, &
    ystates1(lid_nh4), ystates1(lid_no3), ystates1(lid_o2), o2_decomp_depth, &
    pot_f_nit_mol_per_sec, pot_co2_hr, this%pot_f_nit, this%pot_f_denit, cascade_matrix)

  !do integration, in each integration, the stoichiometric matrix is kept as constant
  !so the reaction rate is a function of state variables. Further, for simplicity,
  !the nitrification and denitrification rates have been assumed as linear function
  !nh4 and no3 in soil.
  call this%arenchyma_gas_transport(this%v1eca_bgc_index, dtime)

  !do the stoichiometric matrix separation
  call pd_decomp(nprimvars, nreactions, cascade_matrix(1:nprimvars, 1:nreactions), &
     cascade_matrixp, cascade_matrixd, bstatus)
  if(bstatus%check_status())return
  !save a copy of nh4 stoichiometry
  this%cascade_matnh4(1:this%v1eca_bgc_index%nom_pools)=cascade_matrix(this%v1eca_bgc_index%lid_nh4,1:this%v1eca_bgc_index%nom_pools)
  time = 0._r8
  yf(:) = ystates1(:)

!  do jj=1,this%v1eca_bgc_index%nom_pools
!   call print_reaction(this%v1eca_bgc_index, cascade_matrix, jj)
!  enddo
!  call this%end_massbal_check('midd runbgc')

!  call this%checksum_cascade(this%v1eca_bgc_index)
  call ode_adapt_ebbks1(this, yf, nprimvars, nstvars, time, dtime, ystates1)

  !if(this%v1eca_bgc_index%debug)call this%checksum_cascade(this%v1eca_bgc_index)
  if(this%use_c14)then
    call this%c14decay(this%v1eca_bgc_index, dtime, ystates1)
  endif

  if(this%batch_mode)call this%sum_tot_store(nstvars, ystates1)
  ystatesf(:) = ystates1(:)
  ystatesf(this%v1eca_bgc_index%lid_pot_co2_hr) = pot_co2_hr*dtime
  ystatesf(this%v1eca_bgc_index%lid_o_scalar) = this%decompkf_eca%o_scalar
!  call this%end_massbal_check('end runbgc')
  end associate
  end subroutine runbgc_v1eca
  !-------------------------------------------------------------------------------
  subroutine c14decay(this, v1eca_bgc_index, dtime, ystates1)

  !apply c14 decay to om pools
  use v1ecaBGCIndexType       , only : v1eca_bgc_index_type

  implicit none
  ! !ARGUMENTS:
  class(v1eca_bgc_type)     , intent(in) :: this
  type(v1eca_bgc_index_type)   , intent(in) :: v1eca_bgc_index
  real(r8), intent(in) :: dtime
  real(r8), intent(inout) :: ystates1(:)

  integer :: jj
  integer :: kc14

  associate(                                &
    litr_beg =>  v1eca_bgc_index%litr_beg, &
    Bm_beg =>  v1eca_bgc_index%Bm_beg    , &
    som_beg =>  v1eca_bgc_index%som_beg  , &
    pom_beg =>  v1eca_bgc_index%pom_beg  , &
    wood_beg =>  v1eca_bgc_index%wood_beg, &
    litr_end =>  v1eca_bgc_index%litr_end, &
    som_end =>  v1eca_bgc_index%som_end  , &
    Bm_end =>  v1eca_bgc_index%Bm_end    , &
    pom_end => v1eca_bgc_index%pom_end   , &
    wood_end =>  v1eca_bgc_index%wood_end, &
    c14_loc=>  v1eca_bgc_index%c14_loc   , &
    nelms => v1eca_bgc_index%nelms         &
  )

  call somc14_decay(litr_beg, litr_end, nelms, c14_loc, this%c14decay_const)

  call somc14_decay(wood_beg, wood_end, nelms, c14_loc, this%c14decay_const)

  call somc14_decay(som_beg, som_end, nelms, c14_loc, this%c14decay_som_const)

  call somc14_decay(pom_beg, pom_end, nelms, c14_loc, this%c14decay_pom_const)

  call somc14_decay(Bm_beg, Bm_end, nelms, c14_loc, this%c14decay_Bm_const)

  end associate
  contains
    subroutine somc14_decay(ebeg, eend, nelms, c14_loc, decay_const)
    implicit none
    integer, intent(in) :: ebeg, eend, nelms
    integer, intent(in) :: c14_loc
    real(r8),intent(in) :: decay_const
    integer :: jj, kc14
    do jj = ebeg, eend, nelms
      kc14=jj-1+c14_loc
      ystates1(kc14) = ystates1(kc14)*exp(- decay_const * dtime)
    enddo

    end subroutine somc14_decay
  end subroutine c14decay
  !-------------------------------------------------------------------------------
  subroutine calc_cascade_matrix(this,v1eca_bgc_index, cascade_matrix, frc_c13, frc_c14)
    !
    ! !DESCRIPTION:
    ! calculate cascade matrix for the decomposition model
    !
    ! !USES:
    use MathfuncMod               , only : safe_div
    use v1ecaBGCIndexType       , only : v1eca_bgc_index_type
    use betr_ctrl                 , only : spinup_state => betr_spinup_state
    implicit none
    ! !ARGUMENTS:
    class(v1eca_bgc_type)     , intent(in) :: this
    type(v1eca_bgc_index_type)   , intent(in) :: v1eca_bgc_index
    real(r8)                      , intent(inout)   :: cascade_matrix(v1eca_bgc_index%nstvars, v1eca_bgc_index%nreactions)
    real(r8)                      , intent(in) :: frc_c13, frc_c14
    ! !LOCAL VARIABLES:
    real(r8) :: ftxt, f1, f2
    integer :: k, reac

    associate(                                                             & !
         lid_autr_rt => v1eca_bgc_index%lid_autr_rt                      , & !
         lid_o2    => v1eca_bgc_index%lid_o2                             , & !
         lid_co2   => v1eca_bgc_index%lid_co2                            , & !
         lid_c13_co2=> v1eca_bgc_index%lid_c13_co2                       , & !
         lid_c14_co2=> v1eca_bgc_index%lid_c14_co2                       , & !
         lid_nh4   => v1eca_bgc_index%lid_nh4                            , & !
         lid_ch4   => v1eca_bgc_index%lid_ch4                            , & !
         lid_ar    => v1eca_bgc_index%lid_ar                             , & !
         lid_no3   => v1eca_bgc_index%lid_no3                            , & !
         lid_n2o   => v1eca_bgc_index%lid_n2o                            , & !
         lid_n2    => v1eca_bgc_index%lid_n2                             , & !
         lid_co2_hr=> v1eca_bgc_index%lid_co2_hr                         , & !
         lid_minn_nh4_immob => v1eca_bgc_index%lid_minn_nh4_immob        , & !
         lid_minn_no3_immob => v1eca_bgc_index%lid_minn_no3_immob        , & !
         lid_n2_paere => v1eca_bgc_index%lid_n2_paere                    , & !
         lid_ch4_paere => v1eca_bgc_index%lid_ch4_paere                  , & !
         lid_n2o_paere => v1eca_bgc_index%lid_n2o_paere                  , & !
         lid_o2_paere => v1eca_bgc_index%lid_o2_paere                    , & !
         lid_ar_paere => v1eca_bgc_index%lid_ar_paere                    , & !
         lid_co2_paere => v1eca_bgc_index%lid_co2_paere                  , & !
         lid_c13_co2_paere => v1eca_bgc_index%lid_c13_co2_paere          , & !
         lid_c14_co2_paere => v1eca_bgc_index%lid_c14_co2_paere          , & !
         lid_minp_soluble => v1eca_bgc_index%lid_minp_soluble            , & !
         lid_minp_sorb => v1eca_bgc_index%lid_minp_sorb                  , & !
         lid_plant_minp => v1eca_bgc_index%lid_plant_minp                , & !
         lid_minp_immob => v1eca_bgc_index%lid_minp_immob                , & !

         lid_autr_rt_reac=> v1eca_bgc_index%lid_autr_rt_reac                 , & !
         lid_no3_den  => v1eca_bgc_index%lid_no3_den                     , & !
         lid_plant_minn_nh4_up_reac=> v1eca_bgc_index%lid_plant_minn_nh4_up_reac , & !
         lid_plant_minn_no3_up_reac=> v1eca_bgc_index%lid_plant_minn_no3_up_reac , & !
         lid_plant_minn_nh4  => v1eca_bgc_index%lid_plant_minn_nh4       , &
         lid_plant_minn_no3  => v1eca_bgc_index%lid_plant_minn_no3       , &
         lid_minp_soluble_to_labile_reac => v1eca_bgc_index%lid_minp_soluble_to_labile_reac      , & !
         lid_plant_minp_up_reac => v1eca_bgc_index%lid_plant_minp_up_reac, & !

         lid_n2_aren_reac => v1eca_bgc_index%lid_n2_aren_reac            , & !
         lid_ch4_aren_reac=> v1eca_bgc_index%lid_ch4_aren_reac           , & !
         lid_n2o_aren_reac=> v1eca_bgc_index%lid_n2o_aren_reac           , & !
         lid_o2_aren_reac => v1eca_bgc_index%lid_o2_aren_reac            , & !
         lid_ar_aren_reac => v1eca_bgc_index%lid_ar_aren_reac            , & !
         lid_co2_aren_reac=> v1eca_bgc_index%lid_co2_aren_reac           , & !
         lid_c13_co2_aren_reac=> v1eca_bgc_index%lid_c13_co2_aren_reac   , & !
         lid_c14_co2_aren_reac=> v1eca_bgc_index%lid_c14_co2_aren_reac     & !
         )

    !higher [nh4] makes lower [no3] competitiveness
    !note all reactions are in the form products - substrates = 0, therefore
    !mass balance is automatically ensured.
    !set up first order reactions

    !---------------------------------------------------------------------------------
    !reaction 10, inorganic P non-equilibrium adsorption
    !P_soluble -> p_secondary
    reac = lid_minp_soluble_to_labile_reac
    cascade_matrix(lid_minp_soluble,  reac) = -1._r8
    cascade_matrix(lid_minp_sorb, reac) = 1._r8

    if(this%plant_ntypes>0)then
      !----------------------------------------------------------------------
      !reaction 12, plant mineral nitrogen nh4 uptake
      reac = lid_plant_minn_nh4_up_reac
      !  nh4  -> plant_nitrogen
      cascade_matrix(lid_nh4, reac)        = -1._r8
      cascade_matrix(lid_plant_minn_nh4, reac) = 1._r8

      !----------------------------------------------------------------------
      !reaction 13, plant mineral nitrogen no3 uptake
      reac = lid_plant_minn_no3_up_reac
      !  no3  -> plant_nitrogen
      cascade_matrix(lid_no3, reac)        = -1._r8
      cascade_matrix(lid_plant_minn_no3, reac) = 1._r8

      !----------------------------------------------------------------------
      !reaction 14, plant mineral phosphorus uptake
      reac = lid_plant_minp_up_reac
      ! p_solution -> plant_p
      cascade_matrix(lid_minp_soluble, reac) = -1._r8
      cascade_matrix(lid_plant_minp, reac) = 1._r8

      !----------------------------------------------------------------------
      !reaction 15, ar + o2 -> co2
      reac = lid_autr_rt_reac
      cascade_matrix(lid_co2, reac) =  1._r8
      cascade_matrix(lid_o2,  reac) = -1._r8

      if(this%use_c13)then
        cascade_matrix(lid_c13_co2, reac) =  1._r8 * frc_c13
      endif

      if(this%use_c14)then
        cascade_matrix(lid_c14_co2, reac) =  1._r8 * frc_c14
      endif
    endif
    !--------------------------------------------------------------------
    !arenchyma transport
    !second primary variables
    reac                               = lid_o2_aren_reac
    cascade_matrix(lid_o2, reac)       = -1._r8
    cascade_matrix(lid_o2_paere, reac) = 1._r8
    if ( spinup_state == 0 ) then
       reac                                = lid_ch4_aren_reac
       cascade_matrix(lid_ch4, reac)       = -1._r8
       cascade_matrix(lid_ch4_paere, reac) = 1._r8

       reac                                = lid_ar_aren_reac
       cascade_matrix(lid_ar, reac)        = -1._r8
       cascade_matrix(lid_ar_paere, reac)  = 1._r8

       reac                                = lid_co2_aren_reac
       cascade_matrix(lid_co2, reac)       = -1._r8
       cascade_matrix(lid_co2_paere, reac) = 1._r8

       if(this%use_c13)then
         reac                                = lid_c13_co2_aren_reac
         cascade_matrix(lid_c13_co2, reac)       = -1._r8
         cascade_matrix(lid_c13_co2_paere, reac) = 1._r8
       endif

       if(this%use_c14)then
         reac                                = lid_c14_co2_aren_reac
         cascade_matrix(lid_c14_co2, reac)       = -1._r8
         cascade_matrix(lid_c14_co2_paere, reac) = 1._r8
       endif

       reac                                = lid_n2o_aren_reac
       cascade_matrix(lid_n2o, reac)       = -1._r8
       cascade_matrix(lid_n2o_paere, reac) = 1._r8

       reac                                = lid_n2_aren_reac
       cascade_matrix(lid_n2, reac)        = -1._r8
       cascade_matrix(lid_n2_paere, reac)  = 1._r8
    endif

  end associate
  end subroutine calc_cascade_matrix
  !--------------------------------------------------------------------
  subroutine init_states(this, v1eca_bgc_index, bgc_forc)

  use v1ecaBGCIndexType       , only : v1eca_bgc_index_type
  use JarBgcForcType            , only : JarBGC_forc_type
  implicit none
  class(v1eca_bgc_type)     , intent(inout) :: this
  type(v1eca_bgc_index_type)  , intent(in) :: v1eca_bgc_index
  type(JarBGC_forc_type)  , intent(in) :: bgc_forc

  integer :: j
  associate(                           &
    lid_n2 => v1eca_bgc_index%lid_n2, &
    lid_o2 => v1eca_bgc_index%lid_o2, &
    lid_co2 => v1eca_bgc_index%lid_co2, &
    lid_c13_co2 => v1eca_bgc_index%lid_c13_co2, &
    lid_c14_co2 => v1eca_bgc_index%lid_c14_co2, &
    lid_n2o => v1eca_bgc_index%lid_n2o, &
    lid_ar => v1eca_bgc_index%lid_ar, &
    lid_ch4 => v1eca_bgc_index%lid_ch4,  &
    lid_co2_hr => v1eca_bgc_index%lid_co2_hr, &
    lid_pot_co2_hr => v1eca_bgc_index%lid_pot_co2_hr, &
    lid_n2o_nit => this%v1eca_bgc_index%lid_n2o_nit,&
    lid_nh4_nit => this%v1eca_bgc_index%lid_nh4_nit, &
    lid_no3_den => this%v1eca_bgc_index%lid_no3_den,  &
    lid_n2o_den => this%v1eca_bgc_index%lid_n2o_den,  &
    lid_minn_nh4_immob=> v1eca_bgc_index%lid_minn_nh4_immob , &
    lid_minn_no3_immob => v1eca_bgc_index%lid_minn_no3_immob, &
    lid_minp_immob => v1eca_bgc_index%lid_minp_immob         &
  )
  this%ystates0(:) = bgc_forc%ystates(:)
  this%ystates0(lid_co2_hr) = 0._r8
  this%ystates0(lid_n2o_nit)= 0._r8
  this%ystates0(lid_no3_den)= 0._r8
  this%ystates0(lid_nh4_nit)= 0._r8
  this%ystates0(lid_n2o_den)= 0._r8
  this%ystates0(lid_pot_co2_hr)=0._r8
  this%ystates0(lid_minn_nh4_immob) =0._r8
  this%ystates0(lid_minn_no3_immob) =0._r8
  this%ystates0(lid_minp_immob) =0._r8
  this%ystates1(:) = this%ystates0(:)

  !set conversion parameters for arenchyma transport
  this%scal_f(lid_n2) = bgc_forc%aren_cond_n2
  this%conc_f(lid_n2) = bgc_forc%conc_atm_n2
  this%conv_f(lid_n2) = 1._r8/bgc_forc%n2_g2b

  this%scal_f(lid_o2) = bgc_forc%aren_cond_o2
  this%conc_f(lid_o2) = bgc_forc%conc_atm_o2
  this%conv_f(lid_o2) = 1._r8/bgc_forc%o2_g2b

  this%scal_f(lid_ar) = bgc_forc%aren_cond_ar
  this%conc_f(lid_ar) = bgc_forc%conc_atm_ar
  this%conv_f(lid_ar) = 1._r8/bgc_forc%ar_g2b

  this%scal_f(lid_co2) = bgc_forc%aren_cond_co2
  this%conc_f(lid_co2) = bgc_forc%conc_atm_co2
  this%conv_f(lid_co2) = 1._r8/bgc_forc%co2_g2b

  if(this%use_c13)then
    this%scal_f(lid_c13_co2) = bgc_forc%aren_cond_co2_c13
    this%conc_f(lid_c13_co2) = bgc_forc%conc_atm_co2_c13
    this%conv_f(lid_c13_co2) = 1._r8/bgc_forc%co2_g2b
  endif

  if(this%use_c14)then
    this%scal_f(lid_c14_co2) = bgc_forc%aren_cond_co2_c14
    this%conc_f(lid_c14_co2) = bgc_forc%conc_atm_co2_c14
  endif

  this%scal_f(lid_ch4) = bgc_forc%aren_cond_ch4
  this%conc_f(lid_ch4) = bgc_forc%conc_atm_ch4
  this%conv_f(lid_ch4) = 1._r8/bgc_forc%ch4_g2b

  this%scal_f(lid_n2o) = bgc_forc%aren_cond_n2o
  this%conc_f(lid_n2o) = bgc_forc%conc_atm_n2o
  this%conv_f(lid_n2o) = 1._r8/bgc_forc%n2o_g2b

  this%plant_ntypes = bgc_forc%plant_ntypes
  this%soilorder = bgc_forc%soilorder

  end associate
  end subroutine init_states
  !--------------------------------------------------------------------
  subroutine bgc_integrate(this, ystate, dtime, time, nprimvars, nstvars, dydt)
  !
  !DESCRIPTION
  !calculate the reaction rates
  !In current implementation, no active adsorption of NH4 is involved.
  !The inorganic phosphorus does the transition from soluble->secondary->occlude
  use SOMStateVarUpdateMod , only : calc_dtrend_som_bgc
  use MathfuncMod          , only : lom_type, safe_div
  use tracer_varcon        , only : patomw
  implicit none
  class(v1eca_bgc_type) , intent(inout) :: this
  integer                   , intent(in) :: nstvars
  integer                   , intent(in) :: nprimvars
  real(r8)                  , intent(in) :: dtime
  real(r8)                  , intent(in) :: time
  real(r8)                  , intent(in) :: ystate(nstvars)
  real(r8)                  , intent(out) :: dydt(nstvars)

  !local variables
  real(r8) :: mic_pot_nn_flx  !potential nitrogen uptake to support decomposition
  real(r8) :: mic_pot_np_flx_up  !potential phosphorus uptake to support decomposition
  real(r8) :: mic_pot_np_flx_mn
  real(r8) :: pot_decomp(this%v1eca_bgc_index%nom_pools)
  real(r8) :: rrates(this%v1eca_bgc_index%nreactions)
  real(r8) :: p_dt(1:this%v1eca_bgc_index%nprimvars)
  real(r8) :: d_dt(1:this%v1eca_bgc_index%nprimvars)
  real(r8) :: dydt1(nstvars)
  real(r8) :: pscal(1:nprimvars)
  real(r8) :: rscal(1:this%v1eca_bgc_index%nreactions)
  real(r8) :: ECA_flx_nh4_plants(this%plant_ntypes)
  real(r8) :: ECA_flx_no3_plants(this%plant_ntypes)
  real(r8) :: ECA_flx_phosphorus_plants(this%plant_ntypes)
  real(r8) :: ECA_factor_minp_msurf
  real(r8) :: ECA_factor_phosphorus_mic
  real(r8) :: ECA_factor_nh4_mic
  real(r8) :: ECA_factor_no3_mic
  real(r8) :: ECA_factor_nitrogen_mic
  real(r8) :: ECA_factor_den
  real(r8) :: ECA_factor_nit
  real(r8) :: minp_soluble
  real(r8) :: minn_nh4_sol
  real(r8) :: minn_no3_sol
  real(r8) :: dminn
  real(r8) :: sol_smin_nh4, sol_smin_no3
  real(r8) :: sol_sminp_soluble
  integer  :: jj, it,reac
  integer, parameter  :: itmax = 10
  type(lom_type) :: lom
  type(betr_status_type) :: bstatus
  logical :: lneg
  real(r8) :: scal
  real(r8) :: adsorb_to_labilep
  real(r8) :: desorb_to_solutionp
  real(r8) :: dsolutionp_dt
  associate(                                                                                      &
    nreactions => this%v1eca_bgc_index%nreactions                                              , &
    nom_pools => this%v1eca_bgc_index%nom_pools                                                , &
    lid_nh4 => this%v1eca_bgc_index%lid_nh4                                                    , &
    lid_no3 => this%v1eca_bgc_index%lid_no3                                                    , &
    lid_co2_hr => this%v1eca_bgc_index%lid_co2_hr                                              , &
    lid_co2_somhr => this%v1eca_bgc_index%lid_co2_somhr                                        , &
    lid_co2_lithr => this%v1eca_bgc_index%lid_co2_lithr                                        , &
    lid_co2_cwdhr => this%v1eca_bgc_index%lid_co2_cwdhr                                        , &
    lid_cum_closs=> this%v1eca_bgc_index%lid_cum_closs                                         , &
    varnames=> this%v1eca_bgc_index%varnames                                                   , &
    som1    => this%v1eca_bgc_index%som1                                                       , &
    som2    => this%v1eca_bgc_index%som2                                                       , &
    som3    => this%v1eca_bgc_index%som3                                                       , &
    lit1    => this%v1eca_bgc_index%lit1                                                       , &
    lit2    => this%v1eca_bgc_index%lit2                                                       , &
    lit3    => this%v1eca_bgc_index%lit3                                                       , &
    cwd     => this%v1eca_bgc_index%cwd                                                        , &
    c_loc   => this%v1eca_bgc_index%c_loc                                                      , &
    n_loc   => this%v1eca_bgc_index%n_loc                                                      , &
    p_loc   => this%v1eca_bgc_index%p_loc                                                      , &
    nelms   => this%v1eca_bgc_index%nelms                                                      , &
    lid_plant_minn_no3_pft=> this%v1eca_bgc_index%lid_plant_minn_no3_pft                       , &
    lid_plant_minn_nh4_pft=> this%v1eca_bgc_index%lid_plant_minn_nh4_pft                       , &
    lid_plant_minp_pft=> this%v1eca_bgc_index%lid_plant_minp_pft                               , &
    lid_plant_minp    => this%v1eca_bgc_index%lid_plant_minp                                   , &
    lid_plant_minn_nh4 => this%v1eca_bgc_index%lid_plant_minn_nh4                              , &
    lid_plant_minn_no3 => this%v1eca_bgc_index%lid_plant_minn_no3                              , &
    lid_minn_nh4_immob => this%v1eca_bgc_index%lid_minn_nh4_immob                              , & !
    lid_minn_no3_immob => this%v1eca_bgc_index%lid_minn_no3_immob                              , & !
    lid_minp_soluble => this%v1eca_bgc_index%lid_minp_soluble                                  , &
    lid_minp_sorb => this%v1eca_bgc_index%lid_minp_sorb                                        , &
    lid_minp_soluble_to_labile_reac=> this%v1eca_bgc_index%lid_minp_soluble_to_labile_reac     , &
    lid_autr_rt_reac => this%v1eca_bgc_index%lid_autr_rt_reac                                  , &
    lid_nh4_nit_reac => this%v1eca_bgc_index%lid_nh4_nit_reac                                  , &
    lid_no3_den_reac => this%v1eca_bgc_index%lid_no3_den_reac                                  , &
    lid_plant_minn_nh4_up_reac => this%v1eca_bgc_index%lid_plant_minn_nh4_up_reac              , &
    lid_plant_minn_no3_up_reac => this%v1eca_bgc_index%lid_plant_minn_no3_up_reac              , &
    lid_plant_minp_up_reac => this%v1eca_bgc_index%lid_plant_minp_up_reac                      , &
    lid_supp_minp => this%v1eca_bgc_index%lid_supp_minp                                        , &
    lid_supp_minn => this%v1eca_bgc_index%lid_supp_minn                                        , &
    cascade_matrix=>  this%cascade_matrix                                                        &
  )

!  print*,'somc',ystate((som1-1)*nelms+c_loc)
!  print*,'kdecay',this%k_decay(1:nom_pools)
  dydt(:) = 0._r8
  rrates(:) = 0._r8
  !calculate reaction rates, because arenchyma transport is
  !done, now only calculate for those that are actively
  !reacting. These include: OM pools, plant nutrient uptake
  !microbial nutrient uptake

  call this%censom%calc_pot_min_np_flx(dtime, this%v1eca_bgc_index,  ystate, this%k_decay,&
    this%cascade_matrix, this%alpha_n, this%alpha_p, pot_decomp, mic_pot_nn_flx, mic_pot_np_flx_up, &
    mic_pot_np_flx_mn)

  !update dsolutionp_dt, in g/m3/s
  dsolutionp_dt = this%competECA%dsolutionp_dt+ (mic_pot_np_flx_mn -mic_pot_np_flx_up)   !

  !---------------------
  !update adsorb_to_labilep and desorb_to_solutionp, eq. (A12) in Zhu et al., 2016
  adsorb_to_labilep = (this%competECA%vmax_minsurf_p * this%competECA%kaff_minp_msurf) / &
                       (this%competECA%kaff_minp_msurf+ystate(lid_minp_soluble))**2._r8 * dsolutionp_dt
  ! sign convention: if adsorb_to_labilep_vr(c,j) < 0, then it's desorption
  if (adsorb_to_labilep >= 0._r8) then
    adsorb_to_labilep = max(min(adsorb_to_labilep, &
       this%competECA%vmax_minsurf_p-this%competECA%dlabp_dt),0.0_r8)
    desorb_to_solutionp = 0.0_r8
  else
    desorb_to_solutionp = min(-1._r8*adsorb_to_labilep, this%competECA%dlabp_dt)
    adsorb_to_labilep = 0.0_r8
  end if
  !---------------------

  !do ECA nutrient scaling
  !
  sol_smin_nh4 = ystate(lid_nh4)/(this%competECA%bd*2.76_r8 + this%competECA%h2osoi_vol)
  sol_smin_no3 = ystate(lid_no3)/this%competECA%h2osoi_vol
!  print*,'connh4',ystate(lid_nh4),sol_smin_nh4,dtime
!  print*,'conno3',ystate(lid_no3),sol_smin_no3
  call this%competECA%run_compet_nitrogen(this%non_limit, sol_smin_nh4, sol_smin_no3,&
     this%plant_ntypes, ECA_factor_nit, &
     ECA_factor_den, ECA_factor_nh4_mic, ECA_factor_no3_mic, &
     ECA_flx_nh4_plants,ECA_flx_no3_plants)
!  print*,'ecafnh4', ECA_factor_nit,ECA_factor_nh4_mic
!  print*,'ecafno3', ECA_factor_den,ECA_factor_no3_mic
  ECA_factor_nitrogen_mic = ECA_factor_nh4_mic + ECA_factor_no3_mic
  sol_sminp_soluble=ystate(lid_minp_soluble)/this%competECA%h2osoi_vol
  call this%competECA%run_compet_phosphorus(this%nop_limit, sol_sminp_soluble,  &
      this%plant_ntypes,  ECA_factor_phosphorus_mic, ECA_factor_minp_msurf,&
      ECA_flx_phosphorus_plants)

  !apply ECA factor to obtain actual reaction rate, decomposition
  !plant, nit, den nutrient uptake,
  do jj = 1, nom_pools
    if(this%alpha_n(jj)>0._r8 .and. (.not. this%non_limit))then
      this%cascade_matrixd(lid_no3,jj) = this%cascade_matnh4(jj) * &
          safe_div(ECA_factor_no3_mic,ECA_factor_nitrogen_mic)
      this%cascade_matrixd(lid_nh4,jj) =  this%cascade_matnh4(jj)-this%cascade_matrixd(lid_no3,jj)

      this%cascade_matrix(lid_nh4,jj) = this%cascade_matrixd(lid_nh4,jj)
      this%cascade_matrix(lid_no3,jj) = this%cascade_matrixd(lid_no3,jj)

      this%cascade_matrix(lid_minn_no3_immob,jj) = - this%cascade_matrix(lid_no3,jj)
      this%cascade_matrix(lid_minn_nh4_immob,jj) = - this%cascade_matrix(lid_minn_nh4_immob,jj)

    endif
!    print*,'alpha',jj,this%alpha_n(jj)
    rrates(jj) = pot_decomp(jj)
  enddo
  rrates(lid_nh4_nit_reac) = this%pot_f_nit  !*ECA_factor_nit
  rrates(lid_no3_den_reac) = this%pot_f_denit !*ECA_factor_den
  rrates(lid_minp_soluble_to_labile_reac) = adsorb_to_labilep !*  ECA_factor_minp_msurf  !calculate from eca competition
!  print*,'ixxrrates',rrates(5:7)
  if(this%plant_ntypes>0)then
    rrates(lid_autr_rt_reac) = this%rt_ar                            !authotrophic respiration
    rrates(lid_plant_minn_no3_up_reac) = sum(ECA_flx_no3_plants)     !calculate by ECA competition
    rrates(lid_plant_minn_nh4_up_reac) = sum(ECA_flx_nh4_plants)     !calculate by ECA competition
    rrates(lid_plant_minp_up_reac) =     sum(ECA_flx_phosphorus_plants) !calculate by ECA competition
  endif

  if(this%plant_ntypes==1)then
    do jj = 1, this%plant_ntypes
      this%cascade_matrix(lid_plant_minn_no3_pft(jj),lid_plant_minn_no3_up_reac) = 1._r8
      this%cascade_matrix(lid_plant_minn_nh4_pft(jj),lid_plant_minn_nh4_up_reac) = 1._r8
      this%cascade_matrix(lid_plant_minp_pft(jj),lid_plant_minp_up_reac) = 1._r8
    enddo
  elseif(this%plant_ntypes>=2)then
    do jj = 1, this%plant_ntypes-1
      this%cascade_matrix(lid_plant_minn_no3_pft(jj),lid_plant_minn_no3_up_reac) = &
           safe_div(ECA_flx_no3_plants(jj),rrates(lid_plant_minn_no3_up_reac))
      this%cascade_matrix(lid_plant_minn_nh4_pft(jj),lid_plant_minn_nh4_up_reac) = &
           safe_div(ECA_flx_nh4_plants(jj),rrates(lid_plant_minn_nh4_up_reac))
      this%cascade_matrix(lid_plant_minp_pft(jj),lid_plant_minp_up_reac) = &
           safe_div(ECA_flx_phosphorus_plants(jj),rrates(lid_plant_minp_up_reac))
    enddo
    jj = this%plant_ntypes
    this%cascade_matrix(lid_plant_minn_no3_pft(jj),lid_plant_minn_no3_up_reac) = &
      1._r8 - sum(this%cascade_matrix(lid_plant_minn_no3_pft(1:jj-1),lid_plant_minn_no3_up_reac))
    this%cascade_matrix(lid_plant_minn_nh4_pft(jj),lid_plant_minn_nh4_up_reac) = &
      1._r8 - sum(this%cascade_matrix(lid_plant_minn_nh4_pft(1:jj-1),lid_plant_minn_nh4_up_reac))
    this%cascade_matrix(lid_plant_minp_pft(jj),lid_plant_minp_up_reac) = &
      1._r8 - sum(this%cascade_matrix(lid_plant_minp_pft(1:jj-1),lid_plant_minp_up_reac))
  endif


  it=0
  rscal=0._r8
  do
    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixp(1:nprimvars, 1:nreactions), rrates, p_dt)

    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixd(1:nprimvars, 1:nreactions), rrates, d_dt)

    !update the state variables
    call lom%calc_state_pscal(nprimvars, dtime, ystate(1:nprimvars), p_dt(1:nprimvars),  d_dt(1:nprimvars), &
        pscal(1:nprimvars), lneg, bstatus)

!    print*,it,'pscal',pscal(1:nprimvars)
    if(lneg .and. it<=itmax)then
      call lom%calc_reaction_rscal(nprimvars, nreactions,  pscal(1:nprimvars), &
        this%cascade_matrixd(1:nprimvars, 1:nreactions),rscal, bstatus)
!       if(any(rscal(5:7)<1._r8))then
!        print*,'rxxscal',rscal(5:7)
!        print*,'soilc1'
!        do jj = 1, nprimvars
!          if(this%cascade_matrixd(jj,5)/=0._r8)write(*,'(I2,X,A,2(X,E20.10))')jj,varnames(jj),this%cascade_matrixd(jj,5),pscal(jj)
!        enddo
!        stop
!       endif
      call lom%apply_reaction_rscal(nreactions, rscal(1:nreactions), rrates(1:nreactions))
    else
      call calc_dtrend_som_bgc(nstvars, nreactions, this%cascade_matrix(1:nstvars, 1:nreactions), &
         rrates(1:nreactions), dydt)
      exit
    endif
    it = it + 1
  enddo
!  print*,'ffrrates',rrates(5:7)
!  print*,'soic1',dydt((som1-1)*nelms+c_loc)
!  print*,'co2',dydt(lid_co2_hr),dydt(lid_co2_somhr),dydt(lid_co2_lithr)
  if(lid_supp_minp>0)then
    !check for mineral phosphorous
    minp_soluble=dydt(lid_minp_soluble) * dtime+ ystate(lid_minp_soluble)
    if(minp_soluble<0._r8)then
      dydt(lid_supp_minp) = -minp_soluble/dtime*(1._r8+1.e-10_r8)
      dydt(lid_minp_soluble) = dydt(lid_minp_soluble)  + dydt(lid_supp_minp)
    endif
  endif
  if(lid_supp_minn>0)then
    minn_nh4_sol = dydt(lid_nh4) * dtime + ystate(lid_nh4)
    minn_no3_sol = dydt(lid_no3) * dtime + ystate(lid_no3)
    if(minn_nh4_sol<0._r8)then
      dminn= -minn_nh4_sol/dtime*(1._r8+1.e-10_r8)
      dydt(lid_nh4) = dydt(lid_nh4) + dminn
      dydt(lid_supp_minn) = dminn
    endif
    if(minn_no3_sol<0._r8)then
      dminn= -minn_no3_sol/dtime*(1._r8+1.e-10_r8)
      dydt(lid_no3) = dydt(lid_no3) + dminn
      dydt(lid_supp_minn) = dydt(lid_supp_minn)+ dminn
    endif
  endif
  if(this%batch_mode)dydt(lid_cum_closs)=dydt(lid_co2_hr)
  !print*,'reac'
  !print*,(jj,this%cascade_matrix(lid_co2_hr, jj),rrates(jj),new_line('A'),jj=1,nreactions)

  end associate
  end subroutine bgc_integrate
  !--------------------------------------------------------------------
  subroutine arenchyma_gas_transport(this, v1eca_bgc_index, dtime)
  use v1ecaBGCIndexType       , only : v1eca_bgc_index_type
  implicit none
  class(v1eca_bgc_type)     , intent(inout) :: this
  type(v1eca_bgc_index_type)  , intent(in) :: v1eca_bgc_index
  real(r8), intent(in) :: dtime

  integer :: j
  real(r8) :: y0
  associate(                             &
    lid_n2 => v1eca_bgc_index%lid_n2,   &
    lid_o2 => v1eca_bgc_index%lid_o2,   &
    lid_co2 => v1eca_bgc_index%lid_co2, &
    lid_c13_co2 => v1eca_bgc_index%lid_c13_co2, &
    lid_c14_co2 => v1eca_bgc_index%lid_c14_co2, &
    lid_n2o => v1eca_bgc_index%lid_n2o, &
    lid_ar => v1eca_bgc_index%lid_ar,   &
    lid_ch4 => v1eca_bgc_index%lid_ch4  &
  )

  j = lid_o2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(v1eca_bgc_index%lid_o2_paere) = this%ystates1(j)-y0

  if( spinup_state == 0)then
    j = lid_n2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(v1eca_bgc_index%lid_n2_paere) = this%ystates1(j)-y0

    j = lid_ar; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(v1eca_bgc_index%lid_ar_paere) = this%ystates1(j)-y0

    j = lid_ch4; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(v1eca_bgc_index%lid_ch4_paere) = this%ystates1(j)-y0

    j = lid_co2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(v1eca_bgc_index%lid_co2_paere) = this%ystates1(j)-y0

    if(this%use_c13)then
      j = lid_c13_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(v1eca_bgc_index%lid_c13_co2_paere) = this%ystates1(j)-y0

    endif

    if(this%use_c14)then
      j = lid_c14_co2; y0=this%ystates1(j)
      call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
      this%ystates1(v1eca_bgc_index%lid_c14_co2_paere) = this%ystates1(j)-y0

    endif
    j = lid_n2o; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(v1eca_bgc_index%lid_n2o_paere) = this%ystates1(j)-y0
  endif
  end associate
  contains
    subroutine exp_ode_int(dt, c1, c2, c3, y0)
    !
    ! DESCRIPTION
    ! solve dy/dt=-c1*(c2*y-c3) using analytic solution
    implicit none
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: c1
    real(r8), intent(in) :: c2
    real(r8), intent(in) :: c3
    real(r8), intent(inout) :: y0

    if(c1>0._r8)then
      y0 = c3/c2+(y0-c3/c2)*exp(-c1/c2*dtime)
    endif

    end subroutine exp_ode_int
  end subroutine arenchyma_gas_transport

  !--------------------------------------------------------------------
  subroutine sumup_cnp_msflx_v1eca(this, ystates1, c_mass, n_mass, p_mass,c_flx,n_flx,p_flx, bstatus)

  use tracer_varcon, only : catomw, natomw, patomw
  implicit none
  class(v1eca_bgc_type)     , intent(in) :: this
  real(r8), intent(in)  :: ystates1(:)
  real(r8), intent(out) :: c_mass, n_mass, p_mass
  real(r8), optional, intent(out) :: c_flx, n_flx, p_flx
  type(betr_status_type), optional, intent(out) :: bstatus
  !local variables

  integer :: kc, kn, kp, jj
  associate(                        &
    c_loc=>  this%v1eca_bgc_index%c_loc,&
    n_loc=>  this%v1eca_bgc_index%n_loc,&
    p_loc=>  this%v1eca_bgc_index%p_loc,&
    lid_n2o_nit => this%v1eca_bgc_index%lid_n2o_nit,&
    lid_no3_den => this%v1eca_bgc_index%lid_no3_den,&
    lit1 =>  this%v1eca_bgc_index%lit1, &
    lit2 =>  this%v1eca_bgc_index%lit2, &
    lit3 =>  this%v1eca_bgc_index%lit3, &
    cwd =>   this%v1eca_bgc_index%cwd, &

    som1 =>  this%v1eca_bgc_index%som1, &
    som2 =>  this%v1eca_bgc_index%som2, &
    som3 =>  this%v1eca_bgc_index%som3, &
    nelms => this%v1eca_bgc_index%nelms, &
    lid_nh4=> this%v1eca_bgc_index%lid_nh4, &
    lid_no3=> this%v1eca_bgc_index%lid_no3, &
    lid_plant_minn_nh4 => this%v1eca_bgc_index%lid_plant_minn_nh4, &
    lid_plant_minn_no3 => this%v1eca_bgc_index%lid_plant_minn_no3, &
    lid_co2_hr => this%v1eca_bgc_index%lid_co2_hr, &
    lid_minp_soluble =>  this%v1eca_bgc_index%lid_minp_soluble,  &
    lid_minp_sorb => this%v1eca_bgc_index%lid_minp_sorb, &
    lid_plant_minp => this%v1eca_bgc_index%lid_plant_minp  &
  )
  if(present(bstatus))then
    call bstatus%reset()
    SHR_ASSERT_ALL((size(ystates1) == this%v1eca_bgc_index%nstvars), errMsg(mod_filename,__LINE__),bstatus)
  endif
  c_mass = 0._r8; n_mass = 0._r8; p_mass = 0._r8;
  if(present(c_flx))c_flx=0._r8
  if(present(n_flx))n_flx=0._r8
  if(present(p_flx))p_flx=0._r8
  print*,'lit',lit1,lit2,lit3,cwd,som1,som2,som3
  call sum_omjj(lit1, c_mass, n_mass, p_mass)
  call sum_omjj(lit2, c_mass, n_mass, p_mass)
  call sum_omjj(lit3, c_mass, n_mass, p_mass)
  call sum_omjj(cwd, c_mass, n_mass, p_mass)
  call sum_omjj(som1, c_mass, n_mass, p_mass)
  call sum_omjj(som2, c_mass, n_mass, p_mass)
  call sum_omjj(som3, c_mass, n_mass, p_mass)

  n_mass = n_mass + ystates1(lid_nh4) + ystates1(lid_no3)

  p_mass = p_mass + ystates1(lid_minp_soluble) + ystates1(lid_minp_sorb)

  c_mass = c_mass * catomw
  n_mass = n_mass * natomw
  p_mass = p_mass * patomw

  if(present(c_flx))then
    c_flx = c_flx + ystates1(lid_co2_hr)
    c_flx = c_flx * catomw
    print*,'co2_hr', ystates1(lid_co2_hr)*catomw
  endif

  if(present(n_flx))then
    n_flx=n_flx + ystates1(lid_plant_minn_nh4) + ystates1(lid_plant_minn_no3) &
     + ystates1(lid_n2o_nit) + ystates1(lid_no3_den)
    n_flx = n_flx * natomw
    print*,'n loss',ystates1(lid_plant_minn_nh4) + ystates1(lid_plant_minn_no3)+ystates1(lid_n2o_nit) + ystates1(lid_no3_den)
  endif

  if(present(p_flx))then
    p_flx=p_flx + ystates1(lid_plant_minp)
    p_flx = p_flx * patomw
    print*,'plos',ystates1(lid_plant_minp)
  endif

  end associate
  contains
    subroutine sum_omjj(jj,c_mass,n_mass,p_mass)

    implicit none
    integer,  intent(in) :: jj
    real(r8), intent(inout):: c_mass, n_mass, p_mass

    integer :: kc, kn, kp

    associate(                             &
      c_loc=>  this%v1eca_bgc_index%c_loc,&
      n_loc=>  this%v1eca_bgc_index%n_loc,&
      p_loc=>  this%v1eca_bgc_index%p_loc,&
      nelms => this%v1eca_bgc_index%nelms &
    )
    kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
    c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)
    end associate
    end subroutine sum_omjj
  end subroutine sumup_cnp_msflx_v1eca

  !--------------------------------------------------------------------
  subroutine sumup_cnp_mass(this, header, c_mass, n_mass, p_mass)
  use tracer_varcon         , only : catomw, natomw, patomw
  implicit none
  class(v1eca_bgc_type)     , intent(in) :: this
  character(len=*), intent(in) :: header
  real(r8), intent(out) :: c_mass, n_mass, p_mass
  !local variables

  integer :: kc, kn, kp, jj
  associate(                        &
    c_loc=>  this%v1eca_bgc_index%c_loc,&
    n_loc=>  this%v1eca_bgc_index%n_loc,&
    p_loc=>  this%v1eca_bgc_index%p_loc,&
    lid_n2o => this%v1eca_bgc_index%lid_n2o,&
    lid_n2 => this%v1eca_bgc_index%lid_n2,&
    lit1 =>  this%v1eca_bgc_index%lit1, &
    lit2 =>  this%v1eca_bgc_index%lit2, &
    lit3 =>  this%v1eca_bgc_index%lit3, &
    cwd =>   this%v1eca_bgc_index%cwd, &
    som1 =>  this%v1eca_bgc_index%som1, &
    som2 =>  this%v1eca_bgc_index%som2, &
    som3 =>  this%v1eca_bgc_index%som3, &
    nelms => this%v1eca_bgc_index%nelms, &
    lid_nh4=> this%v1eca_bgc_index%lid_nh4, &
    lid_no3=> this%v1eca_bgc_index%lid_no3, &
    lid_plant_minn_nh4 => this%v1eca_bgc_index%lid_plant_minn_nh4, &
    lid_plant_minn_no3 => this%v1eca_bgc_index%lid_plant_minn_no3, &
    lid_co2 => this%v1eca_bgc_index%lid_co2, &
    lid_minp_soluble =>  this%v1eca_bgc_index%lid_minp_soluble,  &
    lid_minp_sorb => this%v1eca_bgc_index%lid_minp_sorb, &
    lid_plant_minp => this%v1eca_bgc_index%lid_plant_minp, &
    ystates1 => this%ystates1  &
  )
  print*,header
  c_mass = 0._r8; n_mass = 0._r8; p_mass = 0._r8;

  jj=lit1;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=lit2;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=lit3;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=cwd;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som1;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som2;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som3;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)


  n_mass = n_mass + ystates1(lid_nh4) + ystates1(lid_no3)

  p_mass = p_mass + ystates1(lid_minp_soluble) + ystates1(lid_minp_sorb)

  c_mass = c_mass * catomw
  n_mass = n_mass * natomw
  p_mass = p_mass * patomw
  if(p_mass>1.e10_r8)then
     print*,'sum cnp bad mass',p_mass
     stop
  endif
  end associate
  end subroutine sumup_cnp_mass

  !-------------------------------------------------------------------------------
  subroutine ode_adapt_ebbks1(me, y0, nprimeq, neq, t, dt, y)
    ! !DESCRIPTION:
    !first order implicit bkks ode integration with the adaptive time stepping
    !This could be used as an example for the implementation of time-adaptive
    !mbbks1.
    ! !NOTE:
    ! this code should only be used for mass positive ODE integration
    use ODEMOD, only : ebbks, get_rerr, get_tscal
    implicit none
    ! !ARGUMENTS:
    class(v1eca_bgc_type),  intent(inout)  :: me
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: t        ! time stamp
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
    real(r8), intent(out) :: y(neq)   ! updated state variable

    ! !LOCAL VARIABLES:
    real(r8) :: yc(neq)    !coarse time stepping solution
    real(r8) :: yf(neq)    !fine time stepping solution
    real(r8) :: ycp(neq)   !temporary variable
    real(r8) :: f(neq)   ! derivative
    real(r8) :: dt2
    real(r8) :: dtr
    real(r8) :: dt05
    real(r8) :: dtmin
    real(r8) :: tt,tt2     !temporary variables
    logical  :: acc
    real(r8) :: rerr, dt_scal, pscal
    integer  :: n, nJ

    dt2=dt
    dtmin=dt/64._r8
    dtr=dt
    tt=0._r8
    !make a copy of the solution at the current time step
    yc=y0
   !get coarse grid solution
    call me%bgc_integrate(yc, dt, tt, nprimeq, neq, f)
    call ebbks(yc, f, nprimeq, neq, dt, y, pscal)
!    print*,'pscal=',pscal
    return
    do
       if(dt2<=dtmin)then
          call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)
          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)
          dtr=dtr-dt2
          tt=tt+dt2
          y=yc
       else
          !get coarse grid solution
          call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)
          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)

          !get fine grid solution
          dt05=dt2*0.5_r8
          call ebbks(y,f,nprimeq, neq,dt05, yf, pscal)
          tt2=tt+dt05
          ycp=yf
          call me%bgc_integrate(ycp, dt05, tt, nprimeq, neq, f)
          call ebbks(ycp,f,nprimeq, neq,dt05,yf,pscal)

          !determine the relative error
          rerr=get_rerr(yc,yf, neq)*exp(1._r8-1._r8/(pscal+1.e-20))

          !determine time scalar factor
          call get_tscal(rerr,dt_scal,acc)

          if(acc)then
             dtr=dtr-dt2
             tt=tt+dt2
             y=yf
          endif
          dt2=dt2*dt_scal
          dt2=min(dt2,dtr)
       endif
       if(abs(dtr/dt)<1.e-4_r8)exit
    enddo
  end subroutine ode_adapt_ebbks1

  !-------------------------------------------------------------------------------
  subroutine init_cold_v1eca(this, nstvars, ystates)
  implicit none
  class(v1eca_bgc_type)     , intent(inout) :: this
  integer , intent(in) :: nstvars
  real(r8), intent(inout) :: ystates(nstvars)
  ystates(:)=0._r8
  if(nstvars>=0)continue
  end subroutine init_cold_v1eca

  !-------------------------------------------------------------------------------
  subroutine sum_tot_store(this, nstvars, ystates)
  implicit none
  class(v1eca_bgc_type)     , intent(inout) :: this
  integer                   , intent(in)    :: nstvars
  real(r8)                  , intent(inout) :: ystates(nstvars)

  associate(                  &
    lid_totstore => this%v1eca_bgc_index%lid_totstore , &
    lit1       => this%v1eca_bgc_index%lit1        , &
    lit2       => this%v1eca_bgc_index%lit2        , &
    lit3       => this%v1eca_bgc_index%lit3        , &
    som1       => this%v1eca_bgc_index%som1        , &
    som2       => this%v1eca_bgc_index%som2        , &
    som3       => this%v1eca_bgc_index%som3        , &
    cwd        => this%v1eca_bgc_index%cwd         , &

    c_loc      => this%v1eca_bgc_index%c_loc       , & !
    nelms      => this%v1eca_bgc_index%nelms         & !
  )
  ystates(lid_totstore) =ystates((lit1-1)*nelms+c_loc) + ystates((lit2-1)*nelms+c_loc) + &
    ystates((lit3-1)*nelms+c_loc) + ystates((cwd-1)*nelms+c_loc) + &
    ystates((som1-1)*nelms+c_loc) + ystates((som2-1)*nelms+c_loc) + &
    ystates((som3-1)*nelms+c_loc)
  end associate
  end subroutine sum_tot_store
  !-------------------------------------------------------------------------------
  subroutine display_index(this)
  implicit none
  class(v1eca_bgc_type)     , intent(inout) :: this
  print*,'no3',this%v1eca_bgc_index%lid_plant_minn_no3_pft
  print*,'nh4', this%v1eca_bgc_index%lid_plant_minn_nh4_pft

  end subroutine display_index

  !-------------------------------------------------------------------------------
  subroutine print_reaction(bgc_index, cascade_matrix, reaction_id)

  use MathfuncMod, only : num2strf,num2str
  use betr_ctrl  , only : biulog
  implicit none

  type(v1eca_bgc_index_type), intent(in) :: bgc_index
  real(r8)                  , intent(in) :: cascade_matrix(bgc_index%nstvars, bgc_index%nreactions)
  integer                   , intent(in) :: reaction_id
  integer :: jj
  character(len=256) :: reaction_str
  character(len=32) :: str_loc
  real(r8) :: cef
  logical :: first
  integer :: kk
  integer :: kk1
  integer :: pos
  if(reaction_id>bgc_index%nreactions)then
     print*,'reaction does not exist'
     return
  endif

  associate(                           &
      varnames => bgc_index%varnames   &
  )
  reaction_str(:)=''
  first=.true.
  pos=0
  do jj = 1, bgc_index%nprimvars
    cef=cascade_matrix(jj,reaction_id)
    if(abs(cef)>1.e-9_r8)then
      write(str_loc,'(F10.4)')cef
      if(first)then
        first=.false.
      else
        if(cef>0._r8)then
          pos=pos+1
          reaction_str(pos:pos)='+'
        endif
      endif

      do kk = 1, len(str_loc)
        if (str_loc(kk:kk)/=' ')then
          kk1=kk
          exit
        endif
      enddo
      do kk = kk1, kk1+10
        if(ischnum(str_loc(kk:kk)))then
          pos=pos+1
          reaction_str(pos:pos)=str_loc(kk:kk)
        else
          exit
        endif
      enddo
      pos=pos+1
      reaction_str(pos:pos)='*'
      do kk=1,len(varnames(jj))
        if(ischnum(varnames(jj)(kk:kk)))then
          pos=pos+1
          reaction_str(pos:pos)=varnames(jj)(kk:kk)
        else
          exit
        endif
      enddo
    endif
  enddo
  reaction_str(pos+1:256)=''
  write(*,'(A,I2,A,A)')'reaction',reaction_id,':',reaction_str
  end associate
  end subroutine print_reaction

  logical function ischnum(a)

  implicit none
  character(len=1), intent(in) :: a

  integer :: jj
  jj = ichar(a)

  ischnum= (jj>=65 .and. jj<=90) .or. (jj>=97 .and. jj<=122) .or. (jj>=48 .and. jj<=57) .or. jj==95 .or. jj==45 .or. jj==46 .or. jj==43
  return
  end function ischnum
end module v1ecaBGCType
