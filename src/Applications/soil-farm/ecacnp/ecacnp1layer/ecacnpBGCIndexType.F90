module ecacnpBGCIndexType

  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_ctrl    , only : spinup_state => betr_spinup_state
  use gBGCIndexType  , only : gbgc_index_type
  use betr_varcon    , only : var_flux_type, var_state_type
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  integer, parameter :: loc_name_len=64

  type, private :: list_t
    character(len=loc_name_len) :: name
    integer :: id
    integer :: itype
    type(list_t), pointer :: next => null()
  end type list_t

  type, public, extends(gbgc_index_type) :: ecacnp_bgc_index_type
     integer           :: nom_pools                              !not include coarse wood debris

     integer           :: nom_tot_elms
     integer           :: lit1, lit1_dek_reac
     integer           :: lit2, lit2_dek_reac
     integer           :: lit3, lit3_dek_reac
     integer           :: som1, som1_dek_reac
     integer           :: som2, som2_dek_reac
     integer           :: som3, som3_dek_reac
     integer           :: cwd,  cwd_dek_reac
     integer           :: lwd,  lwd_dek_reac
     integer           :: fwd,  fwd_dek_reac
     integer           :: litr_beg, litr_end  !litr group
     integer           :: wood_beg, wood_end  !wood group
     integer           :: som_beg,  som_end   !som group
     integer           :: dom_beg,  dom_end   !dom group
     integer           :: Bm_beg,  Bm_end   !dom group
     integer           :: pom_beg, pom_end
     integer           :: c_loc
     integer           :: n_loc
     integer           :: p_loc
     integer           :: c13_loc = 0
     integer           :: c14_loc = 0
     integer           :: lid_supp_minp                          !supplementary mineral P for spinup purpose
     integer           :: lid_supp_minn
     integer           :: nelms                                  !number of chemical elements in an om pool
                                                                 !reactive primary variables
     integer           :: lid_nh4, lid_nh4_nit_reac              !local position of nh4 in the state variable vector
     integer           :: lid_no3, lid_no3_den_reac              !local position of no3 in the state variable vector
     integer           :: lid_plant_minn_nh4, lid_plant_minn_nh4_up_reac !local position of plant uptake of mineral nitrogen NH4 in the state variable vector
     integer           :: lid_plant_minn_no3, lid_plant_minn_no3_up_reac !
     integer           :: lid_plant_minp, lid_plant_minp_up_reac !local position of plant uptake of mineral P in the state variable vector
     integer           :: lid_minp_soluble, lid_minp_soluble_to_secp_reac    !conversation of adsorbed into secondary phase
     integer           :: lid_minp_secondary,lid_minp_secondary_to_sol_occ_reac   !local position of secondary P in the state variable vector

     integer           :: lid_minp_occlude      !local position of occluded P in the state variable vector

     integer           :: lid_autr_rt, lid_autr_rt_reac             !root autotrophic respiration

                                                                 !non reactive primary variables
     integer           :: lid_ar, lid_ar_aren_reac               !local position of ar in the state variable vector
     integer           :: lid_ch4, lid_ch4_aren_reac             !nonreactive primary variables

                                                                 !secondary variables
     integer           :: lid_o2,  lid_o2_aren_reac              !local position of o2 in the state variable vector
     integer           :: lid_co2, lid_co2_aren_reac             !local position of co2 in the state variable vector
     integer           :: lid_n2,  lid_n2_aren_reac
     integer           :: lid_n2o, lid_n2o_aren_reac
                                                                 !diagnostic variables
     integer           :: lid_n2o_nit                            !n2o production from nitrification, used to for mass balance book keeping
     integer           :: lid_co2_hr                             !co2 production from heterotrophic respiration
     integer           :: lid_c13_co2, lid_c13_co2_aren_reac
     integer           :: lid_c14_co2, lid_c14_co2_aren_reac
     integer           :: lid_no3_den                            !no3 consumption due to denitrification
     integer           :: lid_minn_nh4_immob                     !net mineral NH4 immobilization for decomposition
     integer           :: lid_minn_no3_immob                     !net mineral NO3 immobilization for decomposition
     integer           :: lid_nh4_nit
     integer           :: lid_minp_secondary_trc
     integer           :: lid_minp_occlude_trc
                                                                 !aerechyma transport, diagnostic efflux
     integer           :: lid_minp_immob                         !net P immobilization by aerobic decomposer

     integer           :: lid_ar_paere
     integer           :: lid_n2_paere
     integer           :: lid_o2_paere
     integer           :: lid_co2_paere
     integer           :: lid_c13_co2_paere
     integer           :: lid_c14_co2_paere
     integer           :: lid_ch4_paere
     integer           :: lid_n2o_paere
     integer           :: nprimvars                              !total number of primary variables
     integer           :: nstvars                                !number of equations for the state variabile vector
     integer           :: nreactions                             !seven decomposition pathways plus nitrification, denitrification and plant immobilization
     integer           :: lid_totinput
     integer           :: lid_totstore
     integer           :: lid_cum_closs
     integer , allocatable :: primvarid(:)
     logical , allocatable :: is_aerobic_reac(:)

     integer, allocatable :: lid_plant_minn_no3_pft(:)
     integer, allocatable :: lid_plant_minn_nh4_pft(:)
     integer, allocatable :: lid_plant_minp_pft(:)
     logical, allocatable :: is_cenpool_som(:)
     logical :: debug
     character(len=loc_name_len), allocatable :: varnames(:)
     character(len=loc_name_len), allocatable :: varunits(:)
     character(len=loc_name_len), allocatable :: ompoolnames(:)
     integer                    , allocatable :: vartypes(:)
   contains
     procedure, public  :: Init
     procedure, public  :: hcopy
     procedure, private :: InitPars
     procedure, private :: InitAllocate
     procedure, private :: set_primvar_reac_ids
     procedure, public  :: display_index
  end type ecacnp_bgc_index_type

  contains
  !-----------------------------------------------------------------------
  subroutine list_init(self, name, id, itype)
  implicit none
  type(list_t), pointer :: self
  character(len=*), intent(in) :: name
  integer, intent(inout) :: id
  integer, optional, intent(in) :: itype
  allocate(self)
  nullify(self%next)
  id=id+1;
  write(self%name,'(A)')trim(name)
  self%id=id
  if(present(itype))then
    self%itype=itype
  else
    self%itype=0
  endif
  end subroutine list_init
  !-----------------------------------------------------------------------
  subroutine list_insert(self, name, id, itype)

  implicit none
  type(list_t), pointer :: self
  character(len=*), intent(in) :: name
  integer, intent(inout) :: id
  integer, optional, intent(in) :: itype
  type(list_t), pointer :: next

  allocate(next)
  id=id+1
  write(next%name,'(A)')trim(name)
  next%id=id
  if(present(itype))then
    next%itype=itype
  else
    next%itype=0
  endif
  next%next=> self
  self => next

  end subroutine list_insert
  !-----------------------------------------------------------------------
  subroutine list_free(self)
  implicit none
  type(list_t), pointer :: self
  type(list_t), pointer :: current
  type(list_t), pointer :: elem

  elem => self
  do while(associated(elem))
    current => elem
    elem => current%next
    deallocate(current)
  enddo
  end subroutine list_free
  !-----------------------------------------------------------------------
  function list_next(self)result(next)

  implicit none
  type(list_t), pointer :: self
  type(list_t), pointer :: next

  next => self%next
  end function list_next

  !-------------------------------------------------------------------------------
  subroutine copy_name(num_names, list_name, outnames)

  implicit none
  integer, intent(in) :: num_names
  type(list_t), pointer :: list_name
  character(len=loc_name_len), intent(out) :: outnames(num_names)

  type(list_t), pointer :: next
  integer :: jj
  next => list_name
  do jj = num_names, 1, -1
    write(outnames(jj),'(A)')trim(next%name)
    next=>list_next(next)
  enddo
  end subroutine copy_name
  !-------------------------------------------------------------------------------
  subroutine copy_name_type(num_names, list_name, vartypes)

  implicit none
  integer, intent(in) :: num_names
  type(list_t), pointer :: list_name
  integer, intent(out) :: vartypes(num_names)

  type(list_t), pointer :: next
  integer :: jj
  next => list_name
  do jj = num_names, 1, -1
    vartypes(jj) = next%itype
    next=>list_next(next)
  enddo
  end subroutine copy_name_type
  !-------------------------------------------------------------------------------
  subroutine list_disp(list)
  implicit none
  type(list_t), pointer :: list
  type(list_t), pointer :: next

  next => list
  do while(associated(next))
    write(*,'(A30,X,A,I0)')trim(next%name),'=',next%id
    next=>list_next(next)
  enddo

  end subroutine list_disp
  !-------------------------------------------------------------------------------
  subroutine add_ompool_name(list_name, list_unit, list_pool, prefix, use_c13, use_c14, do_init, vid,uid,pid)
  implicit none
  type(list_t), pointer :: list_name
  type(list_t), pointer :: list_unit
  type(list_t), pointer :: list_pool
  character(len=*), intent(in) :: prefix
  logical, intent(in) :: use_c13, use_c14
  logical, intent(in) :: do_init
  integer, intent(inout) :: vid
  integer, intent(inout) :: uid
  integer, intent(inout) :: pid
  if(do_init)then
    call list_init(list_name, trim(prefix)//'_c',vid, itype=var_state_type)
    call list_init(list_unit, 'mol C m-3',uid)
    call list_init(list_pool, trim(prefix),pid)
  else
    call list_insert(list_name, trim(prefix)//'_c',vid, itype=var_state_type)
    call list_insert(list_unit, 'mol C m-3',uid)
    call list_insert(list_pool, trim(prefix),pid)
  endif
  call list_insert(list_name, trim(prefix)//'_n',vid, itype=var_state_type)
  call list_insert(list_unit, 'mol N m-3',uid)
  call list_insert(list_name, trim(prefix)//'_p',vid, itype=var_state_type)
  call list_insert(list_unit, 'mol P m-3',uid)
  if(use_c13)then
    vid=vid+1;call list_insert(list_name, trim(prefix)//'_c13',vid, itype=var_state_type)
    vid=vid+1;call list_insert(list_unit, 'mol C13 m-3',uid)
  endif
  if(use_c14)then
    vid=vid+1;call list_insert(list_name, trim(prefix)//'_c14',vid, itype=var_state_type)
    vid=vid+1;call list_insert(list_unit, 'mol C14 m-3',uid)
  endif
  end subroutine add_ompool_name

  !-------------------------------------------------------------------------------
  subroutine Init(this, use_c13, use_c14, non_limit, nop_limit, maxpft, batch_mode)
    !
    ! DESCRIPTION:
    ! Initialize ecacnp_bgc type
    ! !USES:
  implicit none
  ! !ARGUMENTS:
  class(ecacnp_bgc_index_type), intent(inout) :: this
  logical, intent(in) :: use_c13
  logical, intent(in) :: use_c14
  logical, intent(in) :: non_limit
  logical, intent(in) :: nop_limit
  integer, optional, intent(in) :: maxpft
  logical, optional, intent(in) :: batch_mode
  ! !LOCAL VARIABLES:
  integer :: maxpft_loc
  logical :: batch_mode_loc
  maxpft_loc = 0; batch_mode_loc=.false.
  this%dom_beg=0; this%dom_end=-1
  if(present(maxpft))maxpft_loc=maxpft
  if(present(batch_mode))batch_mode_loc=batch_mode
  call this%InitPars(maxpft_loc, use_c14, use_c13, non_limit, nop_limit, batch_mode_loc)

  call this%InitAllocate()

  call this%set_primvar_reac_ids(maxpft_loc)

  this%debug = .false.
  end subroutine Init
  !-------------------------------------------------------------------------------

  subroutine InitPars(this, maxpft, use_c14, use_c13, non_limit, nop_limit, batch_mode)
    !
    ! !DESCRIPTION:
    !  describe the layout of the stoichiometric matrix for the reactions
    !           r{1} r{2} r{3} r{4} ... r{n}
    ! s{1}
    ! s{2}
    ! s{3}
    ! s{4}
    ! ...
    ! s{n}
    ! s{n+1}  nonreactive primary variables
    ! s{n+2}
    ! ...
    ! s{m}
    ! s{m+1} diagnostic variables
    ! s{p}
    ! each reaction is associated with a primary species, the secondary species follows after primary species
    ! for the century model, the primary species are seven om pools and nh4, no3 and plant nitrogen
    !
    ! !USES:
    use MathfuncMod   , only : addone, countelm
    use betr_utils    , only : num2str
    use betr_constants, only : betr_string_length_long
    implicit none
    class(ecacnp_bgc_index_type) :: this
    integer, intent(in) :: maxpft
    logical, intent(in) :: use_c13
    logical, intent(in) :: use_c14
    logical, intent(in) :: non_limit
    logical, intent(in) :: nop_limit
    logical, intent(in) :: batch_mode
    ! !LOCAL VARIABLES:
    integer :: itemp
    integer :: ireac   !counter of reactions
    integer :: itemp0, itemp1
    integer :: ielem
    integer :: vid,uid,pid
    integer :: jj
    type(list_t), pointer :: list_name => null()
    type(list_t), pointer :: list_unit => null()
    type(list_t), pointer :: list_pool => null()
    type(list_t), pointer :: list_react=> null()
    character(len=loc_name_len) :: postfix

    this%lid_supp_minp = -1
    this%lid_supp_minn = -1

    itemp = 0; itemp0=0
    ireac = 0
    ielem= 0
    vid = 0;uid=0;pid=0
    this%c_loc = addone(ielem)
    this%n_loc = addone(ielem)
    this%p_loc = addone(ielem)
    if(use_c13)then
      this%c13_loc = addone(ielem)
    endif
    if(use_c14)then
      this%c14_loc = addone(ielem)
    endif
    this%nelms = ielem   !carbon and nitrogen

    !litter group
    this%litr_beg=1
    this%lit1 = addone(itemp); this%lit1_dek_reac = addone(ireac); call list_init(list_react, 'lit1_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'lit1', use_c13, use_c14, do_init=.true., vid=vid,uid=uid,pid=pid)
    this%lit2 = addone(itemp); this%lit2_dek_reac = addone(ireac); call list_insert(list_react, 'lit2_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'lit2', use_c13, use_c14, do_init=.false.,vid=vid,uid=uid,pid=pid)
    this%lit3 = addone(itemp); this%lit3_dek_reac = addone(ireac); call list_insert(list_react, 'lit3_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'lit3', use_c13, use_c14, do_init=.false.,vid=vid,uid=uid,pid=pid)
    this%litr_end=this%litr_beg-1+(this%lit3-this%lit1+1)*this%nelms

    !woody group
    this%wood_beg=this%litr_end+1
    this%cwd  = addone(itemp); this%cwd_dek_reac  = addone(ireac); call list_insert(list_react, 'cwd_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'cwd', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%lwd  = addone(itemp); this%lwd_dek_reac  = addone(ireac); call list_insert(list_react, 'lwd_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'lwd', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%fwd  = addone(itemp); this%fwd_dek_reac  = addone(ireac); call list_insert(list_react, 'fwd_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'fwd', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%wood_end=this%wood_beg-1+(this%fwd-this%cwd+1)*this%nelms

    !microbial biomass group
    this%Bm_beg=this%wood_end+1
    this%som1 = addone(itemp); this%som1_dek_reac = addone(ireac); call list_insert(list_react, 'som1_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'SOM1_MB', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%Bm_end=this%Bm_beg-1+this%nelms

    this%pom_beg=this%Bm_end+1
    this%som2 = addone(itemp); this%som2_dek_reac = addone(ireac); call list_insert(list_react, 'som2_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'SOM2_POM', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%pom_end=this%pom_beg-1+this%nelms

    this%som_beg=this%pom_end+1
    this%som3 = addone(itemp); this%som3_dek_reac = addone(ireac); call list_insert(list_react, 'som3_dek_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'som3', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%som_end=this%som_beg-1+this%nelms

    this%nom_pools = (countelm(this%litr_beg, this%litr_end)+&
       countelm(this%wood_beg,this%wood_end) + &
       countelm(this%som_beg,this%som_end) + &
       countelm(this%Bm_beg,this%Bm_end) + &
       countelm(this%pom_beg,this%pom_end))/this%nelms   !include coarse wood debris

    allocate(this%is_cenpool_som(this%nom_pools)); this%is_cenpool_som(:)=.false.
    this%is_cenpool_som(this%som1)=.true.
    this%is_cenpool_som(this%som2)=.true.
    this%is_cenpool_som(this%som3)=.true.

    itemp               = this%nom_pools*this%nelms

    this%nom_tot_elms    = itemp

    this%lid_minp_secondary = addone(itemp); this%lid_minp_secondary_to_sol_occ_reac=addone(ireac)
    call list_insert(list_react, 'minp_secondary_to_sol_occ_reac', itemp0)
    call list_insert(list_name, 'minp_secondary',vid, itype=var_state_type)
    call list_insert(list_unit, 'mol P m-3',uid)

    this%lid_minp_occlude = addone(itemp);
    call list_insert(list_name, 'minp_occlude',vid, itype=var_state_type)
    call list_insert(list_unit, 'mol P m-3',uid)
    if(maxpft>0)then
      this%lid_plant_minn_nh4_up_reac = addone(ireac); call list_insert(list_react, 'plant_minn_nh4_up_reac', itemp0)
      this%lid_plant_minn_no3_up_reac = addone(ireac); call list_insert(list_react, 'plant_minn_no3_up_reac', itemp0)
      this%lid_plant_minp_up_reac = addone(ireac); call list_insert(list_react, 'plant_minp_up_reac', itemp0)
      this%lid_autr_rt_reac = addone(ireac); call list_insert(list_react, 'autr_rt_reac', itemp0)
    endif
    !non-reactive primary variables
    this%lid_ch4        = addone(itemp);call list_insert(list_name, 'ch4',vid, itype=var_state_type); call list_insert(list_unit, 'mol C m-3',uid)
    this%lid_ar         = addone(itemp);call list_insert(list_name, 'ar',vid, itype=var_state_type); call list_insert(list_unit, 'mol m-3',uid)

    !second primary variables
    this%lid_o2         = addone(itemp);call list_insert(list_name, 'o2',vid, itype=var_state_type); call list_insert(list_unit, 'mol m-3',uid)

    this%lid_co2        = addone(itemp);call list_insert(list_name, 'co2',vid, itype=var_state_type);call list_insert(list_unit,'mol m-3',uid)

    if(use_c13)then
      this%lid_c13_co2 = addone(itemp);call list_insert(list_name, 'co2_c13',vid, itype=var_state_type); call list_insert(list_unit, 'mol C13 m-3',uid)
    endif
    if(use_c14)then
      this%lid_c14_co2 = addone(itemp);call list_insert(list_name, 'co2_c14',vid, itype=var_state_type); call list_insert(list_unit, 'mol C14 m-3',uid)
    endif

    this%lid_n2o        = addone(itemp);call list_insert(list_name, 'n2o',vid, itype=var_state_type); call list_insert(list_unit,'mol N2O m-3',uid)
    this%lid_n2         = addone(itemp);call list_insert(list_name, 'n2',vid, itype=var_state_type); call list_insert(list_unit, 'mol N2 m-3',uid)

    this%nprimvars      = itemp

    if(non_limit)then
      !when N is unlimited, nh4 and no3 are not primary variables
      if(nop_limit .or. spinup_state /= 0)then
        this%lid_nh4        = addone(itemp)
        this%lid_no3        = addone(itemp)
        call list_insert(list_name, 'nh4',vid, itype=var_state_type); call list_insert(list_unit, 'mol N m-3',uid)
        call list_insert(list_name, 'no3',vid, itype=var_state_type); call list_insert(list_unit, 'mol N m-3',uid)
        this%lid_minp_soluble=addone(itemp);
        call list_insert(list_name, 'minp_soluble',vid, itype=var_state_type); call list_insert(list_unit, 'mol P m-3',uid)
      else
        this%lid_minp_soluble=addone(itemp);
        call list_insert(list_name, 'minp_soluble',vid, itype=var_state_type); call list_insert(list_unit, 'mol P m-3',uid)
        this%nprimvars      = this%nprimvars + 1     !primary state variables 14 + 6
        this%lid_nh4        = addone(itemp)
        this%lid_no3        = addone(itemp)
        call list_insert(list_name, 'nh4',vid, itype=var_state_type); call list_insert(list_unit, 'mol N m-3',uid)
        call list_insert(list_name, 'no3',vid, itype=var_state_type); call list_insert(list_unit, 'mol N m-3',uid)
      endif
    else
      this%lid_nh4        = addone(itemp)
      this%lid_no3        = addone(itemp)
      call list_insert(list_name, 'nh4',vid, itype=var_state_type); call list_insert(list_unit, 'mol N m-3',uid)
      call list_insert(list_name, 'no3',vid, itype=var_state_type); call list_insert(list_unit, 'mol N m-3',uid)
      this%nprimvars = this%nprimvars + 2
      this%lid_minp_soluble=addone(itemp);
      call list_insert(list_name, 'minp_soluble',vid, itype=var_state_type); call list_insert(list_unit, 'mol P m-3',uid)
      if(.not. (nop_limit .or. spinup_state /= 0))then
        this%nprimvars = this%nprimvars + 1
      endif
    endif
    this%lid_nh4_nit_reac = addone(ireac); call list_insert(list_react, 'nh4_nit_reac', itemp0)       !this is also used to indicate the nitrification reaction
    this%lid_no3_den_reac = addone(ireac); call list_insert(list_react, 'no3_den_reac', itemp0)       !this is also used to indicate the denitrification reaction

    this%lid_minp_soluble_to_secp_reac = addone(ireac); call list_insert(list_react, 'minp_soluble_to_secp_reac', itemp0)

    if(non_limit)then
      this%lid_supp_minn=addone(itemp)
      call list_insert(list_name, 'supp_minn',vid, itype=var_state_type); call list_insert(list_unit, 'mol N m-3',uid)
    endif
    if(spinup_state /= 0 .or. nop_limit)then
      this%lid_supp_minp=addone(itemp)
      call list_insert(list_name, 'supp_minp_soluble',vid, itype=var_state_type); call list_insert(list_unit, 'mol P m-3',uid)
    endif
    !diagnostic variables
    if(maxpft>0)then
      this%lid_plant_minn_nh4 = addone(itemp)  !this is used to indicate plant mineral nitrogen uptake
      call list_insert(list_name, 'plant_minn_nh4',vid, itype=var_flux_type); call list_insert(list_unit, 'mol P m-3 s-1',uid)

      this%lid_plant_minn_no3 = addone(itemp)  !this is used to indicate plant mineral nitrogen uptake
      call list_insert(list_name, 'plant_minn_no3',vid, itype=var_flux_type); call list_insert(list_unit, 'mol N mol-3 s-1',uid)

      this%lid_plant_minp = addone(itemp)
      call list_insert(list_name, 'plant_minp',vid, itype=var_flux_type); call list_insert(list_unit, 'mol P m-3 s-1',uid)

      this%lid_autr_rt      = addone(itemp)           !this is used to indicate plant autotrophic root respiration
      call list_insert(list_name, 'autr_rt',vid, itype=var_flux_type); call list_insert(list_unit,'mol C m-3 s-1',uid)
    endif
    this%lid_n2o_nit  = addone(itemp);
    call list_insert(list_name, 'n2o_nit',vid, itype=var_flux_type); call list_insert(list_unit, 'mol N2O m-3 s-1',uid)

    this%lid_co2_hr   = addone(itemp);
    call list_insert(list_name, 'co2_hr',vid, itype=var_flux_type); call list_insert(list_unit,'mol C m-3 s-1',uid)

    this%lid_no3_den  = addone(itemp);
    call list_insert(list_name, 'no3_den',vid, itype=var_flux_type); call list_insert(list_unit, 'mol N m-3 s-1',uid)

    this%lid_minn_nh4_immob = addone(itemp);
    call list_insert(list_name, 'minn_nh4_immob',vid, itype=var_flux_type); call list_insert(list_unit, 'mol N m-3 s-1',uid)

    this%lid_minn_no3_immob = addone(itemp);
    call list_insert(list_name, 'minn_no3_immob',vid, itype=var_flux_type); call list_insert(list_unit, 'mol N m-3 s-1',uid)

    this%lid_minp_immob = addone(itemp);
    call list_insert(list_name, 'minp_immob',vid, itype=var_flux_type); call list_insert(list_unit, 'mol P m-3 s-1',uid)

    this%lid_nh4_nit        = addone(itemp);
    call list_insert(list_name, 'nh4_nit',vid, itype=var_flux_type); call list_insert(list_unit,'mol N m-3 s-1',uid)

    !aerechyma transport
    this%lid_o2_paere   = addone(itemp);
    this%lid_o2_aren_reac  = addone(ireac); call list_insert(list_react, 'o2_aren_reac', itemp0)

    call list_insert(list_name, 'o2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)
    if ( spinup_state == 0 ) then
       this%lid_ar_paere   = addone(itemp);
       this%lid_ar_aren_reac  = addone(ireac);  call list_insert(list_react, 'ar_aren_reac', itemp0)   !
       call list_insert(list_name, 'ar_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)

       this%lid_n2_paere   = addone(itemp);
       this%lid_n2_aren_reac  = addone(ireac); call list_insert(list_react, 'n2_aren_reac', itemp0)      !
       call list_insert(list_name, 'n2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol N2 m-3 s-1',uid)

       this%lid_co2_paere  = addone(itemp);
       this%lid_co2_aren_reac = addone(ireac); call list_insert(list_react, 'co2_aren_reac', itemp0)      !
       call list_insert(list_name, 'co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol C m-3 s-1',uid)

       if(use_c13)then
         this%lid_c13_co2_paere  = addone(itemp)
         this%lid_c13_co2_aren_reac = addone(ireac); call list_insert(list_react, 'c13_co2_aren_reac', itemp0)      !
         call list_insert(list_name, 'c13_co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol C13 m-3 s-1',uid)
       endif

       if(use_c14)then
         this%lid_c14_co2_paere  = addone(itemp);
         this%lid_c14_co2_aren_reac = addone(ireac); call list_insert(list_react, 'c14_co2_aren_reac', itemp0)      !
         call list_insert(list_name, 'c14_co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol C14 m-3 s-1',uid)
       endif

       this%lid_ch4_paere  = addone(itemp);
       this%lid_ch4_aren_reac = addone(ireac); call list_insert(list_react, 'ch4_aren_reac', itemp0)   !
       call list_insert(list_name, 'ch4_paere',vid, itype=var_flux_type); call list_insert(list_unit, 'mol C m-3 s-1',uid)

       this%lid_n2o_paere  = addone(itemp);
       this%lid_n2o_aren_reac = addone(ireac); call list_insert(list_react, 'n2o_aren_reac', itemp0)   !
       call list_insert(list_name, 'n2o_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol  m-3 s-1',uid)
    endif

    if(maxpft>0)then
      allocate(this%lid_plant_minn_no3_pft(maxpft));
      allocate(this%lid_plant_minn_nh4_pft(maxpft));
      allocate(this%lid_plant_minp_pft(maxpft));
      do jj = 1, maxpft
        this%lid_plant_minn_no3_pft(jj) = addone(itemp)
        postfix = num2str(jj,'(I2.2)')
        call list_insert(list_name, 'plant_minn_no3_pft_'//trim(postfix),vid, itype=var_flux_type)
        call list_insert(list_unit, 'mol N m-3 s-1',uid)
        this%lid_plant_minn_nh4_pft(jj) = addone(itemp)
        call list_insert(list_name, 'plant_minn_nh4_pft_'//trim(postfix),vid, itype=var_flux_type)
        call list_insert(list_unit,'mol N m-3 s-1',uid)
        this%lid_plant_minp_pft(jj)     = addone(itemp)
        call list_insert(list_name, 'plant_minp_pft_'//trim(postfix),vid, itype=var_flux_type)
        call list_insert(list_unit,'mol P m-3 s-1',uid)
      enddo
    endif
    if(batch_mode)then
      this%lid_totinput  = addone(itemp);call list_insert(list_name, 'c_tot_input',vid, itype=var_state_type); call list_insert(list_unit,'mol m-3',uid)
      this%lid_totstore  = addone(itemp);call list_insert(list_name, 'c_tot_store',vid, itype=var_state_type); call list_insert(list_unit,'mol m-3',uid)
      this%lid_cum_closs  = addone(itemp);call list_insert(list_name, 'c_loss',vid, itype=var_state_type); call list_insert(list_unit,'mol m-3',uid)
    endif
    this%nstvars          = itemp          !totally 14+32 state variables
    this%nreactions = ireac            !seven decomposition pathways plus root auto respiration, nitrification, denitrification and plant immobilization

    allocate(this%primvarid(ireac)); this%primvarid(:) = -1
    allocate(this%is_aerobic_reac(ireac)); this%is_aerobic_reac(:)=.false.

    allocate(this%vartypes(this%nstvars))
    allocate(this%varnames(this%nstvars))
    allocate(this%varunits(this%nstvars))
    allocate(this%ompoolnames(this%nom_pools))

    call copy_name(this%nstvars, list_name, this%varnames(1:this%nstvars))
    call copy_name(this%nstvars, list_unit, this%varunits(1:this%nstvars))
    call copy_name(this%nom_pools, list_pool, this%ompoolnames(1:this%nom_pools))
    call copy_name_type(this%nstvars, list_name, this%vartypes(1:this%nstvars))
!    call list_disp(list_name);call list_disp(list_pool);call list_disp(list_unit)
!    call list_disp(list_react)

    call list_free(list_name)
    call list_free(list_pool)
    call list_free(list_unit)
    call list_free(list_react)
  end subroutine InitPars
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this)
    !
    ! !DESCRIPTION:
    ! memory allocation for the data type specified by this
    !
  implicit none
    ! !ARGUMENTS:
  class(ecacnp_bgc_index_type), intent(inout) :: this


  end subroutine InitAllocate
  !-------------------------------------------------------------------------------
  subroutine set_primvar_reac_ids(this,maxpft_loc)

  implicit none
  class(ecacnp_bgc_index_type), intent(inout)  :: this
  integer, intent(in) :: maxpft_loc
  integer :: reac

  associate(                                      &
    lit1  => this%lit1                       ,    &
    lit2  => this%lit2                       ,    &
    lit3  => this%lit3                       ,    &
    cwd  => this%cwd                         ,    &
    lwd  => this%lwd                         ,    &
    fwd  => this%fwd                         ,    &
    som1  => this%som1                       ,    &
    som2  => this%som2                       ,    &
    som3  => this%som3                       ,    &
    nelms => this%nelms                      ,    &
    c_loc => this%c_loc                      ,    &
    n_loc => this%n_loc                      ,    &
    p_loc => this%p_loc                      ,    &
    lid_nh4=> this%lid_nh4                   ,    &
    lid_no3 => this%lid_no3                  ,    &
    lid_o2 => this%lid_o2                    ,    &
    lid_minp_soluble=> this%lid_minp_soluble ,    &
    lid_minp_secondary => this%lid_minp_secondary &
  )
  !reaction1, lit1 -> s1
  reac=this%lit1_dek_reac;     this%primvarid(reac) = (lit1-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.
  !print*,'1',reac,size(this%primvarid)

  !reaction 2, lit2 -> s1
  reac =this%lit2_dek_reac;   this%primvarid(reac) = (lit2-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.
  !print*,'2',reac,size(this%primvarid)

  !reaction 3, lit3->s2
  reac =this%lit3_dek_reac; this%primvarid(reac) = (lit3-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.
  !print*,'3',reac,size(this%primvarid)

  !reaction 4, SOM1 -> f1*SOM2 + f2*SOm3
  reac =this%som1_dek_reac; this%primvarid(reac) = (som1-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.
  !print*,'4',reac,size(this%primvarid)

  !reaction 5, som2->som1, som3
  reac =this%som2_dek_reac;  this%primvarid(reac) = (som2-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.
  !print*,'5',reac,size(this%primvarid)

  !reaction 6, s3-> s1
  reac = this%som3_dek_reac;  this%primvarid(reac) = (som3-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.
  !print*,'6',reac,size(this%primvarid)

  !reaction 7, cwd -> lit1 + lit2
  reac = this%cwd_dek_reac; this%primvarid(reac) = (cwd-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.
  !print*,'7',reac,size(this%primvarid)

  reac = this%fwd_dek_reac; this%primvarid(reac) = (fwd-1)*nelms+c_loc
  !print*,'8',reac,size(this%primvarid)

  reac = this%lwd_dek_reac; this%primvarid(reac) = (lwd-1)*nelms+c_loc
  !print*,'9',reac,size(this%primvarid)

  !reaction 8, nitrification
  reac = this%lid_nh4_nit_reac; this%primvarid(reac) = this%lid_nh4
  !x is_aerobic_reac(reac) = .true.
  !print*,'10',reac,size(this%primvarid)

  !reaction 9, denitrification
  reac = this%lid_no3_den_reac; this%primvarid(reac) = lid_no3
  !print*,'11',reac,size(this%primvarid)

  !reaction 10, inorganic P non-equilibrium adsorption
  !P_solution -> p_secondary
  reac = this%lid_minp_soluble_to_secp_reac; this%primvarid(reac) = lid_minp_soluble
  !print*,'12',reac,size(this%primvarid)

  !reaction 11, inorganic P non-equilibrium desorption
  ! p_secondary -> P_solution
  reac = this%lid_minp_secondary_to_sol_occ_reac; this%primvarid(reac) = lid_minp_secondary
  !print*,'13',reac,size(this%primvarid)

  if(maxpft_loc>0)then
    !reaction 12, plant mineral nitrogen nh4 uptake
    reac = this%lid_plant_minn_nh4_up_reac; this%primvarid(reac) = lid_nh4
    !print*,'14',reac,size(this%primvarid)

    !reaction 13, plant mineral nitrogen no3 uptake
    reac = this%lid_plant_minn_no3_up_reac; this%primvarid(reac) = lid_no3
    !print*,'15',reac,size(this%primvarid)

    !reaction 14, plant mineral phosphorus uptake
    reac = this%lid_plant_minp_up_reac; this%primvarid(reac) = lid_minp_soluble
    !print*,'16',reac,size(this%primvarid)
  endif
  !reaction 15, autotrophic respiration, ar + o2 -> co2
  !x reac = lid_ar_rt_reac; is_aerobic_reac(reac) = .true.

  !reaction 15, o2 transport through arenchyma
  reac = this%lid_o2_aren_reac; this%primvarid(reac) = lid_o2
  !print*,'17',reac,size(this%primvarid)
  end associate

  end subroutine set_primvar_reac_ids

!-------------------------------------------------------------------------------
  subroutine display_index(this)

  implicit none
  class(ecacnp_bgc_index_type) :: this
  print*,'lit1=', this%lit1
  print*,'lit2=', this%lit2
  print*,'lit3=', this%lit3
  print*,'som1=', this%som1
  print*,'som2=', this%som2
  print*,'som3=', this%som3
  print*,'cwd=',  this%cwd
  print*,'lwd=',  this%lwd
  print*,'fwd=',  this%fwd
  print*,'litrgrp=', this%litr_beg, this%litr_end  !litr group
  print*,'woodgrp=', this%wood_beg, this%wood_end  !wood group
  print*,'somgrp =', this%som_beg,  this%som_end   !som group
  print*,'domgrp =', this%dom_beg,  this%dom_end   !dom group
  print*,'Bmgrp  =', this%Bm_beg,   this%Bm_end   !dom group
  print*,'pomgrp =', this%pom_beg,  this%pom_end
  print*,'supminp=', this%lid_supp_minp                          !supplementary mineral P for spinup purpose
  print*,'supminn=', this%lid_supp_minn
  print*,'lid_nh4=', this%lid_nh4,this%lid_nh4_nit_reac              !local position of nh4 in the state variable vector
  print*,'lid_no3=', this%lid_no3,this%lid_no3_den_reac              !local position of no3 in the state variable vector
  print*,'pltmnh4=', this%lid_plant_minn_nh4, this%lid_plant_minn_nh4_up_reac !local position of plant uptake of mineral nitrogen NH4 in the state variable vector
  print*,'pltmno3=', this%lid_plant_minn_no3, this%lid_plant_minn_no3_up_reac !
  print*,'pltmmnp=', this%lid_plant_minp, this%lid_plant_minp_up_reac !local position of plant uptake of mineral P in the state variable vector
  print*,'mnpsol =', this%lid_minp_soluble, this%lid_minp_soluble_to_secp_reac    !conversation of adsorbed into secondary phase
  print*,'mnp2nd =', this%lid_minp_secondary,this%lid_minp_secondary_to_sol_occ_reac   !local position of secondary P in the state variable vector
  print*,'mnpocl =', this%lid_minp_occlude      !local position of occluded P in the state variable vector
  print*,'atr_rt =', this%lid_autr_rt, this%lid_autr_rt_reac             !root autotrophic respiration
  print*,'ar     =', this%lid_ar, this%lid_ar_aren_reac               !local position of ar in the state variable vector
  print*,'ch4    =', this%lid_ch4, this%lid_ch4_aren_reac             !nonreactive primary variables
  print*,' o2    =', this%lid_o2,  this%lid_o2_aren_reac              !local position of o2 in the state variable vector
  print*,' co2   =', this%lid_co2, this%lid_co2_aren_reac             !local position of co2 in the state variable vector
  print*,'  n2   =', this%lid_n2,  this%lid_n2_aren_reac
  print*,'  n2o  =', this%lid_n2o, this%lid_n2o_aren_reac
  print*,'n2o_nit=', this%lid_n2o_nit                            !n2o production from nitrification, used to for mass balance book keeping
  print*,'co2_hr =', this%lid_co2_hr                             !co2 production from heterotrophic respiration
  print*,'c13_co2=', this%lid_c13_co2, this%lid_c13_co2_aren_reac
  print*,'c14_co2=', this%lid_c14_co2, this%lid_c14_co2_aren_reac
  print*,'no3_den=', this%lid_no3_den                            !no3 consumption due to denitrification
  print*,'nh4imob=', this%lid_minn_nh4_immob                     !net mineral NH4 immobilization for decomposition
  print*,'no3imob=', this%lid_minn_no3_immob                     !net mineral NO3 immobilization for decomposition
  print*,'nh4_nit=', this%lid_nh4_nit
  print*,'mnp2nd =', this%lid_minp_secondary_trc
  print*,'mnpocl =', this%lid_minp_occlude_trc
  print*,'mnpimob=', this%lid_minp_immob                         !net P immobilization by aerobic decomposer
  print*,'arpare =', this%lid_ar_paere
  print*,'n2pare =', this%lid_n2_paere
  print*,'o2pare =', this%lid_o2_paere
  print*,'co2pare=', this%lid_co2_paere
  print*,'c13o2pr=', this%lid_c13_co2_paere
  print*,'c14o2pr=', this%lid_c14_co2_paere
  print*,'ch4pare=', this%lid_ch4_paere
  print*,'n2opare=', this%lid_n2o_paere
  print*,'pno3pft=', this%lid_plant_minn_no3_pft(:)
  print*,'pnh4pft=', this%lid_plant_minn_nh4_pft(:)
  print*,'minppft=', this%lid_plant_minp_pft(:)
  end subroutine display_index
!-------------------------------------------------------------------------------
  subroutine hcopy(this, that)
  implicit none
  class(ecacnp_bgc_index_type) :: this
  class(ecacnp_bgc_index_type), intent(in) :: that

  integer :: jj
  integer :: maxpft

  this%nom_pools = that%nom_pools                              !not include coarse wood debris

  this%nom_tot_elms = that%nom_tot_elms
  this%lit1 = that%lit1; this%lit1_dek_reac = that%lit1_dek_reac
  this%lit2 = that%lit2; this%lit2_dek_reac = that%lit2_dek_reac
  this%lit3 = that%lit3; this%lit3_dek_reac = that%lit3_dek_reac
  this%som1 = that%som1; this%som1_dek_reac = that%som1_dek_reac
  this%som2 = that%som2; this%som2_dek_reac = that%som2_dek_reac
  this%som3 = that%som3; this%som3_dek_reac = that%som3_dek_reac
  this%cwd  = that%cwd;  this%cwd_dek_reac  = that%cwd_dek_reac
  this%lwd  = that%lwd;  this%lwd_dek_reac = that%lwd_dek_reac
  this%fwd  = that%fwd;  this%fwd_dek_reac = that%fwd_dek_reac
  this%litr_beg = that%litr_beg; this%litr_end = that%litr_end  !litr group
  this%wood_beg = that%wood_beg; this%wood_end = that%wood_end  !wood group
  this%som_beg  = that%som_beg;  this%som_end = that%som_end   !som group
  this%dom_beg  = that%dom_beg;  this%dom_end = that%dom_end   !dom group
  this%Bm_beg   = that%Bm_beg;  this%Bm_end = that%Bm_end   !dom group
  this%pom_beg  = that%pom_beg; this%pom_end= that%pom_end
  this%c_loc = that%c_loc
  this%n_loc = that%n_loc
  this%p_loc = that%p_loc
  this%c13_loc = that%c13_loc
  this%c14_loc = that%c14_loc
  this%lid_supp_minp = that%lid_supp_minp                          !supplementary mineral P for spinup purpose
  this%lid_supp_minn = that%lid_supp_minn
  this%nelms = that%nelms                                  !number of chemical elements in an om pool
                                                                 !reactive primary variables
  this%lid_nh4 = that%lid_nh4; this%lid_nh4_nit_reac = that%lid_nh4_nit_reac             !local position of nh4 in the state variable vector
  this%lid_no3 = that%lid_no3; this%lid_no3_den_reac = that%lid_no3_den_reac              !local position of no3 in the state variable vector
  this%lid_plant_minn_nh4 = that%lid_plant_minn_nh4
  this%lid_plant_minn_nh4_up_reac = that%lid_plant_minn_nh4_up_reac !local position of plant uptake of mineral nitrogen NH4 in the state variable vector
  this%lid_plant_minn_no3 = that%lid_plant_minn_no3
  this%lid_plant_minn_no3_up_reac = that%lid_plant_minn_no3_up_reac !
  this%lid_plant_minp = that%lid_plant_minp
  this%lid_plant_minp_up_reac = that%lid_plant_minp_up_reac !local position of plant uptake of mineral P in the state variable vector
  this%lid_minp_soluble = that%lid_minp_soluble
  this%lid_minp_soluble_to_secp_reac = that%lid_minp_soluble_to_secp_reac    !conversation of adsorbed into secondary phase
  this%lid_minp_secondary = that%lid_minp_secondary
  this%lid_minp_secondary_to_sol_occ_reac = that%lid_minp_secondary_to_sol_occ_reac   !local position of secondary P in the state variable vector

  this%lid_minp_occlude = that%lid_minp_occlude      !local position of occluded P in the state variable vector

  this%lid_autr_rt = that%lid_autr_rt; this%lid_autr_rt_reac = that%lid_autr_rt_reac             !root autotrophic respiration

                                                                 !non reactive primary variables
  this%lid_ar = that%lid_ar; this%lid_ar_aren_reac = that%lid_ar_aren_reac               !local position of ar in the state variable vector
  this%lid_ch4= that%lid_ch4; this%lid_ch4_aren_reac = that%lid_ch4_aren_reac             !nonreactive primary variables

                                                                 !secondary variables
  this%lid_o2 = that%lid_o2;  this%lid_o2_aren_reac = that%lid_o2_aren_reac              !local position of o2 in the state variable vector
  this%lid_co2 = that%lid_co2; this%lid_co2_aren_reac = that%lid_co2_aren_reac             !local position of co2 in the state variable vector
  this%lid_n2 = that%lid_n2;  this%lid_n2_aren_reac = that%lid_n2_aren_reac
  this%lid_n2o= that%lid_n2o; this%lid_n2o_aren_reac = that%lid_n2o_aren_reac
                                                                 !diagnostic variables
  this%lid_n2o_nit = that%lid_n2o_nit                            !n2o production from nitrification, used to for mass balance book keeping
  this%lid_co2_hr = that%lid_co2_hr                             !co2 production from heterotrophic respiration
  this%lid_c13_co2 = that%lid_c13_co2;
  this%lid_c13_co2_aren_reac = that%lid_c13_co2_aren_reac
  this%lid_c14_co2 = that%lid_c14_co2;
  this%lid_c14_co2_aren_reac = that%lid_c14_co2_aren_reac
  this%lid_no3_den = that%lid_no3_den                            !no3 consumption due to denitrification
  this%lid_minn_nh4_immob = that%lid_minn_nh4_immob                     !net mineral NH4 immobilization for decomposition
  this%lid_minn_no3_immob = that%lid_minn_no3_immob                     !net mineral NO3 immobilization for decomposition
  this%lid_nh4_nit = that%lid_nh4_nit
  this%lid_minp_secondary_trc = that%lid_minp_secondary_trc
  this%lid_minp_occlude_trc = that%lid_minp_occlude_trc
                                                                 !aerechyma transport, diagnostic efflux
  this%lid_minp_immob = that%lid_minp_immob                         !net P immobilization by aerobic decomposer

  this%lid_ar_paere = that%lid_ar_paere
  this%lid_n2_paere = that%lid_n2_paere
  this%lid_o2_paere = that%lid_o2_paere
  this%lid_co2_paere = that%lid_co2_paere
  this%lid_c13_co2_paere = that%lid_c13_co2_paere
  this%lid_c14_co2_paere = that%lid_c14_co2_paere
  this%lid_ch4_paere = that%lid_ch4_paere
  this%lid_n2o_paere = that%lid_n2o_paere
  this%nprimvars = that%nprimvars                              !total number of primary variables
  this%nstvars = that%nstvars                                !number of equations for the state variabile vector
  this%nreactions = that%nreactions                             !seven decomposition pathways plus nitrification, denitrification and plant immobilization
  this%lid_totinput = that%lid_totinput
  this%lid_totstore = that%lid_totstore
  this%lid_cum_closs = that%lid_cum_closs

  this%debug = that%debug
  if(allocated(that%lid_plant_minn_no3_pft))then
    maxpft=size(that%lid_plant_minn_no3_pft)
  else
    maxpft=0
  endif
  if(maxpft>0)then
     allocate(this%lid_plant_minn_no3_pft(maxpft));
     allocate(this%lid_plant_minn_nh4_pft(maxpft));
     allocate(this%lid_plant_minp_pft(maxpft));
     do jj = 1, maxpft
       this%lid_plant_minn_no3_pft(jj) = that%lid_plant_minn_no3_pft(jj);
       this%lid_plant_minn_nh4_pft(jj) = that%lid_plant_minn_nh4_pft(jj);
       this%lid_plant_minp_pft(jj)     = that%lid_plant_minp_pft(jj);
     enddo
  endif

  if(this%nreactions>0)then
    allocate(this%primvarid(this%nreactions))
    allocate(this%is_aerobic_reac(this%nreactions));
    do jj = 1, this%nreactions
      this%primvarid(jj) = that%primvarid(jj)
      this%is_aerobic_reac(jj) = that%is_aerobic_reac(jj)
    enddo
  endif
  if(this%nstvars>0)then
    allocate(this%vartypes(this%nstvars))
    allocate(this%varnames(this%nstvars))
    allocate(this%varunits(this%nstvars))
    do jj = 1, this%nstvars
      this%vartypes(jj) = that%vartypes(jj)
      this%varnames(jj) = that%varnames(jj)
      this%varunits(jj) = that%varunits(jj)
    enddo
  endif
  if(this%nom_pools>0)then
    allocate(this%ompoolnames(this%nom_pools))
    allocate(this%is_cenpool_som(this%nom_pools))
    do jj = 1, this%nom_pools
      this%ompoolnames(jj) = that%ompoolnames(jj)
      this%is_cenpool_som(jj)=that%is_cenpool_som(jj)
    enddo
  endif
  end subroutine hcopy
end module ecacnpBGCIndexType
