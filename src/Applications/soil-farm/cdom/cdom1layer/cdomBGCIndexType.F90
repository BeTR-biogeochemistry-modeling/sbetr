module cdomBGCIndexType

  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_ctrl      , only : spinup_state => betr_spinup_state
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

  type, public, extends(gbgc_index_type) :: cdom_bgc_index_type
     integer           :: nom_pools                              !not include coarse wood debris

     integer           :: nom_tot_elms
     integer           :: lmet, lmet_dek_reac
     integer           :: lcel, lcel_dek_reac
     integer           :: llig, llig_dek_reac
     integer           :: mic, mic_dek_reac
     integer           :: pom, pom_dek_reac
     integer           :: humus, humus_dek_reac
     integer           :: cwd,  cwd_dek_reac
     integer           :: lwd,  lwd_dek_reac
     integer           :: fwd,  fwd_dek_reac
     integer           :: dom,  dom_dek_reac
     integer           :: litr_beg, litr_end          !litr group
     integer           :: wood_beg, wood_end          !wood group
     integer           :: pom_beg,  pom_end           !som group
     integer           :: dom_beg,  dom_end           !dom group
     integer           :: micbiom_beg,  micbiom_end   !dom group
     integer           :: humus_beg, humus_end
     integer           :: lid_supp_minp
     integer           :: lid_supp_minn
     integer           :: nlit, nwood
     integer           :: c_loc
     integer           :: n_loc
     integer           :: p_loc
     integer           :: e_loc
     integer           :: c13_loc = 0
     integer           :: c14_loc = 0
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

     integer , pointer :: primvarid(:)   => null()
     logical , pointer :: is_aerobic_reac(:)=> null()

     integer, pointer :: lid_plant_minn_no3_pft(:)=> null()
     integer, pointer :: lid_plant_minn_nh4_pft(:)=> null()
     integer, pointer :: lid_plant_minp_pft(:)=> null()
     logical, pointer :: is_ompool_som(:) => null()
     logical, pointer :: is_dom_pool(:) => null()
     logical :: debug
     character(len=loc_name_len), allocatable :: varnames(:)
     character(len=loc_name_len), allocatable :: varunits(:)
     character(len=loc_name_len), allocatable :: ompoolnames(:)
     integer                    , allocatable :: vartypes(:)
   contains
     procedure, public  :: Init
     procedure, private :: InitPars
     procedure, private :: InitAllocate
     procedure, private :: set_primvar_reac_ids
  end type cdom_bgc_index_type

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
  subroutine add_ompool_name(list_name, list_unit, list_pool,prefix, use_c13, use_c14, do_init, vid,uid,pid)
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
    call list_insert(list_name, trim(prefix)//'_c',vid)
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
  subroutine Init(this, use_c13, use_c14, non_limit, nop_limit, maxpft)
    !
    ! DESCRIPTION:
    ! Initialize centurybgc type
    ! !USES:
  implicit none
  ! !ARGUMENTS:
  class(cdom_bgc_index_type), intent(inout) :: this
  logical, intent(in) :: use_c13
  logical, intent(in) :: use_c14
  logical, intent(in) :: non_limit
  logical, intent(in) :: nop_limit
  integer, optional, intent(in) :: maxpft

  ! !LOCAL VARIABLES:
  integer :: maxpft_loc

  maxpft_loc = 0
  this%dom_beg=0; this%dom_end=-1
  if(present(maxpft))maxpft_loc=maxpft
  call this%InitPars(maxpft_loc, use_c14, use_c13, non_limit, nop_limit)

  call this%InitAllocate()

  call this%set_primvar_reac_ids(maxpft_loc)

  this%debug = .false.
  end subroutine Init
  !-------------------------------------------------------------------------------

  subroutine InitPars(this, maxpft, use_c14, use_c13, non_limit, nop_limit)
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
    class(cdom_bgc_index_type) :: this
    integer, intent(in) :: maxpft
    logical, intent(in) :: use_c13
    logical, intent(in) :: use_c14
    logical, intent(in) :: non_limit
    logical, intent(in) :: nop_limit
    ! !LOCAL VARIABLES:
    integer :: itemp
    integer :: ireac   !counter of reactions
    integer :: itemp1
    integer :: ielem
    integer :: vid,uid,pid
    integer :: jj
    type(list_t), pointer :: list_name => null()
    type(list_t), pointer :: list_unit => null()
    type(list_t), pointer :: list_pool => null()
    character(len=loc_name_len) :: postfix

    this%lid_supp_minp = -1
    this%lid_supp_minn = -1
    itemp = 0
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
    this%e_loc = ielem + 1

    !litter group
    this%litr_beg=1; this%nlit=0
    this%lmet = addone(itemp); this%lmet_dek_reac = addone(ireac); this%nlit=this%nlit+1
    call add_ompool_name(list_name, list_unit, list_pool,'Lmet', use_c13, use_c14, do_init=.true., vid=vid,uid=uid,pid=pid)
    this%lcel = addone(itemp); this%lcel_dek_reac = addone(ireac); this%nlit=this%nlit+1
    call add_ompool_name(list_name, list_unit, list_pool,'Lcel', use_c13, use_c14, do_init=.false.,vid=vid,uid=uid,pid=pid)
    this%llig = addone(itemp); this%llig_dek_reac = addone(ireac); this%nlit=this%nlit+1
    call add_ompool_name(list_name, list_unit, list_pool,'Llig', use_c13, use_c14, do_init=.false.,vid=vid,uid=uid,pid=pid)
    this%litr_end=this%litr_beg-1+this%nelms*this%nlit

    !woody group
    this%wood_beg=this%litr_end+1; this%nwood=0
    this%cwd  = addone(itemp); this%cwd_dek_reac  = addone(ireac); this%nwood=this%nwood+1
    call add_ompool_name(list_name, list_unit, list_pool,'CWD', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%lwd  = addone(itemp); this%lwd_dek_reac  = addone(ireac); this%nwood=this%nwood+1
    call add_ompool_name(list_name, list_unit, list_pool,'LWD', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%fwd  = addone(itemp); this%fwd_dek_reac  = addone(ireac); this%nwood=this%nwood+1
    call add_ompool_name(list_name, list_unit, list_pool,'FWD', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%wood_end=this%wood_beg-1+this%nelms*this%nwood

    !microbial biomass group
    this%micbiom_beg=this%wood_end+1
    this%mic = addone(itemp); this%mic_dek_reac = addone(ireac)
    call add_ompool_name(list_name, list_unit, list_pool,'Micbiom', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%micbiom_end=this%micbiom_beg-1+this%nelms

    this%pom_beg=this%micbiom_end+1
    this%pom = addone(itemp); this%pom_dek_reac = addone(ireac)
    call add_ompool_name(list_name, list_unit, list_pool,'POM', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%pom_end=this%pom_beg-1+this%nelms

    this%humus_beg=this%pom_end+1
    this%humus = addone(itemp); this%humus_dek_reac = addone(ireac)
    call add_ompool_name(list_name, list_unit, list_pool,'Humus', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%humus_end=this%humus_beg-1+this%nelms

    !dom group
    this%dom_beg=this%humus_end+1;
    this%dom = addone(itemp); this%dom_dek_reac = addone(ireac)
    call add_ompool_name(list_name, list_unit, list_pool,'DOM', use_c13, use_c14, do_init=.false., vid=vid, uid=uid, pid=pid)
    call list_insert(list_name, 'DOM_e',vid)
    call list_insert(list_unit, 'mol e m-3',uid)
    this%dom_end=this%dom_beg-1+this%nelms + 1

    this%nom_pools = (countelm(this%litr_beg, this%litr_end)+ &
       countelm(this%wood_beg, this%wood_end)               + &
       countelm(this%humus_beg, this%humus_end)             + &
       countelm(this%micbiom_beg, this%micbiom_end)         + &
       countelm(this%pom_beg, this%pom_end)                 + &
       countelm(this%dom_beg, this%dom_end)-1)/this%nelms   !include coarse wood debris

    allocate(this%is_ompool_som(this%nom_pools)); this%is_ompool_som(:)=.false.
    allocate(this%is_dom_pool(this%nom_pools)); this%is_dom_pool(:) = .false.
    this%is_ompool_som(this%mic)    =.true.
    this%is_ompool_som(this%pom)    =.true.
    this%is_ompool_som(this%humus)  =.true.
    this%is_ompool_som(this%dom)    =.true.
    this%is_dom_pool(this%dom)      =.true.

    itemp               = this%nom_pools*this%nelms + 1

    this%nom_tot_elms    = itemp

    this%lid_minp_secondary = addone(itemp); this%lid_minp_secondary_to_sol_occ_reac=addone(ireac)
    call list_insert(list_name, 'minp_secondary',vid, itype=var_state_type)
    call list_insert(list_unit, 'mol P m-3',uid)

    this%lid_minp_occlude = addone(itemp);
    call list_insert(list_name, 'minp_occlude',vid, itype=var_state_type)
    call list_insert(list_unit, 'mol P m-3',uid)
    if(maxpft>0)then
      this%lid_plant_minn_nh4_up_reac = addone(ireac)
      this%lid_plant_minn_no3_up_reac = addone(ireac)
      this%lid_plant_minp_up_reac = addone(ireac)
    endif
    this%lid_autr_rt_reac = addone(ireac)

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
      if(nop_limit  .or. spinup_state /= 0)then
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
      if(.not. (nop_limit  .or. spinup_state /= 0))then
        this%nprimvars = this%nprimvars + 1
      endif
    endif
    this%lid_nh4_nit_reac = addone(ireac)       !this is also used to indicate the nitrification reaction
    this%lid_no3_den_reac = addone(ireac)       !this is also used to indicate the denitrification reaction

    this%lid_minp_soluble_to_secp_reac = addone(ireac)

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
    this%lid_o2_aren_reac  = addone(ireac)

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
    call list_insert(list_name, 'o2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)
    if ( spinup_state == 0 ) then
       this%lid_ar_paere   = addone(itemp);  this%lid_ar_aren_reac  = addone(ireac)   !
       call list_insert(list_name, 'ar_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)

       this%lid_n2_paere   = addone(itemp);  this%lid_n2_aren_reac  = addone(ireac)   !
       call list_insert(list_name, 'n2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol N2 m-3 s-1',uid)

       this%lid_co2_paere  = addone(itemp);  this%lid_co2_aren_reac = addone(ireac)   !
       call list_insert(list_name, 'co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol C m-3 s-1',uid)

       if(use_c13)then
         this%lid_c13_co2_paere  = addone(itemp);  this%lid_c13_co2_aren_reac = addone(ireac)   !
         call list_insert(list_name, 'c13_co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol C13 m-3 s-1',uid)
       endif

       if(use_c14)then
         this%lid_c14_co2_paere  = addone(itemp);  this%lid_c14_co2_aren_reac = addone(ireac)   !
         call list_insert(list_name, 'c14_co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol C14 m-3 s-1',uid)
       endif

       this%lid_ch4_paere  = addone(itemp);  this%lid_ch4_aren_reac = addone(ireac)   !
       call list_insert(list_name, 'ch4_paere',vid, itype=var_flux_type); call list_insert(list_unit, 'mol C m-3 s-1',uid)

       this%lid_n2o_paere  = addone(itemp);  this%lid_n2o_aren_reac = addone(ireac)   !
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

    this%nstvars          = itemp          !totally 14+32 state variables
    this%nreactions = ireac            !seven decomposition pathways plus root auto respiration, nitrification, denitrification and plant immobilization

    allocate(this%primvarid(ireac)); this%primvarid(:) = -1
    allocate(this%is_aerobic_reac(ireac)); this%is_aerobic_reac(:)=.false.

    allocate(this%varnames(this%nstvars))
    allocate(this%varunits(this%nstvars))
    allocate(this%ompoolnames(this%nom_pools))
    allocate(this%vartypes(this%nstvars))

    call copy_name(this%nstvars, list_name, this%varnames(1:this%nstvars))
    call copy_name(this%nstvars, list_unit, this%varunits(1:this%nstvars))
    call copy_name(this%nom_pools, list_pool, this%ompoolnames(1:this%nom_pools))
    call copy_name_type(this%nstvars, list_name, this%vartypes(1:this%nstvars))
!    call list_disp(list_name);call list_disp(list_pool);call list_disp(list_unit)
!    pause
    call list_free(list_name)
    call list_free(list_pool)
    call list_free(list_unit)
  end subroutine InitPars
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this)
    !
    ! !DESCRIPTION:
    ! memory allocation for the data type specified by this
    !
  implicit none
    ! !ARGUMENTS:
  class(cdom_bgc_index_type), intent(inout) :: this



  end subroutine InitAllocate
  !-------------------------------------------------------------------------------
  subroutine set_primvar_reac_ids(this,maxpft_loc)

  implicit none
  class(cdom_bgc_index_type), intent(inout)  :: this
  integer, intent(in) :: maxpft_loc
  integer :: reac

  associate(                                      &
    lmet  => this%lmet                          , &
    lcel  => this%lcel                          , &
    llig  => this%llig                          , &
    cwd  => this%cwd                            , &
    lwd  => this%lwd                            , &
    fwd  => this%fwd                            , &
    mic  => this%mic                            , &
    pom  => this%pom                            , &
    humus  => this%humus                        , &
    dom  => this%dom                            , &
    nelms => this%nelms                         , &
    c_loc => this%c_loc                         , &
    n_loc => this%n_loc                         , &
    p_loc => this%p_loc                         , &
    lid_nh4=> this%lid_nh4                      , &
    lid_no3 => this%lid_no3                     , &
    lid_o2 => this%lid_o2                       , &
    lid_minp_soluble=> this%lid_minp_soluble    , &
    lid_minp_secondary => this%lid_minp_secondary &
  )
  !reaction 1, lmet -> dom
  reac=this%lmet_dek_reac;     this%primvarid(reac) = (lmet-1)*nelms + c_loc

  !reaction 2, lcel -> dom
  reac =this%lcel_dek_reac;    this%primvarid(reac) = (lcel-1)*nelms + c_loc

  !reaction 3, pom->dom
  reac = this%pom_dek_reac; this%primvarid(reac) = (pom-1)*nelms + c_loc

  !reaction 4, humus -> dom
  reac =this%humus_dek_reac;  this%primvarid(reac) = (humus-1)*nelms + c_loc

  !reaction 5, llig->pom
  reac =this%llig_dek_reac; this%primvarid(reac) = (llig-1)*nelms + c_loc

  !reaction 6, dom -> mic
  reac = this%dom_dek_reac; this%primvarid(reac) = (dom-1)*nelms + c_loc

  !reaction 7, mic -> f1*pom + f2*humus
  reac =this%mic_dek_reac; this%primvarid(reac) = (mic-1)*nelms + c_loc

  !reaction 8, cwd -> lmet + lcel
  reac = this%cwd_dek_reac; this%primvarid(reac) = (cwd-1)*nelms+c_loc

  !reaction 9, fwd -> lmet + lcel
  reac = this%fwd_dek_reac; this%primvarid(reac) = (fwd-1)*nelms+c_loc

  !reaction 10, lwd -> lmet + lcel
  reac = this%lwd_dek_reac; this%primvarid(reac) = (lwd-1)*nelms+c_loc

  !reaction 11, nitrification
  reac = this%lid_nh4_nit_reac; this%primvarid(reac) = this%lid_nh4

  !reaction 12, denitrification
  reac = this%lid_no3_den_reac; this%primvarid(reac) = lid_no3

  !reaction 13, inorganic P non-equilibrium adsorption
  !P_solution -> p_secondary
  reac = this%lid_minp_soluble_to_secp_reac; this%primvarid(reac) = lid_minp_soluble

  !reaction 14, inorganic P non-equilibrium desorption
  ! p_secondary -> P_solution
  reac = this%lid_minp_secondary_to_sol_occ_reac; this%primvarid(reac) = lid_minp_secondary
  if(maxpft_loc>0)then
    !reaction 15, plant mineral nitrogen nh4 uptake
    reac = this%lid_plant_minn_nh4_up_reac; this%primvarid(reac) = lid_nh4

    !reaction 16, plant mineral nitrogen no3 uptake
    reac = this%lid_plant_minn_no3_up_reac; this%primvarid(reac) = lid_no3

    !reaction 17, plant mineral phosphorus uptake
    reac = this%lid_plant_minp_up_reac; this%primvarid(reac) = lid_minp_soluble
  endif
  !reaction 18, autotrophic respiration, ar + o2 -> co2

  !reaction 19, o2 transport through arenchyma
  reac = this%lid_o2_aren_reac; this%primvarid(reac) = lid_o2
  end associate
  end subroutine set_primvar_reac_ids
end module cdomBGCIndexType
