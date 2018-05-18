module simicBGCIndexType

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

  type, public, extends(gbgc_index_type) :: simic_bgc_index_type
     integer           :: lit1, lit1_depoly_reac                  !c
     integer           :: lit2, lit2_depoly_reac
     integer           :: lit3, lit3_depoly_reac
     integer           :: cwd , cwd_depoly_reac
     integer           :: lid_micbl, micbl_mort_reac
     integer           :: lid_doc, doc_uptake_reac, doc_sorb_reac
     integer           :: lid_micbd, micbd_depoly_reac
     integer           :: lid_n2
     integer           :: lid_o2   , o2_resp_reac
     integer           :: lid_ar
     integer           :: lid_co2
     integer           :: lid_c13_co2
     integer           :: lid_c14_co2
     integer           :: lid_ch4
     integer           :: lid_co2_hr       !diagnostic variables
     integer           :: lid_doc_e
     integer           :: lid_o2_paere
     integer           :: lid_n2_paere
     integer           :: lid_ar_paere
     integer           :: lid_co2_paere
     integer           :: lid_c13_co2_paere
     integer           :: lid_c14_co2_paere
     integer           :: lid_ch4_paere
     integer           :: litr_beg, litr_end  !litr group
     integer           :: wood_beg, wood_end  !wood group
     integer           :: dom_beg,  dom_end   !dom group
     integer           :: Bm_beg,  Bm_end     !Bm group
     integer           :: pom_beg, pom_end
     integer           :: lid_pom, lid_pom_e, pom_desorb_reac
     integer           :: nelms
     integer           :: c_loc
     integer           :: c13_loc
     integer           :: c14_loc
     integer           :: e_loc
     integer           :: nom_tot_elms
     integer           :: nom_pools
     integer           :: nprimvars        !total number of primary variables
     integer           :: nstvars          !number of equations for the state variabile vector
     integer           :: nreactions
     integer           :: lid_totinput
     integer           :: lid_totstore
     integer           :: lid_cum_closs
     integer , pointer :: primvarid(:)   => null()
     logical :: debug
     character(len=loc_name_len), allocatable :: varnames(:)
     character(len=loc_name_len), allocatable :: varunits(:)
     character(len=loc_name_len), allocatable :: ompoolnames(:)
     integer, allocatable :: vartypes(:)
   contains
     procedure, public  :: Init
     procedure, private :: InitPars
     procedure, private :: InitAllocate
     procedure, private :: set_primvar_reac_ids
  end type simic_bgc_index_type

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
    write(*,'(A30,A,I0)')trim(next%name),' =',next%id
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
    ! Initialize simic type
    ! !USES:
  implicit none
  ! !ARGUMENTS:
  class(simic_bgc_index_type), intent(inout) :: this
  logical, intent(in) :: use_c13
  logical, intent(in) :: use_c14
  logical, intent(in) :: non_limit
  logical, intent(in) :: nop_limit
  integer, optional, intent(in) :: maxpft
  logical, optional, intent(in) :: batch_mode
  ! !LOCAL VARIABLES:
  integer :: maxpft_loc
  logical :: batch_mode_loc = .false.
  maxpft_loc = 0
  this%dom_beg=0; this%dom_end=-1
  if(present(maxpft))maxpft_loc=maxpft
  if(present(batch_mode))batch_mode_loc=batch_mode
  call this%InitPars(maxpft_loc, use_c14, use_c13, non_limit, nop_limit, batch_mode_loc)

  call this%InitAllocate()

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
    class(simic_bgc_index_type) :: this
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

    if(maxpft>=0)continue

    itemp = 0; itemp0=0
    ireac = 0
    ielem= 0
    vid = 0;uid=0;pid=0
    this%c13_loc=0; this%c14_loc=0
    this%c_loc = addone(ielem)

    if(use_c13)then
      this%c13_loc= addone(ielem)
    endif
    if(use_c14)then
      this%c14_loc=addone(ielem)
    endif
    this%nelms = ielem
    this%e_loc = ielem + 1

    !litter group
    this%litr_beg=1
    this%lit1 = addone(itemp);this%lit1_depoly_reac = addone(ireac); call list_init(list_react, 'lit1_depoly_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'lit1', use_c13, use_c14, do_init=.true., vid=vid,uid=uid,pid=pid)
    this%lit2 = addone(itemp);this%lit2_depoly_reac = addone(ireac); call list_insert(list_react, 'lit2_depoly_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'lit2', use_c13, use_c14, do_init=.false.,vid=vid,uid=uid,pid=pid)
    this%lit3 = addone(itemp);this%lit3_depoly_reac = addone(ireac); call list_insert(list_react, 'lit3_depoly_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'lit3', use_c13, use_c14, do_init=.false.,vid=vid,uid=uid,pid=pid)
    this%litr_end = this%litr_beg -1 + (this%lit3-this%lit1+1)*this%nelms

    !woody group
    this%wood_beg=this%litr_end+1
    this%cwd  = addone(itemp);this%cwd_depoly_reac  = addone(ireac); call list_insert(list_react, 'cwd_depoly_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'cwd', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%wood_end=this%wood_beg-1+(this%cwd-this%cwd+1)*this%nelms

    !microbial biomass group
    this%Bm_beg=this%wood_end+1
    this%lid_micbl = addone(itemp); this%micbl_mort_reac  = addone(ireac); call list_insert(list_react, 'micbl_mort_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'MB_live', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%lid_micbd = addone(itemp); this%micbd_depoly_reac = addone(ireac); call list_insert(list_react, 'micbd_depoly_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'MB_dead', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%Bm_end=this%Bm_beg-1+(this%lid_micbd-this%lid_micbl+1)*this%nelms

    this%dom_beg = this%Bm_end + 1
    this%lid_doc = addone(itemp);this%doc_uptake_reac  = addone(ireac); call list_insert(list_react, 'doc_uptake_reac', itemp0)
    this%doc_sorb_reac  = addone(ireac); call list_insert(list_react, 'doc_sorb_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'DOC', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%lid_doc_e = addone(itemp)
    call list_insert(list_name, 'DOC_e',vid)
    call list_insert(list_unit, 'mol e m-3',uid)
    this%dom_end = this%dom_beg - 1 + (this%lid_doc-this%lid_doc+1)*(this%nelms+1)

    this%pom_beg = this%dom_end + 1
    this%lid_pom = addone(itemp);this%pom_desorb_reac  = addone(ireac); call list_insert(list_react, 'pom_desorb_reac', itemp0)
    call add_ompool_name(list_name, list_unit, list_pool,'POM', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
    this%lid_pom_e = addone(itemp)
    call list_insert(list_name, 'POM_e',vid)
    call list_insert(list_unit, 'mol e m-3',uid)
    this%pom_end = this%pom_beg - 1 + (this%lid_pom-this%lid_pom+1)*(this%nelms+1)

    this%nom_pools = (countelm(this%litr_beg, this%litr_end)+&
       countelm(this%wood_beg,this%wood_end) + &
       countelm(this%Bm_beg,this%Bm_end))/this%nelms  + &
       (countelm(this%dom_beg,this%dom_end) + &
       countelm(this%pom_beg, this%pom_end))/(this%nelms+1)   !include coarse wood debris

    itemp               = countelm(this%litr_beg, this%litr_end) + countelm(this%wood_beg,this%wood_end) + &
        countelm(this%Bm_beg,this%Bm_end)  + countelm(this%dom_beg,this%dom_end) + &
        countelm(this%pom_beg,this%pom_end)

    this%nom_tot_elms    = itemp

    !non-reactive primary variables
    this%lid_ar         = addone(itemp);call list_insert(list_name, 'ar',vid, itype=var_state_type); call list_insert(list_unit, 'mol m-3',uid)

    !second primary variables
    this%lid_o2         = addone(itemp);call list_insert(list_name, 'o2',vid, itype=var_state_type); call list_insert(list_unit, 'mol m-3',uid)
    this%o2_resp_reac   = addone(ireac); call list_insert(list_react, 'o2_resp_reac', itemp0)
    this%lid_co2        = addone(itemp);call list_insert(list_name, 'co2',vid, itype=var_state_type);call list_insert(list_unit,'mol m-3',uid)

    if(use_c13)then
      this%lid_c13_co2  = addone(itemp);call list_insert(list_name, 'c13_co2',vid, itype=var_state_type);call list_insert(list_unit,'mol m-3',uid)
    endif
    if(use_c14)then
      this%lid_c14_co2  = addone(itemp);call list_insert(list_name, 'c14_co2',vid, itype=var_state_type);call list_insert(list_unit,'mol m-3',uid)
    endif

    this%lid_n2         = addone(itemp);call list_insert(list_name, 'n2',vid, itype=var_state_type); call list_insert(list_unit, 'mol N2 m-3',uid)

    this%lid_ch4        = addone(itemp);call list_insert(list_name, 'ch4',vid, itype=var_state_type); call list_insert(list_unit, 'mol ch4 m-3',uid)

    this%nprimvars      = itemp

    this%lid_co2_hr     = addone(itemp);call list_insert(list_name, 'co2_hr',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)

    this%lid_o2_paere  = addone(itemp);call list_insert(list_name, 'o2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)

    this%lid_n2_paere  = addone(itemp);call list_insert(list_name, 'n2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)

    this%lid_ar_paere  = addone(itemp);call list_insert(list_name, 'ar_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)

    this%lid_ch4_paere  = addone(itemp);call list_insert(list_name, 'ch4_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)

    this%lid_co2_paere  = addone(itemp);call list_insert(list_name, 'co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)
    if(use_c13)then
      this%lid_c13_co2_paere  = addone(itemp);call list_insert(list_name, 'c13_co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)
    endif
    if(use_c14)then
      this%lid_co2_paere  = addone(itemp);call list_insert(list_name, 'c14_co2_paere',vid, itype=var_flux_type); call list_insert(list_unit,'mol m-3 s-1',uid)
    endif
    if(batch_mode)then
      this%lid_totinput  = addone(itemp);call list_insert(list_name, 'c_tot_input',vid, itype=var_state_type); call list_insert(list_unit,'mol m-3',uid)
      this%lid_totstore  = addone(itemp);call list_insert(list_name, 'c_tot_store',vid, itype=var_state_type); call list_insert(list_unit,'mol m-3',uid)
      this%lid_cum_closs  = addone(itemp);call list_insert(list_name, 'c_loss',vid, itype=var_state_type); call list_insert(list_unit,'mol m-3',uid)
    endif
    this%nstvars          = itemp          !totally 14+32 state variables
    this%nreactions = ireac            !

    allocate(this%primvarid(ireac)); this%primvarid(:) = -1
    allocate(this%vartypes(this%nstvars))
    allocate(this%varnames(this%nstvars))
    allocate(this%varunits(this%nstvars))
    allocate(this%ompoolnames(this%nom_pools))

    call copy_name(this%nstvars, list_name, this%varnames(1:this%nstvars))
    call copy_name(this%nstvars, list_unit, this%varunits(1:this%nstvars))
    call copy_name(this%nom_pools, list_pool, this%ompoolnames(1:this%nom_pools))
    call copy_name_type(this%nstvars, list_name, this%vartypes(1:this%nstvars))
!    call list_disp(list_name); call list_disp(list_pool);call list_disp(list_unit);
!    print*,'nprimvars=',this%nprimvars
!    call list_disp(list_react)

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
  class(simic_bgc_index_type), intent(inout) :: this

  if (this%dummy_compiler_warning) continue
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------
  subroutine set_primvar_reac_ids(this)

  implicit none
  class(simic_bgc_index_type), intent(inout) :: this

  integer :: reac

  reac=this%lit1_depoly_reac;   this%primvarid(reac) = this%lit1
  reac=this%lit2_depoly_reac;   this%primvarid(reac) = this%lit2
  reac=this%lit3_depoly_reac;   this%primvarid(reac) = this%lit3
  reac=this%cwd_depoly_reac;    this%primvarid(reac) = this%cwd
  reac=this%micbd_depoly_reac;  this%primvarid(reac) = this%lid_micbd
  reac=this%micbl_mort_reac;    this%primvarid(reac) = this%lid_micbl
  reac=this%doc_uptake_reac;    this%primvarid(reac) = this%lid_doc
  reac=this%doc_sorb_reac  ;    this%primvarid(reac) = this%lid_doc
  reac=this%o2_resp_reac;       this%primvarid(reac) = this%lid_o2
  reac=this%pom_desorb_reac;    this%primvarid(reac) = this%lid_pom

  end subroutine set_primvar_reac_ids
end module simicBGCIndexType
