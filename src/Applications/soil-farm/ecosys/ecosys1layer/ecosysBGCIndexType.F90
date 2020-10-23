module ecosysBGCIndexType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_ctrl      , only : spinup_state => betr_spinup_state
  use gBGCIndexType  , only : gbgc_index_type
  use betr_varcon    , only : var_flux_type, var_state_type
implicit none
private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

integer, parameter :: loc_name_len=64

  type, public, extends(gbgc_index_type) :: ecosys_bgc_index_type
    integer           :: lit1, lit1_depoly_reac
    integer           :: lit2, lit2_depoly_reac
    integer           :: lit3, lit3_depoly_reac
    integer           :: cwd , cwd_depoly_reac
    integer           :: lid_n2
    integer           :: lid_o2   , o2_resp_reac
    integer           :: lid_ar
    integer           :: lid_co2
    integer           :: lid_c13_co2
    integer           :: lid_c14_co2
    integer           :: lid_ch4
    integer           :: lid_dom, dom_uptake_reac
    integer           :: lid_cue, micbd_depoly_reac
    integer           :: lid_micbl, lid_micbd, micbl_mort_reac
  !diagnostic variables
    integer           :: lid_co2_hr
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
    integer           :: pom_beg,  pom_end   !pom group
    integer           :: Bm_beg,  Bm_end     !microbial group
    integer           :: doc_uptake_reac
    integer           :: nelms
    integer           :: c_loc, n_loc, p_loc
    integer           :: c13_loc, c14_loc
    integer           :: e_loc
    integer           :: nom_tot_elms
    integer           :: nom_pools
    integer           :: nprimvars        !total number of primary variables
    integer           :: nstvars          !number of equations for the state variabile vector
    integer           :: nreactions
    integer , pointer :: primvarid(:)   => null()
    logical           :: debug
    character(len=loc_name_len), allocatable :: varnames(:)
    character(len=loc_name_len), allocatable :: varunits(:)
    character(len=loc_name_len), allocatable :: ompoolnames(:)
    integer, allocatable :: vartypes(:)
  contains
    procedure, public  :: Init
    procedure, private :: InitPars
    procedure, private :: InitAllocate
    procedure, private :: set_primvar_reac_ids
  end type ecosys_bgc_index_type
contains
!-----------------------------------------------------------------------
  subroutine add_ompool_name(list_name, list_unit, list_pool, prefix, use_c13, use_c14, do_init, vid,uid,pid)
  !
  !DESCRIPTION
  !add organic matter pools to the list
  use listMod       , only : list_s, list_init, list_insert
  implicit none
  type(list_s), pointer :: list_name
  type(list_s), pointer :: list_unit
  type(list_s), pointer :: list_pool
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
  ! Initialize ecosys_bgc_index_type
  ! !USES:
  implicit none
  ! !ARGUMENTS:
  class(ecosys_bgc_index_type), intent(inout) :: this
  logical, intent(in) :: use_c13
  logical, intent(in) :: use_c14
  logical, intent(in) :: non_limit
  logical, intent(in) :: nop_limit
  integer, optional, intent(in) :: maxpft
  logical, optional, intent(in) :: batch_mode
  ! !LOCAL VARIABLES:
  integer :: maxpft_loc
  logical :: batch_mode_loc
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
  !
  ! !USES:
  use MathfuncMod   , only : addone, countelm
  use betr_utils    , only : num2str
  use betr_constants, only : betr_string_length_long
  use listMod       , only : list_s, copy_name, list_init, list_insert, list_free, copy_name_type
  implicit none
  class(ecosys_bgc_index_type) :: this
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
  type(list_s), pointer :: list_name => null()
  type(list_s), pointer :: list_unit => null()
  type(list_s), pointer :: list_pool => null()
  type(list_s), pointer :: list_react=> null()
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
  !modify the code below accordingly
  !this%lid_micbl = addone(itemp); this%micbl_mort_reac  = addone(ireac); call list_insert(list_react, 'micbl_mort_reac', itemp0)
  !call add_ompool_name(list_name, list_unit, list_pool,'MB_live', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
  !this%lid_micbd = addone(itemp); this%micbd_depoly_reac = addone(ireac); call list_insert(list_react, 'micbd_depoly_reac', itemp0)
  !call add_ompool_name(list_name, list_unit, list_pool,'MB_dead', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
  !this%Bm_end=this%Bm_beg-1+(this%lid_micbd-this%lid_micbl+1)*this%nelms
  this%Bm_end=this%Bm_beg
  !DOM, only one pool is defined at this moment
  this%dom_beg = this%Bm_end + 1
  this%lid_dom = addone(itemp);this%dom_uptake_reac  = addone(ireac); call list_insert(list_react, 'dom_uptake_reac', itemp0)
  call add_ompool_name(list_name, list_unit, list_pool,'DOC', use_c13, use_c14, do_init=.false., vid=vid,uid=uid,pid=pid)
  this%lid_cue = this%lid_dom + 1
  call list_insert(list_name, 'DOM_e',vid)
  call list_insert(list_unit, 'mol e m-3',uid)
  this%dom_end = this%dom_beg - 1 + (this%lid_dom-this%lid_dom+1)*(this%nelms+1)

  this%nom_pools = (countelm(this%litr_beg, this%litr_end)+&
     countelm(this%wood_beg,this%wood_end) + &
     countelm(this%Bm_beg,this%Bm_end))/this%nelms + &
     countelm(this%dom_beg,this%dom_end)/(this%nelms+1)

  itemp         = countelm(this%litr_beg, this%litr_end)+&
     countelm(this%wood_beg,this%wood_end) + &
     countelm(this%Bm_beg,this%Bm_end) + &
     countelm(this%dom_beg,this%dom_end)
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
  this%nstvars          = itemp
  this%nreactions = ireac

  allocate(this%primvarid(ireac)); this%primvarid(:) = -1
  allocate(this%vartypes(this%nstvars))
  allocate(this%varnames(this%nstvars))
  allocate(this%varunits(this%nstvars))
  allocate(this%ompoolnames(this%nom_pools))

  call copy_name(this%nstvars, list_name, this%varnames(1:this%nstvars))
  call copy_name(this%nstvars, list_unit, this%varunits(1:this%nstvars))
  call copy_name(this%nom_pools, list_pool, this%ompoolnames(1:this%nom_pools))
  call copy_name_type(this%nstvars, list_name, this%vartypes(1:this%nstvars))
  !call list_disp(list_name); call list_disp(list_pool);call list_disp(list_unit)
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
  class(ecosys_bgc_index_type), intent(inout) :: this

  if (this%dummy_compiler_warning) continue
  end subroutine InitAllocate

 !-------------------------------------------------------------------------------
 subroutine set_primvar_reac_ids(this)
 !
 !DESCRIPTION
 !set primary variable for each reaction
 implicit none
 class(ecosys_bgc_index_type), intent(inout) :: this

 !local variables
 integer :: reac

 reac=this%lit1_depoly_reac;   this%primvarid(reac) = this%lit1
 reac=this%lit2_depoly_reac;   this%primvarid(reac) = this%lit2
 reac=this%lit3_depoly_reac;   this%primvarid(reac) = this%lit3
 reac=this%cwd_depoly_reac;    this%primvarid(reac) = this%cwd
 reac=this%micbd_depoly_reac;  this%primvarid(reac) = this%lid_micbd
 reac=this%micbl_mort_reac;    this%primvarid(reac) = this%lid_micbl
 reac=this%dom_uptake_reac;    this%primvarid(reac) = this%lid_dom
 reac=this%o2_resp_reac;       this%primvarid(reac) = this%lid_o2

 end subroutine set_primvar_reac_ids
end module ecosysBGCIndexType

