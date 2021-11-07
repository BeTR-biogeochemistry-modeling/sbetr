module ecosysBGCIndexType
  !DESCRIPTION
  !Soil bgc built based on the ecosys (nitro.f) formulation
  !aqueous chemistry (elctrolyte chemistry and precipitation/dissolution) is done elsewhere

  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_ctrl    , only : spinup_state => betr_spinup_state
  use gBGCIndexType  , only : gbgc_index_type
  use betr_varcon    , only : var_flux_type, var_state_type
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  integer, parameter :: loc_name_len=64
  integer, parameter :: loc_name_len2=128
  character(len=loc_name_len), parameter :: scomp(5)=&
       (/'woody_lit    ','non_woody_lit','manure_lit   ','POM          ','humus        '/)

  type, private :: list_t
    character(len=loc_name_len) :: name  !name of the variable
    integer :: id                        !id of the variable in the index vector
    integer :: itype                     !type of the variable: flux for state
    type(list_t), pointer :: next => null()
  end type list_t

  type, public, extends(gbgc_index_type) :: ecosys_bgc_index_type

     integer :: n_om_complxes           !number of litter om-microbial complexes
     integer :: n_kinetc_compont        !number of kinetic components for each litter complex
     integer :: n_micbr_compont         !number of microbial biomass components for residual
     integer :: n_micbl_compont         !number of microbial biomass components for living biomass
     integer :: n_micbl_guilds          !number of guilds for each microbes
     integer :: n_micbl_fgrup           !number of functional groups associated with om-microbial complex
     integer :: gid_obli_aerobes        !group id for obligate aerobic bacteria
     integer :: gid_nh3_oxidizers       !group id for ammonia oxidizing bacteria
     integer :: gid_facul_anaerobes     !group id for facultative anaerobic bacteria
     integer :: gid_no2_oxidizers       !group id for nitrite oxidizing bacteria
     integer :: gid_fungi               !group id for saprotrophic fungi
     integer :: gid_methanotrophs       !group id for methanotrophs
     integer :: gid_anero_fermenters    !group id for anaerobic fermenters
     integer :: gid_aceto_methanogens   !group id for acetotrohpic methanogens
     integer :: gid_hydrog_methanogens  !group id for hydrogenic methanogens
     integer :: gid_aerob_dizotrophs    !group id for aerobic dizotrophic N2 fixers
     integer :: gid_anero_dizotrophs    !group id for anaerobic dizotrophic N2 fixers
     integer, allocatable :: lid_omcomplx_osc(:,:)        !id for om-microbial complex C
     integer, allocatable :: lid_omcomplx_osc_c13(:,:)    !id for om-microbial complex C13
     integer, allocatable :: lid_omcomplx_osc_c14(:,:)    !id for om-microbial complex C14
     integer, allocatable :: lid_omcomplx_osn(:,:)        !id for om-microbial complex N
     integer, allocatable :: lid_omcomplx_osp(:,:)        !id for om-microbial complex P
     integer, allocatable :: lid_omcomplx_osa(:,:)        !id for om-microbial complex colonized
     integer, allocatable :: lid_mbresidu_orc(:,:)        !id for microbial residual C
     integer, allocatable :: lid_mbresidu_orc_c13(:,:)    !id for microbial residual C
     integer, allocatable :: lid_mbresidu_orc_c14(:,:)    !id for microbial residual C
     integer, allocatable :: lid_mbresidu_orn(:,:)        !id for microbial residual N
     integer, allocatable :: lid_mbresidu_orp(:,:)        !id for microbial residual P
     integer, allocatable :: lid_omsorbed_ohc(:)          !id for adsorbed C
     integer, allocatable :: lid_omsorbed_ohc_c13(:)      !id for adsorbed C13
     integer, allocatable :: lid_omsorbed_ohc_c14(:)      !id for adsorbed C14
     integer, allocatable :: lid_omsorbed_ohn(:)          !id for adsorbed N
     integer, allocatable :: lid_omsorbed_ohp(:)          !id for adsorbed P
     integer, allocatable :: lid_micbioms_omc(:,:,:,:)    !id for living microbial C
     integer, allocatable :: lid_micbioms_omc_c13(:,:,:,:)    !id for living microbial C
     integer, allocatable :: lid_micbioms_omc_c14(:,:,:,:)    !id for living microbial C
     integer, allocatable :: lid_micbioms_omn(:,:,:,:)    !id for living microbial N
     integer, allocatable :: lid_micbioms_omp(:,:,:,:)    !id for living microbial P
     integer, allocatable :: lid_dom_oqc(:)               !id for DOC
     integer, allocatable :: lid_dom_oqc_c13(:)               !id for DOC
     integer, allocatable :: lid_dom_oqc_c14(:)               !id for DOC
     integer, allocatable :: lid_dom_oqn(:)               !id for DON
     integer, allocatable :: lid_dom_oqp(:)               !id for DOP
     integer, allocatable :: lid_acetate_oqa(:)           !id for acetate
     integer :: lid_h2, lid_h2_paere
     integer :: lid_o2, lid_o2_paere
     integer :: lid_ch4,lid_ch4_paere
     integer :: lid_ar, lid_ar_paere
     integer :: lid_c13_ch4, lid_c14_ch4, lid_c13_ch4_paere, lid_c14_ch4_paere
     integer :: lid_co2x, lid_co2x_paere
     integer :: lid_c13_co2x, lid_c13_co2x_paere, lid_c14_co2x, lid_c14_co2x_paere
     integer :: lid_nh3x, lid_nh3x_paere
     integer :: lid_n2, lid_n2_paere
     integer :: lid_n2o, lid_n2o_paere
     integer :: lid_hno2, lid_hno3
     integer :: lid_h0po4, lid_h1po4, lid_h2po4, lid_h3po4
     integer :: lid_co2_hr, lid_n2o_den, lid_n2o_nit

     integer :: lid_supp_minp, lid_supp_minn
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
     integer :: lid_plant_minn_no3, lid_plant_minn_nh4
     integer :: lid_plant_minp
     logical :: debug
     character(len=loc_name_len), allocatable :: varnames(:)
     character(len=loc_name_len), allocatable :: varunits(:)
     integer                    , allocatable :: vartypes(:)
   contains
     procedure, public  :: Init
     procedure, public  :: hcopy
     procedure, private :: InitPars
     procedure, private :: InitAllocate
     procedure, private :: set_primvar_reac_ids
     procedure, public  :: display_index
     procedure, private :: add_omcomplx_names
     procedure, private :: add_complx1d_names
     procedure, private :: add_micbrcomplx_names
     procedure, private :: add_micblcomplx_names
  end type ecosys_bgc_index_type

  contains
  !-----------------------------------------------------------------------
  subroutine list_init(self, name, id, itype)
  !
  ! Initialize the list
  implicit none
  type(list_t)     , pointer       :: self
  character(len=*) , intent(in)    :: name
  integer          , intent(inout) :: id
  integer, optional, intent(in)    :: itype

  allocate(self)
  nullify(self%next)
  id=id+1
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
  subroutine copy_var_type(num_names, list_name, vartypes)

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
  end subroutine copy_var_type
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
  subroutine add_omcomplx_names(this, list_name, list_unit, var, vid, uid)
  !DESCRIPTION
  !add variable name and unit for organic matter complex
  implicit none
  class(ecosys_bgc_index_type) :: this
  type(list_t), pointer :: list_name
  type(list_t), pointer :: list_unit
  character(len=*), intent(in) :: var
  integer, intent(inout) :: vid
  integer, intent(inout) :: uid
  character(len=loc_name_len), parameter :: kcomp(4)=&
       (/'protein   ','carbohydro','cellulose ','lignin    '/)
  integer :: j1, j2

  select case (var)
  case ('C')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_kinetc_compont
    call list_insert(list_name,'C_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol C m-3',uid)
    enddo
    enddo
  case ('C13')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_kinetc_compont
    call list_insert(list_name,'C13_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol C13 m-3',uid)
    enddo
    enddo
  case ('C14')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_kinetc_compont
    call list_insert(list_name,'C14_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol C14 m-3',uid)
    enddo
    enddo
  case ('N')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_kinetc_compont
    call list_insert(list_name,'N_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol N m-3',uid)
    enddo
    enddo
  case ('P')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_kinetc_compont
    call list_insert(list_name,'P_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol P m-3',uid)
    enddo
    enddo
  case default
  !do nothing
  end select
  end subroutine add_omcomplx_names
  !-------------------------------------------------------------------------------
  subroutine add_micbrcomplx_names(this, list_name, list_unit, var, vid, uid)
  !DESCRIPTION
  !add variable name and unit for microbial residual
  implicit none
  class(ecosys_bgc_index_type) :: this
  type(list_t), pointer :: list_name
  type(list_t), pointer :: list_unit
  character(len=*), intent(in) :: var
  integer, intent(inout) :: vid
  integer, intent(inout) :: uid
  character(len=loc_name_len), parameter :: kcomp(2)=&
       (/'labile      ','recalcitrant'/)
  integer :: j1, j2
  select case (var)
  case ('C')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_micbr_compont
    call list_insert(list_name,'C_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid, itype=var_state_type)
    call list_insert(list_unit,'mol C m-3',uid)
    enddo
    enddo
  case ('C13')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_micbr_compont
    call list_insert(list_name,'C13_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol C13 m-3',uid)
    enddo
    enddo
  case ('C14')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_micbr_compont
    call list_insert(list_name,'C14_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol C14 m-3',uid)
    enddo
    enddo
  case ('N')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_micbr_compont
    call list_insert(list_name,'N_'//trim(kcomp(j1))//'_cplx_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol N m-3',uid)
    enddo
    enddo
  case ('P')
    do j2=1, this%n_om_complxes
    do j1=1, this%n_micbr_compont
    call list_insert(list_name,'P_'//trim(kcomp(j1))//'_'//trim(scomp(j2)),vid,itype=var_state_type)
    call list_insert(list_unit,'mol P m-3',uid)
    enddo
    enddo
  case default
  !do nothing
  end select
  end subroutine add_micbrcomplx_names
  !-------------------------------------------------------------------------------
  subroutine add_complx1d_names(this,list_name, list_unit, varname, varunit, vid, uid)
  !DESCRIPTION
  !add variables to 1d complex
  class(ecosys_bgc_index_type) :: this
  type(list_t), pointer :: list_name
  type(list_t), pointer :: list_unit
  character(len=*), intent(in) :: varname
  character(len=*), intent(in) :: varunit
  integer, intent(inout) :: vid
  integer, intent(inout) :: uid
  integer :: j1

  do j1=1, this%n_om_complxes
    call list_insert(list_name,trim(varname)//'_cplx_'//trim(scomp(j1)),vid,itype=var_state_type)
    call list_insert(list_unit,trim(varunit),uid)
  enddo
  end subroutine add_complx1d_names
  !-------------------------------------------------------------------------------
  subroutine add_micblcomplx_names(this, list_name, list_unit, var, vid, uid)
  !DESCRIPTION
  !add variables related to living microbes
  implicit none
  class(ecosys_bgc_index_type) :: this
  type(list_t), pointer :: list_name
  type(list_t), pointer :: list_unit
  character(len=*), intent(in) :: var
  integer, intent(inout) :: vid
  integer, intent(inout) :: uid
  !local variables
  integer :: j1, j2, j3, j4
  character(len=loc_name_len) :: comp
  character(len=loc_name_len), parameter :: fgrp1(7) = &
        (/'obligate aerobe        ','facultative anaerobe   ', &
          'fungi                  ','anaerobic fermenter    ', &
          'acetotrophic methanogen','aerobic diazotroph     ', &
          'anaerobic diazotroph   '/)
  character(len=loc_name_len), parameter :: fgrp2(7) = &
          (/'ammonia oxidizer       ','nitrite oxidizer       ', &
            'methanotroph           ','none                   ', &
            'hydrotrophic methanogen','none                   ', &
            'none                   '/)
  character(len=loc_name_len) :: biomcomp(3)=(/'labile      ','recalcitrant','reserve     '/)
  character(len=loc_name_len) :: gldname
  character(len=loc_name_len2) :: varname
  select case (var)
  case ('C')
   do j4 = 1, this%n_om_complxes
   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='C_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_'//trim(scomp(j4))
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol C m-3',uid)
   enddo
   enddo
   enddo
   enddo

   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='C_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_abstract'
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol C m-3',uid)
   enddo
   enddo
   enddo

  case ('C13')
   do j4 = 1, this%n_om_complxes
   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='C13_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_'//trim(scomp(j4))
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol C13 m-3',uid)
   enddo
   enddo
   enddo
   enddo

   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='C13_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_abstract'
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol C13 m-3',uid)
   enddo
   enddo
   enddo

  case ('C14')
   do j4 = 1, this%n_om_complxes
   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='C14_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_'//trim(scomp(j4))
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol C14 m-3',uid)
   enddo
   enddo
   enddo
   enddo

   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='C14_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_abstract'
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol C14 m-3',uid)
   enddo
   enddo
   enddo

  case ('N')
   do j4 = 1, this%n_om_complxes
   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='N_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_'//trim(scomp(j4))
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol N m-3',uid)
   enddo
   enddo
   enddo
   enddo

   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='N_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_abstract'
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol N m-3',uid)
   enddo
   enddo
   enddo

  case ('P')
   do j4 = 1, this%n_om_complxes
   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='P_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_'//trim(scomp(j4))
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol P m-3',uid)
   enddo
   enddo
   enddo
   enddo

   do j3 = 1, this%n_micbl_fgrup
   do j2 = 1, this%n_micbl_guilds
   write(gldname,'(A,I2.2)')'guild-',j2
   do j1 = 1, this%n_micbl_compont
   varname='P_'//trim(biomcomp(j1))//'_'//trim(gldname)&
      //'_'//trim(fgrp1(j3))//'_cplx_abstract'
   call list_insert(list_name, trim(varname),vid, itype=var_state_type)
   call list_insert(list_unit,'mol P m-3',uid)
   enddo
   enddo
   enddo

  case default
  end select
  end subroutine add_micblcomplx_names
  !-------------------------------------------------------------------------------
  subroutine Init(this, use_c13, use_c14, non_limit, nop_limit, maxpft, batch_mode)
    !
    ! DESCRIPTION:
    ! Initialize ecosys_bgc type
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


  maxpft_loc = 0; batch_mode_loc=.false.
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
    use MathfuncMod   , only : addone, countelm, assign_ids
    use betr_utils    , only : num2str
    use betr_constants, only : betr_string_length_long
    implicit none
    class(ecosys_bgc_index_type) :: this
    integer, intent(in) :: maxpft
    logical, intent(in) :: use_c13
    logical, intent(in) :: use_c14
    logical, intent(in) :: non_limit        !switch to turn off nitrogen limitation
    logical, intent(in) :: nop_limit        !switch to turn off phosphorus limitation
    logical, intent(in) :: batch_mode       !switch to run the model in batch mode
    ! !LOCAL VARIABLES:
    integer :: itemp
    integer :: ireac   !counter of reactions
    integer :: itemp0, itemp1
    integer :: ielem
    integer :: vid,uid,pid
    integer :: jj, j1, j2
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

    this%n_om_complxes = 5
    this%n_kinetc_compont = 4
    this%n_micbr_compont = 2
    this%n_micbl_compont = 3
    this%n_micbl_fgrup = 7
    this%n_micbl_guilds = 1
    this%gid_obli_aerobes = 1
    this%gid_nh3_oxidizers = 1
    this%gid_facul_anaerobes = 2
    this%gid_no2_oxidizers = 2
    this%gid_fungi = 3
    this%gid_methanotrophs = 3
    this%gid_anero_fermenters = 4
    this%gid_aceto_methanogens = 5
    this%gid_hydrog_methanogens = 5
    this%gid_aerob_dizotrophs = 6
    this%gid_anero_dizotrophs = 7

    itemp=0
    this%lid_h2=addone(itemp);
    call list_init(list_name, "H2",vid, itype=var_state_type)
    call list_init(list_unit, 'mol H2 m-3',uid)
    this%lid_o2=addone(itemp);
    call list_insert(list_name, 'O2',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol O2 m-3',uid)
    this%lid_ar=addone(itemp);
    call list_insert(list_name, 'Ar',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol Ar m-3',uid)
    this%lid_nh3x=addone(itemp);
    call list_insert(list_name, 'NH3x',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol NH3x m-3',uid)
    this%lid_n2=addone(itemp)
    call list_insert(list_name, 'N2',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol N2 m-3',uid)
    this%lid_n2o=addone(itemp)
    call list_insert(list_name, 'N2O',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol N2O m-3',uid)
    this%lid_hno2=addone(itemp)
    call list_insert(list_name, 'HNO2',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol HNO2 m-3',uid)
    this%lid_hno3=addone(itemp)
    call list_insert(list_name, 'HNO3',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol HNO3 m-3',uid)
    this%lid_h0po4=addone(itemp)
    call list_insert(list_name, 'H0PO4',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol H0PO4 m-3',uid)
    this%lid_h1po4=addone(itemp)
    call list_insert(list_name, 'H1PO4',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol H1PO4 m-3',uid)
    this%lid_h2po4=addone(itemp)
    call list_insert(list_name, 'H2PO4',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol H2PO4 m-3',uid)
    this%lid_h3po4=addone(itemp)
    call list_insert(list_name, 'H3PO4',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol H3PO4 m-3',uid)
    this%lid_ch4=addone(itemp)
    call list_insert(list_name, 'CH4',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol CH4 m-3',uid)
    this%lid_co2x=addone(itemp)
    call list_insert(list_name, 'CO2x',vid, itype=var_state_type);
    call list_insert(list_unit, 'mol CO2x m-3',uid)
    if(use_c13)then
      this%lid_c13_ch4=addone(itemp)
      call list_insert(list_name, '13CH4',vid, itype=var_state_type);
      call list_insert(list_unit, 'mol 13CH4 m-3',uid)
      this%lid_c13_co2x=addone(itemp)
      call list_insert(list_name, '13CO2x',vid, itype=var_state_type);
      call list_insert(list_unit, 'mol 13CO2x m-3',uid)
    endif

    if(use_c14)then
      this%lid_c14_ch4=addone(itemp)
      call list_insert(list_name, '14CH4',vid, itype=var_state_type);
      call list_insert(list_unit, 'mol 14CH4 m-3',uid)
      this%lid_c14_co2x=addone(itemp)
      call list_insert(list_name, '14CO2x',vid, itype=var_state_type);
      call list_insert(list_unit, 'mol 14CO2x m-3',uid)
    endif

    allocate(this%lid_omcomplx_osc(this%n_kinetc_compont,this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_omcomplx_osc, itemp)
    call this%add_omcomplx_names(list_name, list_unit, 'C', vid=vid,uid=uid)
    if(use_c13)then
      allocate(this%lid_omcomplx_osc_c13(this%n_kinetc_compont,this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_omcomplx_osc_c13, itemp)
      call this%add_omcomplx_names(list_name, list_unit, 'C13', vid=vid,uid=uid)
    endif

    if(use_c14)then
      allocate(this%lid_omcomplx_osc_c14(this%n_kinetc_compont,this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_omcomplx_osc_c14, itemp)
      call this%add_omcomplx_names(list_name, list_unit, 'C14', vid=vid,uid=uid)
    endif
    allocate(this%lid_omcomplx_osn(this%n_kinetc_compont,this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_omcomplx_osn, itemp)
    call this%add_omcomplx_names(list_name, list_unit, 'N', vid=vid,uid=uid)

    allocate(this%lid_omcomplx_osp(this%n_kinetc_compont,this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_omcomplx_osp, itemp)
    call this%add_omcomplx_names(list_name, list_unit, 'P', vid=vid,uid=uid)

    allocate(this%lid_omcomplx_osa(this%n_kinetc_compont,this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_omcomplx_osa, itemp)
    call this%add_omcomplx_names(list_name, list_unit, 'A', vid=vid,uid=uid)

    allocate(this%lid_mbresidu_orc(this%n_micbr_compont, this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_mbresidu_orc, itemp)
    call this%add_micbrcomplx_names(list_name, list_unit, 'C', vid=vid,uid=uid)
    if(use_c13)then
      allocate(this%lid_mbresidu_orc_c13(this%n_micbr_compont, this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_mbresidu_orc_c13, itemp)
      call this%add_micbrcomplx_names(list_name, list_unit, 'C13', vid=vid,uid=uid)
    endif
    if(use_c14)then
      allocate(this%lid_mbresidu_orc_c14(this%n_micbr_compont, this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_mbresidu_orc_c14, itemp)
      call this%add_micbrcomplx_names(list_name, list_unit, 'C14', vid=vid,uid=uid)
    endif

    allocate(this%lid_mbresidu_orn(this%n_micbr_compont, this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_mbresidu_orn, itemp)
    call this%add_micbrcomplx_names(list_name, list_unit, 'N', vid=vid,uid=uid)

    allocate(this%lid_mbresidu_orp(this%n_micbr_compont, this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_mbresidu_orp, itemp)
    call this%add_micbrcomplx_names(list_name, list_unit, 'P', vid=vid,uid=uid)

    allocate(this%lid_omsorbed_ohc(this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_omsorbed_ohc, itemp)
    call this%add_complx1d_names(list_name, list_unit, 'sorbed C', 'mol C m-3',vid=vid,uid=uid)
    if(use_c13)then
      allocate(this%lid_omsorbed_ohc_c13(this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_omsorbed_ohc_c13, itemp)
      call this%add_complx1d_names(list_name, list_unit, 'sorbed C13','mol C13 m-3', vid=vid,uid=uid)
    endif
    if(use_c14)then
      allocate(this%lid_omsorbed_ohc_c14(this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_omsorbed_ohc_c14, itemp)
      call this%add_complx1d_names(list_name, list_unit, 'sorbed C14', 'mol C14 m-3',vid=vid,uid=uid)
    endif

    allocate(this%lid_omsorbed_ohn(this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_omsorbed_ohn, itemp)
    call this%add_complx1d_names(list_name, list_unit, 'sorbed N', 'mol N m-3',vid=vid,uid=uid)

    allocate(this%lid_omsorbed_ohp(this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_omsorbed_ohp, itemp)
    call this%add_complx1d_names(list_name, list_unit, 'sorbed P', 'mol P m-3',vid=vid,uid=uid)

    allocate(this%lid_micbioms_omc(this%n_micbl_compont, this%n_micbl_guilds, this%n_micbl_fgrup, this%n_om_complxes+1))
    itemp=itemp+assign_ids(this%lid_micbioms_omc, itemp)
    call this%add_micblcomplx_names(list_name, list_unit, 'C', vid=vid,uid=uid)
    if(use_c13)then
      allocate(this%lid_micbioms_omc_c13(this%n_micbl_compont, this%n_micbl_guilds, this%n_micbl_fgrup, this%n_om_complxes+1))
      itemp=itemp+assign_ids(this%lid_micbioms_omc_c13, itemp)
      call this%add_micblcomplx_names(list_name, list_unit, 'C13', vid=vid,uid=uid)
    endif
    if(use_c14)then
      allocate(this%lid_micbioms_omc_c14(this%n_micbl_compont, this%n_micbl_guilds, this%n_micbl_fgrup, this%n_om_complxes+1))
      itemp=itemp+assign_ids(this%lid_micbioms_omc_c14, itemp)
      call this%add_micblcomplx_names(list_name, list_unit, 'C14', vid=vid,uid=uid)
    endif

    allocate(this%lid_micbioms_omn(this%n_micbl_compont, this%n_micbl_guilds, this%n_micbl_fgrup, this%n_om_complxes+1))
    itemp=itemp+assign_ids(this%lid_micbioms_omn, itemp)
    call this%add_micblcomplx_names(list_name, list_unit, 'N', vid=vid,uid=uid)

    allocate(this%lid_micbioms_omp(this%n_micbl_compont, this%n_micbl_guilds, this%n_micbl_fgrup, this%n_om_complxes+1))
    itemp=itemp+assign_ids(this%lid_micbioms_omp, itemp)
    call this%add_micblcomplx_names(list_name, list_unit, 'P', vid=vid,uid=uid)

    allocate(this%lid_dom_oqc(this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_dom_oqc, itemp)
    call this%add_complx1d_names(list_name,list_unit, 'DOC','mol C m-3',vid=vid,uid=uid)
    if(use_c13)then
      allocate(this%lid_dom_oqc_c13(this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_dom_oqc_c13, itemp)
      call this%add_complx1d_names(list_name,list_unit, 'DOC_C13','mol C13 m-3',vid=vid,uid=uid)
    endif
    if(use_c14)then
      allocate(this%lid_dom_oqc_c14(this%n_om_complxes))
      itemp=itemp+assign_ids(this%lid_dom_oqc_c14, itemp)
      call this%add_complx1d_names(list_name,list_unit, 'DOC_C14','mol C14 m-3',vid=vid,uid=uid)
    endif

    allocate(this%lid_dom_oqn(this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_dom_oqn, itemp)
    call this%add_complx1d_names(list_name,list_unit, 'DON','mol N m-3',vid=vid,uid=uid)

    allocate(this%lid_dom_oqp(this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_dom_oqp, itemp)
    call this%add_complx1d_names(list_name,list_unit, 'DOP','mol P m-3',vid=vid,uid=uid)

    allocate(this%lid_acetate_oqa(this%n_om_complxes))
    itemp=itemp+assign_ids(this%lid_acetate_oqa, itemp)
    call this%add_complx1d_names(list_name,list_unit, 'acetate','mol m-3',vid=vid,uid=uid)

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


    allocate(this%vartypes(this%nstvars))
    allocate(this%varnames(this%nstvars))
    allocate(this%varunits(this%nstvars))

    call copy_name(this%nstvars, list_name, this%varnames(1:this%nstvars))
    call copy_name(this%nstvars, list_unit, this%varunits(1:this%nstvars))
    call copy_var_type(this%nstvars, list_name, this%vartypes(1:this%nstvars))
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
  class(ecosys_bgc_index_type), intent(inout) :: this


  end subroutine InitAllocate
  !-------------------------------------------------------------------------------
  subroutine set_primvar_reac_ids(this,maxpft_loc)

  implicit none
  class(ecosys_bgc_index_type), intent(inout)  :: this
  integer, intent(in) :: maxpft_loc
  integer :: reac

  associate(                                      &
    lid_o2 => this%lid_o2                         &
  )

  if(maxpft_loc>0)then
  endif

  end associate

  end subroutine set_primvar_reac_ids

!-------------------------------------------------------------------------------
  subroutine display_index(this)

  implicit none
  class(ecosys_bgc_index_type) :: this

  end subroutine display_index
!-------------------------------------------------------------------------------
  subroutine hcopy(this, that)
  !DESCRIPTION
  !copy indices from that to this
  implicit none
  class(ecosys_bgc_index_type) :: this
  class(ecosys_bgc_index_type), intent(in) :: that

  integer :: jj
  integer :: maxpft

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

  end subroutine hcopy
end module ecosysBGCIndexType
