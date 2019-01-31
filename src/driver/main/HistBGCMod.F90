module HistbgcMod

  use histMod  , only : hist_var_str_len, hist_unit_str_len
implicit none


  type, public :: hist_bgc_type

    integer :: nhistvars
    character(len=hist_var_str_len), allocatable :: varl(:)
    character(len=hist_unit_str_len), allocatable :: unitl(:)
    integer                         , allocatable :: vartypes(:)
  contains
    procedure, public :: getvarllen
    procedure, public :: Init
  end type hist_bgc_type

contains

  !-----------------------------------------------------------------------

  subroutine Init(this, reaction_method)
  use listMod, only : list_s, copy_name, list_init, list_insert, list_free, copy_name_type
  use betr_varcon, only : var_flux_type, var_state_type
  use betr_constants           , only : stdout                  ! added for output testing remarks. -zlyu 01/27/2019

  implicit none
  class(hist_bgc_type), intent(inout) :: this
  character(len=*), intent(in) :: reaction_method

  type(list_s), pointer :: hist_list_var
  type(list_s), pointer :: hist_list_unit

  integer :: id1,id2
  id1 = 0; id2 = 0
  if(index(trim(reaction_method),'ecacnp')/=0 .or. &
     index(trim(reaction_method),'keca')/=0 .or. &
     index(trim(reaction_method),'ch4soil')/=0)then

    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'in if case ecacnp, before all calls'
    write(stdout, *) '***************************'
    ! end of the testing 

    call list_init(hist_list_var,'hr',id1, itype=var_flux_type)          ; call list_init(hist_list_unit,'gC m-2 s-1',id2)
    call list_insert(hist_list_var,'f_n2o_nit',id1, itype=var_flux_type) ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'f_denit',id1, itype=var_flux_type)   ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'f_nit',id1, itype=var_flux_type)     ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'co2_soil_flx',id1,itype=var_flux_type); call list_insert(hist_list_unit,'gC m-2 s-1',id2)
    call list_insert(hist_list_var,'nh3_soil_flx',id1,itype=var_flux_type); call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'cwdc',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totlitc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totsomc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totlitc_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totsomc_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'cwdn',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totlitn',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totsomn',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totlitn_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totsomn_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'cwdp',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totlitp',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totsomp',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totlitp_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totsomp_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'smin_nh4',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'smin_no3',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'smin_psol',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'som1c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som2c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som3c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som1n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som2n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som3n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som1p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'som2p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'som3p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    this%nhistvars=id1
  elseif(index(trim(reaction_method),'cdom')/=0)then

    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'in if case cdom, before all calls'
    write(stdout, *) '***************************'
    ! end of the testing 

    call list_init(hist_list_var,'hr',id1, itype=var_flux_type)          ; call list_init(hist_list_unit,'gC m-2 s-1',id2)
    call list_insert(hist_list_var,'f_n2o_nit',id1, itype=var_flux_type) ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'f_denit',id1, itype=var_flux_type)   ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'f_nit',id1, itype=var_flux_type)     ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'cwdc',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totlitc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totsomc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totlitc_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totsomc_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'cwdn',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totlitn',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totsomn',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totlitn_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totsomn_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'cwdp',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totlitp',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totsomp',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totlitp_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totsomp_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'smin_nh4',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'smin_no3',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som1c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som2c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som3c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som1n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som2n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som3n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som1p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gPN m-2',id2)
    call list_insert(hist_list_var,'som2p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'som3p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'domc',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'domn',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'domp',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    this%nhistvars=id1
  elseif(index(trim(reaction_method),'simic')/=0)then

    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'in if case simic, before all calls'
    write(stdout, *) '***************************'
    ! end of the testing 

    call list_init(hist_list_var,'hr',id1, itype=var_flux_type)          ; call list_init(hist_list_unit,'gC m-2 s-1',id2)
    call list_insert(hist_list_var,'cwdc',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totlitc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'DOC',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'micb_live',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'micb_dead',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'POM_C',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totsomc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    this%nhistvars=id1
  elseif(index(trim(reaction_method),'summs')/=0 )then
    ! adding reaction_method 'summs'     -zlyu
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'in if case summs, before all calls'
    write(stdout, *) '***************************'
    ! end of the testing 

    call list_init(hist_list_var,'hr',id1, itype=var_flux_type)          ; call list_init(hist_list_unit,'gC m-2 s-1',id2)

    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'in if case summs, after list_init'
    write(stdout, *) '***************************'
    ! end of the testing 

    call list_insert(hist_list_var,'f_n2o_nit',id1, itype=var_flux_type) ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'f_denit',id1, itype=var_flux_type)   ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'f_nit',id1, itype=var_flux_type)     ; call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'co2_soil_flx',id1,itype=var_flux_type); call list_insert(hist_list_unit,'gC m-2 s-1',id2)
    call list_insert(hist_list_var,'nh3_soil_flx',id1,itype=var_flux_type); call list_insert(hist_list_unit,'gN m-2 s-1',id2)
    call list_insert(hist_list_var,'cwdc',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totlitc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totsomc',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totlitc_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'totsomc_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'cwdn',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totlitn',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totsomn',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totlitn_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'totsomn_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'cwdp',id1, itype=var_state_type)      ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totlitp',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totsomp',id1, itype=var_state_type)   ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totlitp_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'totsomp_1m',id1, itype=var_state_type); call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'smin_nh4',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'smin_no3',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'smin_psol',id1, itype=var_state_type)  ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'som1c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som2c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som3c',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gC m-2',id2)
    call list_insert(hist_list_var,'som1n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som2n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som3n',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gN m-2',id2)
    call list_insert(hist_list_var,'som1p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'som2p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    call list_insert(hist_list_var,'som3p',id1, itype=var_state_type)     ; call list_insert(hist_list_unit,'gP m-2',id2)
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'in if case summs, after list_insert'
    write(stdout, *) '***************************'
    ! end of the testing 

    this%nhistvars=id1
    ! end of adding for reaction_method 'summs'
  endif

    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after if cases'
    write(stdout, *) '***************************'
    ! end of the testing 

  allocate(this%varl(this%nhistvars)); 
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after this%varl'
    write(stdout, *) '***************************'
    ! end of the testing 
  allocate(this%unitl(this%nhistvars))
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after this%unitl'
    write(stdout, *) '***************************'
    ! end of the testing 
  allocate(this%vartypes(this%nhistvars))
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after this%vartypes'
    write(stdout, *) '***************************'
    ! end of the testing 
  call copy_name(this%nhistvars, hist_list_var, this%varl)
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after copy_name_1'
    write(stdout, *) '***************************'
    ! end of the testing 
  call copy_name_type(this%nhistvars, hist_list_var, this%vartypes)
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after copy_name_type'
    write(stdout, *) '***************************'
    ! end of the testing 
  call copy_name(this%nhistvars, hist_list_unit, this%unitl)
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after copy_name_2'
    write(stdout, *) '***************************'
    ! end of the testing 
  call list_free(hist_list_var); call list_free(hist_list_unit)
    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'after list_free'
    write(stdout, *) '***************************'
    ! end of the testing 
  end subroutine Init
  !-----------------------------------------------------------------------
  function getvarllen(this)result(ans)

  use betr_constants           , only : stdout                  ! added for output testing remarks. -zlyu 01/27/2019

  implicit none
  class(hist_bgc_type), intent(inout) :: this
    integer :: ans

    ! testing only, where the run crushed        -zlyu   01/27/2019    
    write(stdout, *) '***************************'
    write(stdout, *) 'inside getvarllen subroutine'
    write(stdout, *) '***************************'
    ! end of the testing 

  ans = this%nhistvars
  end function getvarllen



end module HistbgcMod
