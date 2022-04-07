module ELMBeTRNLMod

  use shr_log_mod , only: errMsg => shr_log_errMsg
  use betr_constants , only : betr_namelist_buffer_size, betr_namelist_buffer_size_ext
  use abortutils , only: endrun
  use betr_constants, only : betr_string_length
implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  character(len=betr_namelist_buffer_size), public :: betr_namelist_buffer
  public :: betr_readNL
  public :: set_betr_cnpbgc
  public :: do_betr_bgc_type
contains


!-------------------------------------------------------------------------------
  function do_betr_bgc_type(type_char)result(ans)
  use betr_ctrl, only : bgc_type
  implicit none
  character(len=*), intent(in) :: type_char


  logical :: ans

  ans = index(bgc_type, type_char)/=0

  end function do_betr_bgc_type



  !-------------------------------------------------------------------------------
  subroutine betr_readNL(NLFilename, use_c13, use_c14, nsoilorder, lbgcalib)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use betr_utils    , only : log2str
    use elm_varctl    , only : iulog, spinup_state
    use tracer_varcon  , only : advection_on, diffusion_on, reaction_on, ebullition_on, reaction_method
    use ApplicationsFactory, only : AppLoadParameters
    use betr_varcon   , only : betr_maxpatch_pft, betr_max_soilorder
    use betr_ctrl     , only : betr_spinup_state,bgc_type
    use tracer_varcon , only : advection_on, diffusion_on, reaction_on, ebullition_on, reaction_method
    use tracer_varcon , only : AA_spinup_on, fix_ip, do_bgc_calibration, bgc_param_file
    use ApplicationsFactory, only : AppInitParameters
    use tracer_varcon , only : use_c13_betr, use_c14_betr
    use BetrStatusType  , only : betr_status_type
    use tracer_varcon, only : lbcalib
    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename              ! Namelist filename
    logical,          intent(in) :: use_c13
    logical,          intent(in) :: use_c14
    integer,          intent(in) :: nsoilorder
    logical, optional,intent(in) :: lbgcalib
                                                            !
                                                            ! !LOCAL VARIABLES:
    integer                      :: ierr                    ! error code
    integer                      :: unitn                   ! unit for namelist file
    character(len=32)            :: subname = 'betr_readNL' ! subroutine name
    type(betr_status_type)       :: bstatus
    !-----------------------------------------------------------------------

    character(len=255):: AppParNLFile
    character(len=1), parameter  :: quote = ''''
    namelist / betr_inparm / reaction_method, &
      advection_on, diffusion_on, reaction_on, ebullition_on, &
      AppParNLFile, AA_spinup_on, fix_ip, do_bgc_calibration, &
      bgc_param_file

    logical :: appfile_on

    if(present(lbgcalib))then
      lbcalib = lbgcalib
    else
      lbcalib = .false.
    endif
    !initialize spinup state
    if(index(bgc_type,'type2_bgc')/=0)then
      betr_spinup_state =spinup_state
    endif
    betr_max_soilorder=nsoilorder
    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    reaction_method = 'mock_run'
    advection_on    = .true.
    diffusion_on    = .true.
    reaction_on     = .true.
    ebullition_on   = .true.
    AA_spinup_on    = .false.
    use_c13_betr    = use_c13
    use_c14_betr    = use_c14
    AppParNLFile    = ''
!    appfile_on      = .false.
    fix_ip          = .false.
    do_bgc_calibration=.false.
    bgc_param_file  = ''
    if ( masterproc )then
       unitn = getavu()
       write(iulog,*) 'Read in betr_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'betr_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, betr_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading betr_inparm namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )
!       if(trim(AppParNLFile)/='')then
!         appfile_on=.true.
!         call LoadFile2String(AppParNLFile, bgc_namelist_buffer)
!       else
!         bgc_namelist_buffer='none'
!       endif
    end if

!    call shr_mpi_bcast(appfile_on, mpicom)
    !pass parameters to all files
!    call shr_mpi_bcast(bgc_namelist_buffer, mpicom)

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(reaction_method, mpicom)
    call shr_mpi_bcast(advection_on, mpicom)
    call shr_mpi_bcast(diffusion_on, mpicom)
    call shr_mpi_bcast(reaction_on, mpicom)
    call shr_mpi_bcast(ebullition_on, mpicom)
    call shr_mpi_bcast(AA_spinup_on, mpicom)
    call shr_mpi_bcast(fix_ip, mpicom)
    call shr_mpi_bcast(do_bgc_calibration, mpicom)
    call shr_mpi_bcast(bgc_param_file, mpicom)

    if(masterproc)then
      write(iulog,*)'&betr_parameters'
      write(iulog,*)'reaction_method=',trim(reaction_method)
      write(iulog,*)'advection_on   =',advection_on
      write(iulog,*)'diffusion_on   =',diffusion_on
      write(iulog,*)'reaction_on    =',reaction_on
      write(iulog,*)'ebullition_on  =',ebullition_on
      write(iulog,*)'AA_spinup_on   =',AA_spinup_on
      write(iulog,*)'fix_ip         =',fix_ip
      write(iulog,*)'do_bgc_calibration=',do_bgc_calibration
    endif

    call AppInitParameters(reaction_method, bstatus)
    if(bstatus%check_status())call endrun(msg=bstatus%print_msg())

  end subroutine betr_readNL

  !-------------------------------------------------------------------------------
  subroutine LoadFile2String(AppParNLFile, bgc_namelist_buffer)

  !
  !DESCRITION
  !turn a namelist file into a big and long string. Now
  !logical data types are not supported.
  implicit none
  character(len=*), intent(in) :: AppParNLFile
  character(len=betr_namelist_buffer_size_ext), intent(out) :: bgc_namelist_buffer

  logical :: exist_file
  character(len=255) :: lineread
  logical :: start
  integer :: ios

  inquire(file=AppParNLFile, exist=exist_file)
  if(.not. exist_file)then
    call endrun(msg='file '//trim(AppParNLFile)//' does not exist '//errmsg(__FILE__, __LINE__))
  endif

  open(unit=999, file=trim(AppParNLFile),status='old',action='read',form='formatted')
  start=.true.
  do
    read(999,'(a)',IOSTAT=ios)lineread
    if(ios>0)then
      call endrun(msg='error reading file '//trim(AppParNLFile)//errmsg(__FILE__, __LINE__))
    elseif(ios==0)then
      if(start)then
        bgc_namelist_buffer=trim(lineread)//new_line('A')
        start=.false.
      else
        bgc_namelist_buffer=trim(bgc_namelist_buffer)//trim(lineread)//new_line('A')
      endif
    else
      exit
    endif
  enddo

  close(999)

  end subroutine LoadFile2String

!-------------------------------------------------------------------------------
    subroutine set_betr_cnpbgc(suplnitro,suplphos, spinup_state)
    !
    !DESCRIPTION
    !set n and p switches of the betr bgc model
    use tracer_varcon, only : is_nitrogen_active, is_phosphorus_active
    implicit none
    character(len=*), intent(in) :: suplnitro
    character(len=*), intent(in) :: suplphos
    integer         , intent(in) :: spinup_state
    integer :: cnpset

    !set
    cnpset=111
    if(trim(suplnitro)=='ALL')cnpset=cnpset-10
    if(trim(suplphos)=='ALL')cnpset=cnpset-1

    select case (cnpset)
    case (100)
      is_nitrogen_active = .false.
      is_phosphorus_active=.false.
    case (110)
      is_nitrogen_active = .true.
      is_phosphorus_active=.false.
    case (101)
      is_nitrogen_active = .false.
      is_phosphorus_active=.true.
    case default
      is_nitrogen_active = .true.
      is_phosphorus_active=.true.
    end select
    !make sure P has full supply during spinup
    !if(spinup_state==1)is_phosphorus_active=.false.
    end subroutine set_betr_cnpbgc
end module ELMBeTRNLMod
