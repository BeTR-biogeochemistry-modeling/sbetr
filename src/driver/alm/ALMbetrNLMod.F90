module ALMBeTRNLMod

  use shr_log_mod , only: errMsg => shr_log_errMsg
  use betr_constants , only : betr_namelist_buffer_size
  use abortutils , only: endrun
implicit none



  character(len=betr_namelist_buffer_size), public :: betr_namelist_buffer
  public :: betr_readNL
contains


  !-------------------------------------------------------------------------------
  subroutine betr_readNL(NLFilename)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use betr_utils    , only : log2str
    use clm_varctl    , only : iulog
    use betr_constants, only : betr_string_length
    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename              ! Namelist filename
                                                            !
                                                            ! !LOCAL VARIABLES:
    integer                      :: ierr                    ! error code
    integer                      :: unitn                   ! unit for namelist file
    character(len=32)            :: subname = 'betr_readNL' ! subroutine name
    !-----------------------------------------------------------------------

    character(len=betr_string_length)      :: reaction_method
    logical                                :: advection_on, diffusion_on, reaction_on, ebullition_on
    character(len=1), parameter  :: quote = ''''
    namelist / betr_inparm / reaction_method, &
      advection_on, diffusion_on, reaction_on, ebullition_on
    logical :: esm_on
    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    reaction_method = ''
    advection_on    = .true.
    diffusion_on    = .true.
    reaction_on     = .true.
    ebullition_on   = .true.
    esm_on          = .true.
    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in betr_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_CanopyHydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, betr_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading betr_inparm namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if
    ! Broadcast namelist variables read in
    call shr_mpi_bcast(reaction_method, mpicom)
    call shr_mpi_bcast(advection_on, mpicom)
    call shr_mpi_bcast(diffusion_on, mpicom)
    call shr_mpi_bcast(reaction_on, mpicom)
    call shr_mpi_bcast(ebullition_on, mpicom)

    write(betr_namelist_buffer,*) '&betr_parameters'//new_line('A'), &
      ' reaction_method='//quote//trim(reaction_method)//quote//new_line('A'), &
      ' esm_on=',trim(log2str(esm_on)),new_line('A'),&
      ' advection_on=',trim(log2str(advection_on)),new_line('A'), &
      ' diffusion_on=',trim(log2str(diffusion_on)),new_line('A'), &
      ' reaction_on=',trim(log2str(reaction_on)),new_line('A'), &
      ' ebullition_on=',trim(log2str(ebullition_on)),new_line('A')//'/'

  end subroutine betr_readNL



end module ALMBeTRNLMod
