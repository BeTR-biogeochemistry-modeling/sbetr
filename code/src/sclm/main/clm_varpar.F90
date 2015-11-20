module clm_varpar
  ! !PUBLIC TYPES:
  implicit none
  save

  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
  integer  :: nlevtrc_soil = 10
  integer  :: nlevsoi = 10
  integer  :: nlevgrnd = 15
  integer  :: maxpatch_glcmec= -1    !some dumb number for stand alone betr


 logical :: crop_prog   = .false.  ! If prognostic crops is turned on

  integer :: ndecomp_pools
  integer :: ndecomp_cascade_transitions
  integer :: i_met_lit, i_cel_lit, i_lig_lit, i_cwd

end module clm_varpar
