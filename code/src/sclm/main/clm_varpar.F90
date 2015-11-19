module clm_varpar
  ! !PUBLIC TYPES:
  implicit none
  save

  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
  integer  :: nlevtrc_soil = 10
  integer  :: nlevsoi = 10
  integer  :: nlevgrnd = 15
  integer  :: maxpatch_glcmec= -1    !some dumb number for stand alone betr
end module clm_varpar
