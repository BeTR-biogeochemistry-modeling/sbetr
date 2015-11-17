module clm_varpar
  ! !PUBLIC TYPES:
  implicit none
  save

  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
  integer  :: nlevtrc_soil = 10
  integer  :: nlevsoi = 10
  integer  :: nlevgrnd = 15
end module clm_varpar