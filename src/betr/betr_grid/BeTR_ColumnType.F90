module BeTR_ColumnType
  !

  ! Column type for data transfer between lsm and BeTR
  ! !PUBLIC TYPES:
  use bshr_kind_mod   , only : r8 => shr_kind_r8
  implicit none
  save
  private

  type, public :: betr_column_type
     ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer , pointer :: landunit (:)   ! index into landunit level quantities
   integer , pointer :: gridcell (:)   ! index into gridcell level quantities

   integer , pointer :: snl      (:)   ! number of snow layers
   real(r8), pointer :: dz       (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: zi       (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
   real(r8), pointer :: z        (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)

  end type betr_column_type

  type(betr_column_type), public :: betr_col  !column data structure (soil/snow/canopy columns)

end module BeTR_ColumnType
