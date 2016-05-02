module BeTR_PatchType
  use bshr_kind_mod, only: r8 => shr_kind_r8
implicit none
  type, public :: betr_patch_type

   real(r8), pointer :: wtcol(:)
   integer,  pointer :: column(:)
   integer,  pointer :: itype(:)
   integer,  pointer :: landunit(:)
   real(r8), pointer :: crop(:)
  end type betr_patch_type

  type(betr_patch_type), public :: betr_pft

end module BeTR_PatchType
