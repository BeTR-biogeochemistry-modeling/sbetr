module BeTR_nitrogenfluxRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  type, public :: betr_nitrogenflux_recv_type
    real(r8), pointer :: f_denit_vr_col(:,:)  => null()
    real(r8), pointer :: f_nit_vr_col(:,:)
    real(r8), pointer :: f_n2o_nit_vr(:,:)

    real(r8), pointer :: smin_no3_to_plant_patch(:)
    real(r8), pointer :: smin_nh4_to_plant_patch(:)
    real(r8), pointer :: fire_decomp_nloss_col(:) => null()

    real(r8), pointer :: f_nit_col(:)
    real(r8), pointer :: f_denit_col(:)
    real(r8), pointer :: f_n2o_nit_col(:)
    real(r8), pointer :: smin_no3_leached_col(:)
    real(r8), pointer :: smin_no3_runoff_col(:)
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
  end type betr_nitrogenflux_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_nitrogenflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_nitrogenflux_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%f_denit_vr_col(begc:endc,lbj:ubj))
  allocate(this%f_nit_vr_col(begc:endc,lbj:ubj))
  allocate(this%f_n2o_nit_vr(begc:endc,lbj:ubj))

  allocate(this%smin_no3_to_plant_patch(begp:endp))
  allocate(this%smin_nh4_to_plant_patch(begp:endp))

  allocate(this%f_nit_col(begc:endc))
  allocate(this%f_denit_col(begc:endc))
  allocate(this%f_n2o_nit_col(begc:endc))
  allocate(this%smin_no3_leached_col(begc:endc))
  allocate(this%smin_no3_runoff_col(begc:endc))
  allocate(this%fire_decomp_nloss_col(begc:endc))
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_nitrogenflux_recv_type)  :: this
  real(r8), intent(in) :: value_column


  this%f_denit_vr_col(:,:) = value_column
  this%f_nit_vr_col(:,:) = value_column
  this%f_n2o_nit_vr(:,:) = value_column

  this%smin_no3_to_plant_patch(:) = value_column
  this%smin_nh4_to_plant_patch(:) = value_column

  this%f_nit_col(:) = value_column
  this%f_denit_col(:) = value_column
  this%f_n2o_nit_col(:) = value_column
  this%smin_no3_leached_col(:) = value_column
  this%smin_no3_runoff_col(:) = value_column
  this%fire_decomp_nloss_col(:) = value_column
  end subroutine reset
end module BeTR_nitrogenfluxRecvType
