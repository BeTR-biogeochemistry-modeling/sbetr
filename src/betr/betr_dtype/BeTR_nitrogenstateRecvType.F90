module BeTR_nitrogenstateRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
  use tracer_varcon  , only : reaction_method
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_nitrogenstate_recv_type
    real(r8), pointer :: cwdn_col(:) => null()
    real(r8), pointer :: totlitn_col(:) => null()
    real(r8), pointer :: totsomn_col(:) => null()
    real(r8), pointer :: sminn_col(:) => null()
    real(r8), pointer :: sminn_nh4_col(:)=> null()
    real(r8), pointer :: sminn_no3_col(:)=> null()
    real(r8), pointer :: totlitn_1m_col(:) => null()
    real(r8), pointer :: totsomn_1m_col(:) => null()
    real(r8), pointer :: cwdn_vr_col(:,:) => null()
    real(r8), pointer :: totlitn_vr_col(:,:) => null()
    real(r8), pointer :: totsomn_vr_col(:,:) => null()
    real(r8), pointer :: sminn_vr_col(:,:) => null()
    real(r8), pointer :: sminn_nh4_vr_col(:,:) => null()
    real(r8), pointer :: sminn_no3_vr_col(:,:) => null()
    real(r8), pointer :: som1n_col(:) => null()
    real(r8), pointer :: som2n_col(:) => null()
    real(r8), pointer :: som3n_col(:) => null()
    real(r8), pointer :: som1n_vr_col(:,:) => null()
    real(r8), pointer :: som2n_vr_col(:,:) => null()
    real(r8), pointer :: som3n_vr_col(:,:) => null()
    real(r8), pointer :: domn_vr_col(:,:) => null()
    real(r8), pointer :: polyn_col(:) => null()
    real(r8), pointer :: monoenzn_col(:) => null()
    real(r8), pointer :: micresn_col(:) => null()
    real(r8), pointer :: polyn_vr_col(:,:) => null()
    real(r8), pointer :: monoenzn_vr_col(:,:) => null()
    real(r8), pointer :: micresn_vr_col(:,:) => null()

  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_nitrogenstate_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_nitrogenstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  implicit none
  class(betr_nitrogenstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%cwdn_col(begc:endc))
  allocate(this%totlitn_col(begc:endc))
  allocate(this%totsomn_col(begc:endc))
  allocate(this%sminn_col(begc:endc))
  allocate(this%sminn_nh4_col(begc:endc))
  allocate(this%sminn_no3_col(begc:endc))
  allocate(this%totlitn_1m_col(begc:endc))
  allocate(this%totsomn_1m_col(begc:endc))
  if(index(reaction_method,'ecacnp')/=0)then
    allocate(this%som1n_col(begc:endc));  this%som1n_col(:) = nan
    allocate(this%som2n_col(begc:endc));  this%som2n_col(:) = nan
    allocate(this%som3n_col(begc:endc));  this%som3n_col(:) = nan
    allocate(this%som1n_vr_col(begc:endc, lbj:ubj));  this%som1n_vr_col(:,:) = nan
    allocate(this%som2n_vr_col(begc:endc, lbj:ubj));  this%som2n_vr_col(:,:) = nan
    allocate(this%som3n_vr_col(begc:endc, lbj:ubj));  this%som3n_vr_col(:,:) = nan
  endif
  if(index(reaction_method,'summs')/=0)then
    allocate(this%polyn_col(begc:endc));  this%polyn_col(:) = nan
    allocate(this%monoenzn_col(begc:endc));  this%monoenzn_col(:) = nan
    allocate(this%micresn_col(begc:endc));   this%micresn_col(:) = nan
    allocate(this%polyn_vr_col(begc:endc, lbj:ubj));  this%polyn_vr_col(:,:) = nan
    allocate(this%monoenzn_vr_col(begc:endc, lbj:ubj));  this%monoenzn_vr_col(:,:) = nan
    allocate(this%micresn_vr_col(begc:endc, lbj:ubj));   this%micresn_vr_col(:,:) = nan
  endif
  allocate(this%cwdn_vr_col(begc:endc,lbj:ubj)); this%cwdn_vr_col(:,:) = nan
  allocate(this%totlitn_vr_col(begc:endc,lbj:ubj)); this%totlitn_vr_col(:,:)=nan
  allocate(this%totsomn_vr_col(begc:endc,lbj:ubj)); this%totsomn_vr_col(:,:) =nan
  allocate(this%sminn_vr_col(begc:endc,lbj:ubj)); this%sminn_vr_col(:,:) =nan
  allocate(this%sminn_nh4_vr_col(begc:endc,lbj:ubj)); this%sminn_nh4_vr_col(:,:) =nan
  allocate(this%sminn_no3_vr_col(begc:endc,lbj:ubj)); this%sminn_no3_vr_col(:,:) =nan
  allocate(this%domn_vr_col(begc:endc, lbj:ubj)); this%domn_vr_col(:,:)=nan
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_nitrogenstate_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%cwdn_vr_col(:,:) = value_column
  this%totlitn_vr_col(:,:) = value_column
  this%totsomn_vr_col(:,:) = value_column
  this%sminn_vr_col(:,:) = value_column
  this%sminn_nh4_vr_col(:,:) = value_column
  this%sminn_no3_vr_col(:,:) = value_column

  if(index(reaction_method,'ecacnp')/=0)then
    this%som1n_vr_col(:,:) = value_column
    this%som2n_vr_col(:,:) = value_column
    this%som3n_vr_col(:,:) = value_column
  endif
  if(index(reaction_method,'summs')/=0)then
    this%polyn_vr_col(:,:) = value_column
    this%monoenzn_vr_col(:,:) = value_column
    this%micresn_vr_col(:,:) = value_column
  endif
  this%domn_vr_col(:,:)=value_column
  end subroutine reset

  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz, zs)
  use tracer_varcon, only : natomw
  implicit none
  class(betr_nitrogenstate_recv_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)
  real(r8), intent(in) :: zs(bounds%begc:bounds%endc,lbj:ubj)


  integer :: c, j

  this%cwdn_col(:) = 0._r8
  this%totlitn_col(:) = 0._r8
  this%totsomn_col(:) = 0._r8
  this%sminn_col(:) = 0._r8
  this%sminn_nh4_col(:) = 0._r8
  this%sminn_no3_col(:) = 0._r8
  this%totlitn_1m_col(:) = 0._r8
  this%totsomn_1m_col(:) = 0._r8

  if(index(reaction_method,'ecacnp')/=0)then
    this%som1n_col(:) = 0._r8
    this%som2n_col(:) = 0._r8
    this%som3n_col(:) = 0._r8
    do j = lbj, ubj
      do c = bounds%begc, bounds%endc
        this%som1n_col(c) =   this%som1n_col(c) + dz(c,j)*this%som1n_vr_col(c,j)
        this%som2n_col(c) =   this%som2n_col(c) + dz(c,j)*this%som2n_vr_col(c,j)
        this%som3n_col(c) =   this%som3n_col(c) + dz(c,j)*this%som3n_vr_col(c,j)
      enddo
    enddo
  endif
  if(index(reaction_method,'summs')/=0)then
    this%polyn_col(:) = 0.0_r8
    this%monoenzn_col(:) = 0.0_r8
    this%micresn_col(:) = 0.0_r8
    do j = lbj, ubj
      do c = bounds%begc, bounds%endc
        this%polyn_col(c) = this%polyn_col(c) + dz(c,j)*this%polyn_vr_col(c,j)
        this%monoenzn_col(c) = this%monoenzn_col(c) + dz(c,j)*this%monoenzn_vr_col(c,j)
        this%micresn_col(c) = this%micresn_col(c) + dz(c,j)*this%micresn_vr_col(c,j)
      enddo
    enddo
  endif
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      this%cwdn_col(c) = this%cwdn_col(c) + dz(c,j) * this%cwdn_vr_col(c,j)
      this%totlitn_col(c) = this%totlitn_col(c) + dz(c,j)*this%totlitn_vr_col(c,j)
      this%totsomn_col(c) = this%totsomn_col(c) + dz(c,j)*this%totsomn_vr_col(c,j)
      this%sminn_col(c) = this%sminn_col(c) + dz(c,j)*this%sminn_vr_col(c,j)
      this%sminn_nh4_col(c) = this%sminn_nh4_col(c) + dz(c,j)*this%sminn_nh4_vr_col(c,j)
      this%sminn_no3_col(c) = this%sminn_no3_col(c) + dz(c,j)*this%sminn_no3_vr_col(c,j)

      if(zs(c,j)<1._r8)then
        if(zs(c,j+1)>1._r8)then
          this%totlitn_1m_col(c) = this%totlitn_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totlitn_vr_col(c,j)
          this%totsomn_1m_col(c) = this%totsomn_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totsomn_vr_col(c,j)
        else
          this%totlitn_1m_col(c) = this%totlitn_1m_col(c) + dz(c,j)*this%totlitn_vr_col(c,j)
          this%totsomn_1m_col(c) = this%totsomn_1m_col(c) + dz(c,j)*this%totsomn_vr_col(c,j)
        endif
      endif
    enddo
  enddo
  end subroutine summary
end module BeTR_nitrogenstateRecvType
