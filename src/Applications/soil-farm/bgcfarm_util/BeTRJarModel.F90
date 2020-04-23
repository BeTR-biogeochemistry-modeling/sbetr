module BeTRJarModel

  use gbetrType              , only : gbetr_type
  use BetrStatusType         , only : betr_status_type
  use BiogeoConType          , only : BiogeoCon_type
  use bshr_kind_mod          , only : r8 => shr_kind_r8
implicit none


  type , extends(gbetr_type), public :: jar_model_type
   character(len=24) :: jarname
  contains
    procedure, public  :: init          => JarModel_init
    procedure, public  :: runbgc        => JarModel_runbgc
    procedure, public  :: UpdateParas   => JarModel_UpdateParas
    procedure, public  :: sumup_cnp_msflx => JarModel_sumup_cnp_msflx
    procedure, public  :: getvarllen    => JarModel_getvarllen
    procedure, public  :: getvarlist    => JarModel_getvarlist
    procedure, public  :: init_cold     => JarModel_init_cold
  end type jar_model_type

  contains

  !-------------------------------------------------------------------------------
  subroutine JarModel_init(this,  biogeo_con,  batch_mode, bstatus)

  implicit none
  class(jar_model_type)      , intent(inout) :: this
  class(BiogeoCon_type)      , intent(in)    :: biogeo_con
  logical                    , intent(in)    :: batch_mode
  type(betr_status_type)     , intent(out)   :: bstatus

  call bstatus%reset()
  end subroutine JarModel_init
  !-------------------------------------------------------------------------------
  subroutine JarModel_getvarlist(this, nstvars, varnames, varunits, vartypes)
  implicit none
  class(jar_model_type)      , intent(inout) :: this
  integer, intent(in) :: nstvars
  character(len=*), intent(out) :: varnames(1:nstvars)
  character(len=*), intent(out) :: varunits(1:nstvars)
  integer         , intent(out) :: vartypes(1:nstvars)

  varnames(:)=''
  varunits(:)=''
  vartypes(:)=0
  if(nstvars>=0)return
  end subroutine JarModel_getvarlist

  !-------------------------------------------------------------------------------
  function JarModel_getvarllen(this)result(ans)

  implicit none
  class(jar_model_type)      , intent(inout) :: this
  integer :: ans

  ans = 0

  end function JarModel_getvarllen
  !-------------------------------------------------------------------------------
  subroutine JarModel_runbgc(this,  is_surflit, dtime, bgc_forc, nstates, ystates0, ystatesf, bstatus)

  use JarBgcForcType            , only : JarBGC_forc_type
  implicit none
  class(jar_model_type)      , intent(inout) :: this
  logical                    , intent(in)    :: is_surflit
  real(r8)                   , intent(in)    :: dtime
  type(JarBGC_forc_type)     , intent(in)    :: bgc_forc
  integer                    , intent(in)    :: nstates
  real(r8)                   , intent(out)   :: ystates0(nstates)
  real(r8)                   , intent(out)   :: ystatesf(nstates)
  type(betr_status_type)     , intent(out)   :: bstatus

  ystates0(:)=0._r8
  ystatesf(:)=0._r8
  call bstatus%reset()
  end subroutine JarModel_runbgc

  !-------------------------------------------------------------------------------
  subroutine JarModel_UpdateParas(this,  biogeo_con, bstatus)
  implicit none
  class(jar_model_type)      , intent(inout) :: this
  class(BiogeoCon_type)       , intent(in)    :: biogeo_con
  type(betr_status_type)     , intent(out)   :: bstatus

  call bstatus%reset()
  end subroutine JarModel_UpdateParas
  !-------------------------------------------------------------------------------

  subroutine JarModel_sumup_cnp_msflx(this, ystates1, c_mass, n_mass, p_mass,c_flx,n_flx,p_flx, bstatus)

  implicit none
  class(jar_model_type), intent(in) :: this
  real(r8), intent(in)  :: ystates1(:)
  real(r8), intent(out) :: c_mass, n_mass, p_mass
  real(r8), optional, intent(out) :: c_flx, n_flx, p_flx
  type(betr_status_type), optional, intent(out)   :: bstatus
  !local variables
  c_mass=0._r8;n_mass=0._r8;p_mass=0._r8
  if(present(c_flx))c_flx=0._r8
  if(present(n_flx))n_flx=0._r8
  if(present(p_flx))p_flx=0._r8

  call bstatus%reset()
  end subroutine JarModel_sumup_cnp_msflx


  !-------------------------------------------------------------------------------
  subroutine JarModel_init_cold(this, nstvars, ystates)

  implicit none
  class(jar_model_type)      , intent(inout) :: this
  integer , intent(in) :: nstvars
  real(r8), intent(inout) :: ystates(nstvars)
  ystates(:)=0._r8
  if(nstvars>=0)continue
  end subroutine JarModel_init_cold
end module BeTRJarModel
