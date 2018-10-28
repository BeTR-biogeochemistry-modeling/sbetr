module fakeBGCType
#include "bshr_assert.h"
  !
  ! !DESCRIPTION:

  ! !USES:
  use bshr_kind_mod             , only : r8 => shr_kind_r8
  use bshr_log_mod              , only : errMsg => shr_log_errMsg
  use betr_varcon               , only : spval => bspval
  use betr_ctrl                 , only : spinup_state => betr_spinup_state
  use gbetrType                 , only : gbetr_type
  use fakeParaType             , only : fake_para_type
  use BetrStatusType            , only : betr_status_type
  use fakeBGCIndexType            , only : fake_bgc_index_type
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: fake_bgc_type
  type(fake_bgc_index_type),     private :: fake_index

    !declare parameters below

  contains
    procedure, public  :: init          => init_fake
    procedure, public  :: runbgc        => runbgc_fake
    procedure, public  :: init_cold     => init_cold_fake
  end type fake_bgc_type
contains

  !-------------------------------------------------------------------------------

  subroutine init_fake(this, batch_mode,  bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(fake_bgc_type) , intent(inout) :: this
  logical, optional, intent(in) :: batch_mode
  type(betr_status_type)      , intent(out) :: bstatus
  !local variables
  character(len=256) :: msg

  end subroutine init_fake
  !-------------------------------------------------------------------------------
  subroutine runbgc_fake(this,  is_surflit, dtime, bgc_forc, nstates, ystates0, ystatesf, bstatus)

  !DESCRIPTION
  !do bgc model integration for one step
  use JarBgcForcType        , only : JarBGC_forc_type
  use BetrStatusType        , only : betr_status_type
  use tracer_varcon         , only : catomw, natomw, patomw
  implicit none
  class(fake_bgc_type)  , intent(inout) :: this
  logical                    , intent(in)    :: is_surflit
  real(r8)                   , intent(in)    :: dtime
  type(JarBGC_forc_type)     , intent(in)    :: bgc_forc
  integer                    , intent(in)    :: nstates
  real(r8)                   , intent(out)   :: ystates0(nstates)
  real(r8)                   , intent(out)   :: ystatesf(nstates)
  type(betr_status_type)     , intent(out)   :: bstatus

  !local variables
  real(r8)               :: time = 0._r8
  character(len=*),parameter :: subname = 'runbgc_fake'

  end subroutine runbgc_fake
  !-------------------------------------------------------------------------------
  subroutine init_cold_fake(this, nstvars, ystates)
  !
  !DESCRPTION
  !do a cold state initialization for batch mode simulation
  implicit none
  class(fake_bgc_type)     , intent(inout) :: this
  integer                   , intent(in)    :: nstvars
  real(r8)                  , intent(inout) :: ystates(nstvars)

  !Initialize necessary state variables below
  end subroutine init_cold_fake
  end module fakeBGCType
