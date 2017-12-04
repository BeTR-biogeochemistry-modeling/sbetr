module BgcSummsDecompType

!
! DESCRIPTIONS
! contains data structure for ECA-century decomposition.

! USES

  use bshr_kind_mod  , only : r8 => shr_kind_r8
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: DecompSumms_type
  real(r8) :: o_scalar       ! fraction by which decomposition is limited by anoxia
  real(r8) :: w_scalar       ! fraction by which decomposition is limited by h2osoi_liqure availability
  real(r8) :: t_scalar       ! fraction by which litter decomposition is limited by temperature
  real(r8) :: depth_scalar   ! depth dependent factor for heteorotrophic respiration

  !parameters
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding

  real(r8) :: vmax_mic
  real(r8) :: vmax_enz
  real(r8) :: kaff_mono_mic
  real(r8) :: kaff_enz_poly
  real(r8) :: mr_mic
  real(r8) :: kappa_mic
  real(r8) :: kaff_mono_msurf
  real(r8) :: kaff_enz_msurf
  real(r8) :: ea_vmax_mic
  real(r8) :: ea_vmax_enz
  real(r8) :: ea_kaff_mono_mic
  real(r8) :: ea_kaff_enz_poly
  real(r8) :: ea_mr_mic
  real(r8) :: ea_kappa_mic
  real(r8) :: ea_kaff_mono_msurf
  real(r8) :: ea_kaff_enz_msurf
  real(r8) :: ref_vmax_mic
  real(r8) :: ref_vmax_enz
  real(r8) :: ref_kaff_mono_mic
  real(r8) :: ref_kaff_enz_poly
  real(r8) :: ref_mr_mic
  real(r8) :: ref_kappa_mic
  real(r8) :: ref_kaff_mono_msurf
  real(r8) :: ref_kaff_enz_msurf

  contains
    procedure, public  :: Init
    procedure, public  :: set_decompk_scalar
    procedure, private :: InitCold
    procedure, private :: InitAllocate
    procedure, private :: initPar
  end type DecompSumms_type

 contains

 !------------------------------------------------------------------------
 subroutine Init(this, bgcsumms_con)

  use BgCConSummsType, only : BgcConSumms_type
  implicit none
  class(DecompSumms_type), intent(inout) :: this
  type(BgcConSumms_type),intent(in) :: bgcsumms_con

  call this%InitAllocate ()

  call this%InitCold ()

  call this%InitPar(bgcsumms_con)

 end subroutine Init
 !------------------------------------------------------------------------
 subroutine InitAllocate(this)
   !
   ! !DESCRIPTION:
   ! Initialize module data structure
   !
   ! !USES:
   use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
   use betr_varcon    , only : spval  => bspval
   !
   ! !ARGUMENTS:
   class(DecompSumms_type), intent(inout) :: this

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this)
    !
    ! !USES:

    implicit none
    !
    ! !ARGUMENTS:
    class(DecompSumms_type), intent(inout) :: this

    this%o_scalar = 1._r8
    this%w_scalar  = 1._r8
    this%t_scalar    = 1._r8
    this%depth_scalar  = 1._r8

  end subroutine initCold

  !-----------------------------------------------------------------------
  subroutine InitPar(this, bgcsumms_con)

  use BgcConSummsType, only : BgcConSumms_type
  implicit none
  class(DecompSumms_type) , intent(inout) :: this
  type(BgcConSumms_type)   , intent(in) :: bgcsumms_con

  ! set "Q10" parameter
  this%Q10 = bgcsumms_con%Q10
  this%froz_q10 = bgcsumms_con%froz_q10
  this%decomp_depth_efolding = bgcsumms_con%decomp_depth_efolding

    this%ea_vmax_mic            = bgcsumms_con%ea_vmax_mic
    this%ea_vmax_enz            = bgcsumms_con%ea_vmax_enz
    this%ea_kaff_mono_mic       = bgcsumms_con%ea_kaff_mono_mic
    this%ea_kaff_enz_poly       = bgcsumms_con%ea_kaff_enz_poly
    this%ea_mr_mic              = bgcsumms_con%ea_mr_mic
    this%ea_kappa_mic           = bgcsumms_con%ea_kappa_mic
    this%ea_kaff_mono_msurf     = bgcsumms_con%ea_kaff_mono_msurf
    this%ea_kaff_enz_msurf      = bgcsumms_con%ea_kaff_enz_msurf
    
    this%ref_vmax_mic           = bgcsumms_con%ref_vmax_mic
    this%ref_vmax_enz           = bgcsumms_con%ref_vmax_enz
    this%ref_kaff_mono_mic      = bgcsumms_con%ref_kaff_mono_mic
    this%ref_kaff_enz_poly      = bgcsumms_con%ref_kaff_enz_poly
    this%ref_mr_mic             = bgcsumms_con%ref_mr_mic
    this%ref_kappa_mic          = bgcsumms_con%ref_kappa_mic
    this%ref_kaff_mono_msurf    = bgcsumms_con%ref_kaff_mono_msurf
    this%ref_kaff_enz_msurf     = bgcsumms_con%ref_kaff_enz_msurf

  end subroutine InitPar
  !-----------------------------------------------------------------------
  subroutine set_decompk_scalar(this, o2b, summs_forc)

  use BgcSummsForcType , only : summseca_forc_type
  use bshr_const_mod     , only : SHR_CONST_TKFRZ
  use BgcSummsMath     , only : interp1
  implicit none
  ! !ARGUMENTS:
  class(DecompSumms_type)     , intent(inout) :: this
  real(r8)                    , intent(in) :: o2b
  type(summseca_forc_type)       , intent(in) :: summs_forc

  ! !LOCAL VARIABLES:
  real(r8), parameter :: normalization_tref = 15._r8 ! reference temperature for normalizaion (degrees C)
  real(r8)            :: tref = 288.15 ! reference temperature (K)
  real(r8)            :: minpsi
  real(r8)            :: maxpsi
  real(r8)            :: normalization_factor
  real(r8)            :: catanf_30
  real(r8)            :: t1
  real(r8)            :: o2w
  real(r8)            :: psi
  !real(r8)            :: depz
  !real(r8)            :: Q10
  !real(r8)            :: froz_q10
  !real(r8)            :: decomp_depth_efolding
  !real(r8)            :: o2_w2b
  !real(r8)            :: soilpsi
  !real(r8)            :: sucsat

        ! Set up parameters to estimate the fraction of active enzymes at the current temperature
        real(r8) :: rgas = 8.31446 ! Universal gas constant (J/K/mol)
        integer :: trangebot = 220
        integer :: trangetop = 340
        integer, dimension(:), allocatable :: temp0    ! Range of temperatures for interpolation <- integers for now, but want reals
        real(r8), dimension(:), allocatable :: deltag0 ! Change in Gibbs free energy for enzymes over temperature range
        real(r8), dimension(:), allocatable :: t_fact0 ! Fraction of active enzymes over temperature range
        integer :: ii                                  ! Array constructor index
        integer :: dimtemp0                            ! Length of temperature range temp0
        real(r8) :: xpar1 = 249.544170969785           ! Number of amino acid residues per enzyme
        real(r8) :: xpar2 = 5341.422691388677          ! Enthalpy change at convergence temperature for enthalpy (J/mol)
        real(r8) :: xpar3 = 5.617549086429             ! Average number of non-polar hydrogen atoms per amino acid residue in enzyme
        real(r8) :: cp                                 ! Heat capacity
        real(r8) :: deltag1                            ! Change in Gibbs free energy for enzymes at reference temperature
        real(r8) :: t_fact1                            ! Fraction of active enzymes at reference temperature
        real(r8) :: deltas_star = 18.1                 ! Entropy change at the convergence temperature for entropy (J/K/mol)
        real(r8) :: ts_star = 385.2                    ! Convergence temperature for entropy (K)
        real(r8) :: th_star = 373.6                    ! Convergence temperature for enthalpy (K)
        real(r8) :: t_fact                             ! Active enzyme fraction at given temperatur
        ! Set up parameters to estimate the change in activity at the current temperature
        real(r8) :: fref0                              ! Modifies non-enzyme and non-equilibrium reactions
        real(r8) :: tinv                               ! Modifies activation energy
        real(r8) :: fref                               ! Modifies non-equilibrium enzymatic reactions
        !real(r8) :: temp                               ! Current temperature

  associate(                                       &
    temp          => summs_forc%temp    ,          &
    depz          => summs_forc%depz    ,          &
    o2_w2b        => summs_forc%o2_w2b  ,          &
    sucsat        => summs_forc%sucsat  ,          & ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
    soilpsi       => summs_forc%soilpsi ,          & ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
    Q10           => this%Q10           ,          &
    froz_q10      => this%froz_q10      ,          &
    decomp_depth_efolding => this%decomp_depth_efolding , &
    ea_vmax_mic   => this%ea_vmax_mic                   , &
    ea_vmax_enz   => this%ea_vmax_enz                    , &
    ea_kaff_mono_mic   => this%ea_kaff_mono_mic   , &
    ea_kaff_enz_poly   => this%ea_kaff_enz_poly    , &
    ea_mr_mic          => this%ea_mr_mic             , &
    ea_kappa_mic       => this%ea_kappa_mic          , &
    ea_kaff_mono_msurf => this%ea_kaff_mono_msurf , &
    ea_kaff_enz_msurf     => this%ea_kaff_enz_msurf      , &
    ref_vmax_mic          => this%ref_vmax_mic          , &        
    vmax_mic              => this%vmax_mic              , &         
    ref_vmax_enz          => this%ref_vmax_enz          , &       
    vmax_enz              => this%vmax_enz              , &       
    ref_kaff_mono_mic     => this%ref_kaff_mono_mic     , &  
    kaff_mono_mic         => this%kaff_mono_mic         , &   
    ref_kaff_enz_poly     => this%ref_kaff_enz_poly     , &  
    kaff_enz_poly         => this%kaff_enz_poly         , &     
    ref_mr_mic            => this%ref_mr_mic            , &   
    mr_mic                => this%mr_mic                , &      
    ref_kappa_mic         => this%ref_kappa_mic         , &    
    kappa_mic             => this%kappa_mic             , &      
    ref_kaff_mono_msurf   => this%ref_kaff_mono_msurf   , &
    kaff_mono_msurf       => this%kaff_mono_msurf       , &
    ref_kaff_enz_msurf    => this%ref_kaff_enz_msurf    , &
    kaff_enz_msurf        => this%kaff_enz_msurf          &     
  )

    catanf_30 = catanf(30._r8)

  !temperature scalar
  this%t_scalar     = 1._r8
  !use Charlie's Q10 based temperature scalar
  if (temp >= SHR_CONST_TKFRZ) then
    this%t_scalar= (Q10**((temp-(SHR_CONST_TKFRZ+25._r8))/10._r8))
  else
    this%t_scalar= (Q10**(-25._r8/10._r8))*(froz_q10**((temp-SHR_CONST_TKFRZ)/10._r8))
  endif

  ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
  normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10**((normalization_tref-25._r8)/10._r8))
  this%t_scalar = this%t_scalar * normalization_factor

  !Update temperature parameters

      temp0 = (/ (ii, ii=trangebot,trangetop) /) ! Generate sequence of temperatures

      ! Fraction of active enzymes using Murphy et al. 1990, see Tang & Riely 2015 Eq.46-48, also see Ratkowsky2005JTB
      cp = -46.+30.*(1.-1.54*(xpar1**(-0.268)))*xpar3 
      deltag0 = xpar2-deltas_star*temp0+cp*(temp0-th_star-temp0*log(temp0/ts_star))
      t_fact0 = 1./(1.+exp(-xpar1*deltag0/(rgas*temp0)))   
      deltag1 = xpar2-deltas_star*Tref+cp*(Tref-th_star-Tref*log(Tref/ts_star))
      t_fact1  = 1./(1.+exp(-xpar1*deltag1/(rgas*Tref)))

      t_fact0=t_fact0/t_fact1 ! Active enzyme fraction in total enzyme vs temperaure
    
      tinv=1./temp-1./tref ! Modifies activation energy
      
      call interp1(temp0, t_fact0, temp, t_fact) ! Interpolate to find fraction of active enzymes at current temperature
      ! This subroutine should return t_fact

      fref=t_fact*(temp/tref) ! Modifies non-equilibrium enzymatic reactions

  !Update parameters
    this%vmax_mic         = ref_vmax_mic *fref*exp(ea_vmax_mic*tinv)
    this%vmax_enz         = ref_vmax_enz *fref*exp(ea_vmax_enz*tinv)
    this%kaff_mono_mic    = ref_kaff_mono_mic *exp(ea_kaff_mono_mic*tinv)
    this%kaff_enz_poly    = ref_kaff_enz_poly *exp(ea_kaff_enz_poly*tinv)
    this%mr_mic           = ref_mr_mic        *exp(ea_mr_mic*tinv)
    this%kappa_mic        = ref_kappa_mic*fref*exp(ea_kappa_mic*tinv)  
    this%kaff_mono_msurf  = ref_kaff_mono_msurf*exp(ea_kaff_mono_msurf*tinv)   
    this%kaff_enz_msurf   = ref_kaff_enz_msurf*exp(ea_kaff_enz_msurf*tinv)

  !h2osoi_liqure scalar, also follows what Charlie has done
  minpsi = -10.0_r8
  this%w_scalar     = 1._r8
  maxpsi = sucsat * (-9.8e-6_r8)   !kg -> MPa
  psi = min(soilpsi,maxpsi)

  ! decomp only if soilpsi is higher than minpsi, some modification is needed for the following
  ! double check the paper by Wilson and Griffin, 1975
  if (psi > minpsi) then
    this%w_scalar = (log(minpsi/psi)/log(minpsi/maxpsi))
  else
    this%w_scalar = 0._r8
  end if

  !oxygen scalar, this is different from what CLM4.5bgc does, I use a M-M formulation to indicate O2 stress
  !and the O2 budget is done on the fly
  o2w = o2b / o2_w2b
  this%o_scalar = o2w/(o2w+0.02_r8)   !the value 0.22 mol O3/m3 is from Arah and Kirk, 2000

  !depth scalar, according to Koven et al. (2013), BG, the depth scalar is needed to resolve the radiocarbon profile
  this%depth_scalar = exp(-depz/decomp_depth_efolding)

  end associate
  end subroutine set_decompk_scalar

  !-----------------------------------------------------------------------
  function catanf(t1)result(ans)
  !DESCRIPTION:
  ! CENTURY T response function
  use bshr_const_mod, only : SHR_CONST_PI
  implicit none
  real(r8), intent(in) :: t1

  real(r8) :: ans
  ans = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

  end function catanf

end module BgcSummsDecompType
