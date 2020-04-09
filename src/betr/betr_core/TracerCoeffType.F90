module TracerCoeffType
  !
  ! DESCRIPTION:
  ! datatype for tracer phase conversion parameters and other scaling parameters
  !
  ! !USES:
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  use BeTR_decompMod      , only : bounds_type  => betr_bounds_type
  use BeTR_landvarconType , only : landvarcon => betr_landvarcon
  use TracerBaseType      , only : tracerbase_type
  !
  ! !PUBLIC TYPES:
  implicit none
#include "bshr_alloc.h"
  private


  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! !PUBLIC DATA:
  !
  !-------------------------------------------------------------------------------
  ! Column tracer phase conversion/transport parameters structure
  !-------------------------------------------------------------------------------
  type, public,  extends(tracerbase_type) ::  TracerCoeff_type
     real(r8), pointer :: aqu2neutralcef_col        (:,:,:)  => null()    !aqueous tracer into neutral aqueous tracer
     real(r8), pointer :: aqu2bulkcef_mobile_col    (:,:,:)  => null()    !coefficient to convert bulk concentration into aqueous phase, (nlevsno+nlevtrc_soil)
     real(r8), pointer :: gas2bulkcef_mobile_col    (:,:,:)  => null()    !coefficient to convert bulk concentration into gaseous phase, (nlevsno+nlevlak+nlevtrc_soil)
     real(r8), pointer :: aqu2equilsolidcef_col     (:,:,:)  => null()    !coefficient to convert solid phase (including ice) into aqueous phase
     real(r8), pointer :: henrycef_col              (:,:,:)  => null()    !henry's law constant
     real(r8), pointer :: bunsencef_col             (:,:,:)  => null()    !bunsen solubility
     real(r8), pointer :: henrycef_snow_col         (:,:,:)  => null()    !henry's law constant
     real(r8), pointer :: bunsencef_snow_col        (:,:,:)  => null()    !bunsen solubility
     real(r8), pointer :: tracer_diffusivity_air_col(:,:)   => null()     !diffusivity in the air
     real(r8), pointer :: aere_cond_col             (:,:)   => null()     !column level aerenchyma conductance (m/s)
     real(r8), pointer :: scal_aere_cond_col        (:,:)  => null()      !column level scaling factor for arenchyma or parenchyma transport
     real(r8), pointer :: rsoi_gas_topsno_col       (:,:)   => null()     !gas diffusivity in top snow layer, this is not used currently
     real(r8), pointer :: diffblkm_topsoi_col       (:,:)   => null()     !bulk molecular diffusivity in top soil layer
     real(r8), pointer :: diffgas_topsoi_col        (:,:)   => null()     !gas diffusivity in top soil layer, this is not used currently
     real(r8), pointer :: hmconductance_col         (:,:,:)  => null()    !geometrically weighted conductances (nlevsno+nlevtrc_soil)
     real(r8), pointer :: bulk_diffus_col           (:,:,:)  => null()
     real(r8), pointer :: aqu_diffus_col            (:,:,:) => null()
     real(r8), pointer :: aqu_diffus0_col           (:,:,:) => null()
     real(r8), pointer :: k_sorbsurf_col            (:,:,:) => null()     !affinity parameter for sorption surface
     real(r8), pointer :: Q_sorbsurf_col            (:,:,:) => null()     !maximum sorption capacity
     real(r8), pointer :: snowres_col               (:,:) => null()
   contains
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, public  :: retrieve_hist
  end type TracerCoeff_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! initialize the data type
    !
    ! !USES:
    use BeTRTracerType, only : BeTRTracer_Type
    implicit none
    ! !ARGUMENTS:
    class(TracerCoeff_type), intent(inout) :: this
    type(bounds_type)    , intent(in) :: bounds
    integer              , intent(in) :: lbj, ubj
    type(BeTRTracer_Type), intent(inout) :: betrtracer_vars

    call this%InitAllocate(bounds, lbj, ubj, betrtracer_vars)
    call this%tracer_base_init()
    call this%InitHistory(betrtracer_vars)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, flag, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! Now it is purposely empty, but will be potentially useful in the future
    ! !USES:
    use BetrTracerType , only : betrtracer_type
    use betr_varcon    , only : spval  => bspval

    !
    ! !ARGUMENTS:
    class(TracerCoeff_type), intent(inout) :: this
    type(bounds_type)    , intent(in)    :: bounds
    character(len=*)     , intent(in)    :: flag                                         ! 'read' or 'write'
    type(BeTRTracer_Type), intent(inout)   :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file

    ! remove compiler warnings for unused dummy args
    if (len(betrtracer_vars%betr_simname) > 0) continue
    if (bounds%begc > 0) continue


  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! Allocate memories for arrays in the datatype
    !
    ! !USES:
    use BeTRTracerType, only : BeTRTracer_Type
    implicit none
    !
    ! !ARGUMENTS:
    class(TracerCoeff_type), intent(inout)  :: this
    type(bounds_type), intent(in)     :: bounds
    integer, intent(in)               :: lbj, ubj
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc, nlevsno
    !---------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc; nlevsno = bounds%nlevsno

    NAN_ALLOC(this%aqu2neutralcef_col(begc:endc, lbj:ubj, 1:betrtracer_vars%ngwmobile_tracer_groups))

    NAN_ALLOC(this%aqu2bulkcef_mobile_col(begc:endc, lbj:ubj, 1:betrtracer_vars%ngwmobile_tracer_groups))

    NAN_ALLOC(this%gas2bulkcef_mobile_col(begc:endc, lbj:ubj, 1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%henrycef_col(begc:endc, lbj:ubj, 1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%bunsencef_col(begc:endc, lbj:ubj, 1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%tracer_diffusivity_air_col(begc:endc, 1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%scal_aere_cond_col(begc:endc, 1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%aere_cond_col(begc:endc,          1:betrtracer_vars%nvolatile_tracer_groups))

!    NAN_ALLOC(this%rsoi_gas_topsno_col(begc:endc, 1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%diffgas_topsoi_col(begc:endc,          1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%hmconductance_col(begc:endc, lbj:ubj, 1:betrtracer_vars%ntracer_groups))

    NAN_ALLOC(this%aqu2equilsolidcef_col(begc:endc, lbj:ubj, 1:betrtracer_vars%nsolid_equil_tracer_groups))

    NAN_ALLOC(this%bulk_diffus_col (begc:endc, lbj:ubj, 1:betrtracer_vars%ntracer_groups))

    NAN_ALLOC(this%aqu_diffus_col (begc:endc, lbj:ubj, 1:betrtracer_vars%ngwmobile_tracer_groups))

    NAN_ALLOC(this%aqu_diffus0_col (begc:endc, lbj:ubj, 1:betrtracer_vars%ngwmobile_tracer_groups))

    NAN_ALLOC(this%k_sorbsurf_col(begc:endc, lbj:ubj, 1:betrtracer_vars%nsolid_equil_tracers))

    NAN_ALLOC(this%Q_sorbsurf_col(begc:endc, lbj:ubj, 1:betrtracer_vars%nsolid_equil_tracers))

    NAN_ALLOC(this%diffblkm_topsoi_col(begc:endc, 1:betrtracer_vars%ngwmobile_tracer_groups))

    NAN_ALLOC(this%snowres_col(begc:endc, 1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%bunsencef_snow_col(begc:endc, -nlevsno+1:0,1:betrtracer_vars%nvolatile_tracer_groups))

    NAN_ALLOC(this%henrycef_snow_col(begc:endc, -nlevsno+1:0,1:betrtracer_vars%nvolatile_tracer_groups))
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! History fields initialization
    !
    ! !USES:
    !use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use betr_varcon    , only : spval => bspval
    use BeTRTracerType , only : BeTRTracer_Type
    use betr_constants  , only : betr_var_name_length
    !
    ! !ARGUMENTS:
    class(TracerCoeff_type), intent(inout) :: this
    type(BeTRTracer_Type), intent(inout) :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: jj, kk, trcid
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays
    integer :: num2d, num1d, it
    character(len=betr_var_name_length) :: tracername
    !use the interface provided from CLM
    associate(                                                                     &
         ntracer_groups       => betrtracer_vars%ntracer_groups                  , &
         tracer_group_memid   => betrtracer_vars%tracer_group_memid              , &
         ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups   , &
         nsolid_equil_tracers => betrtracer_vars%nsolid_equil_tracers            , &
         is_volatile          => betrtracer_vars%is_volatile                     , &
         volatilegroupid      => betrtracer_vars%volatilegroupid                   &
         )

     num2d = 0; num1d= 0
     do it = 1, 2
       do jj = 1, ntracer_groups
         trcid = tracer_group_memid(jj,1)
         tracername =  betrtracer_vars%get_tracername(trcid)
         if(jj <= ngwmobile_tracer_groups)then

            if(is_volatile(trcid))then
               kk = volatilegroupid(jj)

               call this%add_hist_var1d (it, num1d, fname='SCAL_ARENCHYMA_'//tracername, units='none',                &
                    avgflag='A', long_name='scaling factor for tracer transport through arenchyma for '//trim(tracername), &
                    default='inactive')

               call this%add_hist_var1d (it, num1d, fname='ARENCHYMA_'//tracername, units='m/s',                         &
                    avgflag='A', long_name='conductance for tracer transport through arenchyma for '//trim(tracername),    &
                    default='inactive')

               call this%add_hist_var1d (it, num1d, fname='CDIFF_TOPSOI_'//tracername, units='none',                     &
                    avgflag='A', long_name='gas diffusivity in top soil layer for '//trim(tracername),                     &
                    default='inactive')

               call this%add_hist_var2d (it, num2d, fname='CGAS2BULK_'//tracername, units='none', type2d='levtrc',          &
                    avgflag='A', long_name='converting factor from gas to bulk phase for '//trim(tracername),              &
                    default='inactive')
            endif

            call this%add_hist_var2d (it, num2d, fname='CAQU2BULK_vr_'//tracername, units='none', type2d='levtrc',            &
                 avgflag='A', long_name='converting factor from aqeous to bulk phase for '//trim(tracername),              &
                 default='inactive')

         endif

         call this%add_hist_var2d (it, num2d, fname='HMCONDC_vr_'//tracername, units='none', type2d='levtrc',            &
              avgflag='A', long_name='bulk conductance for '//trim(tracername),                           &
              default='inactive')
      enddo
      if(it==1)call this%alloc_hist_list(num1d, num2d)
      num2d = 0; num1d= 0
    enddo

    end associate
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    !  do cold initialization
    !
    ! !USES:
    use betr_varcon , only : spval  => bspval
    !
    ! !ARGUMENTS:
    class(TracerCoeff_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l       ! index

    !-----------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
      this%aqu2neutralcef_col        (c,:,:) = 0._r8
      this%aqu2bulkcef_mobile_col    (c,:,:) = 0._r8
      this%henrycef_col              (c,:,:) = 0._r8
      this%bunsencef_col             (c,:,:) = 0._r8
      this%hmconductance_col         (c,:,:) = 0._r8
      if(size(this%gas2bulkcef_mobile_col,3)>0)    this%gas2bulkcef_mobile_col    (c,:,:) = 0._r8
      if(size(this%tracer_diffusivity_air_col,2)>0)this%tracer_diffusivity_air_col(c,:)   = 0._r8
      if(size(this%scal_aere_cond_col,2)>0)        this%scal_aere_cond_col        (c,:)   = 0._r8
      if(size(this%aere_cond_col,2)>0)             this%aere_cond_col             (c,:)   = 0._r8
!      if(size(this%rsoi_gas_topsno_col,2)>0)       this%rsoi_gas_topsno_col       (c,:)   = 0._r8
      if(size(this%diffgas_topsoi_col,2)>0)        this%diffgas_topsoi_col   (c,:)   = 0._r8

      if(size(this%k_sorbsurf_col,3)>0)            this%k_sorbsurf_col (c,:,:) = 1._r8
      if(size(this%Q_sorbsurf_col,3)>0)            this%Q_sorbsurf_col (c,:,:) = 1._r8
      if(size(this%diffblkm_topsoi_col,2)>0)       this%diffblkm_topsoi_col(c,:)   = 0._r8
      if(size(this%snowres_col,2)>0) this%snowres_col(c,:) = 0._r8
    enddo

  end subroutine InitCold


  !----------------------------------------------------------------
  subroutine retrieve_hist(this, bounds, lbj, ubj, state_2d, state_1d, betrtracer_vars)
  !
  !DESCRIPTION
  !retrieve data for history output
  use MathfuncMod, only : addone
  use BeTRTracerType , only : BeTRTracer_Type
  implicit none
  class(TracerCoeff_type), intent(inout) :: this
  type(bounds_type)    , intent(in)  :: bounds
  integer, intent(in) :: lbj, ubj
  real(r8), intent(inout) :: state_2d(bounds%begc:bounds%endc, lbj:ubj,1:this%num_hist2d)
  real(r8), intent(inout) :: state_1d(bounds%begc:bounds%endc, 1:this%num_hist1d)
  type(BeTRTracer_Type)  , intent(in)  :: betrtracer_vars
  integer :: begc, endc
  integer :: jj, kk, trcid
  integer :: idtemp1d, idtemp2d

  associate(                                                       &
         ntracer_groups       => betrtracer_vars%ntracer_groups                  , &
         tracer_group_memid   => betrtracer_vars%tracer_group_memid              , &
         ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups   , &
         nsolid_equil_tracers => betrtracer_vars%nsolid_equil_tracers            , &
         is_volatile          => betrtracer_vars%is_volatile                     , &
         volatilegroupid      => betrtracer_vars%volatilegroupid                 , &
         tracernames          => betrtracer_vars%tracernames                       &
   )
  begc = bounds%begc; endc=bounds%endc
  idtemp1d = 0; idtemp2d = 0

      do jj = 1, ntracer_groups
         trcid = tracer_group_memid(jj,1)
         if(jj <= ngwmobile_tracer_groups)then

            if(is_volatile(trcid))then
               kk = volatilegroupid(jj)
               state_1d(:,addone(idtemp1d)) = this%scal_aere_cond_col(begc:endc, kk)

               state_1d(:,addone(idtemp1d)) = this%aere_cond_col(begc:endc, kk)

               state_1d(:,addone(idtemp1d)) = this%diffgas_topsoi_col(begc:endc, kk)

               state_2d(:,:,addone(idtemp2d)) = this%gas2bulkcef_mobile_col(:,:,kk)
            endif
            state_2d(:,:,addone(idtemp2d)) = this%aqu2bulkcef_mobile_col(:,:,jj)
         endif
         state_2d(:,:,addone(idtemp2d)) =  this%hmconductance_col(:,:,jj)
      enddo
  end associate
  end subroutine retrieve_hist
end module TracerCoeffType
