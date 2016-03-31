module BetrType

#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  !  subroutines for betr application
  !
  !  !USES:
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
  use BeTR_decompMod      , only : bounds_type  => betr_bounds_type
  use betr_ctrl           , only : iulog  => biulog
  use betr_constants      , only : betr_string_length
  use tracer_varcon       , only : nlevsoi  => betr_nlevsoi
  use tracer_varcon        , only : nlevsno => betr_nlevsno
  use betr_varcon           , only : spval => bspval
  use BeTR_PatchType        , only : pft => betr_pft
  use BeTR_ColumnType       , only : col => betr_col
  use BeTR_LandunitType     , only : lun => betr_lun

  use BGCReactionsMod, only : bgc_reaction_type
  use PlantSoilBGCMod, only : plant_soilbgc_type
  use BeTRTracerType, only : betrtracer_type
  use TracerCoeffType, only : TracerCoeff_type
  use TracerFluxType, only : TracerFlux_type
  use TracerStateType, only : TracerState_type
  use tracerboundarycondType, only : tracerboundarycond_type
  use BeTR_aerocondType, only : betr_aerecond_type
  use BeTR_CarbonFluxType, only : betr_carbonflux_type
  use BeTR_WaterstateType, only : betr_waterstate_type
  use BeTR_CNStateType, only : betr_cnstate_type
  use betr_time_manager   , only : get_nstep
  use EcophysConType, only : ecophyscon_type
  use BetrBGCMod
  implicit none

  private

  character(len=*), parameter :: filename = __FILE__

  type, public :: betr_type
     ! namelist control variables
     character(len=betr_string_length) :: reaction_method
     logical :: diffusion_on = .true.
     logical :: advection_on = .true.
     logical :: reaction_on  = .true.

     ! internal types
     class(bgc_reaction_type), allocatable, public :: bgc_reaction
     class(plant_soilbgc_type), allocatable, public :: plant_soilbgc

     type(BeTRtracer_type), public :: tracers
     type(TracerCoeff_type), public :: tracercoeffs
     type(TracerFlux_type), public :: tracerfluxes
     type(TracerState_type), public :: tracerstates
     type(tracerboundarycond_type), public :: tracerboundaryconds
     type(betr_cnstate_type), public :: cnstates
     type(betr_aerecond_type), public :: aereconds
     type(betr_carbonflux_type), public :: carbonfluxes
     type(betr_waterstate_type), public :: waterstate
     type(betr_cnstate_type), public :: cnstate

     real(r8), private, pointer :: h2osoi_liq_copy(:,:)
     real(r8), private, pointer :: h2osoi_ice_copy(:,:)

     ! FIXME(bja, 201603) replace LSM specific types!

   contains
     procedure, public :: Init
     procedure, public :: step_without_drainage
     procedure, public :: step_with_drainage
     procedure, public :: calc_dew_sub_flux
     procedure, public :: betr_lsm_flux_state_sendback
     procedure, public :: pre_diagnose_soilcol_water_flux
     procedure, public :: diagnose_advect_water_flux
     procedure, public :: diagnose_drainage_water_flux
     procedure, public :: diagnose_dtracer_freeze_thaw
     procedure, public :: check_mass_err
     procedure, public :: Enter_tracer_LayerAdjustment
     procedure, public :: Exit_tracer_LayerAdjustment
     procedure, public :: tracer_snowcapping
     procedure, public :: tracer_DivideSnowLayers
     procedure, public :: tracer_CombineSnowLayers
     procedure, private :: ReadNamelist

  end type betr_type

contains

!-------------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, bounds, waterstate, cnstate, ecophyscon)

    use babortutils, only : endrun
    use bshr_log_mod, only : errMsg => shr_log_errMsg
    use BeTR_decompMod, only : betr_bounds_type
    use betr_constants, only : betr_namelist_buffer_size
    use ReactionsFactory, only : create_bgc_reaction_type, &
         create_plant_soilbgc_type
    use TransportMod          , only : init_transportmod
    use TracerParamsMod       , only : tracer_param_init

    implicit none

    class(betr_type), intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    type(betr_bounds_type), intent(in) :: bounds
    type(betr_waterstate_type), intent(in) :: waterstate
    type(betr_cnstate_type), intent(in) :: cnstate
    type(ecophyscon_type), intent(in), optional :: ecophyscon

    type(ecophyscon_type) :: junk

    integer :: lbj, ubj

    if (present(ecophyscon)) then
       write(*,*) 'ERROR: ecophyscon not implemented in BeTR class '
       call endrun(msg=errMsg(filename,__LINE__))
    end if

    lbj = bounds%lbj;  ubj = bounds%ubj

    this%waterstate%h2osoi_liq_col => waterstate%h2osoi_liq_col
    this%waterstate%h2osoi_ice_col => waterstate%h2osoi_ice_col
    this%cnstate%isoilorder  => cnstate%isoilorder

    call this%ReadNamelist(namelist_buffer)

    allocate(this%bgc_reaction, &
         source=create_bgc_reaction_type(this%reaction_method))
    allocate(this%plant_soilbgc, &
         source=create_plant_soilbgc_type(this%reaction_method))

    call this%tracers%init_scalars()

    call this%bgc_reaction%Init_betrbgc(bounds, lbj, ubj, this%tracers)

    call this%aereconds%Init(bounds)

    call init_transportmod()

    call this%tracerstates%Init(bounds, lbj, ubj, this%tracers)

    call this%tracerfluxes%Init(bounds,  lbj, ubj, this%tracers)

    call this%tracercoeffs%Init(bounds, lbj, ubj, this%tracers)

    call this%tracerboundaryconds%Init(bounds, this%tracers)

    !inside Init_plant_soilbgc, specific plant soil bgc coupler data type will be created
    call this%plant_soilbgc%Init_plant_soilbgc(bounds, lbj, ubj)

    !initialize state variable
    call this%bgc_reaction%initCold(bounds,  this%tracers, this%waterstate, this%tracerstates)

    !initialize boundary condition type
    call this%bgc_reaction%init_boundary_condition_type(bounds, this%tracers, this%tracerboundaryconds)

    !initialize the betr parameterization module
    call tracer_param_init(bounds)

    ! FIXME(bja, 201603) ecophyscon is not currently being used, so we
    ! are explicitly passing initialized junk data.
    call this%bgc_reaction%init_betr_lsm_bgc_coupler(&
         bounds, this%plant_soilbgc, &
         this%tracers, this%tracerstates, this%cnstates, &
         junk)

   allocate(this%h2osoi_liq_copy(bounds%begc:bounds%endc, 1:nlevsoi));  this%h2osoi_liq_copy(:, :) = spval
   allocate(this%h2osoi_ice_copy(bounds%begc:bounds%endc, 1:nlevsoi));  this%h2osoi_ice_copy(:, :) = spval

  end subroutine Init

!-------------------------------------------------------------------------------
  subroutine ReadNamelist(this, namelist_buffer)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use betr_ctrl     , only : iulog => biulog
    use babortutils      , only : endrun
    use bshr_log_mod     , only : errMsg => shr_log_errMsg

    use betr_constants, only : stdout, betr_string_length_long, betr_namelist_buffer_size

    implicit none
    ! !ARGUMENTS:
    class(betr_type), intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer

    !
    ! !LOCAL VARIABLES:
    integer :: nml_error
    character(len=*), parameter :: subname = 'ReadNamelist'
    character(len=betr_string_length) :: reaction_method
    character(len=betr_string_length_long) :: ioerror_msg
    logical :: advection_on, diffusion_on, reaction_on

    !-----------------------------------------------------------------------

    namelist / betr_parameters / &
         reaction_method, &
         advection_on, diffusion_on, reaction_on

    reaction_method = ''
    advection_on = .true.
    diffusion_on = .true.
    reaction_on = .true.

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=betr_parameters, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading betr_parameters namelist "//errmsg(filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr bgc type :'
       write(stdout, *)
       write(stdout, *) ' betr_parameters namelist settings :'
       write(stdout, *)
       write(stdout, betr_parameters)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    this%reaction_method = reaction_method

    this%advection_on = advection_on
    this%diffusion_on = diffusion_on
    this%reaction_on = reaction_on

  end subroutine ReadNamelist

  !-------------------------------------------------------------------------------
  subroutine step_without_drainage(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,         &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, waterflux_vars, temperature_vars,&
       chemstate_vars, cnstate_vars, canopystate_vars,  carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! run betr code one time step forward, without drainage calculation

    ! !USES:
    use betr_time_manager             , only         : get_step_size
    use tracerfluxType               , only          : tracerflux_type
    use tracerstatetype              , only          : tracerstate_type
    use tracercoeffType              , only          : tracercoeff_type
    use TracerBoundaryCondType       , only          : TracerBoundaryCond_type
    use BetrTracerType               , only          : betrtracer_type
    use TracerParamsMod              , only          : set_phase_convert_coeff, set_multi_phase_diffusion, calc_tracer_infiltration
    use TracerParamsMod              , only          : get_zwt, calc_aerecond, betr_annualupdate
    use BeTR_SoilStateType           , only          : betr_soilstate_type
    use BeTR_WaterStateType          , only          : BeTR_Waterstate_Type
    use BeTR_TemperatureType         , only          : betr_temperature_type
    use BeTR_ChemStateType           , only          : betr_chemstate_type
    use BeTR_WaterfluxType           , only          : betr_waterflux_type
    use BGCReactionsMod              , only          : bgc_reaction_type
    use BeTR_atm2lndType             , only          : betr_atm2lnd_type
    use BeTR_SoilHydrologyType       , only          : betr_soilhydrology_type
    use PlantSoilBGCMod              , only          : plant_soilbgc_type
    use BeTR_CNStateType             , only          : betr_cnstate_type
    use BeTR_CarbonFluxType          , only          : betr_carbonflux_type
    use BeTR_aerocondType            , only          : betr_aerecond_type
    use betr_ctrl                    , only          : is_active_betr_bgc
    use BeTR_CanopyStateType         , only          : betr_canopystate_type

    !
    ! !ARGUMENTS :
    class(betr_type), intent(inout) :: this
    type(bounds_type)                , intent(in)    :: bounds                     ! bounds
    integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    integer                          , intent(in)    :: num_soilp
    integer                          , intent(in)    :: filter_soilp(:)            ! pft filter
    type(betr_waterstate_type)       , intent(in)    :: waterstate_vars            ! water state variables
    type(betr_soilstate_type)        , intent(in)    :: soilstate_vars             ! column physics variable
    type(betr_temperature_type)      , intent(in)    :: temperature_vars           ! energy state variable
    type(betr_chemstate_type)        , intent(in)    :: chemstate_vars
    type(betr_atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(betr_soilhydrology_type)    , intent(in)    :: soilhydrology_vars
    type(betr_cnstate_type)          , intent(inout) :: cnstate_vars
    type(betr_canopystate_type)      , intent(in)    :: canopystate_vars
    type(betr_carbonflux_type)       , intent(in)    :: carbonflux_vars
    type(betr_waterflux_type)        , intent(inout) :: waterflux_vars


    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'run_betr_one_step_without_drainage'
    real(r8)           :: dtime2, dtime
    real(r8)           :: Rfactor(bounds%begc:bounds%endc, bounds%lbj:bounds%ubj,1:this%tracers%ngwmobile_tracers) !retardation factor
    integer            :: j
    integer            :: lbj, ubj

    lbj = bounds%lbj; ubj = bounds%ubj

    dtime = get_step_size()

    !initialize extra parameters
    dtime2 = dtime * 0.5_r8


    call stage_tracer_transport(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, atm2lnd_vars, &
         carbonflux_vars, soilstate_vars, waterstate_vars, waterflux_vars, temperature_vars, soilhydrology_vars, &
         chemstate_vars, this%aereconds, canopystate_vars, this%tracers, this%tracercoeffs, &
         this%tracerboundaryconds, this%tracerfluxes, this%bgc_reaction, Rfactor, this%advection_on)

    call surface_tracer_hydropath_update(bounds, num_soilc, filter_soilc, &
       waterstate_vars, waterflux_vars, soilhydrology_vars, this%tracers, this%tracerstates, &
       this%tracercoeffs,  this%tracerfluxes)

    call this%bgc_reaction%calc_bgc_reaction(bounds, lbj, ubj, &
         num_soilc,                                       &
         filter_soilc,                                    &
         num_soilp,                                       &
         filter_soilp,                                    &
         this%tracerboundaryconds%jtops_col,              &
         dtime,                                           &
         this%tracers,                                    &
         this%tracercoeffs,                               &
         cnstate_vars,                                    &
         this%tracerstates,                               &
         this%tracerfluxes,                               &
         this%tracerboundaryconds,                        &
         this%plant_soilbgc)

    call tracer_gws_transport(bounds, num_soilc, filter_soilc, Rfactor, waterstate_vars, &
      waterflux_vars, this%tracers, this%tracerboundaryconds, this%tracercoeffs, &
      this%tracerstates, this%tracerfluxes, this%bgc_reaction, this%advection_on, this%diffusion_on)

    call calc_ebullition(bounds, 1, ubj,                                 &
         this%tracerboundaryconds%jtops_col,                             &
         num_soilc,                                                      &
         filter_soilc,                                                   &
         atm2lnd_vars%forc_pbot_downscaled_col,                          &
         col%zi(bounds%begc:bounds%endc, 0:ubj),                         &
         col%dz(bounds%begc:bounds%endc, 1:ubj),                         &
         dtime,                                                          &
         soilhydrology_vars%fracice_col(bounds%begc:bounds%endc, 1:ubj), &
         soilhydrology_vars%zwts_col(bounds%begc:bounds%endc),           &
         this%tracers,                                                   &
         this%tracercoeffs,                                              &
         this%tracerstates,                                              &
         this%tracerfluxes%tracer_flx_ebu_col(bounds%begc:bounds%endc, 1:this%tracers%nvolatile_tracers))

    if (is_active_betr_bgc) then
       !update nitrogen storage pool
       call this%plant_soilbgc%plant_soilbgc_summary(bounds, lbj, ubj, num_soilc,                   &
            filter_soilc,                                                                              &
            col%dz(bounds%begc:bounds%endc,1:ubj),                                                     &
            this%tracers, this%tracerfluxes)
    endif

  end subroutine step_without_drainage


  !--------------------------------------------------------------------------------
  subroutine step_with_drainage(this, bounds,  num_soilc, filter_soilc, jtops, waterflux_vars)
    !
    ! !DESCRIPTION:
    ! do tracer update due to drainage
    !
    ! !USES:
    use tracerfluxType        , only : tracerflux_type
    use tracerstatetype       , only : tracerstate_type
    use tracercoeffType       , only : tracercoeff_type
    use MathfuncMod           , only : safe_div
    use BeTR_WaterFluxType    , only : betr_waterflux_type
    ! !ARGUMENTS:
    class(betr_type), intent(inout) :: this
    type(bounds_type),        intent(in)    :: bounds
    integer,                  intent(in)    :: num_soilc                          ! number of columns in column filter_soilc
    integer,                  intent(in)    :: filter_soilc(:)                    ! column filter_soilc
    integer,                  intent(in)    :: jtops(bounds%begc: )
    type(betr_waterflux_type)    , intent(in)    :: waterflux_vars

    ! !LOCAL VARIABLES:
    real(r8) :: aqucon
    integer  :: fc, c, j, k
    integer  :: lbj, ubj


    SHR_ASSERT_ALL((ubound(jtops)         == (/bounds%endc/))      , errMsg(filename,__LINE__))

    associate(                                                                   & !
         ngwmobile_tracers        => this%tracers%ngwmobile_tracers         , & !
         groupid                  => this%tracers%groupid                   , & !
         is_h2o                   => this%tracers%is_h2o                    , & !
         is_advective             => this%tracers%is_advective              , & !
         aqu2bulkcef_mobile       => this%tracercoeffs%aqu2bulkcef_mobile_col   , & !
         tracer_conc_mobile       => this%tracerstates%tracer_conc_mobile_col   , & !
         tracer_conc_grndwater    => this%tracerstates%tracer_conc_grndwater_col, & !
         dz                       => col%dz                                    , & !
         tracer_flx_drain         => this%tracerfluxes%tracer_flx_drain_col      , & !
         qflx_drain_vr            => waterflux_vars%qflx_drain_vr_col          , & ! Output  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/step) (to river +)
         qflx_totdrain            => waterflux_vars%qflx_totdrain_col            & ! Output  : [real(r8) (:,:) ]  (m H2o/step)

         )
      lbj = bounds%lbj; ubj = bounds%ubj

      do j = lbj, ubj
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            if(j>=jtops(c))then
               do k = 1, ngwmobile_tracers
                  !obtain aqueous concentration
                  if(.not. is_advective(k))cycle
                  if(qflx_drain_vr(c,j) > 0._r8)then
                    aqucon = safe_div(tracer_conc_mobile(c,j,k),aqu2bulkcef_mobile(c,j,groupid(k)))
                  else
                    !when drainage is negative, tracer comes from groundwater
                    aqucon = tracer_conc_grndwater(c,k)
                  endif
                  !when drainage is negative, assume the flux is magically coming from other groundwater sources

                  tracer_flx_drain(c,k)     = tracer_flx_drain(c,k)  + aqucon * qflx_drain_vr(c,j)
                  tracer_conc_mobile(c,j,k) =  tracer_conc_mobile(c,j,k) - aqucon * qflx_drain_vr(c,j)/dz(c,j)
                  if(tracer_conc_mobile(c,j,k)<0._r8)then
                     tracer_flx_drain(c,k) = tracer_flx_drain(c,k)+tracer_conc_mobile(c,j,k)*dz(c,j)
                     tracer_conc_mobile(c,j,k)=0._r8
                  endif
               enddo
            endif
         enddo
      enddo

      !diagnose gas pressure
      call diagnose_gas_pressure(bounds, lbj, ubj, num_soilc, filter_soilc, &
           this%tracers, this%tracercoeffs, this%tracerstates)

    end associate
  end subroutine step_with_drainage


  !--------------------------------------------------------------------------------
  subroutine calc_dew_sub_flux(this, bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars, betrtracer_vars, tracerflux_vars, tracerstate_vars)
    !
    ! DESCRIPTION:
    ! calculate water flux from dew formation, and sublimation
    ! !USES:
    use betr_time_manager      , only : get_step_size
    use BeTR_WaterfluxType    , only : betr_waterflux_type
    use BeTR_WaterstateType   , only : betr_waterstate_type
    use tracerfluxType        , only : tracerflux_type
    use tracerstatetype       , only : tracerstate_type
    use betr_varcon           , only : spval => bspval
    use tracer_varcon         , only : nlevtrc_soil  => betr_nlevtrc_soil
    use BeTR_landvarconType   , only : landvarcon  => betr_landvarcon

    ! !ARGUMENTS:
    class(betr_type), intent(inout) :: this
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_hydrologyc             ! number of column soil points in column filter_soilc
    integer                   , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
    type(betr_waterstate_type)     , intent(in)    :: waterstate_vars
    type(betr_waterflux_type)      , intent(in)    :: waterflux_vars
    class(betrtracer_type)     , intent(in)    :: betrtracer_vars            ! betr configuration information
    type(tracerflux_type)     , intent(inout) :: tracerflux_vars            ! tracer flux
    type(tracerstate_type)    , intent(inout) :: tracerstate_vars           ! tracer state variables data structure

    ! !LOCAL VARIABLES:
    real(r8) :: dtime
    real(r8) :: tot1, totz, tot0, tot2, tot3, tot4
    integer :: fc, c, j, l, ll, ll2

    ! remove compiler warnings about unused dummy args
    if (this%advection_on) continue
    if (bounds%begc > 0) continue

   associate(                                                               & !
         snl                =>    col%snl                                 ,  & ! Input:  [integer  (:)   ]  number of snow layers
         dz                 =>    col%dz                                  ,  & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         h2osoi_ice         =>    waterstate_vars%h2osoi_ice_col          ,  & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq         =>    waterstate_vars%h2osoi_liq_col          ,  & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)
         frac_h2osfc        =>    waterstate_vars%frac_h2osfc_col         ,  & ! Input:  [real(r8) (:)   ]
         qflx_dew_grnd      =>    waterflux_vars%qflx_dew_grnd_col        ,  & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_dew_snow      =>    waterflux_vars%qflx_dew_snow_col        ,  & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s
         qflx_sub_snow_vol  =>    waterflux_vars%qflx_sub_snow_vol_col    ,  & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s)
         qflx_snow2topsoi   =>    waterflux_vars%qflx_snow2topsoi_col     ,  & ! Output: [real(r8) (:)   ]
         qflx_h2osfc2topsoi =>    waterflux_vars%qflx_h2osfc2topsoi_col   ,  & !
         tracer_flx_dew_grnd=>    tracerflux_vars%tracer_flx_dew_grnd_col ,  & !
         tracer_flx_dew_snow=>    tracerflux_vars%tracer_flx_dew_snow_col ,  & !
         tracer_flx_sub_snow=>    tracerflux_vars%tracer_flx_sub_snow_col ,  & !
         tracer_conc_mobile =>    tracerstate_vars%tracer_conc_mobile_col ,  & !
         tracer_conc_frozen =>    tracerstate_vars%tracer_conc_frozen_col ,  & !
         is_h2o             =>    betrtracer_vars%is_h2o                  ,  & !
         tracernames        =>    betrtracer_vars%tracernames             ,  &
         frozenid           =>    betrtracer_vars%frozenid                ,  &
         clandunit          =>    col%landunit                             , & ! Input:  [integer  (:)   ]  columns's landunit
         ltype              =>    lun%itype                                , & ! Input:  [integer  (:)   ]  landunit type
         ngwmobile_tracers  =>    betrtracer_vars%ngwmobile_tracers          &
         )


         !initialize the following fluxes to zero
      do fc = 1, num_hydrologyc
        c = filter_soilc_hydrologyc(fc)
        qflx_snow2topsoi(c) = 0._r8
        qflx_h2osfc2topsoi(c) = 0._r8
      enddo

      !do tracer update from dew and sublimation
      dtime = get_step_size()

      do j = 1, ngwmobile_tracers
         !now only do water isotope tracer
         if(.not. is_h2o(j))cycle

         do fc = 1, num_hydrologyc
            c = filter_soilc_hydrologyc(fc)
            l = clandunit(c)
            if (ltype(l)/= landvarcon%istsoil .and. ltype(l)/= landvarcon%istcrop)cycle
            if(snl(c)+1>=1)then
               tracer_flx_dew_grnd(c, j) = (1._r8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtime
               tracer_flx_dew_snow(c, j) = (1._r8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtime
               tracer_flx_sub_snow(c,j) = tracer_conc_frozen(c,1,frozenid(j)) * dz(c,1) * qflx_sub_snow_vol(c)
            else
               tracer_flx_dew_grnd(c, j) = 0._r8
               tracer_flx_dew_snow(c, j) = 0._r8
               tracer_flx_sub_snow(c, j) = 0._r8
            endif
         enddo
      enddo

      !apply those fluxes
      do j = 1, ngwmobile_tracers
         if(.not. is_h2o(j))cycle
         do fc = 1, num_hydrologyc
            c = filter_soilc_hydrologyc(fc)
            l = clandunit(c)
            if (ltype(l) /= landvarcon%istsoil .and. ltype(l) /= landvarcon%istcrop)cycle
            tracer_conc_mobile(c,1,j) = tracer_conc_mobile(c,1,j) + tracer_flx_dew_grnd(c, j)/dz(c,1)
            tracer_conc_frozen(c,1,frozenid(j)) = tracer_conc_frozen(c,1,frozenid(j)) + &
                 (tracer_flx_dew_snow(c, j)-tracer_flx_sub_snow(c,j))/dz(c,1)
         enddo
      enddo
    end associate

  end subroutine calc_dew_sub_flux


  !--------------------------------------------------------------------------------

  subroutine check_mass_err(this, c, trcid, ubj, dz, betrtracer_vars, tracerstate_vars, tracerflux_vars)

  !
  !temporary mass balance check
  use tracerfluxType        , only : tracerflux_type
  use tracerstatetype       , only : tracerstate_type

  implicit none

  class(betr_type), intent(inout) :: this
  integer, intent(in) :: ubj
  integer, intent(in) :: c, trcid
  real(r8), intent(in):: dz(1:ubj)
  class(betrtracer_type)   , intent(in) :: betrtracer_vars  ! tracer info data structure
  type(tracerflux_type)    , intent(in) :: tracerflux_vars  ! tracer flux
  class(tracerstate_type)  , intent(in) :: tracerstate_vars ! tracer state variables data structure
  real(r8) :: totmass, err

  ! remove compiler warnings about unused dummy args
  if (this%advection_on) continue

  associate(                                                                            &
       beg_tracer_molarmass      => tracerstate_vars%beg_tracer_molarmass_col         , &
       tracer_flx_netpro         => tracerflux_vars%tracer_flx_netpro_col             , &
       tracer_flx_netphyloss     => tracerflux_vars%tracer_flx_netphyloss_col         , &
       ntracers                  => betrtracer_vars%ntracers                          , &
       is_adsorb                 => betrtracer_vars%is_adsorb                         , &
       adsorbid                  => betrtracer_vars%adsorbid                          , &
       is_frozen                 => betrtracer_vars%is_frozen                         , &
       frozenid                  => betrtracer_vars%frozenid                            &
       )

  totmass=tracerstate_vars%int_mass_mobile_col(1,ubj,c,trcid,dz(1:ubj))

  if(is_frozen(trcid))then
     totmass = totmass + &
          tracerstate_vars%int_mass_frozen_col(1,ubj,c,frozenid(trcid),dz(1:ubj))
  endif
  call tracerflux_vars%flux_summary(c, betrtracer_vars)
  err=beg_tracer_molarmass(c,trcid)-totmass  &
       + tracer_flx_netpro(c,trcid)-tracer_flx_netphyloss(c,trcid)
  print*,'err',err
  end associate
  end subroutine check_mass_err


   !------------------------------------------------------------------------
   subroutine diagnose_dtracer_freeze_thaw(this, bounds, num_nolakec, filter_nolakec, &
     waterstate_vars)
   !
   ! DESCRIPTION
   ! aqueous tracer partition based on freeze-thaw
   !
   ! USES
   !
   use BeTR_LandunitType     , only : lun => betr_lun
   use BeTR_WaterStateType   , only : betr_waterstate_type
   use tracerstatetype       , only : tracerstate_type
   use BeTRTracerType        , only : betrtracer_type
   use BeTR_landvarconType   , only : landvarcon => betr_landvarcon

   !
   ! Arguments
   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)      , intent(in)    :: bounds
   integer                , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points
   type(betr_waterstate_type)  , intent(in)    :: waterstate_vars


   integer :: j, fc, c, l, k
   real(r8) :: thaw_frac, freeze_frac
   real(r8) :: dtracer

   if (bounds%begc > 0) continue

   associate(                                                   &
     frozenid         => this%tracers%frozenid        , &
     is_frozen        => this%tracers%is_frozen       , &
     ngwmobile_tracers=> this%tracers%ngwmobile_tracers,&
     h2osoi_liq       => waterstate_vars%h2osoi_liq_col  , &
     h2osoi_ice       => waterstate_vars%h2osoi_ice_col  , &
     tracer_conc_mobile=> this%tracerstates%tracer_conc_mobile_col , &
     tracer_conc_frozen=> this%tracerstates%tracer_conc_frozen_col   &

   )
   ! the tracer concentration change between frozen and liquid
   !  phases due to freeze and thaw only occurs to non-volatile tracers
   ! and water tracers

   do j = 1, nlevsoi
     do fc = 1, num_nolakec
       c =  filter_nolakec(fc)
       l = col%landunit(c)
       if(lun%itype(l) == landvarcon%istsoil)then
         do k = 1, ngwmobile_tracers
           !if it is a frozenable tracer, do it
           if(is_frozen(k))then
             if(h2osoi_liq(c,j) > this%h2osoi_liq_copy(c,j))then
               !thaw, solid to aqueous
               thaw_frac = min(1._r8-h2osoi_ice(c,j)/this%h2osoi_ice_copy(c,j),1._r8)
               dtracer = tracer_conc_frozen(c,j,frozenid(k)) * thaw_frac
               tracer_conc_frozen(c,j,frozenid(k)) = tracer_conc_frozen(c,j,frozenid(k)) - dtracer
               tracer_conc_mobile(c,j, k) = tracer_conc_mobile(c,j, k) + dtracer
             else
               !freeze, aqueous to solid
               freeze_frac = min(1._r8 - h2osoi_liq(c,j)/this%h2osoi_liq_copy(c,j),1._r8)
               dtracer = tracer_conc_mobile(c,j,k) * freeze_frac          !some modifier are needed to account for change in solubility
               tracer_conc_frozen(c,j,frozenid(k)) = tracer_conc_frozen(c,j,frozenid(k)) + dtracer
               tracer_conc_mobile(c,j, k) = tracer_conc_mobile(c,j, k) - dtracer
             endif
           endif
         enddo
       endif
     enddo
   enddo

   end associate
   end subroutine diagnose_dtracer_freeze_thaw

  !------------------------------------------------------------------------
  subroutine betr_lsm_flux_state_sendback(this, bounds,  num_soilc, filter_soilc )

   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)      , intent(in)    :: bounds
   integer                , intent(in)    :: num_soilc                        ! number of column non-lake points in column filter
   integer                , intent(in)    :: filter_soilc(:)                  ! column filter for non-lake points

   ! remove compiler warnings about unused dummy args
   if (this%advection_on) continue
   if (bounds%begc > 0) continue
   if (num_soilc > 0) continue
   if (size(filter_soilc) > 0) continue

!x   call this%bgc_reaction%betr_lsm_flux_state_sendback(bounds, &
!x       num_soilc, filter_soilc,                    &
!x       this%tracerstate_vars, this%tracerflux_vars,  this%betrtracer_vars)
  end subroutine betr_lsm_flux_state_sendback


   !------------------------------------------------------------------------
   subroutine pre_diagnose_soilcol_water_flux(this, bounds, num_nolakec, filter_nolakec, waterstate_vars)
   !
   ! DESCRIPTION
   ! pre diagnose advective water fluxes at different soil interfaces

   use BeTR_WaterStateType        , only : betr_waterstate_type
   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)      , intent(in)    :: bounds
   integer                , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points
   type(betr_waterstate_type)  , intent(in)    :: waterstate_vars
   integer :: j, fc, c, l

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue

   do j = 1, nlevsoi
     do fc = 1, num_nolakec
       c =  filter_nolakec(fc)
       this%h2osoi_liq_copy(c,j) = waterstate_vars%h2osoi_liq_col(c,j)
       this%h2osoi_ice_copy(c,j) = waterstate_vars%h2osoi_ice_col(c,j)
     enddo
   enddo
   end subroutine pre_diagnose_soilcol_water_flux

   !------------------------------------------------------------------------
   subroutine diagnose_advect_water_flux(this, bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
        waterstate_vars, qflx_bot, waterflux_vars)
   !
   ! DESCRIPTION
   ! diagnose advective water fluxes between different soil layers
   !
   use BeTR_WaterStateType  , only : betr_waterstate_type
   use BeTR_WaterFluxType   , only : betr_waterflux_type
   use betr_time_manager     , only : get_step_size
   use betr_varcon          , only : denh2o => bdenh2o
   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
   integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
   integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
   type(betr_waterflux_type)    , intent(inout) :: waterflux_vars
   type(betr_waterstate_type)   , intent(in)     :: waterstate_vars
   real(r8)                , intent(in)    :: qflx_bot(bounds%begc: )  ! mm H2O/s water exchange rate between soil col and aquifer

   !local variables
   integer :: j, fc, c
   real(r8):: dtime
   real(r8):: diff
   real(r8):: infl_tmp
   real(r8):: scal
   logical :: tf

   if (num_urbanc > 0) continue
   if (size(filter_urbanc) > 0) continue

   SHR_ASSERT_ALL((ubound(qflx_bot)    == (/bounds%endc/))         , errMsg(filename,__LINE__))

   associate(                                                             & !
     h2osoi_liq          =>    waterstate_vars%h2osoi_liq_col           , &
     h2osoi_ice          =>    waterstate_vars%h2osoi_ice_col           , &
     qflx_rootsoi        =>    waterflux_vars%qflx_rootsoi_col          , & ! Iput  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/s) (+ = to atm)
     qflx_adv            =>    waterflux_vars%qflx_adv_col              , & ! Output: [real(r8) (:,:) ]  water flux at interfaces       (m H2O/s) (- = to atm)
     qflx_gross_infl_soil=>    waterflux_vars%qflx_gross_infl_soil_col  , & ! Output: [real(r8) (:)] gross infiltration (mm H2O/s)
     qflx_infl           =>    waterflux_vars%qflx_infl_col             , & ! Output: [real(r8) (:)] infiltration, mm H2O/s
     qflx_gross_evap_soil=>    waterflux_vars%qflx_gross_evap_soil_col    & ! Output: [real(r8) (:)] gross evaporation (mm H2O/s)
   )

   ! get time step
   dtime = get_step_size()
   !start from the bottom layer, because the water exchange between vadose zone soil and aquifer and plant root is known
   !the water flux at uppper surface can be inferred using the mass balance approach
   do j = nlevsoi, 1, -1
     do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if(j==nlevsoi)then
         qflx_adv(c,j) = qflx_bot(c) * 1.e-3_r8                                 ! m/s
       else
         qflx_adv(c,j) = 1.e-3_r8 * (h2osoi_liq(c,j+1)-this%h2osoi_liq_copy(c,j+1))/dtime + qflx_adv(c,j+1) + qflx_rootsoi(c,j+1)
       endif
     enddo
   enddo

   ! correct gross infiltration and gross evaporation
   ! (h2osoi_liq(c,1)-h2osoi_liq_copy(c,1))/dtime=qflx_infl-q_out-qflx_rootsoi
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)

     !obtain the corrected infiltration
     qflx_infl(c) = (h2osoi_liq(c,1)-this%h2osoi_liq_copy(c,1))/dtime + (qflx_rootsoi(c,1)+qflx_adv(c,1))*1.e3_r8
     !the predicted net infiltration
     infl_tmp=qflx_gross_infl_soil(c)-qflx_gross_evap_soil(c)
     diff=qflx_infl(c)-infl_tmp
     if(abs(diff)/=0._r8)then
       if(infl_tmp==0._r8)then
         if(diff>0._r8)then
           !the corrected infiltration > net infiltration
           qflx_gross_infl_soil(c)=qflx_gross_infl_soil(c) + diff
         else
           qflx_gross_evap_soil(c)=qflx_gross_evap_soil(c)-diff
         endif
       else
         scal = (1._r8+diff/infl_tmp)
         qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c) * scal
         qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c) * scal
         if(qflx_gross_evap_soil(c)<0._r8)then
           !no negative evaporation allowed
           qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c)-qflx_gross_evap_soil(c)
           qflx_gross_evap_soil(c) = 0._r8
         endif
         if(qflx_gross_infl_soil(c)<0._r8)then
           qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c)-qflx_gross_infl_soil(c)
           qflx_gross_infl_soil(c) = 0._r8
         endif
       endif
     endif
     qflx_adv(c,0) = qflx_gross_infl_soil(c) *1.e-3_r8  !surface infiltration
   enddo

   end associate
   end subroutine diagnose_advect_water_flux

   !------------------------------------------------------------------------
   subroutine diagnose_drainage_water_flux(this, bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
        waterstate_vars,  waterflux_vars)
   !
   ! DESCRIPTION
   ! diagnose advective water fluxes between different soil layers
   !

   use BeTR_WaterFluxType   , only : betr_waterflux_type
   use BeTR_WaterStateType  , only : betr_waterstate_type
   use betr_varcon          , only : denh2o => bdenh2o
   use betr_time_manager     , only : get_step_size

   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
   integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
   integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
   type(betr_waterflux_type)    , intent(inout) :: waterflux_vars
   type(betr_waterstate_type)   , intent(in)    :: waterstate_vars

   !local variables
   integer :: j, fc, c
   real(r8):: dtime

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue
   if (num_urbanc > 0) continue
   if (size(filter_urbanc) > 0) continue

   associate(                                                           & !
     h2osoi_liq         =>    waterstate_vars%h2osoi_liq_col          , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
     h2osoi_ice         =>    waterstate_vars%h2osoi_ice_col          , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
     qflx_drain_vr      =>    waterflux_vars%qflx_drain_vr_col        , & ! Output  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/step) (to river +)
     qflx_totdrain      =>    waterflux_vars%qflx_totdrain_col          & ! Output  : [real(r8) (:,:) ]  (m H2o/step)
   )

   ! get time step
   dtime = get_step_size()
   !start from the bottom layer, because the water exchange between vadose zone soil and aquifer and plant root is known
   !the water flux at uppper surface can be inferred using the mass balance approach
   !also, the flux through transpiration has been already calcualted in the no drainage code, therefore no
   !the change in water content is solely due to subsurface draiange.
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     qflx_totdrain(c) = 0._r8
   enddo
   do j = nlevsoi, 1, -1
     do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qflx_drain_vr(c,j) = this%h2osoi_liq_copy(c,j)-h2osoi_liq(c,j)   !kg/m2/step
       if(j==1)then
         !the following line is due to the ice check that is added in drainage calculation
         qflx_drain_vr(c,j) = qflx_drain_vr(c,j) + (this%h2osoi_ice_copy(c,j)-h2osoi_ice(c,j))
       endif
       !the following line will allow negative drainage
       qflx_drain_vr(c,j) =qflx_drain_vr(c,j) * 1.e-3_r8               !kg/m2/(kg/m3)/step = m/step

       qflx_totdrain(c) = qflx_totdrain(c) + qflx_drain_vr(c,j)
     enddo
   enddo

   end associate
   end subroutine diagnose_drainage_water_flux

   !------------------------------------------------------------------------

   subroutine tracer_DivideSnowLayers(this, bounds, num_snowc, filter_snowc, divide_matrix)
   !
   ! DESCRIPTIONS
   ! divide tracers in snow
   !
   ! USES
   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8), intent(in) :: divide_matrix(bounds%begc: , 1: , 1: )

   !temporary variable
   real(r8) :: tracer_copy(bounds%begc:bounds%endc, 1:nlevsno)
   integer  :: jj, jj1

   SHR_ASSERT_ALL((ubound(divide_matrix,1)         == bounds%endc)      , errMsg(filename,__LINE__))
   SHR_ASSERT_ALL((ubound(divide_matrix,2)         == nlevsno)      , errMsg(filename,__LINE__))
   SHR_ASSERT_ALL((ubound(divide_matrix,3)         == nlevsno)      , errMsg(filename,__LINE__))

   associate(                             &
      tracer_conc_frozen_col => this%tracerstates%tracer_conc_frozen_col , &
      tracer_conc_mobile_col => this%tracerstates%tracer_conc_mobile_col , &
      tracer_conc_solid_passive_col => this%tracerstates%tracer_conc_solid_passive_col, &
      ntracers         => this%tracers%ntracers, &
      ngwmobile_tracers=> this%tracers%ngwmobile_tracers , &
      is_frozen => this%tracers%is_frozen &
   )

   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
        !make a copy
        call tracer_copy_a2b_div(bounds, num_snowc, filter_snowc, col%snl, &
           tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:0,jj),tracer_copy)
        !remapping
        call tracer_col_mapping_div(bounds, num_snowc, filter_snowc, col%snl, divide_matrix, &
           tracer_copy, tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj))

        if(is_frozen(jj))then
          !make a copy
          call tracer_copy_a2b_div(bounds, num_snowc, filter_snowc, col%snl, &
            tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:0,jj),tracer_copy)
          !remapping
          call tracer_col_mapping_div(bounds, num_snowc, filter_snowc, col%snl, divide_matrix, &
             tracer_copy, tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:0,jj))
        endif
     else
       !x!currently, no solid tracer is allowed to mix with snow
       !! jj1=jj-ngwmobile_tracers
       !x!tracer_conc_solid_passive_col
       !x!copy
       !x!remapping
     endif
   enddo
   end associate
   end subroutine tracer_DivideSnowLayers
   !------------------------------------------------------------------------

   subroutine tracer_CombineSnowLayers(this, bounds, num_snowc, filter_snowc, combine_matrix)
   !
   ! DESCRIPTIONS
   ! combine tracers in snow
   !
   !! USES
   use tracer_varcon, only : nlevsno => betr_nlevsno
   !!
   ! ARGUMENTS
   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8), intent(in) :: combine_matrix(bounds%begc: ,-nlevsno+1: ,-nlevsno+1: )

   !temporary variables
   real(r8) :: tracer_copy(bounds%begc:bounds%endc,-nlevsno+1:1)
   integer  :: jj, jj1

   SHR_ASSERT_ALL((ubound(combine_matrix,1)         == bounds%endc)  , errMsg(filename,__LINE__))
   SHR_ASSERT_ALL((ubound(combine_matrix,2)         == 1)      , errMsg(filename,__LINE__))
   SHR_ASSERT_ALL((ubound(combine_matrix,3)         == 1)      , errMsg(filename,__LINE__))

   associate(                             &
      tracer_conc_frozen_col => this%tracerstates%tracer_conc_frozen_col , &
      tracer_conc_mobile_col => this%tracerstates%tracer_conc_mobile_col , &
      tracer_conc_solid_passive_col => this%tracerstates%tracer_conc_solid_passive_col, &
      ntracers         => this%tracers%ntracers, &
      ngwmobile_tracers=> this%tracers%ngwmobile_tracers , &
      is_frozen => this%tracers%is_frozen &
   )


   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
        !make a copy
        call tracer_copy_a2b_comb(bounds, num_snowc, filter_snowc, col%snl, &
           tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj),tracer_copy)
        !remapping
        call tracer_col_mapping_comb(bounds, num_snowc, filter_snowc, col%snl, combine_matrix, &
           tracer_copy, tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj))
        if(is_frozen(jj))then
          !make a copy
          call tracer_copy_a2b_comb(bounds, num_snowc, filter_snowc, col%snl, &
            tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj),tracer_copy)
          !remapping
          call tracer_col_mapping_comb(bounds, num_snowc, filter_snowc, col%snl, combine_matrix, &
             tracer_copy, tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj))
        endif
     else
       !x!currently, no solid tracer is allowed to mix with snow
       !! jj1=jj-ngwmobile_tracers
       !x!tracer_conc_solid_passive_col
       !x!copy
       !x!remapping
     endif
   enddo

   end associate
   end subroutine tracer_CombineSnowLayers

   !------------------------------------------------------------------------
   subroutine Enter_tracer_LayerAdjustment(this, bounds, num_snowc, filter_snowc )
   !
   !! DESCRIPTION
   ! prepare tracer for entering mass adjustment in response to snow dynamics

   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_snowc(:) ! column filter for soil points

   !temporary variables
   integer :: fc, c, j, jj

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue

   associate(                             &
      tracer_conc_frozen_col => this%tracerstates%tracer_conc_frozen_col , &
      tracer_conc_mobile_col => this%tracerstates%tracer_conc_mobile_col , &
      tracer_conc_solid_passive_col => this%tracerstates%tracer_conc_solid_passive_col, &
      ntracers         => this%tracers%ntracers, &
      ngwmobile_tracers=> this%tracers%ngwmobile_tracers , &
      is_frozen => this%tracers%is_frozen , &
      snl => col%snl &
   )
   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
       do j = -nlevsno + 1, 1
         do fc = 1, num_snowc
           c = filter_snowc(fc)
           if(j >= snl(c) + 1)then
             tracer_conc_mobile_col(c,j,jj) = tracer_conc_mobile_col(c,j,jj) /col%dz(c,j)
             if(is_frozen(jj))then
               tracer_conc_frozen_col(c,j,jj) = tracer_conc_frozen_col(c,j,jj)/col%dz(c,j)
             endif
           endif
         enddo
       enddo
     else
       !no passive tracer adjustment at this moment, can be used for aerosols or mineral particles
     endif
   enddo
   end associate
   end subroutine Enter_tracer_LayerAdjustment

   !------------------------------------------------------------------------
   subroutine Exit_tracer_LayerAdjustment(this, bounds, num_snowc, filter_snowc )
   !! DESCRIPTION
   ! prepare tracer for exit mass adjustment in response to snow dynamics
   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_snowc(:) ! column filter for soil points

   !temporary variables
   integer :: fc, c, j, jj

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue

   associate(                             &
      tracer_conc_frozen_col => this%tracerstates%tracer_conc_frozen_col , &
      tracer_conc_mobile_col => this%tracerstates%tracer_conc_mobile_col , &
      tracer_conc_solid_passive_col => this%tracerstates%tracer_conc_solid_passive_col, &
      ntracers         => this%tracers%ntracers, &
      ngwmobile_tracers=> this%tracers%ngwmobile_tracers , &
      is_frozen => this%tracers%is_frozen , &
      snl => col%snl &
   )
   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
       do j = -nlevsno + 1, 1
         do fc = 1, num_snowc
           c = filter_snowc(fc)
           if(j >= snl(c) + 1)then
             tracer_conc_mobile_col(c,j,jj) = tracer_conc_mobile_col(c,j,jj) *col%dz(c,j)
             if(is_frozen(jj))then
               tracer_conc_frozen_col(c,j,jj) = tracer_conc_frozen_col(c,j,jj)*col%dz(c,j)
             endif
           endif
         enddo
       enddo
     else
       !no passive tracer adjustment at this moment, can be used for aerosols or mineral particles
     endif
   enddo
   end associate
   end subroutine Exit_tracer_LayerAdjustment
   !------------------------------------------------------------------------

   subroutine tracer_snowcapping(this, bounds, num_snowc, filter_snowc, waterflux_vars)

   !do tracer update due to snow capping
   use BeTR_WaterfluxType           , only          : betr_waterflux_type
   implicit none
   class(betr_type), intent(inout) :: this
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_snowc(:) ! column filter for soil points
   type(betr_waterflux_type), intent(in)   :: waterflux_vars

   ! remove compiler warnings about unused dummy args
   if (this%advection_on) continue
   if (bounds%begc > 0) continue
   if (num_snowc > 0) continue
   if (size(filter_snowc) > 0) continue
   if (size(waterflux_vars%qflx_drain_vr_col) > 0) continue

  end subroutine tracer_snowcapping

end module BetrType
