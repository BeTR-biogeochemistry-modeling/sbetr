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
  use clm_time_manager   , only : get_nstep
  use EcophysConType, only : ecophyscon_type
  use BetrBGCMod
  implicit none

  private

  character(len=*), parameter :: filename = __FILE__

  type, public :: betr_type
     character(len=betr_string_length) :: reaction_method

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

     ! FIXME(bja, 201603) replace LSM specific types!

   contains
     procedure, public :: Init
     procedure, public :: step_without_drainage
     procedure, public :: step_with_drainage
     procedure, public :: calc_dew_sub_flux
     procedure, public :: check_mass_err
     procedure, private :: ReadNamelist
  end type betr_type

contains

!-------------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, bounds, waterstate, cnstate, ecophyscon)

    use abortutils, only : endrun
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

    lbj = bounds%lbj
    ubj = bounds%ubj

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

  end subroutine Init

!-------------------------------------------------------------------------------
  subroutine ReadNamelist(this, namelist_buffer)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use clm_varctl   , only : iulog
    use abortutils      , only : endrun
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


    !-----------------------------------------------------------------------

    namelist / betr_parameters / reaction_method

    reaction_method = ''

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

  end subroutine ReadNamelist

  !-------------------------------------------------------------------------------
  subroutine step_without_drainage(this, bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp,         &
       atm2lnd_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, waterflux_vars, temperature_vars,&
       chemstate_vars, cnstate_vars, canopystate_vars,  carbonflux_vars, betrtracer_vars, &
       bgc_reaction, betr_aerecond_vars, tracerboundarycond_vars, tracercoeff_vars, tracerstate_vars, &
       tracerflux_vars, plant_soilbgc_coupler)
    !
    ! !DESCRIPTION:
    ! run betr code one time step forward, without drainage calculation

    ! !USES:
    use clm_time_manager             , only          : get_step_size
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
    use BeTR_ColumnType              , only          : col => betr_col
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
    integer                          , intent(in)    :: lbj, ubj                   ! lower and upper bounds, make sure they are > 0

    type(betr_waterstate_type)       , intent(in)    :: waterstate_vars            ! water state variables
    type(betr_soilstate_type)        , intent(in)    :: soilstate_vars             ! column physics variable
    type(betr_temperature_type)      , intent(in)    :: temperature_vars           ! energy state variable
    type(betr_chemstate_type)        , intent(in)    :: chemstate_vars
    class(betrtracer_type)           , intent(in)    :: betrtracer_vars            ! betr configuration information
    class(bgc_reaction_type)         , intent(in)    :: bgc_reaction
    type(betr_atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(betr_soilhydrology_type)    , intent(in)    :: soilhydrology_vars
    type(betr_cnstate_type)          , intent(inout) :: cnstate_vars
    type(betr_canopystate_type)      , intent(in)    :: canopystate_vars
    type(betr_carbonflux_type)       , intent(in)    :: carbonflux_vars
    type(betr_aerecond_type)         , intent(inout) :: betr_aerecond_vars
    type(betr_waterflux_type)        , intent(inout) :: waterflux_vars
    type(tracerboundarycond_type)    , intent(inout) :: tracerboundarycond_vars
    type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    class(plant_soilbgc_type)        , intent(inout) :: plant_soilbgc_coupler !

    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'run_betr_one_step_without_drainage'
    real(r8)           :: dtime2, dtime
    real(r8)           :: Rfactor(bounds%begc:bounds%endc, lbj:ubj,1:betrtracer_vars%ngwmobile_tracers) !retardation factor
    integer            :: j


    dtime = get_step_size()

    !initialize extra parameters
    dtime2 = dtime * 0.5_r8


    call stage_tracer_transport(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         carbonflux_vars, soilstate_vars, waterstate_vars, waterflux_vars, temperature_vars, soilhydrology_vars, &
         chemstate_vars, betr_aerecond_vars, canopystate_vars, betrtracer_vars, tracercoeff_vars, &
         tracerboundarycond_vars, tracerflux_vars, bgc_reaction, Rfactor)

    call surface_tracer_hydropath_update(bounds, num_soilc, filter_soilc, &
       waterstate_vars, waterflux_vars, soilhydrology_vars, betrtracer_vars, tracerstate_vars, &
       tracercoeff_vars,  tracerflux_vars)

    call bgc_reaction%calc_bgc_reaction(bounds, lbj, ubj, &
         num_soilc,                                       &
         filter_soilc,                                    &
         num_soilp,                                       &
         filter_soilp,                                    &
         tracerboundarycond_vars%jtops_col,               &
         dtime,                                           &
         betrtracer_vars,                                 &
         tracercoeff_vars,                                &
         cnstate_vars,                                    &
         tracerstate_vars,                                &
         tracerflux_vars,                                 &
         tracerboundarycond_vars,                         &
         plant_soilbgc_coupler)

    call tracer_gws_transport(bounds, num_soilc, filter_soilc, Rfactor, waterstate_vars, &
      waterflux_vars, betrtracer_vars, tracerboundarycond_vars, tracercoeff_vars, &
      tracerstate_vars, tracerflux_vars, bgc_reaction)

    call calc_ebullition(bounds, 1, ubj,                                 &
         tracerboundarycond_vars%jtops_col,                              &
         num_soilc,                                                      &
         filter_soilc,                                                   &
         atm2lnd_vars%forc_pbot_downscaled_col,                          &
         col%zi(bounds%begc:bounds%endc, 0:ubj),                         &
         col%dz(bounds%begc:bounds%endc, 1:ubj),                         &
         dtime,                                                          &
         soilhydrology_vars%fracice_col(bounds%begc:bounds%endc, 1:ubj), &
         soilhydrology_vars%zwts_col(bounds%begc:bounds%endc),           &
         betrtracer_vars,                                                &
         tracercoeff_vars,                                               &
         tracerstate_vars,                                               &
         tracerflux_vars%tracer_flx_ebu_col(bounds%begc:bounds%endc, 1:betrtracer_vars%nvolatile_tracers))

    if (is_active_betr_bgc) then
       !update nitrogen storage pool
       call plant_soilbgc_coupler%plant_soilbgc_summary(bounds, lbj, ubj, num_soilc,                   &
            filter_soilc,                                                                              &
            col%dz(bounds%begc:bounds%endc,1:ubj),                                                     &
            betrtracer_vars, tracerflux_vars)
    endif

  end subroutine step_without_drainage


  !--------------------------------------------------------------------------------
  subroutine step_with_drainage(this, bounds, lbj, ubj, num_soilc, filter_soilc, &
       jtops, waterflux_vars, betrtracer_vars, tracercoeff_vars, tracerstate_vars,  tracerflux_vars)
    !
    ! !DESCRIPTION:
    ! do tracer update due to drainage
    !
    ! !USES:
    use tracerfluxType        , only : tracerflux_type
    use tracerstatetype       , only : tracerstate_type
    use tracercoeffType       , only : tracercoeff_type
    use BeTR_ColumnType       , only : col => betr_col
    use MathfuncMod           , only : safe_div
    use BeTR_WaterFluxType    , only : betr_waterflux_type
    ! !ARGUMENTS:
    class(betr_type), intent(inout) :: this
    type(bounds_type),        intent(in)    :: bounds
    integer,                  intent(in)    :: lbj, ubj
    integer,                  intent(in)    :: num_soilc                          ! number of columns in column filter_soilc
    integer,                  intent(in)    :: filter_soilc(:)                    ! column filter_soilc
    integer,                  intent(in)    :: jtops(bounds%begc: )
    type(betr_waterflux_type)    , intent(in)    :: waterflux_vars
    class(betrtracer_type),    intent(in)    :: betrtracer_vars                    ! betr configuration information
    type(tracercoeff_type),   intent(in)    :: tracercoeff_vars                   ! tracer phase conversion coefficients
    type(tracerflux_type),    intent(inout) :: tracerflux_vars
    type(tracerstate_type),   intent(inout) :: tracerstate_vars                   ! tracer state variables data structure

    ! !LOCAL VARIABLES:
    real(r8) :: aqucon
    integer  :: fc, c, j, k


    SHR_ASSERT_ALL((ubound(jtops)         == (/bounds%endc/))      , errMsg(filename,__LINE__))

    associate(                                                                   & !
         ngwmobile_tracers        => betrtracer_vars%ngwmobile_tracers         , & !
         groupid                  => betrtracer_vars%groupid                   , & !
         is_h2o                   => betrtracer_vars%is_h2o                    , & !
         is_advective             => betrtracer_vars%is_advective              , & !
         aqu2bulkcef_mobile       => tracercoeff_vars%aqu2bulkcef_mobile_col   , & !
         tracer_conc_mobile       => tracerstate_vars%tracer_conc_mobile_col   , & !
         tracer_conc_grndwater    => tracerstate_vars%tracer_conc_grndwater_col, & !
         dz                       => col%dz                                    , & !
         tracer_flx_drain         => tracerflux_vars%tracer_flx_drain_col      , & !
         qflx_drain_vr            => waterflux_vars%qflx_drain_vr_col          , & ! Output  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/step) (to river +)
         qflx_totdrain            => waterflux_vars%qflx_totdrain_col            & ! Output  : [real(r8) (:,:) ]  (m H2o/step)

         )
      if(get_nstep()==8)return
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
           betrtracer_vars, tracercoeff_vars, tracerstate_vars)

    end associate
  end subroutine step_with_drainage


  !--------------------------------------------------------------------------------
  subroutine calc_dew_sub_flux(this, bounds, num_hydrologyc, filter_soilc_hydrologyc, &
       waterstate_vars, waterflux_vars, betrtracer_vars, tracerflux_vars, tracerstate_vars)
    !
    ! DESCRIPTION:
    ! calculate water flux from dew formation, and sublimation
    ! !USES:
    use clm_time_manager      , only : get_step_size
    use BeTR_ColumnType       , only : col => betr_col
    use BeTR_LandunitType     , only : lun => betr_lun
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
end module BetrType
