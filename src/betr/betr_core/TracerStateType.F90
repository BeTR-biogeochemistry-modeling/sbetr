module TracerStateType
  !
  ! !DESCRIPTION:
  !  data type for state variables used in betr
  !
  ! !USES:
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  use BeTR_decompMod      , only : bounds_type  => betr_bounds_type
  use BeTR_LandunitType   , only : lun => betr_lun
  use BeTR_ColumnType     , only : col => betr_col
  use betr_ctrl           , only : iulog => biulog
  use betr_varcon         , only : spval => bspval, ispval => bispval
  use BeTR_landvarconType , only : landvarcon => betr_landvarcon
  use MathfuncMod         , only : dot_sum
  use TracerBaseType      , only : tracerbase_type
  !
  ! !PUBLIC TYPES:
  implicit none

  private
  !
  ! !PUBLIC DATA:
  !

  type, public, extends(tracerbase_type) ::  TracerState_type
     ! Column tracer state variables
     real(r8), pointer :: tracer_conc_surfwater_col     (:,:)      !tracer concentration in the hydraulic head
     real(r8), pointer :: tracer_conc_aquifer_col       (:,:)      !tracer concentration in the flux to aquifer
     real(r8), pointer :: tracer_conc_grndwater_col     (:,:)      !tracer concentration in the flux to groundwater, include lateral drainage and discharge to aquifer
     real(r8), pointer :: tracer_col_molarmass_col      (:,:)      !for error tracking, column tracer mass
     real(r8), pointer :: tracer_conc_atm_col           (:,:)      !colum volatile tracer in the atmosphere
     real(r8), pointer :: tracer_P_gas_col              (:,:)      !total gas pressure at different depth, sum of different gas species.
     real(r8), pointer :: tracer_P_gas_frac_col         (:,:,:)    !fraction of the volatile species in the overall pressure
     real(r8), pointer :: tracer_soi_molarmass_col      (:,:)      !vertically integrated tracer content (mol tracer/m2), only in the soil
     real(r8), pointer :: tracer_conc_mobile_col        (:,:,:)    !tracer concentration in each layer (mol/m3) (snow/ponding water + soil)
     real(r8), pointer :: tracer_conc_solid_equil_col   (:,:,:)    !tracer concentration in adsorbed/solid phase for each layer (mol/m3) (soil), which is in equilibrium with mobile phase
     real(r8), pointer :: tracer_conc_solid_passive_col (:,:,:)    !tracer concentration in passive solid phase, which is not in equilibrium with mobile phase. e.g. polymers, or protected monomers, or ice
     real(r8), pointer :: tracer_conc_frozen_col        (:,:,:)    !place holder, tracer concentration in frozen layer for unsaturated part for nonvolatile species
     !real(r8), pointer :: tracer_conc_bubble_col        (:,:,:)    !place holder, a bubble pool to track the lake ebullition in freeze-thaw period, [col, levels, tracer]
     real(r8), pointer :: beg_tracer_molarmass_col      (:,:)      !column integrated tracer mass
     real(r8), pointer :: end_tracer_molarmass_col      (:,:)      !column integrated tracer mass
     real(r8), pointer :: errtracer_col                 (:,:)      !column mass balance error

   contains
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: Reset
     procedure, public  :: int_mass_mobile_col
     procedure, public  :: int_mass_frozen_col
     procedure, public  :: int_mass_adsorb_col
     procedure, public  :: int_mass_solid_col
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, public  :: retrieve_hist
  end type TracerState_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! initialize the data type

    ! !USES:
    use BeTRTracerType, only : BeTRTracer_Type

    implicit none
    ! !ARGUMENTS:
    class(TracerState_type)           :: this
    type(bounds_type)    , intent(in) :: bounds
    integer              , intent(in) :: lbj, ubj
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars

    call this%InitAllocate(bounds, lbj, ubj, betrtracer_vars)
    call this%tracer_base_init()
    call this%InitHistory(bounds, betrtracer_vars)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! allocate memory for arraies in the data type

    ! !USES:
    use BeTRTracerType, only : BeTRTracer_Type
    implicit none
    !
    ! !ARGUMENTS:
    class(TracerState_type)           :: this
    type(bounds_type), intent(in)     :: bounds
    integer, intent(in)               :: lbj, ubj
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: ngwmobile_tracers, ntracers
    integer :: nsolid_equil_tracers
    integer :: nsolid_passive_tracers
    integer :: nvolatile_tracers
    integer :: nfrozen_tracers
    !---------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    ngwmobile_tracers      = betrtracer_vars%ngwmobile_tracers
    ntracers               =  betrtracer_vars%ntracers
    nsolid_equil_tracers   = betrtracer_vars%nsolid_equil_tracers
    nsolid_passive_tracers = betrtracer_vars%nsolid_passive_tracers
    nvolatile_tracers      = betrtracer_vars%nvolatile_tracers
    nfrozen_tracers        = betrtracer_vars%nfrozen_tracers
    allocate(this%tracer_P_gas_col              (begc:endc, lbj:ubj))             ; this%tracer_P_gas_col         (:,:) = nan
    allocate(this%tracer_conc_surfwater_col     (begc:endc, 1:ngwmobile_tracers)) ; this%tracer_conc_surfwater_col(:,:) = nan
    allocate(this%tracer_conc_aquifer_col       (begc:endc, 1:ngwmobile_tracers)) ; this%tracer_conc_aquifer_col  (:,:) = nan
    allocate(this%tracer_conc_grndwater_col     (begc:endc, 1:ngwmobile_tracers)) ; this%tracer_conc_grndwater_col(:,:) = nan
    allocate(this%tracer_col_molarmass_col      (begc:endc, 1:ntracers))          ; this%tracer_col_molarmass_col (:,:) = nan
    allocate(this%tracer_soi_molarmass_col      (begc:endc, 1:ntracers))          ; this%tracer_soi_molarmass_col (:,:) = nan
    allocate(this%errtracer_col                 (begc:endc, 1:ntracers))          ; this%errtracer_col            (:,:) = nan
    allocate(this%tracer_conc_atm_col           (begc:endc, 1:nvolatile_tracers))
    this%tracer_conc_atm_col      (:,:) = nan
    allocate(this%tracer_conc_mobile_col        (begc:endc, lbj:ubj, 1:ngwmobile_tracers))
    this%tracer_conc_mobile_col       (:,:,:) =  nan

    allocate(this%tracer_conc_solid_equil_col   (begc:endc, lbj:ubj, 1:nsolid_equil_tracers))
    this%tracer_conc_solid_equil_col  (:,:,:) = nan

    allocate(this%tracer_conc_solid_passive_col (begc:endc, lbj:ubj, 1:nsolid_passive_tracers))
    this%tracer_conc_solid_passive_col(:,:,:) = nan

    allocate(this%tracer_P_gas_frac_col         (begc:endc, lbj:ubj, 1:nvolatile_tracers))
    this%tracer_P_gas_frac_col        (:,:,:) = nan
    allocate(this%tracer_conc_frozen_col        (begc:endc, lbj:ubj, 1:nfrozen_tracers))
    this%tracer_conc_frozen_col (:,:,:) = nan
    allocate(this%beg_tracer_molarmass_col      (begc:endc, 1:ntracers))
    this%beg_tracer_molarmass_col(:,:) = nan
    allocate(this%end_tracer_molarmass_col      (begc:endc, 1:ntracers))
    this%end_tracer_molarmass_col(:,:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! History fields initialization
    !
    ! !USES:
    use BeTRTracerType, only: BeTRTracer_Type
    !
    ! !ARGUMENTS:
    class(TracerState_type)           :: this
    type(bounds_type)    , intent(in) :: bounds
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: jj, kk

    associate(                                                       &
         ntracers          =>  betrtracer_vars%ntracers            , &
         ngwmobile_tracers =>  betrtracer_vars%ngwmobile_tracers   , &
         is_volatile       =>  betrtracer_vars%is_volatile         , &
         is_isotope        =>  betrtracer_vars%is_isotope          , &
         is_h2o            =>  betrtracer_vars%is_h2o              , &
         volatileid        =>  betrtracer_vars%volatileid          , &
         tracernames       =>  betrtracer_vars%tracernames         , &
         is_frozen         =>  betrtracer_vars%is_frozen           , &
         frozenid          =>  betrtracer_vars%frozenid              &
         )

      this%num_hist1d=0
      this%num_hist2d=0
      call this%add_hist_var2d(fname='TRACER_P_GAS', units='Pa', type2d='levtrc',  &
           avgflag='A', long_name='total gas pressure')

      do jj = 1, ntracers
         if(jj<= ngwmobile_tracers)then

            call this%add_hist_var1d (fname=trim(tracernames(jj))//'_TRACER_CONC_SURFWATER', units='mol m-3', &
                 avgflag='A', long_name='head concentration for tracer '//trim(tracernames(jj)), &
                 default='inactive')

            call this%add_hist_var1d (fname=trim(tracernames(jj))//'_TRACER_CONC_AQUIFER', units='mol m-3', &
                 avgflag='A', long_name='quifier concentration for tracer '//trim(tracernames(jj)), &
                 default='inactive')

            call this%add_hist_var1d (fname=trim(tracernames(jj))//'_TRACER_CONC_GRNDWATER', units='mol m-3', &
                 avgflag='A', long_name='groundwater concentration for tracer '//trim(tracernames(jj)), &
                 default='inactive')

            call this%add_hist_var2d (fname=trim(tracernames(jj))//'_TRACER_CONC_MOBILE', units='mol m-3', type2d='levtrc',  &
                 avgflag='A', long_name='gw-mobile phase for tracer '//trim(tracernames(jj)))

            if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then
               call this%add_hist_var2d (fname=trim(tracernames(jj))//'_TRACER_P_GAS_FRAC', units='none', type2d='levtrc',  &
                    avgflag='A', long_name='fraction of gas phase contributed by '//trim(tracernames(jj)))
            endif

            if(is_frozen(jj))then
               call this%add_hist_var2d (fname=trim(tracernames(jj))//'_TRACER_CONC_FROZEN', units='mol m-3', type2d='levtrc',  &
                    avgflag='A', long_name='frozen phase for tracer '//trim(tracernames(jj)))
            endif
         else
            kk = jj - ngwmobile_tracers
            call this%add_hist_var2d (fname=trim(tracernames(jj))//'TRACER_CONC_SOLID_PASSIVE', units='mol m-3', type2d='levtrc',  &
                 avgflag='A', long_name='passive solid phase for tracer '//trim(tracernames(jj)), &
                 default='inactive')
         endif
         call this%add_hist_var1d (fname=trim(tracernames(jj))//'_TRCER_SOI_MOLAMASS', units='mol m-2', &
              avgflag='A', long_name='total molar mass in soil for '//trim(tracernames(jj)), &
              default='inactive')

         call this%add_hist_var1d (fname=trim(tracernames(jj))//'_TRCER_COL_MOLAMASS', units='mol m-2', &
              avgflag='A', long_name='total molar mass in the column (soi+snow) for '//trim(tracernames(jj)), &
              default='inactive')

      enddo

    end associate
  end subroutine InitHistory

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use betr_ctrl      , only : iulog  => biulog
    use BeTRTracerType , only : BeTRTracer_Type
    use restUtilMod    , only : restartvar
    use ncdio_pio      , only : file_desc_t
    use ncdio_pio      , only : ncd_double
    !
    implicit none
    ! !ARGUMENTS:
    class(TracerState_type)              :: this
    type(bounds_type)    , intent(in)    :: bounds
    class(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
    character(len=*)     , intent(in)    :: flag                                         ! 'read' or 'write'
    type(BeTRTracer_Type), intent(in)    :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: j,c,jj,kk ! indices
    logical :: readvar      ! determine if variable is on initial file
    real(r8), pointer :: ptr1d(:)
    real(r8), pointer :: ptr2d(:,:)

    ! remove unused dummy arg compiler warning
    if (bounds%begc > 0) continue

    associate(                                                       &
         ntracers          =>  betrtracer_vars%ntracers            , &
         ngwmobile_tracers =>  betrtracer_vars%ngwmobile_tracers   , &
         is_adsorb         =>  betrtracer_vars%is_adsorb           , &
         adsorbid          =>  betrtracer_vars%adsorbid            , &
         tracernames       =>  betrtracer_vars%tracernames         , &
         is_frozen         =>  betrtracer_vars%is_frozen           , &
         frozenid          =>  betrtracer_vars%frozenid              &
         )

      do jj = 1, ntracers
         if(jj<= ngwmobile_tracers)then

            ptr1d => this%tracer_conc_aquifer_col(:, jj)
            call restartvar(ncid=ncid, flag=flag, varname=trim(tracernames(jj))//'_TRACER_CONC_AQUIFER', &
                 xtype=ncd_double,  dim1name='column', long_name='',  units='', &
                 interpinic_flag='interp' , readvar=readvar, data=ptr1d)

            ptr2d => this%tracer_conc_mobile_col(:, :, jj)
            call restartvar(ncid=ncid, flag=flag, varname=trim(tracernames(jj))//'_TRACER_CONC_MOIBLE', xtype=ncd_double,  &
                 dim1name='column',dim2name='levtrc', switchdim=.true., &
                 long_name='',  units='', fill_value=spval, &
                 interpinic_flag='interp', readvar=readvar, data=ptr2d)

            if(is_adsorb(jj))then
               ptr2d => this%tracer_conc_solid_equil_col(:, :, adsorbid(jj))
               call restartvar(ncid=ncid, flag=flag, varname=trim(tracernames(jj))//'_TRACER_CONC_SOLID_EQUIL', xtype=ncd_double,  &
                    dim1name='column',dim2name='levtrc', switchdim=.true., &
                    long_name='',  units='', fill_value=spval, &
                    interpinic_flag='interp', readvar=readvar, data=ptr2d)
            endif
            if(is_frozen(jj))then
              ptr2d => this%tracer_conc_frozen_col(:, :, frozenid(jj))
              call restartvar(ncid=ncid, flag=flag, varname=trim(tracernames(jj))//'_TRACER_CONC_FROZEN', xtype=ncd_double,  &
                   dim1name='column',dim2name='levtrc', switchdim=.true., &
                   long_name='',  units='', fill_value=spval, &
                   interpinic_flag='interp', readvar=readvar, data=ptr2d)
            endif
         else
            kk = jj - ngwmobile_tracers
            ptr2d => this%tracer_conc_solid_passive_col(:, :, kk)
            call restartvar(ncid=ncid, flag=flag, varname=trim(tracernames(jj))//'TRACER_CONC_SOLID_PASSIVE', xtype=ncd_double,  &
                 dim1name='column',dim2name='levtrc', switchdim=.true., &
                 long_name='',  units='', &
                 interpinic_flag='interp', readvar=readvar, data=ptr2d)
         endif

      enddo
    end associate
  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Reset(this, column)
    !
    ! !DESCRIPTION:
    !  reset state variables
    !
    ! !ARGUMENTS:
    class(TracerState_type)             :: this
    integer               , intent(in)  :: column     ! column index
    !-----------------------------------------------------------------------

    ! remove unused dummy arg compiler warnings
    if (size(this%tracer_conc_surfwater_col) > 0) continue
    if (column > 0) continue

  end subroutine Reset

  !-----------------------------------------------------------------------
  function int_mass_mobile_col(this, lbj, ubj, c, j, dz)result(int_mass)
  !DESCRIPTION
  !integrate mobile tracer mass, gas+aqueous
  class(TracerState_type) :: this
  integer, intent(in)     :: lbj, ubj
  integer, intent(in)     :: c, j
  real(r8), intent(in)    :: dz(lbj:ubj)
  real(r8)                :: int_mass

  int_mass = dot_sum(this%tracer_conc_mobile_col(c,lbj:ubj,j), dz)

  end function int_mass_mobile_col

  !-----------------------------------------------------------------------
  function int_mass_frozen_col(this, lbj, ubj, c, j, dz)result(int_mass)
  !DESCRIPTION
  !integrate frozen tracer mass
  class(TracerState_type) :: this
  integer, intent(in)     :: lbj, ubj
  integer, intent(in)     :: c, j
  real(r8), intent(in)    :: dz(lbj:ubj)
  real(r8)                :: int_mass

  int_mass = dot_sum(this%tracer_conc_frozen_col(c,lbj:ubj,j), dz)

  end function int_mass_frozen_col


  !-----------------------------------------------------------------------
  function int_mass_adsorb_col(this, lbj, ubj, c, j, dz)result(int_mass)
  !DESCRIPTION
  !integrate adsorbed tracer mass
  class(TracerState_type) :: this
  integer, intent(in)     :: lbj, ubj
  integer, intent(in)     :: c, j
  real(r8), intent(in)    :: dz(lbj:ubj)
  real(r8)                :: int_mass

  int_mass = dot_sum(this%tracer_conc_solid_equil_col(c,lbj:ubj,j), dz)

  end function int_mass_adsorb_col

  !-----------------------------------------------------------------------
  function int_mass_solid_col(this, lbj, ubj, c, j, dz)result(int_mass)
  !DESCRIPTION
  !integrate solid tracer mass
  class(TracerState_type) :: this
  integer, intent(in)     :: lbj, ubj
  integer, intent(in)     :: c, j
  real(r8), intent(in)    :: dz(lbj:ubj)
  real(r8)                :: int_mass

  int_mass = dot_sum(this%tracer_conc_solid_passive_col(c,lbj:ubj,j), dz)

  end function int_mass_solid_col

  !----------------------------------------------------------------
  subroutine retrieve_hist(this, bounds, lbj, ubj, state_2d, state_1d, betrtracer_vars)
  !
  !DESCRIPTION
  !retrieve data for history output
  use MathfuncMod, only : addone
  use BeTRTracerType , only : BeTRTracer_Type
  implicit none
  class(TracerState_type) :: this
  type(bounds_type)    , intent(in)  :: bounds
  integer, intent(in) :: lbj, ubj
  real(r8), intent(inout) :: state_2d(bounds%begc:bounds%endc, lbj:ubj,1:this%num_hist2d)
  real(r8), intent(inout) :: state_1d(bounds%begc:bounds%endc, 1:this%num_hist1d)
  type(BeTRTracer_Type)  , intent(in)  :: betrtracer_vars
  integer :: begc, endc
  integer :: jj, kk
  integer :: idtemp1d, idtemp2d

  associate(                                                       &
         ntracers          =>  betrtracer_vars%ntracers            , &
         ngwmobile_tracers =>  betrtracer_vars%ngwmobile_tracers   , &
         is_volatile       =>  betrtracer_vars%is_volatile         , &
         is_isotope        =>  betrtracer_vars%is_isotope          , &
         is_h2o            =>  betrtracer_vars%is_h2o              , &
         volatileid        =>  betrtracer_vars%volatileid          , &
         tracernames       =>  betrtracer_vars%tracernames         , &
         is_frozen         =>  betrtracer_vars%is_frozen           , &
         frozenid          =>  betrtracer_vars%frozenid              &
   )
  begc = bounds%begc; endc=bounds%endc
  idtemp1d = 0; idtemp2d = 0
  state_2d(begc:endc, lbj:ubj, addone(idtemp2d))= this%tracer_P_gas_col(begc:endc, lbj:ubj)

  do jj = 1, ntracers
    if(jj<= ngwmobile_tracers)then

      state_1d(begc:endc,addone(idtemp1d)) = this%tracer_conc_surfwater_col(begc:endc,jj)

      state_1d(begc:endc, addone(idtemp1d)) =this%tracer_conc_aquifer_col(begc:endc, jj)

      state_1d(begc:endc, addone(idtemp1d))=this%tracer_conc_grndwater_col(begc:endc, jj)

      state_2d(begc:endc, lbj:ubj, addone(idtemp2d))= this%tracer_conc_mobile_col(begc:endc, lbj:ubj, jj)

      if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then
        state_2d(begc:endc, lbj:ubj, addone(idtemp2d))= this%tracer_P_gas_frac_col(begc:endc,lbj:ubj, volatileid(jj))
      endif

      if(is_frozen(jj))then
        state_2d(begc:endc, lbj:ubj, addone(idtemp2d)) = this%tracer_conc_frozen_col(begc:endc,lbj:ubj, frozenid(jj))
      endif
    else
      kk = jj - ngwmobile_tracers
      state_2d(begc:endc, lbj:ubj, addone(idtemp2d)) = this%tracer_conc_solid_passive_col(begc:endc, lbj:ubj, kk)
    endif

    state_1d(begc:endc, addone(idtemp1d))=this%tracer_soi_molarmass_col(begc:endc, jj)

    state_1d(begc:endc, addone(idtemp1d))= this%tracer_col_molarmass_col(begc:endc, jj)

  enddo

  end associate
  end subroutine retrieve_hist
end module TracerStateType
