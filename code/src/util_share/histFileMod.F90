module histFileMod


   use shr_kind_mod   , only : r8 => shr_kind_r8
   use decompMod      , only : bounds_type
implicit none

  integer , private, parameter :: no_snow_MIN = 1                 ! minimum valid value for this flag
  integer , public , parameter :: no_snow_normal = 1              ! normal treatment, which should be used for most fields (use spval when snow layer not present)
  integer , public , parameter :: no_snow_zero = 2                ! average in a 0 value for times when the snow layer isn't present
  integer , private, parameter :: no_snow_MAX = 2                 ! maximum valid value for this flag
  integer , private, parameter :: no_snow_unset = no_snow_MIN - 1 ! flag specifying that field is NOT a multi-layer snow field

contains
  subroutine hist_addfld1d (fname, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_pft, ptr_lnd, &
                        ptr_atm, p2c_scale_type, c2l_scale_type, &
                        l2g_scale_type, set_lake, set_nolake, set_urb, set_nourb, &
                        set_noglcmec, set_spec, default)
    !
    ! !DESCRIPTION:
    ! Initialize a single level history field. The pointer, ptrhist,
    ! is a pointer to the clmtype array that the history buffer will use.
    ! The value of type1d passed to masterlist\_add\_fld determines which of the
    ! 1d type of the output and the beginning and ending indices the history
    ! buffer field). Default history contents for given field on all tapes
    ! are set by calling [masterlist\_make\_active] for the appropriate tape.
    ! After the masterlist is built, routine [htapes\_build] is called for an
    ! initial or branch run to initialize the actual history tapes.
    !
    ! !USES:

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)           :: fname          ! field name
    character(len=*), intent(in)           :: units          ! units of field
    character(len=1), intent(in)           :: avgflag        ! time averaging flag
    character(len=*), intent(in)           :: long_name      ! long name of field
    character(len=*), optional, intent(in) :: type1d_out     ! output type (from clmtype)
    real(r8)        , optional, pointer    :: ptr_gcell(:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_pft(:)     ! pointer to pft array
    real(r8)        , optional, pointer    :: ptr_lnd(:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_atm(:)     ! pointer to atm array
    real(r8)        , optional, intent(in) :: set_lake       ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake     ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb        ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb      ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_noglcmec   ! value to set non-glacier_mec to
    real(r8)        , optional, intent(in) :: set_spec       ! value to set special to
    character(len=*), optional, intent(in) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,g                 ! indices
    integer :: hpindex                 ! history buffer pointer index
    character(len=8) :: l_type1d       ! 1d data type
    character(len=8) :: l_type1d_out   ! 1d output type
    character(len=8) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
    type(bounds_type):: bounds         ! boudns 
    character(len=16):: l_default      ! local version of 'default'
    character(len=*),parameter :: subname = 'hist_addfld1d'
    
    
    
  end subroutine hist_addfld1d
    
    
  subroutine hist_addfld2d (fname, type2d, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_pft, ptr_lnd, ptr_atm, &
                        p2c_scale_type, c2l_scale_type, l2g_scale_type, &
                        set_lake, set_nolake, set_urb, set_nourb, set_spec, &
                        no_snow_behavior, default)
    !
    ! !DESCRIPTION:
    ! Initialize a single level history field. The pointer, ptrhist,
    ! is a pointer to the clmtype array that the history buffer will use.
    ! The value of type1d passed to masterlist\_add\_fld determines which of the
    ! 1d type of the output and the beginning and ending indices the history
    ! buffer field). Default history contents for given field on all tapes
    ! are set by calling [masterlist\_make\_active] for the appropriatae tape.
    ! After the masterlist is built, routine [htapes\_build] is called for an
    ! initial or branch run to initialize the actual history tapes.
    !
    ! !USES:

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname                      ! field name
    character(len=*), intent(in) :: type2d                     ! 2d output type
    character(len=*), intent(in) :: units                      ! units of field
    character(len=1), intent(in) :: avgflag                    ! time averaging flag
    character(len=*), intent(in) :: long_name                  ! long name of field
    character(len=*), optional, intent(in) :: type1d_out       ! output type (from clmtype)
    real(r8)        , optional, pointer    :: ptr_atm(:,:)     ! pointer to atm array
    real(r8)        , optional, pointer    :: ptr_lnd(:,:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_gcell(:,:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:,:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:,:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_pft(:,:)     ! pointer to pft array
    real(r8)        , optional, intent(in) :: set_lake         ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake       ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb          ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb        ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_spec         ! value to set special to
    integer         , optional, intent(in) :: no_snow_behavior ! if a multi-layer snow field, behavior to use for absent snow layers (should be one of the public no_snow_* parameters defined above)
    character(len=*), optional, intent(in) :: p2c_scale_type   ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type   ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type   ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default          ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,g                 ! indices
    integer :: num2d                   ! size of second dimension (e.g. number of vertical levels)
    integer :: hpindex                 ! history buffer index
    character(len=8) :: l_type1d         ! 1d data type
    character(len=8) :: l_type1d_out     ! 1d output type
    character(len=8) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
    type(bounds_type):: bounds         ! bounds 
    character(len=16):: l_default      ! local version of 'default'
    character(len=*),parameter :: subname = 'hist_addfld2d'
    
    
  end subroutine hist_addfld2d
    
end module histFileMod    
