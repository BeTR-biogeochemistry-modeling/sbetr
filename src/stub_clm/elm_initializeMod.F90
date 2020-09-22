module elm_initializeMod

  implicit none
  save
  public

  contains

  subroutine initialize(bounds)
    !
    ! !DESCRIPTION:
    ! CLM initialization
  use decompMod       , only : bounds_type
  use elm_instMod     , only : elm_inst
  implicit none
  type(bounds_type), intent(in) :: bounds

  call elm_inst(bounds)

  end subroutine initialize
end module elm_initializeMod
