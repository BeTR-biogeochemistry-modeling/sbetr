program main
!DESCRIPTION
!test the math subroutines
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use InterpolationMod, only: cmass_interp
  use BetrStatusType   , only : betr_status_type
implicit none

  integer :: nx, ny, nxi
  real(r8), allocatable :: x(:)
  real(r8), allocatable :: Y(:,:)
  real(r8), allocatable :: xi(:)
  real(r8), allocatable :: Yi(:,:)

  type(betr_status_type) :: bstatus
  integer :: j1, j2

  nx = 10; ny = 3; nxi=5

  allocate(x(0:nx))
  allocate(Y(0:nx,ny))
  allocate(xi(0:nxi))
  allocate(Yi(0:nxi,ny))


  do j1 = 0, nx
    x(j1)=j1*0.1_r8
  enddo

  do j1 = 0, nxi
    xi(j1)=x(j1)+0.5_r8
  enddo

  do j2=1,ny
    do j1 = 0, nx
      y(j1,j2)=j2
    enddo
  enddo

  call cmass_interp(nx+1, x(0:nx), ny, Y(0:nx,1:ny), nxi+1, xi(0:nxi), YI(0:nxi,1:ny), bstatus)

  do j1 = 0, nxi
    print*,(yi(j1,j2),j2=1,ny)
  enddo

  deallocate(x)
  deallocate(y)
  deallocate(xi)
  deallocate(yi)

end program main
