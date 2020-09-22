subroutine spm_test


  use bshr_kind_mod, only : r8 => shr_kind_r8
  use SparseMatMod, only : spm_list_type, sparseMat_type, spm_axpy
  use SparseMatMod, only : spm_list_init, spm_list_insert, spm_list_to_mat
  use SparseMatMod, only : spm_print, set_spm_elm

  implicit none
  class(sparseMat_type), pointer :: spm
  type(spm_list_type), pointer :: spm_list
  integer :: nelms
  real(r8) :: x(5),y(5)
  integer :: errinfo

  call spm_list_init(spm_list, 1._r8, 1, 1, nelms)

  call spm_list_insert(spm_list, -1._r8, 2, 1, nelms)

  call spm_list_insert(spm_list, -3._r8, 4, 1, nelms)

  call spm_list_insert(spm_list, -2._r8, 1, 2, nelms)

  call spm_list_insert(spm_list, 5._r8, 2, 2, nelms)

  call spm_list_insert(spm_list, 4._r8, 3, 3, nelms)

  call spm_list_insert(spm_list, 6._r8, 4, 3, nelms)

  call spm_list_insert(spm_list, 4._r8, 5, 3, nelms)

  call spm_list_insert(spm_list, -4._r8, 1, 4, nelms)

  call spm_list_insert(spm_list, 2._r8, 3, 4, nelms)

  call spm_list_insert(spm_list, 7._r8, 4, 4, nelms)

  call spm_list_insert(spm_list, 8._r8, 2, 5, nelms)

  call spm_list_insert(spm_list, -5._r8, 5, 5, nelms)

  call spm_list_to_mat(spm_list, spm, nelms, 5)

  print*,'spm_disp1'
  call spm_print(spm)

  x=(/1._r8,1._r8,1._r8,1._r8,1._r8/)
  y=0._r8
  call spm_axpy(5, 5, 1._r8, x, spm, y, errinfo)
  print*,'errinfo',errinfo
  print*,y
  call set_spm_elm(spm, 3, 4, 100._r8)
  print*,'spm_disp2'
  call spm_print(spm)
end subroutine spm_test


program main
!DESCRIPTION
!test the math subroutines
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use InterpolationMod, only: cmass_interp
  use BetrStatusType   , only : betr_status_type
  use MathfuncMod, only : polyval
implicit none

  integer :: nx, ny, nxi
  real(r8), allocatable :: x(:)
  real(r8), allocatable :: Y(:,:)
  real(r8), allocatable :: xi(:)
  real(r8), allocatable :: Yi(:,:)

  type(betr_status_type) :: bstatus
  integer :: j1, j2
  real(r8) :: a(3)

  a=(/3._r8,2._r8,1._r8/)

  print*,polyval(a,5._r8),polyval(a,7._r8),polyval(a,9._r8)

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
  call spm_test()
end program main
