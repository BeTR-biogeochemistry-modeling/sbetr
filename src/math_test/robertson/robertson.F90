  program main
  !test the robertson model
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  use robertsonMod, only : get_robertson_rates, get_robertson_matrix
  use MathfuncMod   , only : pd_decomp, flux_correction_fullm
  use BetrStatusType, only : betr_status_type
  use LinearAlgebraMod, only : sparse_gemv, taxpy
  use SparseMatMod, only : sparseMat_type, flux_correction_spm, pd_decomp_spm
  use SparseMatMod, only : spm_axpy, spm_pack, spm_print
  implicit none

  real(r8) :: robmat(3,3)
  real(r8) :: robmatd(3,3)
  real(r8) :: robmatp(3,3)
  real(r8) :: ystates(3)
  integer :: nstep
  type(betr_status_type) :: bstatus
  real(r8), parameter :: epsil = 1.e-8_r8
  real(r8) :: dtime, time
  real(r8) :: rrates(3)
  real(r8) :: dydt(3)
  real(r8) :: ystatef(3)
  class(sparseMat_type), pointer :: spmf
  class(sparseMat_type), pointer :: spmd
  class(sparseMat_type), pointer :: spmp
  integer :: errinfo
  integer :: nstepf
  nstepf=10000
  call bstatus%reset()

  call get_robertson_matrix(robmat)

  call pd_decomp(3, 3, robmat(1:3,1:3), robmatp, robmatd, bstatus)

  !set initial condition
  ystates(1:3)=(/1._r8,epsil,epsil/)
  dtime=0.01
  time=0.
  nstep=0
  do
    !evolve one time step
    call get_robertson_rates(ystates, rrates)

    call flux_correction_fullm(3, 3, robmatp, robmatd,&
      dtime, ystates, rrates, bstatus)

    dydt = 0._r8
    call sparse_gemv('N',3, 3, robmat, &
        3, rrates, 3, dydt)

    call taxpy(3, dtime, dydt, 1, ystates, 1)

    nstep=nstep+1
    if(nstep>=nstepf)exit
  enddo
  ystatef=ystates

  !create the sparse matrix
  call spm_pack(robmat, spmf)

  call pd_decomp_spm(spmf, spmd, spmp)

  ystates(1:3)=(/1._r8,epsil,epsil/)
  dtime=0.01
  time=0.
  nstep=0

  do
    !evolve one time step
    call get_robertson_rates(ystates, rrates)

    call flux_correction_spm(3, 3, spmp, spmd, dtime, ystates, rrates, bstatus)

    dydt=0._r8
    call spm_axpy(3, 3, 1._r8, rrates, spmf, dydt, errinfo)

    call taxpy(3, dtime, dydt, 1, ystates, 1)

    nstep=nstep+1
    print*, ystates
    if(nstep>=nstepf)exit
  enddo
!  print*, ystatef

  end program main
