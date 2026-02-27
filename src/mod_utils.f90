module mod_utils
    use mpi
    use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
    use diezdecomp_core, only: diezdecomp_props_transp, diezdecomp_transp_execute, diezdecomp_stream_default
#if defined(_OPENACC)
    use openacc, only: acc_handle_kind
#else
    use, intrinsic :: iso_fortran_env, only: acc_handle_kind => int64
#endif
    implicit none

    contains

    subroutine read_column_r8(A, n, n_prev, buffer, my_file, icol, irank)
      implicit none
      real(dp)           :: A(0:*), buffer(0:*), temp_r0, temp_r1, temp_r2
      character(len=*) :: my_file
      integer          :: n, n_prev, nb, icol, k, irank, mpi_ierr
      if ((icol<0).or.(icol>3)) error stop 'Unrecognized icol value'
      if (irank==0) then 
        open(19,file=my_file)
            read(19,*) nb
            do k=0,nb-1
                if (icol==0) read(19,*)                            buffer(k)
                if (icol==1) read(19,*) temp_r0,                   buffer(k)
                if (icol==2) read(19,*) temp_r0, temp_r1,          buffer(k)
                if (icol==3) read(19,*) temp_r0, temp_r1, temp_r2, buffer(k)
            end do
        close(19)
      end if
      call MPI_BCAST(nb    , 1  , MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(buffer, nb , MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
      if ((n_prev+n-1).ge.(nb)) error stop 'ERROR!: mismatch in (n_prev, n, nb)'
      A(0:n-1) = buffer(n_prev:(n_prev+n-1))
    end subroutine


    subroutine set_val_slice(A,orig_k0,orig_k1,val_orig)
      implicit none
      real(dp)  :: A(0:,0:,0:), val_orig, val
      integer :: n0,n1,i,j,k,k0,k1,k_min,k_max,orig_k0,orig_k1
      val =  val_orig
      k0  =  orig_k0
      k1  =  orig_k1
      if (k0<0) k0 = k0 + size(A,3)
      if (k1<0) k1 = k1 + size(A,3)
      n0     =  size(A,1)
      n1     =  size(A,2)
      k_min  =  min(k0,k1)
      k_max  =  max(k0,k1)
      !$acc wait
      !$acc parallel loop  collapse(3) default(present)
      do     k = k_min, k_max
        do   j = 0, n1-1
          do i = 0, n0-1
              A(i,j,k) =  val_orig
          end do
        end do
      end do 
      !$acc wait
    end subroutine
    
  function quick_pfix(i0,t,n) result(i)
      integer, intent (in) :: i0,t,n 
      integer              :: i
      i = i0
      if (abs((i0-n)-t) < abs(i-t)) i = i0-n
      if (abs((i0+n)-t) < abs(i-t)) i = i0+n
  end function

  subroutine cumsum_int(A)
    integer :: A(0:), total, i
    total = 0
    do i=0,size(A,1)-1
      total = total + A(i)
      A(i)  = total
    end do 
  end subroutine

  subroutine quick_binary_scalar_reduction(A,b,n,mode_op, utils_red_scalar)
    implicit none
    real(dp)             :: A(0:*), b
    real(dp),contiguous  :: utils_red_scalar(0:)
    integer              :: n,mode_op
    !$acc wait
    call batched_binary_reduction_inner(A,utils_red_scalar,n,1,mode_op)
    !$acc wait
    !$acc update self(utils_red_scalar)
    !$acc wait
    b = utils_red_scalar(0)
  end subroutine

  subroutine batched_binary_reduction_inner(A,B,n,batches,mode_op)
    real(dp)   :: A(0:*), B(0:*)
    integer  :: n,batches,i,j,d,K,mode_op
    d = 1
    !$acc wait
    do while ((d<n).and.((n-d-1)>=0))
        K  = 2*d
        if (mode_op == 0) then 
          !$acc parallel loop collapse(2) private(i,j) default(present)
          do   j=0,batches-1
            do i=0,(n-d-1),K
              A(j*n+i) = A(j*n+i) + A(j*n+i+d)
            end do
          end do
        else if (mode_op == 1) then 
          !$acc parallel loop collapse(2) private(i,j) default(present)
          do   j=0,batches-1
            do i=0,(n-d-1),K
              A(j*n+i) = max(A(j*n+i), A(j*n+i+d))
            end do
          end do
        else 
          if (mode_op.ne.(-1)) error stop 'mode_op.ne.(-1)' 
          !$acc parallel loop collapse(2) private(i,j) default(present)
          do   j=0,batches-1
            do i=0,(n-d-1),K
              A(j*n+i) = min(A(j*n+i), A(j*n+i+d))
            end do
          end do
        end if
        d = 2*d
    end do
    !$acc parallel loop collapse(1) private(j) default(present)
    do j=0,batches-1
      B(j) = A(j*n)
    end do
    !$acc wait
  end subroutine
  
  subroutine diezdecomp_transp_execute_generic_buf_assumed(this, p_in, p_out, work, stream)
    implicit none
    type(diezdecomp_props_transp)      :: this
    real(dp), target                   :: p_in(0:*), p_out(0:*)
    real(dp), target                   :: work(0:)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    stream_internal  =  diezdecomp_stream_default
    if (present(stream)) stream_internal = stream
    call diezdecomp_transp_execute(this, p_in, p_out, stream_internal, work)
    if (.not.present(stream)) then
      !$acc wait(stream_internal)
    end if
  end subroutine

  subroutine batched_binary_reduction(A,B,n,batches,mode_sum,nxz_glob_norm_avg,&
                                      obj_fwd,obj_bwd,n_tr,batches_tr,B_tr,B_tr_1d,work_tr)
    type(diezdecomp_props_transp) :: obj_fwd,obj_bwd
    real(dp)  ::  A(0:*), B(0:*), B_tr(0:*), B_tr_1d(0:*), work_tr(0:), temp_real
    integer   ::  n,batches,n_tr,batches_tr,i,j, nxz_glob_norm_avg, mode_op
    logical   ::  mode_sum
    mode_op = 0
    if (.not.mode_sum) mode_op = 1
    call batched_binary_reduction_inner(A,B,n,batches,mode_op)
    !$acc wait
    call diezdecomp_transp_execute_generic_buf_assumed(obj_fwd,B,B_tr,work_tr)
    !$acc wait
    call batched_binary_reduction_inner(B_tr,B_tr_1d,n_tr,batches_tr,mode_op)
    ! propagate 1d
    !$acc parallel loop collapse(2) private(i,j) default(present)
    do   j=0,batches_tr-1
      do i=0,(n_tr-1)
        B_tr(j*n_tr+i) = B_tr_1d(j)
      end do
    end do 
    !$acc wait
    call diezdecomp_transp_execute_generic_buf_assumed(obj_bwd,B_tr,B,work_tr)
    !$acc wait
    if (nxz_glob_norm_avg>0) then 
      if (.not.mode_sum) error stop '.not.mode_sum'
      temp_real  =  1.d0/nxz_glob_norm_avg
      !$acc parallel loop collapse(1) private(j) default(present)
      do   j=0,batches-1
        B(j) = B(j)*temp_real
      end do
    end if 
  end subroutine

    ! ---------------------- read ibm data ----------------------
  subroutine ibm_read_coeffs(fname_ibm_coeffs, all_ijk, all_c, n_ghosts, n_fluids, &
                             nx_glob, nz_glob, ny_glob, only_size, irank_mpi, nproc_mpi)
    character(len=*)  ::  fname_ibm_coeffs
    integer           ::  all_ijk(0:,0:,0:)
    real(dp)          ::  all_c(0:,0:)
    integer           ::  k,j, nb, mpi_ierr, n_ghosts, n_fluids, nx_glob, nz_glob, ny_glob, &
                          irank_mpi, nproc_mpi
    logical           ::  only_size
    if (irank_mpi==0) then
      open(19,file=fname_ibm_coeffs, action='read', status='old', form='unformatted', access='stream')
        read(19) nb
        if (.not.only_size) then 
          do k=0,nb-1
            read(19)             all_ijk(:,0,k), &
                     all_c(1,k), all_ijk(:,1,k), &
                     all_c(2,k), all_ijk(:,2,k), all_c(0,k)
            do j=0,2
              all_ijk(0,j,k) = modulo(all_ijk(0,j,k),nx_glob)
              all_ijk(1,j,k) = modulo(all_ijk(1,j,k),nz_glob)
              all_ijk(2,j,k) = modulo(all_ijk(2,j,k),ny_glob)
            end do
          end do
        end if
      close(19)
    end if
    call MPI_BCAST(nb          , 1  , MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
    if (.not.only_size) then 
      call MPI_BCAST(all_ijk, 9*nb, MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(all_c  , 3*nb, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
    end if 
    n_ghosts = nb
    n_fluids = 2
  end subroutine

  subroutine ibm_fill_slabs_2d(A, k_first, my_file, A_full, nx_loc, nz_loc, di_loc, dj_loc, nx_glob, nz_glob, irank_mpi, nproc_mpi)
    implicit none
    character(len=*) ::  my_file
    integer          ::  A(0:,0:), A_full(0:,0:), irank_mpi, nproc_mpi, i,j, k_first, &
                         nx_loc, nz_loc, di_loc, dj_loc, nx_glob, nz_glob, mpi_ierr
    if (irank_mpi==0) then 
      open(19,file=my_file, action='read', status='old', form='unformatted', access='stream')
        block 
          integer :: ni, nj
          read(19)   ni, nj
          if (ni.ne.nx_glob) error stop 'ni.ne.nx_glob'
          if (nj.ne.nz_glob) error stop 'nj.ne.nz_glob'
        end block
        do   j=0,nz_glob-1
          do i=0,nx_glob-1
            read(19)  A_full(i,j)
          end do
        end do 
      close(19)
    end if
    call MPI_BCAST(A_full , nx_glob*nz_glob, MPI_INTEGER, 0, mpi_comm_world, mpi_ierr)
    A      = A_full(di_loc:(di_loc+nx_loc-1), dj_loc:(dj_loc+nz_loc-1))
    A      = A - k_first
    A_full = -1
  end subroutine

end module

! def batched_binary_reduction(A,B,n,batches,mode_sum):
!     d = 1
!     while ((d<n) and ((n-d-1)>=0)):
!         K  = 2*d
!         if (mode_sum): 
!             for j in range(batches):
!                 for i in range(0,(n-d),K):
!                     A[j*n+i] = A[j*n+i] + A[j*n+i+d]
!         else:
!             for j in range(batches):
!                 for i in range(0,(n-d),K):
!                     A[j*n+i] = max(A[j*n+i], A[j*n+i+d])
!         d = 2*d
!     for j in range(batches):
!         B[j] = A[j*n]
! def test_binary_reduction():
!     import numpy as np
!     import random 
!     np.random.seed(0)
!     random.seed(0)
!     n_trials = 10000
!     for _ in range(n_trials):
!         n0, n1, n2 = [random.randint(1,5) for _ in range(3)]
!         A    = (5*(np.random.rand(n0,n1,n2)-0.5))
!         ref  = [A[:,:,k].sum() for k in range(n2)]
!         result = [None for _ in range(n2)]
!         batched_binary_reduction([A[i,j,k] for k in range(n2) for j in range(n1) for i in range(n0)], 
!                                  result, n0*n1, n2, True)
!         assert np.fabs(np.array(ref)-np.array(result)).max() < 1e-12
!     print('@test_binary_reduction: all tests passed!')
! test_binary_reduction()
