
program main
    use mpi
    use openacc
    use poisson_solver_multigpu, only: init_spec_poisson_multigpu, run_spec_poisson_multigpu, get_lo_n_bounds_1d
    use diezdecomp_api_generic
    use diezDecomp_api_ibm
    use, intrinsic :: iso_c_binding, only: C_INT, c_intptr_t, C_PTR, C_LOC
    use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
    use, intrinsic :: iso_c_binding
#if ((!defined(_USE_FFTW))&&(!defined(_USE_HIPFFT)))
    use cudafor 
    use cufft
#endif
    implicit none
#if defined(_USE_FFTW)
  include "fftw3.f03"
#endif
    integer                    ::  mpi_ierr,irank_mpi,nproc_mpi, ndev, n_ibm, nhalo, nskip, n_ghosts, n_fluids, mydev
    real(dp)    , allocatable  ::  A_input(:,:,:), A_loc(:,:,:),B_output(:,:,:), buffer_ibm(:), &
                                   A_ibm(:,:,:), A_ibm_loc(:,:,:), A_loc_padded(:,:,:), &
                                   P_band_abc_glob(:,:), P_band_abc_00_glob(:,:), P_ref(:,:,:), temp_arr_n_xzy(:,:,:), &
                                   all_c(:,:)
    integer, allocatable       ::  all_ijk(:,:,:)
    type(diezDecomp_ibm_type)  ::  obj_ibm_A
    real(dp)                   ::  inv_dx, inv_dz, max_error_poi, max_error_ibm
    character(len=26)          ::  fname_info
    character(len=29)          ::  fname_A_input
    character(len=30)          ::  fname_B_output
    character(len=44)          ::  fname_ibm_coeffs_A
    character(len=27)          ::  fname_A_ibm
    character(len=31)          ::  fname_err
#if !defined(_USE_FFTW)
    integer(acc_device_kind) ::dev_type
#endif
    ! -------------------------- Spectral Poisson solver arrays --------------------------
    include "../../src/poisson_solver_mgpu_dtdma/declare_poisson_solver_multigpu.f90"
    ! ----------------------------------------------------

    call mpi_init     (mpi_ierr)
    call MPI_COMM_RANK(mpi_comm_world,irank_mpi,mpi_ierr)
    call MPI_COMM_SIZE(mpi_comm_world,nproc_mpi,mpi_ierr)

#if defined(_USE_HIPFFT)
    dev_type  =  acc_get_device_type()
    ndev      = 8
    mydev     =  mod(irank_mpi,ndev)
    call acc_set_device_num(mydev,dev_type)
    call acc_init(dev_type)
#else
#if !defined(_USE_FFTW)
    dev_type  =  acc_get_device_type()
    mpi_ierr  =  cudaGetDeviceCount(ndev)
    mydev     =  mod(irank_mpi,ndev)
    mpi_ierr = cudaSetDevice(mydev)
    call acc_set_device_num(mydev,dev_type)
    call acc_init(dev_type)
#endif
#endif

    block
      character(len=10)  ::  trial_folder
      call get_command_argument(1, trial_folder)
      write(fname_info        ,'(A7,A10,A9)')  'trials/',trial_folder, '/info.dat'
      write(fname_A_input     ,'(A7,A10,A12)') 'trials/',trial_folder, '/A_input.dat'
      write(fname_B_output    ,'(A7,A10,A13)') 'trials/',trial_folder, '/B_output.dat'
      write(fname_A_ibm       ,'(A7,A10,A10)') 'trials/',trial_folder, '/A_ibm.dat'
      write(fname_ibm_coeffs_A,'(A7,A10,A27)') 'trials/',trial_folder, '/geom_data/ibm_coeffs_A.dat'
      write(fname_err         ,'(A7,A10,A14)') 'trials/',trial_folder, '/max_error.dat'
    end block

    if (irank_mpi == 0) then
      block
        integer :: nproc_check
        real(dp) :: dx, dz
        open(54,file=fname_info)
          read(54,*) nproc_check
          read(54,*) nx_glob, nz_glob, ny_glob, mpi_divs_z, mpi_divs_y, dx, dz, n_ibm, nhalo, nskip
          if (nproc_check.ne.nproc_mpi) error stop '(nproc_check.ne.nproc_mpi)'
          allocate(P_band_abc_glob   (0:ny_glob-1,0:2), &
                   P_band_abc_00_glob(0:ny_glob-1,0:2))
          inv_dx  =  1.d0/dx
          inv_dz  =  1.d0/dz
          block
            integer :: i
            do i=0,ny_glob-1
              read(54,*) P_band_abc_glob(i,:)
            end do
            do i=0,ny_glob-1
              read(54,*) P_band_abc_00_glob(i,:)
            end do
          end block
        close(54)
      end block
    end if
     call MPI_BCAST(nx_glob    ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(nz_glob    ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(ny_glob    ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(mpi_divs_y ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(mpi_divs_z ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(inv_dx     ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(inv_dz     ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(n_ibm      ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(nhalo      ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(nskip      ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)

     if (irank_mpi > 0) then
       allocate(P_band_abc_glob   (0:ny_glob-1,0:2), &
                P_band_abc_00_glob(0:ny_glob-1,0:2))
     end if

     call MPI_BCAST(P_band_abc_glob    ,3*ny_glob,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(P_band_abc_00_glob ,3*ny_glob,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)

    allocate(A_input (0:nx_glob-1,0:nz_glob-1,0:ny_glob-1), &
             B_output(0:nx_glob-1,0:nz_glob-1,0:ny_glob-1), &
             A_ibm   (0:nx_glob-1,0:nz_glob-1,0:ny_glob-1))

    if (irank_mpi == 0) then
      block
        integer :: nx_check, nz_check, ny_check,i,j,k
        open(19,file=fname_A_input, action='read', status='old', form='unformatted', access='stream')
          read(19) nx_check, nz_check, ny_check
          if (nx_glob.ne.nx_check) error stop 'nx_glob.ne.nx_check'
          if (nz_glob.ne.nz_check) error stop 'nz_glob.ne.nz_check'
          if (ny_glob.ne.ny_check) error stop 'ny_glob.ne.ny_check'
          do     k=0,ny_glob-1
            do   j=0,nz_glob-1
              do i=0,nx_glob-1
                read(19) A_input(i,j,k)
              end do
            end do
          end do
        close(19)
      end block

      block
        integer :: nx_check, nz_check, ny_check,i,j,k
        open(19,file=fname_B_output, action='read', status='old', form='unformatted', access='stream')
          read(19) nx_check, nz_check, ny_check
          if (nx_glob.ne.nx_check) error stop 'nx_glob.ne.nx_check'
          if (nz_glob.ne.nz_check) error stop 'nz_glob.ne.nz_check'
          if (ny_glob.ne.ny_check) error stop 'ny_glob.ne.ny_check'
          do     k=0,ny_glob-1
            do   j=0,nz_glob-1
              do i=0,nx_glob-1
                read(19) B_output(i,j,k)
              end do
            end do
          end do
        close(19)
      end block

      block
        integer :: nx_check, nz_check, ny_check,i,j,k
        open(19,file=fname_A_ibm, action='read', status='old', form='unformatted', access='stream')
          read(19) nx_check, nz_check, ny_check
          if (nx_glob.ne.nx_check) error stop 'nx_glob.ne.nx_check'
          if (nz_glob.ne.nz_check) error stop 'nz_glob.ne.nz_check'
          if (ny_glob.ne.ny_check) error stop 'ny_glob.ne.ny_check'
          do     k=0,ny_glob-1
            do   j=0,nz_glob-1
              do i=0,nx_glob-1
                read(19) A_ibm(i,j,k)
              end do
            end do
          end do
        close(19)
      end block

    end if

     call MPI_BCAST(A_input  ,nx_glob*nz_glob*ny_glob, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(B_output ,nx_glob*nz_glob*ny_glob, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
     call MPI_BCAST(A_ibm    ,nx_glob*nz_glob*ny_glob, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)

    include "../../src/poisson_solver_mgpu_dtdma/alloc_poisson_solver_multigpu.f90"
    allocate(P_ref,mold=P)
    allocate(temp_arr_n_xzy,mold=P)
    !$acc enter data create(temp_arr_n_xzy)
    allocate(A_loc,mold=P)
    allocate(A_ibm_loc,mold=P)

    block
      integer ::  i0,j0,k0,i1,j1,k1
      i0           =  0
      j0           =  lo_z_blocks_mpi(mpi_pos_z)
      k0           =  lo_y_blocks_mpi(mpi_pos_y)
      i1           =  nx_glob                       - 1
      j1           =  j0 + nz_blocks_mpi(mpi_pos_z) - 1
      k1           =  k0 + ny_blocks_mpi(mpi_pos_y) - 1
      P_ref        =  B_output(i0:i1, j0:j1, k0:k1)
      A_loc        =  A_input(i0:i1, j0:j1, k0:k1)
      A_ibm_loc    =  A_ibm  (i0:i1, j0:j1, k0:k1)
      P_band_a     =  P_band_abc_glob(k0:k1,0)
      P_band_b     =  P_band_abc_glob(k0:k1,1)
      P_band_c     =  P_band_abc_glob(k0:k1,2)
      P_band_a_00  =  P_band_abc_00_glob(k0:k1,0)
      P_band_b_00  =  P_band_abc_00_glob(k0:k1,1)
      P_band_c_00  =  P_band_abc_00_glob(k0:k1,2)
    end block

    call init_spec_poisson_multigpu(inv_poi_B, poi_A, poi_C, poi_slabs_tips, &
                                    ts_band_a, ts_band_b, ts_cp, a_x, a_z, poi_is_00x, poi_is_00z, &
                                    P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00, P, &
                                    spec_x, spec_z, P_z, &
                                    nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                                    nx_transp, nz_transp_py, nz_loc, ny_loc, inv_dx, inv_dz, &
                                    plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z, &
                                    tr_poi_yz, buffer_transp, lo_xt_blocks_mpi, temp_arr_n_xzy)
    P = A_loc ! fftw3 overwrites "P"
    !$acc update device(P)
    call run_spec_poisson_multigpu(inv_poi_B, poi_A, poi_C, poi_slabs_tips, &
                                   ts_band_a, ts_band_b, ts_cp, a_x, a_z, poi_is_00x, poi_is_00z, &
                                   P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00, P, &
                                   spec_x, spec_z, P_z, trispec, &
                                   nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                                   nx_transp, nz_transp_py, nx_loc, nz_loc, ny_loc, &
                                   plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z, &
                                   tr_poi_xy, tr_poi_yx, tr_poi_yz, tr_poi_zy, buffer_transp)

    !$acc update host(P)
    !$acc wait
    associate(arr => P)
      block
        real(dp) :: temp_real, ref_sum
        integer  :: temp_count, full_count
        temp_real  = sum(arr)
        temp_count = size(arr,1)*size(arr,2)*size(arr,3)
        call MPI_Barrier(mpi_comm_world, mpi_ierr)
        call MPI_Allreduce(temp_real , ref_sum, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        call MPI_Allreduce(temp_count, full_count, 1, MPI_INT, mpi_sum, mpi_comm_world, mpi_ierr)
        call MPI_Barrier(mpi_comm_world, mpi_ierr)
        arr     = arr - ref_sum/(full_count+0.d0)
      end block
    end associate

    associate(arr => P_ref)
      block
        real(dp) :: temp_real, ref_sum
        integer  :: temp_count, full_count
        temp_real  = sum(arr)
        temp_count = size(arr,1)*size(arr,2)*size(arr,3)
        call MPI_Barrier(mpi_comm_world, mpi_ierr)
        call MPI_Allreduce(temp_real , ref_sum, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        call MPI_Allreduce(temp_count, full_count, 1, MPI_INT, mpi_sum, mpi_comm_world, mpi_ierr)
        call MPI_Barrier(mpi_comm_world, mpi_ierr)
        arr     = arr - ref_sum/(full_count+0.d0)
      end block
    end associate

    block
      real(dp) :: temp_real
      temp_real  =  maxval(abs(P - P_ref))
      call MPI_Barrier(mpi_comm_world, mpi_ierr)
      call MPI_Allreduce(temp_real, max_error_poi, 1, MPI_DOUBLE_PRECISION, mpi_max, mpi_comm_world, mpi_ierr)
      call MPI_Barrier(mpi_comm_world, mpi_ierr)
    end block
    if (max_error_poi>1e-10) then
      write(6,*) 'ERROR!! ', max_error_poi;flush(6)
      error stop 'max_error_poi>1e-10'
    end if

    ! ---------------------- read ibm ----------------------
    block
      ! subroutine diezdecomp_ibm_read_ibm_coeffs(this, fname_ibm_coeffs_A, irank_mpi)
      integer  ::  k,j, nb, mpi_ierr
      if (irank_mpi==0) then
        open(19,file=fname_ibm_coeffs_A, action='read', status='old', form='unformatted', access='stream')
            read(19) nb
            allocate(all_ijk(0:2,0:2,0:nb-1))
            allocate(all_c  (    0:2,0:nb-1))
            do k=0,nb-1
              read(19)                  all_ijk(:,0,k), &
                       all_c(1,k), all_ijk(:,1,k), &
                       all_c(2,k), all_ijk(:,2,k), all_c(0,k)
              do j=0,2
                all_ijk(0,j,k) = modulo(all_ijk(0,j,k),nx_glob)
                all_ijk(1,j,k) = modulo(all_ijk(1,j,k),nz_glob)
                all_ijk(2,j,k) = modulo(all_ijk(2,j,k),ny_glob)
              end do
            end do
        close(19)
      end if
      call MPI_BCAST(nb          , 1  , MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      if (irank_mpi>0) then
        allocate(all_ijk(0:2,0:2,0:nb-1))
        allocate(all_c  (    0:2,0:nb-1))
      end if
      call MPI_BCAST(all_ijk, 9*nb, MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(all_c  , 3*nb, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
      ! end subroutine diezdecomp_ibm_read_ibm_coeffs
      n_ghosts = nb
      n_fluids = 2
    end block

    ! ---------- ibm ----------
    block
      integer  ::  A_shape(0:2),nside,i0,j0,k0
      nside    =  nhalo+nskip
      A_shape  =  (/ size(A_loc,1)+2*nside, size(A_loc,2)+2*nside, size(A_loc,3)+2*nside /)
      allocate(A_loc_padded(0:A_shape(0)-1, 0:A_shape(1)-1, 0:A_shape(2)-1))
      block
        integer :: i,j,k,i2,j2,k2
        i0           =  0
        j0           =  lo_z_blocks_mpi(mpi_pos_z)
        k0           =  lo_y_blocks_mpi(mpi_pos_y)
        do     i = nskip,(A_shape(0)-1-nskip)
          do   j = nskip,(A_shape(1)-1-nskip)
            do k = nskip,(A_shape(2)-1-nskip)
              i2                   = modulo(i0 + i-nside, nx_glob)
              j2                   = modulo(j0 + j-nside, nz_glob)
              k2                   = modulo(k0 + k-nside, ny_glob)
              A_loc_padded(i,j,k)  =  A_input(i2,j2,k2)
            end do
          end do
        end do
      end block
      block
        integer :: aux_lo_00(0:0)
        integer :: aux_nx_00(0:0)
        integer :: pos_mpi_ref(0:2)
        aux_lo_00   = 0
        aux_nx_00   = nx_glob
        pos_mpi_ref = (/ 0, mpi_pos_z, mpi_pos_y /)
        call diezdecomp_ibm_init(obj_ibm_A, all_ijk, all_c, n_ghosts, n_fluids, A_shape, &
                                 aux_lo_00, lo_z_blocks_mpi, lo_y_blocks_mpi, &
                                 aux_nx_00,   nz_blocks_mpi,   ny_blocks_mpi, pos_mpi_ref, &
                                 nx_glob, nz_glob, ny_glob, nhalo, nskip, irank_mpi, nproc_mpi)
        allocate(buffer_ibm(0:obj_ibm_A%wsize_buf_mpi-1))
        !$acc enter data copyin(A_loc_padded,buffer_ibm)
        !$acc wait
        call diezdecomp_ibm_exec(obj_ibm_A, A_loc_padded, buffer_ibm)
        !$acc wait
        !$acc exit data copyout(A_loc_padded)
      end block
      block
        real(dp) :: temp_real
        temp_real  =  maxval(abs(A_loc_padded(nside:(A_shape(0)-1-nside),&
                                              nside:(A_shape(1)-1-nside),&
                                              nside:(A_shape(2)-1-nside) ) - A_ibm_loc))
        call MPI_Barrier(mpi_comm_world, mpi_ierr)
        call MPI_Allreduce(temp_real, max_error_ibm, 1, MPI_DOUBLE_PRECISION, mpi_max, mpi_comm_world, mpi_ierr)
        call MPI_Barrier(mpi_comm_world, mpi_ierr)
      end block
    end block

    if (irank_mpi ==0) then
      open(54,file=fname_err)
        write(54,*) 'nx_glob        =  ', nx_glob
        write(54,*) 'nz_glob        =  ', nz_glob
        write(54,*) 'ny_glob        =  ', ny_glob
        write(54,*) 'mpi_divs_z     =  ', mpi_divs_z
        write(54,*) 'mpi_divs_y     =  ', mpi_divs_y
        write(54,*) 'nproc_mpi      =  ', nproc_mpi
        write(54,*) 'max_error_poi  =  ', max_error_poi
        write(54,*) 'n_ibm          =  ', n_ibm
        write(54,*) 'nhalo          =  ', nhalo
        write(54,*) 'nskip          =  ', nskip
        write(54,*) 'max_error_ibm  =  ', max_error_ibm
        flush(54)
      close(54)
        write(6,*) '------------ Results ------------'
        write(6,*) 'fname_err      =  ', fname_err
        write(6,*) 'nx_glob        =  ', nx_glob
        write(6,*) 'nz_glob        =  ', nz_glob
        write(6,*) 'ny_glob        =  ', ny_glob
        write(6,*) 'mpi_divs_z     =  ', mpi_divs_z
        write(6,*) 'mpi_divs_y     =  ', mpi_divs_y
        write(6,*) 'nproc_mpi      =  ', nproc_mpi
        write(6,*) 'max_error_poi  =  ', max_error_poi
        write(6,*) 'n_ibm          =  ', n_ibm
        write(6,*) 'nhalo          =  ', nhalo
        write(6,*) 'nskip          =  ', nskip
        write(6,*) 'max_error_ibm  =  ', max_error_ibm
        flush(6)
    endif

    include "../../src/poisson_solver_mgpu_dtdma/dealloc_poisson_solver_multigpu.f90"

    call MPI_Barrier(mpi_comm_world, mpi_ierr)
    call MPI_finalize(mpi_ierr)
  end program
