module poisson_solver_multigpu
  use  openacc
  use  mpi
  use, intrinsic :: iso_c_binding, only: C_INT, c_intptr_t, C_PTR, C_LOC
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  use mod_utils
#if defined(_USE_HIPFFT)
  use hipfort
  use hipfort_check 
  use hipfort_hipfft
#else
#if defined(_USE_FFTW)
  use, intrinsic :: iso_c_binding
#else
  use  cufft
  use  cudafor
#endif
#endif
  use  diezdecomp_api_generic
  implicit none

  type planinfo_hipfft 
    ! Notes:
    !  1) This is generic data structure. It can be compiled without defined(_USE_HIPFFT).
    !  2) Technically, "batch_max" should be defined as the maximum size of: integer(c_int).
    type(C_PTR)  ::  plan_0, plan_1
    integer      ::  nreps_0, nreps_1, batch_0, batch_1, dj_in, dj_out
    logical      ::  is_fwd 
    integer      ::  batch_max =   2**16
  end type

  contains
  subroutine init_spec_poisson_multigpu(inv_poi_B, poi_A, poi_C, poi_slabs_tips, &
                                        ts_band_a, ts_band_b, ts_cp, a_x, a_z, poi_is_00x, poi_is_00z, &
                                        P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00, P, &
                                        spec_x, spec_z, P_z, &
                                        nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                                        nx_transp, nz_transp_py, nz_loc, ny_loc, inv_dx, inv_dz, &
                                        plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z, utils_red_scalar, &
                                        tr_poi_yz, buffer_transp, lo_xt_blocks_mpi)
    implicit none
    integer     ::  i,j,k, mpi_ierr, &
                    nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                    nx_transp, nz_transp_py, nz_loc, ny_loc, lo_xt_blocks_mpi(0:), poi_is_00x(0:), poi_is_00z(0:)
    real(dp)    ::  const_pi, inv_dx, inv_dz, min_abs_div_bacp, &
                    inv_poi_B(0:,0:,0:), poi_A(0:,0:,0:), poi_C(0:,0:,0:), poi_slabs_tips(0:,0:,0:), &
                    ts_band_a(0:,0:,0:), ts_band_b(0:,0:,0:), ts_cp(0:,0:,0:), a_x(0:), a_z(0:), &
                    P_band_a(0:), P_band_b(0:), P_band_c(0:), P_band_a_00(0:), P_band_b_00(0:), P_band_c_00(0:), P(0:,0:,0:), &
                    buffer_transp(0:), &
                    spec_x(0:,0:,0:), spec_z(0:,0:,0:), P_z(0:,0:,0:), &
                    band_a, band_b, band_c
    real(dp), contiguous :: utils_red_scalar(0:)
    type(diezdecomp_props_transp) :: tr_poi_yz
#if defined(_USE_FFTW)
    type(C_PTR)           ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#else 
#if defined(_USE_HIPFFT)
    type(planinfo_hipfft) ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#else
    integer               ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#endif
#endif
#if defined(_USE_FFTW)
    include "fftw3.f03"
     integer(C_FFTW_R2R_KIND), dimension(1) :: kind_fwd, kind_bwd
     kind_fwd = FFTW_R2HC
     kind_bwd = FFTW_HC2R
#endif
  block 
    integer ::  ndims_x(1), ndims_out_x(1), ndims_z(1), ndims_out_z(1)
    ndims_x      =  (/ nx_glob /)
    ndims_z      =  (/ nz_glob /)
  ! ----------------------- general setup -----------------------

    const_pi             =  4.D0*DATAN(1.D0)

  ! ------------------- BEGIN: fft plans -------------------
#if defined(_USE_HIPFFT)
    call init_hipfft_fwd_bwd(plan_fwd_x, plan_bwd_x, nx_glob, n_fft_x, nz_loc*ny_loc   )
    call init_hipfft_fwd_bwd(plan_fwd_z, plan_bwd_z, nz_glob, n_fft_z, nx_transp*ny_loc)
#else
#if defined(_USE_FFTW)
    ndims_out_x  =  (/ n_fft_x /)
    ndims_out_z  =  (/ n_fft_z /)
    ! call my_fftw_plan_many_dft(plan_fwd_x, ndims_x, nz_loc*ny_loc   , P , nx_glob, spec_x, ndims_out_x, n_fft_x, .true. )
    ! call my_fftw_plan_many_dft(plan_bwd_x, ndims_x, nz_loc*ny_loc   , P , nx_glob, spec_x, ndims_out_x, n_fft_x, .false.)
    ! call my_fftw_plan_many_dft(plan_fwd_z, ndims_z, nx_transp*ny_loc, P , nz_glob, spec_z, ndims_out_z, n_fft_z, .true. )
    ! call my_fftw_plan_many_dft(plan_bwd_z, ndims_z, nx_transp*ny_loc, P , nz_glob, spec_z, ndims_out_z, n_fft_z, .false.)
    plan_fwd_x  =  fftw_plan_many_r2r(1     , ndims_x    , nz_loc*ny_loc , &
                                      P     , ndims_x    , 1             , nx_glob , &
                                      spec_x, ndims_out_x, 1             , n_fft_x , &
                                      kind_fwd, FFTW_MEASURE)
    plan_bwd_x  =  fftw_plan_many_r2r(1     , ndims_x    , nz_loc*ny_loc , &
                                      spec_x, ndims_out_x, 1             , n_fft_x , &
                                      P     , ndims_x    , 1             , nx_glob , &
                                      kind_bwd, FFTW_MEASURE)
    plan_fwd_z  =  fftw_plan_many_r2r(1     , ndims_z    , nx_transp*ny_loc , &
                                      P_z   , ndims_z    , 1                , nz_glob , &
                                      spec_z, ndims_out_z, 1                , n_fft_z , &
                                      kind_fwd, FFTW_MEASURE)
    plan_bwd_z  =  fftw_plan_many_r2r(1     , ndims_z    , nx_transp*ny_loc , &
                                      spec_z, ndims_out_z, 1                , n_fft_z , &
                                      P_z   , ndims_z    , 1                , nz_glob , &
                                      kind_bwd, FFTW_MEASURE)
#else
    ndims_out_x  =  (/ n_fft_x/2 /)
    ndims_out_z  =  (/ n_fft_z/2 /)
    block 
      integer :: ierr
      ierr = ierr + cufftPlanMany(plan_fwd_x , 1, ndims_x, &
                                  ndims_x    , 1, nx_glob, &
                                  ndims_out_x, 1, n_fft_x/2, &
                                  CUFFT_D2Z, nz_loc*ny_loc)
      ierr = ierr + cufftPlanMany(plan_bwd_x , 1, ndims_x, &
                                  ndims_out_x, 1, n_fft_x/2, &
                                  ndims_x    , 1, nx_glob, &
                                  CUFFT_Z2D, nz_loc*ny_loc)
      ierr = ierr + cufftPlanMany(plan_fwd_z , 1, ndims_z, &
                                  ndims_z    , 1, nz_glob, &
                                  ndims_out_z, 1, n_fft_z/2, &
                                  CUFFT_D2Z, nx_transp*ny_loc)
      ierr = ierr + cufftPlanMany(plan_bwd_z , 1, ndims_z, &
                                  ndims_out_z, 1, n_fft_z/2, &
                                  ndims_z    , 1, nz_glob, &
                                  CUFFT_Z2D, nx_transp*ny_loc)
    end block
#endif
#endif
  end block
  ! ------------------- END: fft plans -------------------
    
  ! ------------------- BEGIN: a_x a_z definition -------------------
  block 
    integer :: di,ii,jj
    di = lo_xt_blocks_mpi(mpi_pos_z)
    !$acc parallel loop  private(ii,jj) collapse(2) default(present) 
    do      i=0,nx_transp-1 
      do  j=0,n_fft_z-1
#if defined(_USE_FFTW)
        !  r(0) r(1) r(2) ... r(n/2) i((n+1)/2-1) ... i(2) i(1)
        if ((i+di)<=(nx_glob/2)) then; ii = (i+di) ;else; ii = max(1,(nx_glob+1)/2-1 - ((i+di)-((nx_glob)/2 + 1))) ;end if
        if ( j    <=(nz_glob/2)) then; jj =  j     ;else; jj = max(1,(nz_glob+1)/2-1 - ( j    -((nz_glob)/2 + 1))) ;end if
#else
        ii = (i+di)/2
        jj = j/2
#endif
        a_x(i)  =  -4*(inv_dx**2)*(sin(ii*const_pi/nx_glob)**2)
        a_z(j)  =  -4*(inv_dz**2)*(sin(jj*const_pi/nz_glob)**2)
        poi_is_00x(i) = 0
        poi_is_00z(j) = 0
        if (ii==0) poi_is_00x(i) = 1
        if (jj==0) poi_is_00z(j) = 1
      end do
    end do
  end block
    ! ------------------- END: a_xz definition -------------------
    
    ! NOTE: P-bands are now computed outside (to support different BC's)
    !$acc update device(P_band_a   , P_band_b   , P_band_c   )
    !$acc update device(P_band_a_00, P_band_b_00, P_band_c_00)

    if (mpi_divs_y==1) then 

      !$acc parallel loop  collapse(2) private(band_a,band_b,band_c) default(present)
      do   i = 0, nx_transp-1
        do j = 0, n_fft_z-1
            band_b = P_band_b_00(0)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_b(0)*(1-(poi_is_00x(i)*poi_is_00z(j)))
            band_c = P_band_c_00(0)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_c(0)*(1-(poi_is_00x(i)*poi_is_00z(j)))
            poi_C(j,i,0)   = band_c/(band_b+(a_x(i)+a_z(j)))
            if (ny_loc>=2) then 
              do k = 1, ny_loc-1
                band_a = P_band_a_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_a(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
                band_b = P_band_b_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_b(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
                band_c = P_band_c_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_c(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
                poi_C(j,i,k) = band_c/(band_b+(a_x(i)+a_z(j)) - band_a*poi_C(j,i,k-1))
              end do
            end if
        end do
      end do
  
      ! begin: min_abs_div_bacp
      min_abs_div_bacp = 1e20
#if !defined(_BIN_REDUCTION)
      !$acc parallel loop  collapse(2) default(present) private(band_a,band_b) reduction(min:min_abs_div_bacp)
      do     i = 0, nx_transp-1
          do j = 0, n_fft_z-1 
            band_b             = P_band_b_00(0)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_b(0)*(1-(poi_is_00x(i)*poi_is_00z(j)))
            min_abs_div_bacp   = min(min_abs_div_bacp, abs(band_b + (a_x(i)+a_z(j))))
            if (ny_loc>=2) then 
              do k = 1,ny_loc-1
                band_a = P_band_a_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_a(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
                band_b = P_band_b_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_b(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
                min_abs_div_bacp = min(min_abs_div_bacp, abs(band_b + (a_x(i)+a_z(j)) - band_a * poi_C(j,i,k-1)))
              end do
            end if
          end do
      end do
#else
      call set_val_slice(spec_z,0,-1, 9.9e20)
      !$acc parallel loop  collapse(2) default(present) private(band_a,band_b) 
      do     i = 0, nx_transp-1
          do j = 0, n_fft_z-1  
            band_b             = P_band_b_00(0)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_b(0)*(1-(poi_is_00x(i)*poi_is_00z(j)))
            spec_z(j,i,0)   = ( abs(band_b + (a_x(i)+a_z(j))))
            if (ny_loc>=2) then 
              do k = 1,ny_loc-1
                band_a = P_band_a_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_a(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
                band_b = P_band_b_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_b(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
                spec_z(j,i,k) = ( abs(band_b + (a_x(i)+a_z(j)) - band_a * poi_C(j,i,k-1)))
              end do
            end if
          end do
      end do
      !$acc wait
      call MPI_BARRIER(mpi_comm_world,mpi_ierr)
      call quick_binary_scalar_reduction(spec_z, min_abs_div_bacp, nx_transp*n_fft_z*ny_loc, -1, utils_red_scalar)
      !$acc wait
      call MPI_BARRIER(mpi_comm_world,mpi_ierr)
#endif
      !$acc wait
      ! only irank=0 is active
      write(6,'(A52,E14.6)') 'INFO: init_spec_poisson_multigpu, min_abs_div_bacp: ', min_abs_div_bacp;flush(6) 
    else
      ! ------------------- BEGIN: build poi_A, poi_B, poi_C -------------------
      !     temporarily, poi_B is stored inside inv_poi_B
      !$acc parallel loop  collapse(3) private(band_a,band_b,band_c) default(present) 
      do      k = 0, ny_loc-1
        do    i = 0, nx_transp-1 
          do  j = 0, n_fft_z-1 
            band_a = P_band_a_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_a(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
            band_b = P_band_b_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_b(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
            band_c = P_band_c_00(k)*(poi_is_00x(i)*poi_is_00z(j)) + P_band_c(k)*(1-(poi_is_00x(i)*poi_is_00z(j)))
            poi_A    (j,i,k) = band_a
            inv_poi_B(j,i,k) = band_b + (a_x(i)+a_z(j))
            poi_C    (j,i,k) = band_c
          end do
        end do
      end do
      if (ny_loc>=3) then
        !$acc parallel loop  collapse(2) default(present) 
        do         i = 0, nx_transp-1 
            do     j = 0, n_fft_z-1
                do k = 2, ny_loc-1,1
                    inv_poi_B(j,i,k) = inv_poi_B(j,i,k) + (-poi_A(j,i,k) / inv_poi_B(j,i,k-1)) * poi_C(j,i,k-1)
                    poi_A    (j,i,k) =                    (-poi_A(j,i,k) / inv_poi_B(j,i,k-1)) * poi_A(j,i,k-1)
                end do
            end do
        end do
        !$acc parallel loop  collapse(2) default(present) 
        do         i = 0, nx_transp-1 
            do     j = 0, n_fft_z-1
                do k = ny_loc-3, 0, -1
                    if (k==0) then
                        inv_poi_B(j,i,k) =  inv_poi_B(j,i,k) + (-poi_C(j,i,k) / inv_poi_B(j,i,k+1)) * poi_A(j,i,k+1)
                    else 
                        poi_A(j,i,k)     =  poi_A(j,i,k)     + (-poi_C(j,i,k) / inv_poi_B(j,i,k+1)) * poi_A(j,i,k+1)
                    end if
                    poi_C(j,i,k)         =                     (-poi_C(j,i,k) / inv_poi_B(j,i,k+1)) * poi_C(j,i,k+1)
                end do
            end do
        end do
      end if
      
      ! ------------------- END: build poi_A, poi_B, poi_C -------------------
      
      ! -------------- ts_band_a --------------
      !$acc parallel loop  collapse(2) default(present) 
       do     i = 0, nx_transp-1 
           do j = 0, n_fft_z-1
               poi_slabs_tips(j,i,0)  =  poi_A(j,i,0       )
               poi_slabs_tips(j,i,1)  =  poi_A(j,i,ny_loc-1)
           end do
       end do

      !$acc wait
      call diezdecomp_transp_execute_generic_buf(tr_poi_yz, poi_slabs_tips, ts_band_a, buffer_transp)
      !$acc wait

      ! -------------- ts_band_b --------------
      !$acc parallel loop  collapse(2) default(present) 
      do     i = 0, nx_transp-1 
          do j = 0, n_fft_z-1
              poi_slabs_tips(j,i,0)  =  inv_poi_B(j,i,0       )
              poi_slabs_tips(j,i,1)  =  inv_poi_B(j,i,ny_loc-1)
          end do
      end do

      !$acc wait
      call diezdecomp_transp_execute_generic_buf(tr_poi_yz, poi_slabs_tips, ts_band_b, buffer_transp)
      !$acc wait
  
      ! -------------- ts_band_c (temporarily inside ts_cp) --------------
      !$acc parallel loop  collapse(2) default(present) 
      do    i =0, nx_transp-1 
         do j =0, n_fft_z-1
           poi_slabs_tips(j,i,0)  =  poi_C(j,i,0       )
           poi_slabs_tips(j,i,1)  =  poi_C(j,i,ny_loc-1)
         end do
      end do
      
      !$acc wait
      call diezdecomp_transp_execute_generic_buf(tr_poi_yz, poi_slabs_tips, ts_cp, buffer_transp)
      !$acc wait
  
      ! -------------- ts_cp (finish) --------------
      !$acc parallel loop  collapse(2) default(present) 
      do   i = 0, nx_transp-1
        do j = 0, nz_transp_py-1
          ts_cp(j,i,0)     = ts_cp(j,i,0)/(ts_band_b(j,i,0))
          do k = 1, (2*mpi_divs_y)-1
              ts_cp(j,i,k) = ts_cp(j,i,k)/(ts_band_b(j,i,k) - ts_band_a(j,i,k)*ts_cp(j,i,k-1))
          end do
        end do
      end do
      
      ! -------------- inv_poi_B (finish) --------------
      !$acc parallel loop  collapse(3) default(present) 
      do         j = 0, n_fft_z-1
          do     i = 0, nx_transp-1
              do k = 0, ny_loc-1
                  inv_poi_B(j,i,k) =  1.d0 / inv_poi_B(j,i,k)
              end do
          end do
      end do
    end if
  end subroutine

  subroutine run_spec_poisson_multigpu(inv_poi_B, poi_A, poi_C, poi_slabs_tips, &
                                       ts_band_a, ts_band_b, ts_cp, a_x, a_z, poi_is_00x, poi_is_00z, &
                                       P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00, P, &
                                       spec_x, spec_z, P_z, trispec, &
                                       nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                                       nx_transp, nz_transp_py, nx_loc, nz_loc, ny_loc, &
                                       plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z, &
                                       tr_poi_xy, tr_poi_yx, tr_poi_yz, tr_poi_zy, buffer_transp)
    implicit none
    integer     ::  i,j,k, mpi_ierr, &
                    nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                    nx_transp, nz_transp_py, nx_loc, nz_loc, ny_loc, poi_is_00x(0:), poi_is_00z(0:), poi_is_00_ji
    real(dp)    ::  inv_poi_B(0:,0:,0:), poi_A(0:,0:,0:), poi_C(0:,0:,0:), &
                    ts_band_a(0:,0:,0:), ts_band_b(0:,0:,0:), ts_cp(0:,0:,0:), a_x(0:), a_z(0:), &
                    P_band_a(0:), P_band_b(0:), P_band_c(0:), P_band_a_00(0:), P_band_b_00(0:), P_band_c_00(0:), &
                    P(0:,0:,0:), buffer_transp(0:), &
                    spec_x(0:,0:,0:), spec_z(0:,0:,0:), P_z(0:,0:,0:), poi_slabs_tips(0:,0:,0:), trispec(0:,0:,0:), &
                    band_a, band_b, band_c, a_xz_ji
    type(diezdecomp_props_transp) :: tr_poi_xy, tr_poi_yx, tr_poi_yz, tr_poi_zy
#if defined(_USE_FFTW)
    type(C_PTR)           ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#else 
#if defined(_USE_HIPFFT)
    type(planinfo_hipfft) ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#else
    integer               ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#endif
#endif

#if defined(_USE_FFTW)
    include "fftw3.f03"
#endif
  !$acc wait
#if defined(_USE_HIPFFT)
    call quick_hipfft_exec(plan_fwd_x, P, spec_x)
#else 
#if defined(_USE_FFTW)
    call fftw_execute_r2r(plan_fwd_x, P, spec_x)
#else
  block 
    integer :: ierr
    !$acc host_data use_device(P,spec_x)
    ierr = ierr + cufftExecD2Z(plan_fwd_x, P, spec_x)
    !$acc end host_data
    ierr = ierr + cudaDeviceSynchronize()
end block
#endif
#endif
  !$acc wait
  call diezdecomp_transp_execute_generic_buf(tr_poi_xy, spec_x, P_z, buffer_transp)
  !$acc wait
#if defined(_USE_HIPFFT)
    call quick_hipfft_exec(plan_fwd_z, P_z, spec_z)
#else 
#if defined(_USE_FFTW)
    call fftw_execute_r2r(plan_fwd_z, P_z, spec_z)
#else
  block 
    integer :: ierr
    !$acc host_data use_device(P_z,spec_z)
    ierr = ierr + cufftExecD2Z(plan_fwd_z, P_z, spec_z)
    !$acc end host_data
    ierr = ierr + cudaDeviceSynchronize()
end block
#endif
#endif
  !$acc wait
    if (mpi_divs_y==1) then 
      !$acc parallel loop  collapse(2) private(band_a,band_b,a_xz_ji,poi_is_00_ji) default(present) 
      do   i = 0, nx_transp-1
        do j = 0, n_fft_z-1 
          a_xz_ji       =  (a_x(i)+a_z(j))
          poi_is_00_ji  =  (poi_is_00x(i)*poi_is_00z(j))
          band_b        =  P_band_b_00(0)*poi_is_00_ji + P_band_b(0)*(1-poi_is_00_ji)
          spec_z(j,i,0) = spec_z(j,i,0)/(band_b+a_xz_ji)
        end do
      end do
      if (ny_loc>=2) then
        !$acc parallel loop  collapse(2) private(band_a,band_b,a_xz_ji,poi_is_00_ji) default(present) 
        do   i = 0, nx_transp-1
          do j = 0, n_fft_z-1 
            a_xz_ji       =  (a_x(i)+a_z(j))
            poi_is_00_ji  =  (poi_is_00x(i)*poi_is_00z(j))
            do k = 1,ny_loc-1,1
              band_a = P_band_a_00(k)*poi_is_00_ji + P_band_a(k)*(1-poi_is_00_ji)
              band_b = P_band_b_00(k)*poi_is_00_ji + P_band_b(k)*(1-poi_is_00_ji)
              spec_z(j,i,k) = (spec_z(j,i,k) - band_a*spec_z(j,i,k-1))/(band_b+a_xz_ji - &
                                                                        band_a*poi_C(j,i,k-1))
            end do
            do k = ny_loc-2,0,-1
              spec_z(j,i,k)  =  spec_z(j,i,k) - spec_z(j,i,k+1) * poi_C(j,i,k)
            end do
          end do
        end do
      end if
      
    else
      ! ----------------------- BEGIN: pre-processing spec_z -----------------------
      if (ny_loc>=3) then
        !$acc parallel loop  collapse(2) private(i,j,k,band_a,band_c,poi_is_00_ji) default(present) 
        do         j = 0, n_fft_z-1
            do     i = 0, nx_transp-1
                poi_is_00_ji  =  (poi_is_00x(i)*poi_is_00z(j))
                do k=2,ny_loc-1,1
                    band_a         = P_band_a_00(k)*poi_is_00_ji + P_band_a(k)*(1-poi_is_00_ji)
                    spec_z(j,i,k)  =  spec_z(j,i,k)  -  band_a * inv_poi_B(j,i,k-1) * spec_z(j,i,k-1)
                end do
                do k=ny_loc-3,0,-1
                    band_c         = P_band_c_00(k)*poi_is_00_ji + P_band_c(k)*(1-poi_is_00_ji)
                    spec_z(j,i,k)  =  spec_z(j,i,k)  -  band_c * inv_poi_B(j,i,k+1) * spec_z(j,i,k+1)
                end do
            end do
        end do
      end if

      ! ----------------------- END: pre-processing spec -----------------------

      ! ----------------------- BEGIN: solve spec_tips -----------------------
      !$acc parallel loop  collapse(2) default(present) 
      do     i =0,nx_transp-1  
          do j =0,n_fft_z-1
              poi_slabs_tips(j,i,0)  =  spec_z(j,i,0   )
              poi_slabs_tips(j,i,1)  =  spec_z(j,i,ny_loc-1)
          end do
      end do
      !$acc wait
      call diezdecomp_transp_execute_generic_buf(tr_poi_yz, poi_slabs_tips, trispec, buffer_transp)
      !$acc wait
  
      ! ----------- BEGIN: tridiagonal solver for trispec -----------
      !$acc parallel loop  collapse(2) default(present) 
      do   i = 0, nx_transp-1 
        do j = 0, nz_transp_py-1 
          trispec(j,i,0) = trispec(j,i,0)/ts_band_b(j,i,0)
          do k = 1,(2*mpi_divs_y)-1
              trispec(j,i,k)  =  (trispec(j,i,k) - ts_band_a(j,i,k)*trispec(j,i,k-1))/(ts_band_b(j,i,k)               - &
                                                                                       ts_band_a(j,i,k)*ts_cp(j,i,k-1))
          end do
          do k = (2*mpi_divs_y)-2,0,-1
              trispec(j,i,k)  =  trispec(j,i,k) - trispec(j,i,k+1) * ts_cp(j,i,k)
          end do
        end do
      end do
      ! ----------- END: tridiagonal solver for trispec -----------
      !$acc wait
      call diezdecomp_transp_execute_generic_buf(tr_poi_zy, trispec, poi_slabs_tips, buffer_transp)
      !$acc wait
  
      !$acc parallel loop  collapse(2) default(present) 
      do     i =0,nx_transp-1 
          do j =0,n_fft_z-1
              spec_z(j,i,0       ) = poi_slabs_tips(j,i,0)
              spec_z(j,i,ny_loc-1) = poi_slabs_tips(j,i,1)
          end do
      end do
      ! ----------------------- END: solve spec_tips -----------------------

      ! ----------------------- BEGIN: rebuild spec -----------------------
      if (ny_loc>=3) then 
        !$acc parallel loop  collapse(3) default(present) 
        do          k = 1, ny_loc-2,1
            do      i = 0, nx_transp-1
                do  j = 0, n_fft_z-1
                    spec_z(j,i,k) = (spec_z(j,i,k)                     - & 
                                     poi_C(j,i,k)*spec_z(j,i,ny_loc-1) - &
                                     poi_A(j,i,k)*spec_z(j,i,0        ) )*inv_poi_B(j,i,k)
                end do
            end do
        end do
      end if
    end if
    ! ----------------------- END: rebuild spec -----------------------
    
    ! ----------------------- BEGIN: irrft2 -----------------------
  !$acc wait
#if defined(_USE_HIPFFT)
    call quick_hipfft_exec(plan_bwd_z, spec_z, P_z)
#else 
#if defined(_USE_FFTW)
    call fftw_execute_r2r(plan_bwd_z, spec_z, P_z)
#else
  block 
    integer :: ierr
    !$acc host_data use_device(P_z,spec_z)
    ierr = ierr + cufftExecZ2D(plan_bwd_z, spec_z, P_z)
    !$acc end host_data
    ierr = ierr + cudaDeviceSynchronize()
end block
#endif
#endif

  !$acc wait
  call diezdecomp_transp_execute_generic_buf(tr_poi_yx, P_z, spec_x, buffer_transp)
  !$acc wait

#if defined(_USE_HIPFFT)
    call quick_hipfft_exec(plan_bwd_x, spec_x, P)
#else 
#if defined(_USE_FFTW)
    call fftw_execute_r2r(plan_bwd_x, spec_x, P)
#else
  block 
    integer :: ierr
    !$acc host_data use_device(P,spec_x)
    ierr = ierr + cufftExecZ2D(plan_bwd_x, spec_x, P)
    !$acc end host_data
    ierr = ierr + cudaDeviceSynchronize()
end block
#endif
#endif
    ! ----------------------- END: irrft2 -----------------------
  !$acc wait
   block 
     real(dp) :: inv_nxnz
     inv_nxnz             =  1.d0/(nx_glob*nz_glob + 0.d0)
     !$acc parallel loop  collapse(3) default(present) 
     do     k = 0, ny_loc-1
       do   j = 0, nz_loc-1
         do i = 0, nx_loc-1 
           P(i,j,k) = P(i,j,k)*inv_nxnz
         end do
       end do
     end do
   end block

    ! ----------------------- END: irrft2 -----------------------
  end subroutine

  subroutine get_lo_n_bounds_1d(mpi_divs, n_glob, lo_blocks_mpi, n_blocks_mpi) 
    integer :: mpi_divs, n_glob, lo_blocks_mpi(0:), n_blocks_mpi(0:)
    block 
      integer :: i,n_base
      n_base        =  n_glob/mpi_divs
      n_blocks_mpi  =  n_base
      do i=0,(n_glob - n_base*mpi_divs)-1
        n_blocks_mpi(mod(i,mpi_divs)) = n_blocks_mpi(mod(i,mpi_divs)) + 1
      end do
    end block 
    if (sum(n_blocks_mpi).ne.n_glob) error stop 'sum(n_blocks_mpi).ne.n_glob'
    block  
      integer :: i,total
      total = 0
      do i=0,mpi_divs-1
        lo_blocks_mpi(i) = total
        total            = total + n_blocks_mpi(i)
      end do 
    end block
  end subroutine

  subroutine init_hipfft_fwd_bwd(plan_fwd, plan_bwd, n_glob, n_fft, batch)
    type(planinfo_hipfft) ::  plan_fwd, plan_bwd
    integer               ::  n_glob, n_fft, batch, ndims(1), ndims_out(1), batch_max, &
                              nreps_0, nreps_1, batch_0, batch_1
    ndims      =  (/ n_glob /)
    ndims_out  =  (/ n_fft/2 /)
    batch_max  =  plan_fwd%batch_max

    if (plan_fwd%batch_max.ne.plan_bwd%batch_max) error stop 'plan_fwd%batch_max.ne.plan_bwd%batch_max'
    if (batch_max<=0) error stop 'batch_max<=0'
    if (batch>batch_max) then 
      nreps_0 = batch / batch_max ; nreps_1 = 1
      batch_0 = batch_max         ; batch_1 = batch - nreps_0*batch_max
      if (batch_1<=0) error stop 'batch_1<=0'
    else
      nreps_0 = 1                 ; nreps_1 = 0
      batch_0 = batch             ; batch_1 = 0
    end if 

    plan_fwd%nreps_0  = nreps_0  ;  plan_bwd%nreps_0  = nreps_0
    plan_fwd%nreps_1  = nreps_1  ;  plan_bwd%nreps_1  = nreps_1
    plan_fwd%dj_in    = n_glob   ;  plan_bwd%dj_in    = n_fft  
    plan_fwd%dj_out   = n_fft    ;  plan_bwd%dj_out   = n_glob 
    plan_fwd%batch_0  = batch_0  ;  plan_bwd%batch_0  = batch_0
    plan_fwd%batch_1  = batch_1  ;  plan_bwd%batch_1  = batch_1
    plan_fwd%is_fwd   = .true.   ;  plan_bwd%is_fwd   = .false.
#if defined(_USE_HIPFFT)
    call hipfftCheck(hipfftPlanMany_(plan_fwd%plan_0 , 1 , c_loc(ndims),&
                                     c_loc(ndims)    , 1 , n_glob      ,&
                                     c_loc(ndims_out), 1 , n_fft/2     ,&
                                     HIPFFT_D2Z          , batch_0     ))
    call hipfftCheck(hipfftPlanMany_(plan_bwd%plan_0 , 1 , c_loc(ndims),&
                                     c_loc(ndims_out), 1 , n_fft/2     ,&
                                     c_loc(ndims)    , 1 , n_glob      ,&
                                     HIPFFT_Z2D          , batch_0     ))
    if (nreps_1>0) then 
      call hipfftCheck(hipfftPlanMany_(plan_fwd%plan_1 , 1 , c_loc(ndims),&
                                       c_loc(ndims)    , 1 , n_glob      ,&
                                       c_loc(ndims_out), 1 , n_fft/2     ,&
                                       HIPFFT_D2Z          , batch_1     ))
      call hipfftCheck(hipfftPlanMany_(plan_bwd%plan_1 , 1 , c_loc(ndims),&
                                       c_loc(ndims_out), 1 , n_fft/2     ,&
                                       c_loc(ndims)    , 1 , n_glob      ,&
                                       HIPFFT_Z2D          , batch_1     ))
    end if 
#else
        error stop '!defined(_USE_HIPFFT)'
#endif
  end subroutine

  subroutine  quick_hipfft_exec(plan, A, B)
    type(planinfo_hipfft) :: plan
    real(dp)              :: A(0:*), B(0:*)
    integer               :: i,j_in,j_out
    !$acc wait
    !$acc host_data use_device(A,B)
    if (plan%nreps_0>0) then 
      do i=0,(plan%nreps_0 - 1)
        j_in  = i * plan%batch_0 * plan%dj_in
        j_out = i * plan%batch_0 * plan%dj_out
#if defined(_USE_HIPFFT)
        if (plan%is_fwd) then 
          call hipfftCheck(hipfftExecD2Z(plan%plan_0,c_loc(A(j_in)),c_loc(B(j_out))))
        else
          call hipfftCheck(hipfftExecZ2D(plan%plan_0,c_loc(A(j_in)),c_loc(B(j_out))))
        end if 
#else
        error stop '!defined(_USE_HIPFFT)'
#endif
      end do
    end if 
    if (plan%nreps_1>0) then 
        if (plan%nreps_1.ne.1) error stop 'plan%nreps_1.ne.1'
        j_in  = plan%nreps_0 * plan%batch_0 * plan%dj_in
        j_out = plan%nreps_0 * plan%batch_0 * plan%dj_out
#if defined(_USE_HIPFFT)
        if (plan%is_fwd) then 
          call hipfftCheck(hipfftExecD2Z(plan%plan_1,c_loc(A(j_in)),c_loc(B(j_out))))
        else
          call hipfftCheck(hipfftExecZ2D(plan%plan_1,c_loc(A(j_in)),c_loc(B(j_out))))
        end if 
#else
        error stop '!defined(_USE_HIPFFT)'
#endif
    end if
    !$acc end host_data
#if defined(_USE_HIPFFT)
    call hipfftCheck(hipDeviceSynchronize())
#else
        error stop '!defined(_USE_HIPFFT)'
#endif
    !$acc wait
  end subroutine
end module poisson_solver_multigpu


