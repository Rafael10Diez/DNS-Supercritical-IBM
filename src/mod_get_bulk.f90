module mod_get_bulk
  use  openacc
  use  mpi
  use  mod_utils
  use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr, c_loc
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none
  contains
  subroutine bulk_compute(U, V, W, T, P_now, P, rho_arr, vol_U_bulk, vol_T_bulk, vol_div_avg, nx_glob, nz_loc, ny_loc, &
                          dy_cell, temp_arr_n_xzy, use_rough_surf, &
                          bulk_U_bot_inds, bulk_U_top_inds, bulk_T_bot_inds, bulk_T_top_inds, bulk_use_rho_weighting, &
                          inv_dx, inv_dz, inv_dy_p, &
                          bulk_last_iter, iters_full, bulk_use_target_Ub, force_Sf_x_pressure, use_incompressible, &
                          dtmax, inv_dtmax, bulk_nprint, snap_iter_start, snap_freq, nstep, &
                          min_delta_U_bulk, max_delta_U_bulk, min_Sf_x_mid, max_Sf_x_mid, prof_1d_avg_N, &
                          !__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]
                          mpi_divs_y, mpi_pos_z, mpi_pos_y, obj_ranks_0zy, prof_1d_temp, utils_red_scalar, &
                          bulk_time_start_print, bulk_time_finish_print, irank_mpi)
    use  ieee_arithmetic
    use diezdecomp_api_generic, only: diezdecomp_parsed_mpi_ranks
    real(dp) ::  U(-1:,-1:,-1:), V(-1:,-1:,-1:), W(-1:,-1:,-1:), T(-1:,-1:,-1:), P_now(-1:,-1:,-1:), P(0:,0:,0:), &
                 rho_arr(-1:,-1:,-1:), &
                 vol_U_bulk, vol_T_bulk, vol_div_avg, dy_cell(0:), temp_arr_n_xzy(0:,0:,0:), &
                 inv_dx, inv_dz, inv_dy_p(0:), prof_1d_temp(0:), &
                 dtmax, inv_dtmax, &
                 min_delta_U_bulk, max_delta_U_bulk, min_Sf_x_mid, max_Sf_x_mid, &
                 bulk_time_start_print, bulk_time_finish_print
    real(dp), contiguous :: utils_red_scalar(0:)
    !__[PYTHON_HOLDER_PROF1D_COMPBULK_FUNC_DECL_VARS]
    type(diezdecomp_parsed_mpi_ranks)  ::  obj_ranks_0zy
    integer  ::  nx_glob, nz_loc, ny_loc, &
                 bulk_U_bot_inds(0:,0:), bulk_U_top_inds(0:,0:), bulk_T_bot_inds(0:,0:), bulk_T_top_inds(0:,0:), &
                 bulk_last_iter, iters_full, &
                 mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                 bulk_nprint, snap_iter_start, snap_freq, nstep, prof_1d_avg_N, irank_mpi
    logical  ::  use_rough_surf, bulk_use_rho_weighting, bulk_use_target_Ub, force_Sf_x_pressure, use_incompressible
    !$acc wait
    bulk_time_finish_print = MPI_WTIME()
    block
      real(dp)  :: U_bulk, T_bulk, P_maxabs, div_avg, div_maxabs, temp_real
      
      call get_U_bulk(U, U_bulk, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, dy_cell, utils_red_scalar, &
                      temp_arr_n_xzy, use_rough_surf, bulk_U_bot_inds, bulk_U_top_inds, bulk_use_rho_weighting)
      call get_T_bulk(T, T_bulk, vol_T_bulk, nx_glob, nz_loc, ny_loc, rho_arr, dy_cell, utils_red_scalar, &
                      temp_arr_n_xzy, use_rough_surf, bulk_T_bot_inds, bulk_T_top_inds, bulk_use_rho_weighting)
    
      block 
        integer :: i, j, k, mpi_ierr
        
        ! ---------------- check if (vol_div_avg) was defined ----------------
        if (vol_div_avg<0) then
          temp_real  =  (nx_glob+0.d0)*(nz_loc+0.d0)*(ny_loc+0.d0)
          call MPI_Allreduce(temp_real, vol_div_avg, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        end if
    
        ! ---------------- process div_avg ----------------
        div_avg    = 0
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(3) default(present) reduction(+:div_avg)
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                div_avg  =  div_avg + (U(i+1, j  , k  ) - U(i,j,k))*inv_dx + &
                                      (W(i  , j+1, k  ) - W(i,j,k))*inv_dz + &
                                      (V(i  , j  , k+1) - V(i,j,k))*inv_dy_p(k)
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
        !$acc parallel loop collapse(3) default(present) 
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                temp_arr_n_xzy(i,j,k)  =  (U(i+1, j  , k  ) - U(i,j,k))*inv_dx + &
                                          (W(i  , j+1, k  ) - W(i,j,k))*inv_dz + &
                                          (V(i  , j  , k+1) - V(i,j,k))*inv_dy_p(k)
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, div_avg, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        !$acc wait
        temp_real  =  div_avg/vol_div_avg
        call MPI_Allreduce(temp_real, div_avg, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
    
        ! ---------------- process div_maxabs ----------------
        div_maxabs = -1e20
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(3) default(present) reduction(max:div_maxabs)
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                div_maxabs  =  max(div_maxabs, abs((U(i+1, j  , k  ) - U(i,j,k))*inv_dx + &
                                                   (W(i  , j+1, k  ) - W(i,j,k))*inv_dz + &
                                                   (V(i  , j  , k+1) - V(i,j,k))*inv_dy_p(k)))
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, -9.9e20)
        !$acc parallel loop collapse(3) default(present) 
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                temp_arr_n_xzy(i,j,k)  =  ( abs((U(i+1, j  , k  ) - U(i,j,k))*inv_dx + &
                                                (W(i  , j+1, k  ) - W(i,j,k))*inv_dz + &
                                                (V(i  , j  , k+1) - V(i,j,k))*inv_dy_p(k)))
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, div_maxabs, nx_glob*nz_loc*ny_loc, 1, utils_red_scalar)
#endif
        !$acc wait
        temp_real = div_maxabs
        call MPI_Allreduce(temp_real, div_maxabs, 1, MPI_DOUBLE_PRECISION, mpi_max, mpi_comm_world, mpi_ierr)
    
        ! ---------------- process P_maxabs ----------------
        P_maxabs = -1e20
 if (use_incompressible) then 
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(3) default(present) reduction(max:P_maxabs)
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                P_maxabs  =  max(P_maxabs, abs(P(i,j,k)))
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, -9.9e20)
        !$acc parallel loop collapse(3) default(present) 
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                temp_arr_n_xzy(i,j,k)  =  abs(P(i,j,k))
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, P_maxabs, nx_glob*nz_loc*ny_loc, 1, utils_red_scalar)
#endif
 else
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(3) default(present) reduction(max:P_maxabs)
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                P_maxabs  =  max(P_maxabs, abs(P_now(i,j,k)))
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, -9.9e20)
        !$acc parallel loop collapse(3) default(present) 
        do      k=0, ny_loc -1
          do    j=0, nz_loc -1
            do  i=0, nx_glob-1
                temp_arr_n_xzy(i,j,k)  =  abs(P_now(i,j,k))
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, P_maxabs, nx_glob*nz_loc*ny_loc, 1, utils_red_scalar)
#endif
 end if
        !$acc wait
        temp_real  =  P_maxabs
        call MPI_Allreduce(temp_real, P_maxabs, 1, MPI_DOUBLE_PRECISION, mpi_max, mpi_comm_world, mpi_ierr)
      end block
    
      if (irank_mpi == 0) then
        block
          real(dp) :: avg_time_iter
          avg_time_iter  = (bulk_time_finish_print-bulk_time_start_print)/(iters_full-bulk_last_iter + 0.d0)
          write(6,"(A37,I9,6E14.5)") 'i t/i Ub Tb |Pmax| div_avg |div_max|:', iters_full, avg_time_iter, &
                                      U_bulk, T_bulk, P_maxabs, div_avg, div_maxabs
          flush(6)
        end block
      end if
    
      if (ieee_is_nan(U_bulk ).or.ieee_is_nan(T_bulk    ).or.ieee_is_nan(P_maxabs).or. &
          ieee_is_nan(div_avg).or.ieee_is_nan(div_maxabs)) then
        write(6,*) 'Error: NaN detected. Solver will stop.'; flush(6)
        open(19,file='NaN_detected.txt')
          write(19,*) 'Error: NaN detected. Solver will stop.'
        close(19)
        error stop 'Error: NaN detected. Solver will stop.'
      end if
    
    end block
    
    bulk_last_iter         =  iters_full
    !$acc wait
    
    block
      integer :: signal_runtime, mpi_ierr
      signal_runtime = 0
      ! ---------------------------------------- Overview File "signal_runtime.txt" ----------------------------------------
      ! 0                         ! signal_runtime (see below)
      ! dtmax                     ! new [real(dp):: dtmax                     ] (for signal_runtime == 1)
      ! bulk_nprint               ! new [integer :: bulk_nprint               ] (for signal_runtime == 1)
      ! snap_iter_start snap_freq ! new [integer :: snap_iter_start, snap_freq] (for signal_runtime == 1)
      !
      ! ! signal_runtime:
      ! !     0    :  continue
      ! !     1    :  change parameters (dtmax, bulk_nprint, snap_iter_start, snap_freq)
      ! !     other:  stop
      if (irank_mpi == 0) then
        open(19,file="signal_runtime.txt")
          read(19,*)   signal_runtime
          if (signal_runtime==1) then
            read(19,*) dtmax
            read(19,*) bulk_nprint
            read(19,*) snap_iter_start, snap_freq
            if ((snap_iter_start<0).or.(snap_freq<0)) then
              write(6,'(A68,2I10)') 'INFO: disabling intermediate snapshots [snap_iter_start,snap_freq]: ',&
              snap_iter_start, snap_freq;flush(6)
              snap_iter_start = nstep+100
              snap_freq       = nstep+100
            end if
          end if
        close(19)
      end if
      call MPI_BCAST(signal_runtime , 1, MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(dtmax          , 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(bulk_nprint    , 1, MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(snap_iter_start, 1, MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(snap_freq      , 1, MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
      inv_dtmax  =  1.d0/dtmax
      if ((signal_runtime<0).or.(signal_runtime>1))   then ! perform action
        error stop 'signal_runtime<0 or signal_runtime>1 called'
      end if
    end block
    
    if (bulk_use_target_Ub) then 
      block
        integer  :: mpi_ierr
        real(dp) :: temp_real
        temp_real = min_delta_U_bulk
        call MPI_Allreduce(temp_real, min_delta_U_bulk, 1, MPI_DOUBLE_PRECISION, mpi_min, mpi_comm_world, mpi_ierr)
        temp_real = max_delta_U_bulk
        call MPI_Allreduce(temp_real, max_delta_U_bulk, 1, MPI_DOUBLE_PRECISION, mpi_max, mpi_comm_world, mpi_ierr)
      end block
      if (irank_mpi == 0) then 
          write(6,"(A24,I9,2E14.5)") 'i min/max(delta_U_bulk):', iters_full, min_delta_U_bulk, max_delta_U_bulk
      end if 
      min_delta_U_bulk =  1e20
      max_delta_U_bulk = -1e20
    end if
    
    if (force_Sf_x_pressure) then
      block
        integer  :: mpi_ierr
        real(dp) :: temp_real
        temp_real = min_Sf_x_mid
        call MPI_Allreduce(temp_real, min_Sf_x_mid, 1, MPI_DOUBLE_PRECISION, mpi_min, mpi_comm_world, mpi_ierr)
        temp_real = max_Sf_x_mid
        call MPI_Allreduce(temp_real, max_Sf_x_mid, 1, MPI_DOUBLE_PRECISION, mpi_max, mpi_comm_world, mpi_ierr)
      end block
      if (irank_mpi == 0) then 
          write(6,"(A20,I9,2E14.5)") 'i min/max(Sf_x_mid):', iters_full, min_Sf_x_mid, max_Sf_x_mid
      end if 
      min_Sf_x_mid     =  1e20
      max_Sf_x_mid     = -1e20
    end if
    
      if (prof_1d_avg_N>0) then
        block 
          integer  :: k
          real(dp) :: temp_real
          temp_real = 1.0d0/(prof_1d_avg_N+0.d0)
          !$acc parallel loop collapse(1) default(present)
          do    k=0, ny_loc-1
    !{x}(k)  =  {x}(k)*temp_real
    !__[PYTHON_HOLDER_PROF1D_COMPBULK_0]
          end do
          !$acc wait
        end block
    !__[PYTHON_HOLDER_PROF1D_COMPBULK_1]
    ! !$acc update host({x})
    ! block 
    !   integer :: k,iloc,mpi_ierr
    !   if (irank_mpi == 0) then
    !     open(19,file=prof_1d_FILE)
    !       do iloc=0,nproc_mpi-1
    !         if (iloc>0) then
    !           call MPI_Recv(prof_1d_temp, ny, MPI_DOUBLE_PRECISION, iloc, 1, mpi_comm_world, mpi_status_ignore, mpi_ierr)
    !         else 
    !           prof_1d_temp = {x}
    !         end if
    !         do  k=0,ny-1
    !           write(19,*) prof_1d_temp(k) 
    !         end do
    !       end do 
    !     close(19)
    !   else
    !     call MPI_Send({x}, ny, MPI_DOUBLE_PRECISION, 0, 1, mpi_comm_world, mpi_ierr)
    !   end if
    !   prof_1d_temp = -1
    ! end block
      end if
      block 
        integer :: k
        !$acc parallel loop collapse(1) default(present)
        do    k=0, ny_loc-1
    !__[PYTHON_HOLDER_PROF1D_COMPBULK_2]
    ! {x}(k)  =  0.
        end do
      end block
      prof_1d_avg_N = 0

    !$acc wait
    bulk_time_start_print  =  MPI_WTIME()
  end subroutine 

  subroutine get_U_bulk(U, U_bulk, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, dy_cell, utils_red_scalar, &
                        temp_arr_n_xzy, use_rough_surf, bulk_U_bot_inds, bulk_U_top_inds, bulk_use_rho_weighting)
    real(dp) :: U(-1:,-1:,-1:), rho_arr(-1:,-1:,-1:), U_bulk, vol_U_bulk, dy_cell(0:), temp_arr_n_xzy(0:,0:,0:)
    real(dp), contiguous :: utils_red_scalar(0:)
    logical  ::  use_rough_surf, bulk_use_rho_weighting
    integer  ::  nx_glob, nz_loc, ny_loc, bulk_U_bot_inds(0:,0:), bulk_U_top_inds(0:,0:)
    block
      integer :: i, j, k, mpi_ierr
        
      ! ---------------- check if (vol_U_bulk) was defined ----------------
      if (vol_U_bulk<0) then
        vol_U_bulk = 0
        if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
           !$acc parallel loop collapse(2) default(present) reduction(+:vol_U_bulk)
           do      j = 0, nz_loc -1
             do    i = 0, nx_glob-1
               do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
                 vol_U_bulk = vol_U_bulk + dy_cell(k)
               end do
             end do
           end do
#else
           call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
           !$acc parallel loop collapse(2) default(present) 
           do      j = 0, nz_loc -1
             do    i = 0, nx_glob-1
               do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
                 temp_arr_n_xzy(i,j,k) = dy_cell(k)
               end do
             end do
           end do
           call quick_binary_scalar_reduction(temp_arr_n_xzy, vol_U_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        else
#if !defined(_BIN_REDUCTION)
           !$acc parallel loop collapse(3) default(present) reduction(+:vol_U_bulk)
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                 vol_U_bulk = vol_U_bulk + dy_cell(k)
               end do
             end do
           end do
#else
           call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
           !$acc parallel loop collapse(3) default(present) 
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                 temp_arr_n_xzy(i,j,k) = dy_cell(k)
               end do
             end do
           end do
           call quick_binary_scalar_reduction(temp_arr_n_xzy, vol_U_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        end if
        !$acc wait
        block
          real(dp)  ::  temp_real
          temp_real  =  vol_U_bulk
          call MPI_Allreduce(temp_real, vol_U_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        end block
        !$acc wait
      end if
    
      if (.not.bulk_use_rho_weighting) then
        ! ---------------- process U_bulk ----------------
        U_bulk = 0
        if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
          !$acc parallel loop collapse(2) default(present) reduction(+:U_bulk)
          do      j = 0, nz_loc -1
            do    i = 0, nx_glob-1
              do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
                U_bulk = U_bulk + U(i,j,k)*dy_cell(k)
              end do
            end do
          end do
#else
          call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
          !$acc parallel loop collapse(2) default(present) 
          do      j = 0, nz_loc -1
            do    i = 0, nx_glob-1
              do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
                temp_arr_n_xzy(i,j,k) = U(i,j,k)*dy_cell(k)
              end do
            end do
          end do
          call quick_binary_scalar_reduction(temp_arr_n_xzy, U_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        else
#if !defined(_BIN_REDUCTION)
          !$acc parallel loop collapse(3) default(present) reduction(+:U_bulk)
          do      k  = 0, ny_loc -1
            do    j  = 0, nz_loc -1
              do  i  = 0, nx_glob-1
                U_bulk = U_bulk + U(i,j,k)*dy_cell(k)
              end do
            end do
          end do
#else
          call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
          !$acc parallel loop collapse(3) default(present) 
          do      k  = 0, ny_loc -1
            do    j  = 0, nz_loc -1
              do  i  = 0, nx_glob-1
                temp_arr_n_xzy(i,j,k) = U(i,j,k)*dy_cell(k)
              end do
            end do
          end do
          call quick_binary_scalar_reduction(temp_arr_n_xzy, U_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        end if
    
        !$acc wait
        block
          real(dp)  ::  temp_real
          temp_real  =  U_bulk/vol_U_bulk
          call MPI_Allreduce(temp_real, U_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        end block
        !$acc wait
      else
        block
          real(dp) :: U_rho__bulk, rho_bulk
    
            ! ---------------- check if (rho_bulk) was defined ----------------
            rho_bulk = 0
            if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
              !$acc parallel loop collapse(2) default(present) reduction(+:rho_bulk)
              do      j = 0, nz_loc -1
                do    i = 0, nx_glob-1
                  do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
                    rho_bulk  =  rho_bulk + (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                  end do
                end do
              end do
#else
              call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
              !$acc parallel loop collapse(2) default(present) 
              do      j = 0, nz_loc -1
                do    i = 0, nx_glob-1
                  do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
                    temp_arr_n_xzy(i,j,k)  =  (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                  end do
                end do
              end do
              call quick_binary_scalar_reduction(temp_arr_n_xzy, rho_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
            else
#if !defined(_BIN_REDUCTION)
              !$acc parallel loop collapse(3) default(present) reduction(+:rho_bulk)
              do      k  = 0, ny_loc -1
                do    j  = 0, nz_loc -1
                  do  i  = 0, nx_glob-1
                    rho_bulk  =  rho_bulk + (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                  end do
                end do
              end do
#else
              call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
              !$acc parallel loop collapse(3) default(present) 
              do      k  = 0, ny_loc -1
                do    j  = 0, nz_loc -1
                  do  i  = 0, nx_glob-1
                    temp_arr_n_xzy(i,j,k)  =  (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                  end do
                end do
              end do
              call quick_binary_scalar_reduction(temp_arr_n_xzy, rho_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
            endif
            !$acc wait
            block
              real(dp) :: temp_real
              temp_real  =  rho_bulk/vol_U_bulk
              call MPI_Allreduce(temp_real, rho_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
            end block
            !$acc wait
    
          ! ---------------- check if (U_rho__bulk) was defined ----------------
          U_rho__bulk = 0
          if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(2) default(present) reduction(+:U_rho__bulk)
            do      j  = 0, nz_loc -1
              do    i  = 0, nx_glob-1
                do  k  = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
                  U_rho__bulk  =  U_rho__bulk + U(i,j,k)*(0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
            !$acc parallel loop collapse(2) default(present) 
            do      j  = 0, nz_loc -1
              do    i  = 0, nx_glob-1
                do  k  = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
     temp_arr_n_xzy(i,j,k) = U(i,j,k)*(0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, U_rho__bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
          else
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(3) default(present) reduction(+:U_rho__bulk)
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                  U_rho__bulk  =  U_rho__bulk + U(i,j,k)*(0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
            !$acc parallel loop collapse(3) default(present) 
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
     temp_arr_n_xzy(i,j,k) = U(i,j,k)*(0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, U_rho__bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
          endif
    
          !$acc wait
          block
            real(dp) :: temp_real
            temp_real  =  U_rho__bulk/vol_U_bulk
            call MPI_Allreduce(temp_real, U_rho__bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
          end block
          !$acc wait
          U_bulk = U_rho__bulk/rho_bulk
        end block
      end if
    end block
  end subroutine

  subroutine get_T_bulk(T, T_bulk, vol_T_bulk, nx_glob, nz_loc, ny_loc, rho_arr, dy_cell, utils_red_scalar, &
                        temp_arr_n_xzy, use_rough_surf, bulk_T_bot_inds, bulk_T_top_inds, bulk_use_rho_weighting)
    real(dp) :: T(-1:,-1:,-1:), rho_arr(-1:,-1:,-1:), T_bulk, vol_T_bulk, dy_cell(0:), temp_arr_n_xzy(0:,0:,0:)
    real(dp), contiguous :: utils_red_scalar(0:)
    logical  ::  use_rough_surf, bulk_use_rho_weighting
    integer  ::  nx_glob, nz_loc, ny_loc, bulk_T_bot_inds(0:,0:), bulk_T_top_inds(0:,0:)
    block 
      integer :: i, j, k, mpi_ierr
    
      ! ---------------- check if (vol_T_bulk) was defined ----------------
      if (vol_T_bulk<0) then
        vol_T_bulk = 0
        if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
          !$acc parallel loop collapse(2) default(present) reduction(+:vol_T_bulk)
          do      j = 0, nz_loc -1
            do    i = 0, nx_glob-1
              do  k = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                vol_T_bulk = vol_T_bulk + dy_cell(k)
              end do
            end do
          end do
#else
          call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
          !$acc parallel loop collapse(2) default(present) 
          do      j = 0, nz_loc -1
            do    i = 0, nx_glob-1
              do  k = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                temp_arr_n_xzy(i,j,k) = dy_cell(k)
              end do
            end do
          end do
          call quick_binary_scalar_reduction(temp_arr_n_xzy, vol_T_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        else
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(3) default(present) reduction(+:vol_T_bulk)
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                  vol_T_bulk = vol_T_bulk + dy_cell(k)
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
            !$acc parallel loop collapse(3) default(present) 
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                  temp_arr_n_xzy(i,j,k) = dy_cell(k)
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, vol_T_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        end if
        !$acc wait
        block
          real(dp)  ::  temp_real
          temp_real  =  vol_T_bulk
          call MPI_Allreduce(temp_real, vol_T_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        end block
        !$acc wait
      end if
    
      if (.not.bulk_use_rho_weighting) then  
        ! ---------------- process T_bulk ----------------
        T_bulk = 0
        if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
          !$acc parallel loop collapse(2) default(present) reduction(+:T_bulk)
          do      j = 0, nz_loc -1
            do    i = 0, nx_glob-1
              do  k = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                T_bulk = T_bulk + T(i,j,k)*dy_cell(k)
              end do
            end do
          end do
#else
          call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
          !$acc parallel loop collapse(2) default(present) 
          do      j = 0, nz_loc -1
            do    i = 0, nx_glob-1
              do  k = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                temp_arr_n_xzy(i,j,k) = T(i,j,k)*dy_cell(k)
              end do
            end do
          end do
          call quick_binary_scalar_reduction(temp_arr_n_xzy, T_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        else
#if !defined(_BIN_REDUCTION)
          !$acc parallel loop collapse(3) default(present) reduction(+:T_bulk)
          do      k  = 0, ny_loc -1
            do    j  = 0, nz_loc -1
              do  i  = 0, nx_glob-1
                T_bulk = T_bulk + T(i,j,k)*dy_cell(k)
              end do
            end do
          end do
#else
          call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
          !$acc parallel loop collapse(3) default(present) 
          do      k  = 0, ny_loc -1
            do    j  = 0, nz_loc -1
              do  i  = 0, nx_glob-1
                temp_arr_n_xzy(i,j,k) = T(i,j,k)*dy_cell(k)
              end do
            end do
          end do
          call quick_binary_scalar_reduction(temp_arr_n_xzy, T_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
        end if
    
        !$acc wait
        block
          real(dp)  ::  temp_real
          temp_real  =  T_bulk/vol_T_bulk
          call MPI_Allreduce(temp_real, T_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        end block
        !$acc wait
      else
        block
          real(dp) :: T_rho__bulk, rho_bulk
    
          ! ---------------- check if (rho_bulk) was defined ----------------
          rho_bulk = 0
          if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(2) default(present) reduction(+:rho_bulk)
            do      j = 0, nz_loc -1
              do    i = 0, nx_glob-1
                do  k = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                  rho_bulk  =  rho_bulk + (rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
            !$acc parallel loop collapse(2) default(present) 
            do      j = 0, nz_loc -1
              do    i = 0, nx_glob-1
                do  k = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                  temp_arr_n_xzy(i,j,k)  =  (rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, rho_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
          else
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(3) default(present) reduction(+:rho_bulk)
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                  rho_bulk  =  rho_bulk + (rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
            !$acc parallel loop collapse(3) default(present) 
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                  temp_arr_n_xzy(i,j,k)  =  (rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, rho_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
          endif
          !$acc wait
          block
            real(dp) :: temp_real
            temp_real  =  rho_bulk/vol_T_bulk
            call MPI_Allreduce(temp_real, rho_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
          end block
          !$acc wait
    
          ! ---------------- check if (T_rho__bulk) was defined ----------------
          T_rho__bulk = 0
          if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(2) default(present) reduction(+:T_rho__bulk)
            do      j  = 0, nz_loc -1
              do    i  = 0, nx_glob-1
                do  k  = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                  T_rho__bulk  =  T_rho__bulk + T(i,j,k)*(rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
            !$acc parallel loop collapse(2) default(present) 
            do      j  = 0, nz_loc -1
              do    i  = 0, nx_glob-1
                do  k  = bulk_T_bot_inds(i,j),bulk_T_top_inds(i,j)
                  temp_arr_n_xzy(i,j,k)  =  T(i,j,k)*(rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, T_rho__bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
          else
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(3) default(present) reduction(+:T_rho__bulk)
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                  T_rho__bulk  =  T_rho__bulk + T(i,j,k)*(rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
            !$acc parallel loop collapse(3) default(present) 
            do      k  = 0, ny_loc -1
              do    j  = 0, nz_loc -1
                do  i  = 0, nx_glob-1
                  temp_arr_n_xzy(i,j,k)  =  T(i,j,k)*(rho_arr(i,j,k))*dy_cell(k)
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, T_rho__bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
          endif
    
          !$acc wait
          block
            real(dp) :: temp_real
            temp_real  =  T_rho__bulk/vol_T_bulk
            call MPI_Allreduce(temp_real, T_rho__bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
          end block
          !$acc wait
          T_bulk = T_rho__bulk/rho_bulk
        end block
      end if
    end block
  end subroutine

  subroutine get_rho_mu_bulk(rho_bulk, mu_bulk, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, mu_arr, dy_cell, &
                             utils_red_scalar, temp_arr_n_xzy, use_rough_surf, bulk_U_bot_inds, bulk_U_top_inds)
    real(dp) :: rho_bulk, mu_bulk, rho_arr(-1:,-1:,-1:), mu_arr(-1:,-1:,-1:), vol_U_bulk, dy_cell(0:), temp_arr_n_xzy(0:,0:,0:)
    real(dp), contiguous :: utils_red_scalar(0:)
    logical  ::  use_rough_surf
    integer  ::  nx_glob, nz_loc, ny_loc, bulk_U_bot_inds(0:,0:), bulk_U_top_inds(0:,0:)

    block
      integer :: i, j, k, mpi_ierr
    
      if (vol_U_bulk<0) error stop 'vol_U_bulk<0'
    
      ! ---------------- check if (rho_bulk) was defined ----------------
      rho_bulk = 0
      if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(2) default(present) reduction(+:rho_bulk)
        do      j = 0, nz_loc -1
          do    i = 0, nx_glob-1
            do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
              rho_bulk  =  rho_bulk + (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
        !$acc parallel loop collapse(2) default(present) 
        do      j = 0, nz_loc -1
          do    i = 0, nx_glob-1
            do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
              temp_arr_n_xzy(i,j,k)  =  (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, rho_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
      else
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(3) default(present) reduction(+:rho_bulk)
        do      k  = 0, ny_loc -1
          do    j  = 0, nz_loc -1
            do  i  = 0, nx_glob-1
              rho_bulk  =  rho_bulk + (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
        !$acc parallel loop collapse(3) default(present) 
        do      k  = 0, ny_loc -1
          do    j  = 0, nz_loc -1
            do  i  = 0, nx_glob-1
              temp_arr_n_xzy(i,j,k)  =  (0.5_dp*rho_arr(i,j,k) + 0.5_dp*rho_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, rho_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
      endif
      !$acc wait
      block
        real(dp) :: temp_real
        temp_real  =  rho_bulk/vol_U_bulk
        call MPI_Allreduce(temp_real, rho_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
      end block
      !$acc wait
    
      ! ---------------- check if (mu_bulk) was defined ----------------
      mu_bulk = 0
      if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(2) default(present) reduction(+:mu_bulk)
        do      j = 0, nz_loc -1
          do    i = 0, nx_glob-1
            do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
              mu_bulk  =  mu_bulk + (0.5_dp*mu_arr(i,j,k) + 0.5_dp*mu_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
        !$acc parallel loop collapse(2) default(present) 
        do      j = 0, nz_loc -1
          do    i = 0, nx_glob-1
            do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
              temp_arr_n_xzy(i,j,k)  =  (0.5_dp*mu_arr(i,j,k) + 0.5_dp*mu_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, mu_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
      else
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(3) default(present) reduction(+:mu_bulk)
        do      k  = 0, ny_loc -1
          do    j  = 0, nz_loc -1
            do  i  = 0, nx_glob-1
              mu_bulk  =  mu_bulk + (0.5_dp*mu_arr(i,j,k) + 0.5_dp*mu_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
        !$acc parallel loop collapse(3) default(present) 
        do      k  = 0, ny_loc -1
          do    j  = 0, nz_loc -1
            do  i  = 0, nx_glob-1
              temp_arr_n_xzy(i,j,k)  =  (0.5_dp*mu_arr(i,j,k) + 0.5_dp*mu_arr(i-1,j,k))*dy_cell(k)
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, mu_bulk, nx_glob*nz_loc*ny_loc, 0, utils_red_scalar)
#endif
      endif
      !$acc wait
      block
        real(dp) :: temp_real
        temp_real  =  mu_bulk/vol_U_bulk
        call MPI_Allreduce(temp_real, mu_bulk, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
      end block
      !$acc wait
    end block
  end subroutine 
  
  subroutine write_prof_1d(x, my_file, mpi_divs_y, mpi_pos_z, mpi_pos_y, obj_ranks_0zy, prof_1d_temp, ny_loc)
    use diezdecomp_api_generic, only: diezdecomp_parsed_mpi_ranks
    implicit none
    real(dp)                           ::  x(0:), prof_1d_temp(0:)
    character(len=*)                   ::  my_file
    integer                            ::  mpi_divs_y, mpi_pos_z, mpi_pos_y, ny_loc
    type(diezdecomp_parsed_mpi_ranks)  ::  obj_ranks_0zy
    !$acc update host(x)
    block 
      integer :: k,iloc,ny_here,mpi_ierr,other
      call mpi_barrier(mpi_comm_world,mpi_ierr)
      if   (mpi_pos_z == 0) then
        if (mpi_pos_y == 0) then
          open(19,file=my_file)
            do iloc=0,mpi_divs_y-1
              if (iloc>0) then
                other = obj_ranks_0zy%mpi_ranks(0,0,iloc)
                call MPI_Recv(ny_here     , 1      , MPI_INT             , other, 1, mpi_comm_world, mpi_status_ignore, mpi_ierr)
                call MPI_Recv(prof_1d_temp, ny_here, MPI_DOUBLE_PRECISION, other, 1, mpi_comm_world, mpi_status_ignore, mpi_ierr)
              else 
                ny_here                   = ny_loc
                prof_1d_temp(0:ny_here-1) = x(0:ny_here-1)
              end if
              do  k=0,ny_here-1
                write(19,*) prof_1d_temp(k) 
              end do
            end do 
          close(19)
        else
          other = obj_ranks_0zy%mpi_ranks(0,0,0)
          call MPI_Send(ny_loc, 1     , MPI_INT             , other, 1, mpi_comm_world, mpi_ierr)
          call MPI_Send(x     , ny_loc, MPI_DOUBLE_PRECISION, other, 1, mpi_comm_world, mpi_ierr)
        end if
      end if
      prof_1d_temp = -1
      call mpi_barrier(mpi_comm_world,mpi_ierr)
    end block
  end subroutine
end module
