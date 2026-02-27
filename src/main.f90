program main
    ! ---------------------------------------- begin: use modules ----------------------------------------
    use  ieee_arithmetic
    use  openacc
    use  mpi
    use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
#if defined(_USE_HIPFFT)
    use  hipfort
    use  hipfort_check 
    use  hipfort_hipfft
    use, intrinsic :: iso_c_binding, only: C_INT, C_INTPTR_T, C_PTR, C_LOC 
#else
#if defined(_USE_FFTW)
    use, intrinsic :: iso_c_binding
#else 
    use  cufft
    use  cudafor
#endif
#endif
    use  diezdecomp_api_generic
    use  diezDecomp_api_ibm
    use  poisson_solver_multigpu, only: init_spec_poisson_multigpu, run_spec_poisson_multigpu, get_lo_n_bounds_1d, &
                                        planinfo_hipfft ! the datatype might be be used, and contains only generic info
    use  mod_avg_apply
    use  mod_cfd_subroutines
    use  mod_get_bulk
    use  mod_io
    use  mod_utils
    use  mod_cfd_ibm
    ! ---------------------------------------- end: use modules ----------------------------------------

    ! ---------------------------------------- begin: declare variables ----------------------------------------
    implicit none
#if defined(_USE_FFTW)
    include 'fftw3.f03'
#endif
    integer  ::  irank_mpi, nproc_mpi, nskip, nhalo, kmin_Vhat
#if !defined(_USE_FFTW)
    integer(acc_device_kind) ::dev_type
#endif
    include "inc_declare_main_params.f90"
    include "poisson_solver_mgpu_dtdma/declare_poisson_solver_multigpu.f90"
    include "ibm_setup/declare_ibm_params.f90"
    ! ---------------------------------------- end: declare variables ----------------------------------------

    block 
        integer:: mpi_ierr, ierr, ndev, mydev,i
        call mpi_init     (mpi_ierr)
        call MPI_COMM_RANK(mpi_comm_world,irank_mpi,mpi_ierr)
        call MPI_COMM_SIZE(mpi_comm_world,nproc_mpi,mpi_ierr)
#if !defined(_USE_FFTW)
        dev_type  =  acc_get_device_type()
#if defined(_USE_HIPFFT)
        if (irank_mpi==0) then 
          block 
            character(len=16) :: aux_str          
            call get_command_argument(3, aux_str)
            read(aux_str,*) ndev
          end block
        end if
        call MPI_BCAST(ndev,1,MPI_INTEGER, 0, mpi_comm_world, mpi_ierr)
        mydev     =  mod(irank_mpi,ndev)
#else 
        mpi_ierr  =  cudaGetDeviceCount(ndev)
        mydev     =  mod(irank_mpi,ndev)
        ierr      =  ierr + cudaSetDevice(mydev)
#endif
        call acc_set_device_num(mydev,dev_type)
        call acc_init(dev_type)
        if (irank_mpi==0) then 
          write(6,*) '--------------- MPI GPU Binding ---------------';flush(6)
        end if
        do i=0,nproc_mpi-1
         call MPI_BARRIER(mpi_comm_world,mpi_ierr)
         if (i==irank_mpi) then 
           write(6,'(A33,4I6)') '    (irank_mpi, nproc_mpi, mydev, ndev): ',irank_mpi,nproc_mpi,mydev,ndev;flush(6)
         end if
         call MPI_BARRIER(mpi_comm_world,mpi_ierr)
        end do
#endif
    end block 
    include "inc_set_read_params.f90"
    include "poisson_solver_mgpu_dtdma/alloc_poisson_solver_multigpu.f90"
    nhalo      =  1
    nskip      =  0
    kmin_Vhat  =  0
    if (mpi_pos_y==0) kmin_Vhat = 1
    include "inc_set_alloc_arrays.f90"

if (use_rough_surf) then
    include "ibm_setup/set_ibm_bands.f90"
end if

    ! --------------- allocate DNS variables ---------------
    include "inc_read_UVWT.f90"

    ! --------------- begin: initial setup variables ---------------
    iters_full = 0

    if (use_rough_surf) then
      call cfd_periodicity(hl, U_hat, V_hat, W_hat, T, work_tr, mpi_pos_y, mpi_divs_y, &
                           wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                           wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
      call  ibm_part_2(U_hat, V_hat, W_hat, T, &
                       rough_U_bot_inds, rough_U_top_inds, rough_W_bot_inds, rough_W_top_inds , &
                       rough_V_bot_inds, rough_V_top_inds, rough_T_bot_inds, rough_T_top_inds , &
                       wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                       wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                       nx_glob      , nz_loc       , ny_loc                      )
      if (use_rough_surf_quady) then 
      call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                           wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                           wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
      call  ibm_part_2(U, V, W, T, &
                       rough_U_bot_inds, rough_U_top_inds, rough_W_bot_inds, rough_W_top_inds , &
                       rough_V_bot_inds, rough_V_top_inds, rough_T_bot_inds, rough_T_top_inds , &
                       wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                       wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                       nx_glob      , nz_loc       , ny_loc                      )
      end if
      call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                           wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                           wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
      call ibm_part_1(obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div, U, V, W, T, div_now, buffer_ibm, &
                      wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                      wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                      mpi_pos_y    , mpi_divs_y   , use_incompressible)
      call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                           wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                           wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
    else 
      call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                           wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                           wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
    end if

    if (.not.use_incompressible) then 
      call update_props(T_ext, T, work_tr, wall_BC_T_bot, wall_BC_T_top, &
                          mu_arr, cond_arr, rho_arr, buoyancy_arr, enth_arr, cp_arr , div_jacobian_arr, &
                          rough_T_bot_inds, rough_T_top_inds, nx_glob, nz_loc, ny_loc, &
                          lo_z_blocks_mpi, lo_y_blocks_mpi, mpi_pos_z, mpi_pos_y, &
                          use_incompressible, hl, use_rough_surf)
      if (force_Sf_x_pressure) then 
        block 
          integer :: i,j,k
          !$acc parallel loop  collapse(3) default(present)
          do     k = -1, ny_loc
            do   j = -1, nz_loc
              do i = -1, nx_glob
                rho_arr_old(i,j,k)  =  rho_arr(i,j,k)
              end do
            end do
          end do
        end block
      end if 
    end if

    if (use_utau_work) then 
      call build_utau_work(U, V, W, Ruw_pad_ddy_a, Ruw_pad_ddy_b, Ruw_pad_ddy_c, &
                             mu_arr, rho_arr, T_ext, utau_now, &
                             inv_dx, inv_dz, inv_dy_p, &
                             wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, &
                             wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, &
                             inv_dx__2, inv_dz__2, work_tr, buffer_ibm, &
                             nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, &
                             hl, obj_ibm_div, use_incompressible, use_utau_work, use_rough_surf)
      block 
        integer :: i,j,k
        !$acc parallel loop  collapse(3) default(present)
        do     k = 0, ny_loc-1
          do   j = 0, nz_loc-1
            do i = 0, nx_glob-1
              utau_old(i,j,k)  =  utau_now(i,j,k)
            end do
          end do
        end do
      end block
    end if

    if (.not.use_incompressible) then 
      call build_rhs_energy(T, cond_arr, div_now, div_jacobian_arr, RHS_energy, &
                              inv_dx2, inv_dz2, inv_dx__2, inv_dz__2, &
                              Ruw_d2dy2_a, Ruw_d2dy2_b, Ruw_d2dy2_c, Ruw_ddy_a, Ruw_ddy_b, Ruw_ddy_c, Sq, &
                              utau_now, utau_old, nx_glob, nz_loc, ny_loc, &
                              rough_T_bot_inds, rough_T_top_inds, &
                              hl, work_tr, buffer_ibm, mpi_pos_y, mpi_divs_y, wall_BC_T_bot, wall_BC_T_top, &
                              use_utau_work, use_incompressible, use_rough_surf, obj_ibm_T)
    end if
    ! --------------- end: initial setup variables ---------------

    ! ------------------------ general checks ------------------------
    if (1.000000000000001 == 1.d0) then
        write(6,*) 'ERROR: -r8 flag disabled. Simulation terminated.'; flush(6)
        error stop 'ERROR: -r8 flag disabled. Simulation terminated.'
    else
        if (irank_mpi ==0) then 
          write(6,*) '\nSUCCESS: -r8 flag enabled. Simulation will run.'; flush(6)
        end if
    end if
    flush(6)

    !     ------------------------ report Ub,Tb ------------------------
    !---------------------------------------- initialize fft solver ----------------------------------------
    call init_spec_poisson_multigpu(inv_poi_B, poi_A, poi_C, poi_slabs_tips, &
                                    ts_band_a, ts_band_b, ts_cp, a_x, a_z, poi_is_00x, poi_is_00z, &
                                    P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00, P, &
                                    spec_x, spec_z, P_z, &
                                    nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                                    nx_transp, nz_transp_py, nz_loc, ny_loc, inv_dx, inv_dz, &
                                    plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z, utils_red_scalar, &
                                    tr_poi_yz, buffer_transp, lo_xt_blocks_mpi)
    !---------------------------------------- header ----------------------------------------
if (irank_mpi==0) then 
  write(6,*) "\n----------------- Iteration Parameters -----------------"
  write(6,*) "nproc_mpi                 : ", nproc_mpi
  write(6,*) "mpi_divs_z                : ", mpi_divs_z
  write(6,*) "mpi_divs_y                : ", mpi_divs_y
  write(6,*) "nstep                     : ", nstep
  write(6,*) "bulk_nprint               : ", bulk_nprint
  write(6,*) "avg_iter_start            : ", avg_iter_start
  write(6,*) "avg_nraw                  : ", avg_nraw      
  write(6,*) "avg_freq                  : ", avg_freq      
  write(6,*) "avg_bin_nbins             : ", avg_bin_nbins 
  write(6,*) "snap_iter_start           : ", snap_iter_start
  write(6,*) "snap_freq                 : ", snap_freq      
  write(6,*) "nu                        : ", nu
  write(6,*) "cond                      : ", cond  
  write(6,*) "Sf_x                      : ", Sf_x  
  write(6,*) "Sf_z                      : ", Sf_z  
  write(6,*) "Sq                        : ", Sq    
  write(6,*) "dtmax                     : ", dtmax
  write(6,*) "dx                        : ", 1.d0/inv_dx
  write(6,*) "dz                        : ", 1.d0/inv_dz
  write(6,*) "nx_glob                   : ", nx_glob
  write(6,*) "nz_glob                   : ", nz_glob
  write(6,*) "ny_glob                   : ", ny_glob
  write(6,*) "use_rough_surf            : ", use_rough_surf
  write(6,*) "use_rough_surf_quady      : ", use_rough_surf_quady
  write(6,*) "use_filtering             : ", use_filtering
  write(6,*) "use_incompressible        : ", use_incompressible
  write(6,*) "use_utau_work             : ", use_utau_work
  write(6,*) "wall_BC_T_bot             : ", wall_BC_T_bot
  write(6,*) "wall_BC_T_top             : ", wall_BC_T_top
  write(6,*) "wall_BC_U_bot             : ", wall_BC_U_bot
  write(6,*) "wall_BC_U_top             : ", wall_BC_U_top
  write(6,*) "wall_BC_V_bot             : ", wall_BC_V_bot
  write(6,*) "wall_BC_V_top             : ", wall_BC_V_top
  write(6,*) "wall_BC_W_bot             : ", wall_BC_W_bot
  write(6,*) "wall_BC_W_top             : ", wall_BC_W_top
  write(6,*) "Gb_x                      : ", Gb_x  
  write(6,*) "Gb_y                      : ", Gb_y  
  write(6,*) "Gb_z                      : ", Gb_z  
  write(6,*) "filtering_steps           : ", filtering_steps 
  write(6,*) "filtering_first           : ", filtering_first 
  write(6,*) "filtering_nstop           : ", filtering_nstop 
  write(6,*) "extrap_T_wall_0           : ", extrap_T_wall_0
  write(6,*) "extrap_T_wall_1           : ", extrap_T_wall_1
  write(6,*) "extrap_T_wall_m1          : ", extrap_T_wall_m1
  write(6,*) "extrap_T_wall_m2          : ", extrap_T_wall_m2
  write(6,*) "bulk_target_Ub            : ", bulk_target_Ub
  write(6,*) "bulk_use_target_Ub        : ", bulk_use_target_Ub
  write(6,*) "bulk_use_rho_weighting    : ", bulk_use_rho_weighting
  write(6,*) "bulk_target_Reb           : ", bulk_target_Reb
  write(6,*) "Ly                        : ", Ly
  write(6,*) "bulk_use_target_Reb       : ", bulk_use_target_Reb
  write(6,*) "force_Sf_x_pressure       : ", force_Sf_x_pressure
  write(6,*) "force_Sf_x_accel          : ", force_Sf_x_accel
  

  block 
    logical :: aux_bool

#if defined(_USE_HIPFFT)
    aux_bool = .true.
#else
    aux_bool = .false.
#endif
  write(6,*) 'defined(_USE_HIPFFT)      : ', aux_bool

#if defined(_BIN_REDUCTION)
    aux_bool = .true.
#else
    aux_bool = .false.
#endif
  write(6,*) 'defined(_BIN_REDUCTION)   : ', aux_bool

#if defined(_USE_FFTW)
    aux_bool = .true.
#else
    aux_bool = .false.
#endif
  write(6,*) 'defined(_USE_FFTW)        : ', aux_bool

  end block
  write(6,*) ' '
  flush(6)
endif

    block 
        real(dp) :: dtmax_old
        dtmax_old  =  dtmax
        inv_dtmax  =  0.
        dtmax      =  0.
        call build_vars(wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                        wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                        extrap_T_wall_0, extrap_T_wall_1, extrap_T_wall_m1, extrap_T_wall_m2, &
                        buffer_ibm, work_tr, dtmax, inv_dtmax, rho_min, Sf_x_mid, &
                        Ruw_ddy_a  , Ruw_ddy_b  , Ruw_ddy_c, Ruw_d2dy2_a, Ruw_d2dy2_b, Ruw_d2dy2_c, &
                        nu, cond, inv_dx, inv_dx2, inv_dx__2, inv_dz, inv_dz2, inv_dz__2, &
                        inv_dy_p, inv_dy_p_pad, inv_dy_u, Sf_x, Sf_z, Sq, &
                        utau_now, utau_old, temp_arr_n_xzy, Ru_old, Rw_old, Rv_old, Rt_old, P, inv_rho_dt_scale_Sf_x, &
                        mu_arr, cond_arr, buoyancy_arr, rho_arr, rho_arr_old, &
                        enth_arr, cp_arr, div_jacobian_arr, &
                        U, V, W, T, T_ext, RHS_energy, &
                        U_hat, V_hat, W_hat, div_now, slab_00k, P_now, P_old, &
                        Gb_x, Gb_y, Gb_z, rho_profile_y, iters_full, &
                        Rv_interp_uw_bottom, Rv_interp_uw_top, &
                        Rv_ddy_a  , Rv_ddy_b  , Rv_ddy_c, Rv_d2dy2_a, Rv_d2dy2_b, Rv_d2dy2_c, &
                        tr_red1d_B, tr_red1d_B_1d, utils_red_scalar, &
                        lo_z_blocks_mpi, lo_y_blocks_mpi, &
                        mpi_pos_z, mpi_pos_y, mpi_divs_y, nx_glob, nz_loc, ny_loc, nz_glob, kmin_vhat, &
                        ny_glob, ny_loc_00k, irank_mpi, nproc_mpi, mpi_divs_z, use_rough_surf_quady, &
                        use_rough_surf, use_incompressible, use_utau_work, force_Sf_x_semilocal, force_Sf_x_pressure, &
                        obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div, &
                        rough_U_bot_inds, rough_U_top_inds, rough_W_bot_inds, rough_W_top_inds, &
                        rough_V_bot_inds, rough_V_top_inds, rough_T_bot_inds, rough_T_top_inds, &
                        Ruw_pad_ddy_a, Ruw_pad_ddy_b, Ruw_pad_ddy_c, &
                        tr_red1d_fwd , tr_red1d_bwd, tr_io_fwd, hl)
        dtmax      =  dtmax_old
        inv_dtmax  =  1.d0/dtmax
    end block

  call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                       wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                       wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
  if (irank_mpi == 0) then 
    write(6,*) '\n-------------------- Begin: iterations --------------------';flush(6)
  end if

      bulk_last_iter         = -1
      !$acc wait 
      bulk_time_start_print  =  MPI_WTIME()
      time_start_total       =  MPI_WTIME()
      call bulk_compute(U, V, W, T, P_now, P, rho_arr, vol_U_bulk, vol_T_bulk, vol_div_avg, nx_glob, nz_loc, ny_loc, &
                        dy_cell, temp_arr_n_xzy, use_rough_surf, &
                        bulk_U_bot_inds, bulk_U_top_inds, bulk_T_bot_inds, bulk_T_top_inds, bulk_use_rho_weighting, &
                        inv_dx, inv_dz, inv_dy_p, &
                        bulk_last_iter, iters_full, bulk_use_target_Ub, force_Sf_x_pressure, use_incompressible, &
                        dtmax, inv_dtmax, bulk_nprint, snap_iter_start, snap_freq, nstep, &
                        min_delta_U_bulk, max_delta_U_bulk, min_Sf_x_mid, max_Sf_x_mid, prof_1d_avg_N, &
                        !__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]
                        mpi_divs_y, mpi_pos_z, mpi_pos_y, obj_ranks_0zy, prof_1d_temp, utils_red_scalar, &
                        bulk_time_start_print, bulk_time_finish_print, irank_mpi)

   !---------------------------------------- explicit gpu/cpu memory measurement ----------------------------------------
    call execute_command_line("nvidia-smi          > nvidia_smi.txt"   , wait=.true.) ! print gpu memory usage, etc.
    call execute_command_line("top -n 1 -o %MEM -b > cpu_top_stats.txt", wait=.true.) ! print cpu stats
   
   !---------------------------------------- main iterations ----------------------------------------
    do iters_full = 1, nstep
      call build_vars(wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                      wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                      extrap_T_wall_0, extrap_T_wall_1, extrap_T_wall_m1, extrap_T_wall_m2, &
                      buffer_ibm, work_tr, dtmax, inv_dtmax, rho_min, Sf_x_mid, &
                      Ruw_ddy_a  , Ruw_ddy_b  , Ruw_ddy_c, Ruw_d2dy2_a, Ruw_d2dy2_b, Ruw_d2dy2_c, &
                      nu, cond, inv_dx, inv_dx2, inv_dx__2, inv_dz, inv_dz2, inv_dz__2, &
                      inv_dy_p, inv_dy_p_pad, inv_dy_u, Sf_x, Sf_z, Sq, &
                      utau_now, utau_old, temp_arr_n_xzy, Ru_old, Rw_old, Rv_old, Rt_old, P, inv_rho_dt_scale_Sf_x, &
                      mu_arr, cond_arr, buoyancy_arr, rho_arr, rho_arr_old, &
                      enth_arr, cp_arr, div_jacobian_arr, &
                      U, V, W, T, T_ext, RHS_energy, &
                      U_hat, V_hat, W_hat, div_now, slab_00k, P_now, P_old, &
                      Gb_x, Gb_y, Gb_z, rho_profile_y, iters_full, &
                      Rv_interp_uw_bottom, Rv_interp_uw_top, &
                      Rv_ddy_a  , Rv_ddy_b  , Rv_ddy_c, Rv_d2dy2_a, Rv_d2dy2_b, Rv_d2dy2_c, &
                      tr_red1d_B, tr_red1d_B_1d, utils_red_scalar, &
                      lo_z_blocks_mpi, lo_y_blocks_mpi, &
                      mpi_pos_z, mpi_pos_y, mpi_divs_y, nx_glob, nz_loc, ny_loc, nz_glob, kmin_vhat, &
                      ny_glob, ny_loc_00k, irank_mpi, nproc_mpi, mpi_divs_z, use_rough_surf_quady, &
                      use_rough_surf, use_incompressible, use_utau_work, force_Sf_x_semilocal, force_Sf_x_pressure, &
                      obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div, &
                      rough_U_bot_inds, rough_U_top_inds, rough_W_bot_inds, rough_W_top_inds, &
                      rough_V_bot_inds, rough_V_top_inds, rough_T_bot_inds, rough_T_top_inds, &
                      Ruw_pad_ddy_a, Ruw_pad_ddy_b, Ruw_pad_ddy_c, &
                      tr_red1d_fwd , tr_red1d_bwd, tr_io_fwd, hl)
      call run_spec_poisson_multigpu(inv_poi_B, poi_A, poi_C, poi_slabs_tips, &
                                     ts_band_a, ts_band_b, ts_cp, a_x, a_z, poi_is_00x, poi_is_00z, &
                                     P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00, P, &
                                     spec_x, spec_z, P_z, trispec, &
                                     nx_glob, nz_glob, n_fft_x, n_fft_z, mpi_divs_z, mpi_divs_y, mpi_pos_z, mpi_pos_y, &
                                     nx_transp, nz_transp_py, nx_loc, nz_loc, ny_loc, &
                                     plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z, &
                                     tr_poi_xy, tr_poi_yx, tr_poi_yz, tr_poi_zy, buffer_transp)
      call update_vars(nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, kmin_Vhat, &
                       P, &
                       P_per, P_now, P_old, work_tr, &
                       U    , V    , W,dtmax, inv_dx, inv_dz, &
                       U_hat, V_hat, W_hat, &
                       rho_arr, mu_arr, &
                       dPdy_ddy_hi, dPdy_ddy_lo, utils_red_scalar, &
                       wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, &
                       wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, &
                       bulk_target_Ub, vol_U_bulk, dy_cell, temp_arr_n_xzy, &
                       bulk_U_bot_inds, bulk_U_top_inds, min_delta_U_bulk, max_delta_U_bulk, Ly, bulk_target_Reb, &
                       bulk_use_rho_weighting, force_Sf_x_accel, force_Sf_x_pressure, bulk_use_target_Reb, &
                       bulk_use_target_Ub, use_rough_surf, init_Sf_x_mid_prev, &
                       inv_rho_dt_scale_Sf_x, Sf_x_mid, Sf_x_mid_prev, min_Sf_x_mid, max_Sf_x_mid, rho_min, &
                       use_incompressible, hl)

if ((mod(iters_full,bulk_nprint)==0).or.(iters_full<=100).or.(iters_full==nstep)) then ! numerically, the first 100 iterations are the most interesting
  call bulk_compute(U, V, W, T, P_now, P, rho_arr, vol_U_bulk, vol_T_bulk, vol_div_avg, nx_glob, nz_loc, ny_loc, &
                    dy_cell, temp_arr_n_xzy, use_rough_surf, &
                    bulk_U_bot_inds, bulk_U_top_inds, bulk_T_bot_inds, bulk_T_top_inds, bulk_use_rho_weighting, &
                    inv_dx, inv_dz, inv_dy_p, &
                    bulk_last_iter, iters_full, bulk_use_target_Ub, force_Sf_x_pressure, use_incompressible, &
                    dtmax, inv_dtmax, bulk_nprint, snap_iter_start, snap_freq, nstep, &
                    min_delta_U_bulk, max_delta_U_bulk, min_Sf_x_mid, max_Sf_x_mid, prof_1d_avg_N, &
                    !__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]
                    mpi_divs_y, mpi_pos_z, mpi_pos_y, obj_ranks_0zy, prof_1d_temp, utils_red_scalar, &
                    bulk_time_start_print, bulk_time_finish_print, irank_mpi)
end if

if (use_filtering) then
  if (((mod(iters_full, filtering_steps) == 0).or.(iters_full==filtering_first)).and.(iters_full < filtering_nstop)) then 
    call filtering_UVWT(iters_full, nx_glob, nz_loc, ny_loc, kmin_Vhat, irank_mpi, temp_arr_n_xzy, U, V, W, T)
    call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                         wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                         wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
    if (.not.use_incompressible) then 
      call update_props(T_ext, T, work_tr, wall_BC_T_bot, wall_BC_T_top, &
                          mu_arr, cond_arr, rho_arr, buoyancy_arr, enth_arr, cp_arr , div_jacobian_arr, &
                          rough_T_bot_inds, rough_T_top_inds, nx_glob, nz_loc, ny_loc, &
                          lo_z_blocks_mpi, lo_y_blocks_mpi, mpi_pos_z, mpi_pos_y, &
                          use_incompressible, hl, use_rough_surf)
      call build_rhs_energy(T, cond_arr, div_now, div_jacobian_arr, RHS_energy, &
                              inv_dx2, inv_dz2, inv_dx__2, inv_dz__2, &
                              Ruw_d2dy2_a, Ruw_d2dy2_b, Ruw_d2dy2_c, Ruw_ddy_a, Ruw_ddy_b, Ruw_ddy_c, Sq, &
                              utau_now, utau_old, nx_glob, nz_loc, ny_loc, &
                              rough_T_bot_inds, rough_T_top_inds, &
                              hl, work_tr, buffer_ibm, mpi_pos_y, mpi_divs_y, wall_BC_T_bot, wall_BC_T_top, &
                              use_utau_work, use_incompressible, use_rough_surf, obj_ibm_T)
    end if
  end if 
endif
    call avg_apply(nx_glob, nz_glob, ny_glob, nz_loc, ny_loc, irank_mpi, nproc_mpi, &
        !__[PYTHON_HOLDER_AVG_VARS]
        !__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]
    avg_bin_nbins, prof_1d_avg_N, avg_N_now, avg_nraw, iters_full, avg_iter_start, avg_freq, avg_inv_nraw, binary_cycle,&
    U,V,W,T,P_now,P_old,P_per,P,slab_00k,&
    mu_arr           , cond_arr         , rho_arr          , buoyancy_arr     , &
    enth_arr         , cp_arr           , div_jacobian_arr , &
    temp_arr_n_xzy, prof_1d_temp, tr_red1d_fwd, tr_red1d_bwd, tr_io_fwd, &
    mpi_divs_z, ny_loc_00k, tr_red1d_B, tr_red1d_B_1d, work_tr)

    if (((iters_full >= snap_iter_start).and.(mod(iters_full-snap_iter_start,snap_freq)==0)).or.(iters_full == nstep)) then
      block
        character(len=36)  ::  my_file
        write(my_file,'(A23,I0.9,A4)') './end/array_end_U_iter_',iters_full,'.dat'
        call cfd_write_arr(my_file, U, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
        write(my_file,'(A23,I0.9,A4)') './end/array_end_V_iter_',iters_full,'.dat'
        call cfd_write_arr(my_file, V, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
        write(my_file,'(A23,I0.9,A4)') './end/array_end_W_iter_',iters_full,'.dat'
        call cfd_write_arr(my_file, W, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
        write(my_file,'(A23,I0.9,A4)') './end/array_end_T_iter_',iters_full,'.dat'
        call cfd_write_arr(my_file, T, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
      end block
      if (use_incompressible) then 
        block
          character(len=36)  ::  my_file
          write(my_file,'(A23,I0.9,A4)') './end/array_end_P_iter_',iters_full,'.dat'
          call cfd_write_arr(my_file, P, 0, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                             irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
        end block
      else
        block
          character(len=40)  ::  my_file
          write(my_file,'(A27,I0.9,A4)') './end/array_end_P_now_iter_',iters_full,'.dat'
          call cfd_write_arr(my_file, P_now, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                             irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
          write(my_file,'(A27,I0.9,A4)') './end/array_end_P_old_iter_',iters_full,'.dat'
          call cfd_write_arr(my_file, P_old, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                             irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
        end block
      end if
    end if
  end do
  block
    integer :: mpi_ierr
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
  end block
  if (irank_mpi==0) then 
    write(6,*) '-------------------- End: iterations --------------------';flush(6)
  end if
  block
    integer :: mpi_ierr
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    call mpi_finalize(mpi_ierr)
  end block
end program
