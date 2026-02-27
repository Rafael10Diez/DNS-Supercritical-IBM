real(dp)   ::  dtmax, inv_dx, inv_dz, nu, cond, Sf_x, Sf_z, Sq, Gb_x, Gb_z, Gb_y, &
               inv_dx__2, inv_dz__2, inv_dx2, inv_dz2, inv_dtmax, vol_U_bulk, vol_T_bulk, vol_div_avg, rho_min, &
               bulk_target_Ub, min_delta_U_bulk, max_delta_U_bulk, bulk_target_Reb, Ly, Sf_x_mid, Sf_x_mid_prev, &
               min_Sf_x_mid, max_Sf_x_mid

logical  :: use_rough_surf, use_filtering, use_incompressible, use_utau_work, use_rough_surf_quady, &
            bulk_use_target_Ub, bulk_use_rho_weighting, bulk_use_target_Reb, force_Sf_x_pressure, force_Sf_x_accel,&
            init_Sf_x_mid_prev, force_Sf_x_semilocal
integer  :: filtering_steps, filtering_first, filtering_nstop, prof_1d_avg_N

real(dp) , allocatable ::  U                    (:,:,:) ,&
                           W                    (:,:,:) ,&
                           V                    (:,:,:) ,&
                           T                    (:,:,:) ,&
                           T_ext                (:,:,:) ,&
                           RHS_energy           (:,:,:) ,&
                           U_hat                (:,:,:) ,&
                           W_hat                (:,:,:) ,&
                           V_hat                (:,:,:) ,&
                           Ru_old               (:,:,:) ,&
                           Rw_old               (:,:,:) ,&
                           Rv_old               (:,:,:) ,&
                           Rt_old               (:,:,:) ,&
                           Ruw_ddy_a            (:)     ,&
                           Ruw_ddy_b            (:)     ,&
                           Ruw_ddy_c            (:)     ,&
                           Ruw_d2dy2_a          (:)     ,&
                           Ruw_d2dy2_b          (:)     ,&
                           Ruw_d2dy2_c          (:)     ,&
                           Ruw_pad_d2dy2_a      (:)     ,&
                           Ruw_pad_d2dy2_b      (:)     ,&
                           Ruw_pad_d2dy2_c      (:)     ,&
                           inv_dy_p             (:)     ,&
                           inv_dy_u             (:)     ,&
                           inv_dy_p_pad         (:)     ,&
                           dPdy_ddy_lo          (:)     ,&
                           dPdy_ddy_hi          (:)     ,&
                           Rv_interp_uw_bottom  (:)     ,&
                           Rv_interp_uw_top     (:)     ,&
                           Rv_ddy_a             (:)     ,&
                           Rv_ddy_b             (:)     ,&
                           Rv_ddy_c             (:)     ,&
                           Ruw_pad_ddy_a        (:)     ,&
                           Ruw_pad_ddy_b        (:)     ,&
                           Ruw_pad_ddy_c        (:)     ,&
                           Rv_d2dy2_a           (:)     ,&
                           Rv_d2dy2_b           (:)     ,&
                           Rv_d2dy2_c           (:)     ,&
                           buffer_1d            (:)     ,&
                           dy_cell              (:)     ,&
                           buffer_io            (:)     ,&
                           P_per                (:,:,:) ,&
                           div_now              (:,:,:) ,&
                           P_now                (:,:,:) ,&
                           P_old                (:,:,:) ,&
                           mu_arr               (:,:,:) ,&
                           cond_arr             (:,:,:) ,&
                           rho_arr              (:,:,:) ,&
                           enth_arr             (:,:,:) ,&
                           buoyancy_arr         (:,:,:) ,&
                           cp_arr               (:,:,:) ,&
                           div_jacobian_arr     (:,:,:) ,&
                           utau_old             (:,:,:) ,&
                           utau_now             (:,:,:) ,&
                           rho_arr_old          (:,:,:) ,&
                           inv_rho_dt_scale_Sf_x(:,:,:) ,&
                           temp_arr_n_xzy       (:,:,:) ,&
                           rho_profile_y        (:)     ,&
                           prof_1d_temp         (:)     ,&
                           utils_red_scalar     (:)
integer              ::  nstep, iters_full

! bulk parameters
integer     ::  bulk_last_iter, bulk_nprint
real(dp)    ::  bulk_time_start_print, time_start_total, bulk_time_finish_print

!__[PYTHON_HOLDER_PROF1D_DECLARE]
! real(dp), allocatable :: {x}(:)

! avg parameters
!__[PYTHON_HOLDER_AVG_DECLARE]
! real(dp), allocatable :: {x}(:,:,:)
! real(dp), allocatable :: binary_{x}(:,:,:,:)

 real(dp)  ::  avg_inv_nraw
 integer   ::  avg_N_now, avg_iter_start, avg_nraw, avg_freq, avg_bin_nbins

 integer , allocatable  ::  binary_cycle(:)

  integer  :: snap_iter_start, snap_freq

  real(dp) :: wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
              wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
              extrap_T_wall_0 , extrap_T_wall_1, &
              extrap_T_wall_m1, extrap_T_wall_m2

type(diezdecomp_parsed_mpi_ranks) ::  obj_ranks_00k  , obj_ranks_0zy
integer                           ::    mpi_pos_00k  , ny_loc_00k
integer, allocatable              ::  lo_00k_blocks_mpi(:), n00k_blocks_mpi(:)

type(diezdecomp_props_transp)     :: tr_red1d_fwd , tr_red1d_bwd, tr_io_fwd, tr_io_bwd
real(dp), allocatable             :: tr_red1d_B(:,:,:), tr_red1d_B_1d(:), work_tr(:), slab_00k(:,:,:)

type(diezdecomp_props_halo)  ::  hl(0:2)
