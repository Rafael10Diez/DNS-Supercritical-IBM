module mod_cfd_subroutines
  use  openacc
  use  mpi
  use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr, c_loc
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  use  mod_utils
  use  mod_io
  implicit none
  contains

  subroutine build_vars(wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
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
    use diezdecomp_api_generic, only: diezdecomp_props_transp, diezdecomp_props_halo, diezdecomp_halos_execute_generic
    use diezDecomp_api_ibm    , only: diezDecomp_ibm_type, diezdecomp_ibm_exec
    use mod_cfd_ibm
    implicit none
    real(dp)              ::  wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                              wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                              extrap_T_wall_0, extrap_T_wall_1, extrap_T_wall_m1, extrap_T_wall_m2, &
                              nu, cond, inv_dx, inv_dx2, inv_dx__2, inv_dz, inv_dz2, inv_dz__2, &
                              Gb_x, Gb_y, Gb_z, rho_min, Sf_x_mid, &
                              Sf_x, Sf_z, Sq, dtmax, inv_dtmax
    real(dp), contiguous  ::  buffer_ibm(0:), work_tr(0:), &
                              Ruw_pad_ddy_a(0:), Ruw_pad_ddy_b(0:), Ruw_pad_ddy_c(0:), &
                              Ruw_ddy_a(0:)  , Ruw_ddy_b(0:)  , Ruw_ddy_c(0:), &
                              Ruw_d2dy2_a(0:), Ruw_d2dy2_b(0:), Ruw_d2dy2_c(0:), &
                              inv_dy_p(0:), inv_dy_p_pad(0:), inv_dy_u(0:), &
                              utau_now             ( 0:, 0:, 0:), &
                              utau_old             ( 0:, 0:, 0:), &
                              temp_arr_n_xzy       ( 0:, 0:, 0:), &
                              Ru_old               ( 0:, 0:, 0:), &
                              Rw_old               ( 0:, 0:, 0:), &
                              Rv_old               ( 0:, 0:, 0:), &
                              Rt_old               ( 0:, 0:, 0:), &
                              P                    ( 0:, 0:, 0:), &
                              inv_rho_dt_scale_Sf_x( 0:, 0:, 0:), &
                              mu_arr               (-1:,-1:,-1:), &
                              cond_arr             (-1:,-1:,-1:), &
                              buoyancy_arr         (-1:,-1:,-1:), &
                              rho_arr              (-1:,-1:,-1:), &
                              rho_arr_old          (-1:,-1:,-1:), &
                              enth_arr             (-1:,-1:,-1:), &
                              cp_arr               (-1:,-1:,-1:), &
                              div_jacobian_arr     (-1:,-1:,-1:), &
                              U                    (-1:,-1:,-1:), &
                              V                    (-1:,-1:,-1:), &
                              W                    (-1:,-1:,-1:), &
                              T                    (-1:,-1:,-1:), &
                              T_ext                (-1:,-1:,-1:), &
                              U_hat                (-1:,-1:,-1:), &
                              V_hat                (-1:,-1:,-1:), &
                              W_hat                (-1:,-1:,-1:), &
                              div_now              (-1:,-1:,-1:), &
                              P_now                (-1:,-1:,-1:), &
                              P_old                (-1:,-1:,-1:), &
                              RHS_energy           ( 0:, 0:, 0:), &
                              slab_00k             ( 0:, 0:, 0:), &
                              rho_profile_y(0:), &
                              Rv_interp_uw_bottom(0:), Rv_interp_uw_top(0:), &
                              Rv_ddy_a(0:)  , Rv_ddy_b(0:)  , Rv_ddy_c(0:), &
                              Rv_d2dy2_a(0:), Rv_d2dy2_b(0:), Rv_d2dy2_c(0:), &
                              tr_red1d_B(0:,0:,0:), tr_red1d_B_1d(0:)
    real(dp), contiguous :: utils_red_scalar(0:)
    integer                        ::  mpi_pos_z, mpi_pos_y, mpi_divs_y, nx_glob, nz_loc, ny_loc, nz_glob, kmin_vhat, &
                                       ny_glob, ny_loc_00k, irank_mpi, nproc_mpi, mpi_divs_z, &
                                       lo_z_blocks_mpi(0:), lo_y_blocks_mpi(0:), iters_full, &
                        rough_U_bot_inds(0:,0:), rough_U_top_inds(0:,0:), rough_W_bot_inds(0:,0:), rough_W_top_inds(0:,0:), &
                        rough_V_bot_inds(0:,0:), rough_V_top_inds(0:,0:), rough_T_bot_inds(0:,0:), rough_T_top_inds(0:,0:)
    logical                        ::  use_rough_surf, use_incompressible, use_utau_work, &
                                       force_Sf_x_semilocal, force_Sf_x_pressure, use_rough_surf_quady
    type(diezDecomp_ibm_type)      ::  obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div
    type(diezdecomp_props_transp)  ::  tr_red1d_fwd , tr_red1d_bwd, tr_io_fwd
    type(diezdecomp_props_halo)    ::  hl(0:2)

    if (use_rough_surf) then
      if (use_rough_surf_quady) then 
        call  ibm_part_2(U, V, W, T, &
                         rough_U_bot_inds, rough_U_top_inds, rough_W_bot_inds, rough_W_top_inds , &
                         rough_V_bot_inds, rough_V_top_inds, rough_T_bot_inds, rough_T_top_inds , &
                         wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                         wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                         nx_glob      , nz_loc       , ny_loc                      )
        call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                             wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                             wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
      end if 
      call ibm_part_1(obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div, U, V, W, T, div_now, buffer_ibm, &
                      wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                      wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                      mpi_pos_y    , mpi_divs_y   , use_incompressible)
    else 
      call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                           wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                           wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
    end if

    ! ------------------------- begin: div_now walls -------------------------
    if (.not.use_incompressible) then
      block
        integer :: i,j
        if (mpi_pos_y == 0) then
          !$acc parallel loop  collapse(2) default(present)
          do   j = -1, nz_loc
            do i = -1, nx_glob
              div_now(i,j,-1) = extrap_T_wall_0 * div_now(i,j,0) + &
                                extrap_T_wall_1 * div_now(i,j,1)
            end do
          end do
        end if
        if (mpi_pos_y == (mpi_divs_y-1)) then
          !$acc parallel loop  collapse(2) default(present)
          do   j = -1, nz_loc
            do i = -1, nx_glob
              div_now(i,j,ny_loc) = extrap_T_wall_m1 * div_now(i,j,ny_loc-1) + &
                                    extrap_T_wall_m2 * div_now(i,j,ny_loc-2)
            end do
          end do
        end if
      end block
    end if
      ! ------------------------- end: div_now walls -------------------------

    block
      integer :: i,j,k

      ! --------------------- build_U_hat ---------------------
      block
        real(dp) :: Ru_now_ijk
        if (use_incompressible) then
          !$acc parallel loop  collapse(3) default(present) private(Ru_now_ijk)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Ru_now_ijk  =  (-(U(i,j,k))*((U(i+1,j,k)-U(i-1,j,k))*inv_dx__2)&
                -(0.25*(V(i,j,k)+V(i-1,j,k)+V(i,j,k+1)+&
                V(i-1,j,k+1)))*(U(i,j,k-1)*Ruw_ddy_a(k)+U(i,j,k)*&
                Ruw_ddy_b(k)+U(i,j,k+1)*Ruw_ddy_c(k))&
                -(0.25*(W(i,j,k)+W(i-1,j,k)+W(i,j+1,k)+&
                W(i-1,j+1,k)))*((U(i,j+1,k)-U(i,j-1,k))*inv_dz__2)+&
                nu*((U(i+1,j,k)-2*U(i,j,k)+U(i-1,j,k))*inv_dx2+&
                U(i,j,k-1)*Ruw_d2dy2_a(k)+U(i,j,k)*&
                Ruw_d2dy2_b(k)+U(i,j,k+1)*Ruw_d2dy2_c(k)+(U(i,j+1,k)&
                -2*U(i,j,k)+U(i,j-1,k))*inv_dz2)+Sf_x)
                U_hat(i,j,k)  =  U(i,j,k) + dtmax * (1.5*Ru_now_ijk - 0.5*Ru_old(i,j,k))
                Ru_old(i,j,k) =  Ru_now_ijk
              end do
            end do
          end do
        else
          !$acc parallel loop  collapse(3) default(present) private(Ru_now_ijk)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Ru_now_ijk  =  (-((U(i,j,k))*(((U(i+1,j,k)-U(i-1,j,k))*inv_dx__2))+(0.25*(W(i,j,k)+&
                               W(i-1,j,k)+W(i,j+1,k)+W(i-1,j+1,k)))*(((U(i,j+1,k)-&
                               U(i,j-1,k))*inv_dz__2))+(0.25*(V(i,j,k)+V(i-1,j,k)+V(i,j,k+1)+&
                               V(i-1,j,k+1)))*((U(i,j,k-1)*Ruw_ddy_a(k)+U(i,j,k)*Ruw_ddy_b(k)+&
                               U(i,j,k+1)*Ruw_ddy_c(k))))+&
                               (((0.5*(((mu_arr(i-1,j,k)))+((mu_arr(i,j,k))))))*(((U(i+1,j,k)-2*U(i,j,k)+&
                               U(i-1,j,k))*inv_dx2)+((U(i,j+1,k)-2*U(i,j,k)+U(i,j-1,k))*inv_dz2)+(&
                               U(i,j,k-1)*Ruw_d2dy2_a(k)+U(i,j,k)*Ruw_d2dy2_b(k)+&
                               U(i,j,k+1)*Ruw_d2dy2_c(k))+0.3333333333333333_dp*((div_now(i,j,k)-div_now(i-1,j,k))*inv_dx))&
                               +((2*(((((U(i+1,j,k)))-&
                               ((U(i-1,j,k))))*inv_dx__2))-0.6666666666666666_dp*(0.5*(((div_now(i-1,j,k)))+&
                               ((div_now(i,j,k))))))*(((((mu_arr(i,j,k)))&
                               -((mu_arr(i-1,j,k))))*inv_dx))+(((((U(i,j,k-1)))*Ruw_ddy_a(k)+&
                               ((U(i,j,k)))*Ruw_ddy_b(k)+((U(i,j,k+1)))*Ruw_ddy_c(k))+(((0.5*(V(i,j,k)&
                               +V(i,j,k+1)))-(0.5*(V(i-1,j,k)+&
                               V(i-1,j,k+1))))*inv_dx)))*(((0.5*(mu_arr(i-1,j,k-1)+mu_arr(i,j,k-1)))*Ruw_ddy_a(k)&
                               +(0.5*(mu_arr(i-1,j,k)+mu_arr(i,j,k)))*Ruw_ddy_b(k)+(0.5*(mu_arr(i-1,j,k+1)+&
                               mu_arr(i,j,k+1)))*Ruw_ddy_c(k)))+((((((U(i,j+1,k)))-&
                               ((U(i,j-1,k))))*inv_dz__2)+(((0.5*(W(i,j,k)+W(i,j+1,k)))-&
                               (0.5*(W(i-1,j,k)+W(i-1,j+1,k))))*inv_dx)))*((((0.5*(mu_arr(i-1,j+1,k)+&
                               mu_arr(i,j+1,k)))-(0.5*(mu_arr(i-1,j-1,k)+mu_arr(i,j-1,k))))*inv_dz__2)))+Sf_x+&
                               (Gb_x*(0.5*(buoyancy_arr(i,j,k)+&
                               buoyancy_arr(i-1,j,k)))))/((0.5*(((rho_arr(i-1,j,k)))+((rho_arr(i,j,k)))))))
                U_hat(i,j,k)  =  U(i,j,k) + dtmax * (1.5*Ru_now_ijk - 0.5*Ru_old(i,j,k))
                Ru_old(i,j,k) =  Ru_now_ijk
              end do
            end do
          end do
          if (force_Sf_x_pressure) then
            if (force_Sf_x_semilocal) then
              !$acc parallel loop  collapse(3) default(present)
              do     k = 0, ny_loc -1
                do   j = 0, nz_loc -1
                  do i = 0, nx_glob-1
                    temp_arr_n_xzy(i,j,k)  =  (1.5*((0.5*(((rho_arr(i-1,j,k)))     + &
                                                          ((rho_arr(i,j,k))))))    - &
                                               0.5*((0.5*(((rho_arr_old(i-1,j,k))) + &
                                                          ((rho_arr_old(i,j,k)))))))
                  end do
                end do
              end do
              call batched_binary_reduction(temp_arr_n_xzy, rho_profile_y, nx_glob*nz_loc, ny_loc, .true., nx_glob*nz_glob,&
              tr_red1d_fwd,tr_red1d_bwd,mpi_divs_z,ny_loc_00k,tr_red1d_B,tr_red1d_B_1d,work_tr)
            end if
            if (abs(Sf_x)>1e-20) error stop 'abs(Sf_x)>1e-20'
            !$acc parallel loop  collapse(3) default(present)
            do     k = 0, ny_loc -1
              do   j = 0, nz_loc -1
                do i = 0, nx_glob-1
                  inv_rho_dt_scale_Sf_x(i,j,k)  =  dtmax * (1.5/((0.5*(((rho_arr(i-1,j,k)))     + &
                                                                       ((rho_arr(i,j,k))))))    - &
                                                            0.5/((0.5*(((rho_arr_old(i-1,j,k))) + &
                                                                       ((rho_arr_old(i,j,k)))))))*rho_profile_y(k)

                  U_hat(i,j,k) = U_hat(i,j,k) + inv_rho_dt_scale_Sf_x(i,j,k)*Sf_x_mid
                end do
              end do
            end do
          end if
        end if
      end block

      ! --------------------- build_V_hat ---------------------
      block
        real(dp) :: Rv_now_ijk
        if (use_incompressible) then
          !$acc parallel loop  collapse(3) default(present) private(Rv_now_ijk)
          do     k = kmin_Vhat, ny_loc-1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Rv_now_ijk = (-(0.5*(U(i,j,k-1)+U(i+1,j,k-1))*Rv_interp_uw_bottom(k)&
                +0.5*(U(i,j,k)+&
                U(i+1,j,k))*Rv_interp_uw_top(k))*((V(i+1,j,k)-&
                V(i-1,j,k))*inv_dx__2)-(V(i,j,k))*(V(i,j,k-1)*&
                Rv_ddy_a(k)+V(i,j,k)*Rv_ddy_b(k)+V(i,j,k+1)*&
                Rv_ddy_c(k))-(0.5*(W(i,j,k-1)+&
                W(i,j+1,k-1))*Rv_interp_uw_bottom(k)+0.5*(W(i,j,k)+&
                W(i,j+1,k))*Rv_interp_uw_top(k))*((V(i,j+1,k)-&
                V(i,j-1,k))*inv_dz__2)+nu*((V(i+1,j,k)-2*V(i,j,k)+&
                V(i-1,j,k))*inv_dx2+V(i,j,k-1)*Rv_d2dy2_a(k)&
                +V(i,j,k)*Rv_d2dy2_b(k)+V(i,j,k+1)*Rv_d2dy2_c(k)&
                +(V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))*inv_dz2))
                V_hat(i,j,k)  =  V(i,j,k) + dtmax * (1.5*Rv_now_ijk - 0.5*Rv_old(i,j,k))
                Rv_old(i,j,k) =  Rv_now_ijk
              end do
            end do
          end do
        else
          !$acc parallel loop  collapse(3) default(present) private(Rv_now_ijk)
          do     k = kmin_Vhat, ny_loc-1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Rv_now_ijk = (-((0.5*(U(i,j,k-1)+U(i+1,j,k-1))*Rv_interp_uw_bottom(k)+0.5*(U(i,j,k)+&
                             U(i+1,j,k))*Rv_interp_uw_top(k))*(((V(i+1,j,k)-V(i-1,j,k))*inv_dx__2))+&
                             (0.5*(W(i,j,k-1)+W(i,j+1,k-1))*Rv_interp_uw_bottom(k)+0.5*(W(i,j,k)+&
                             W(i,j+1,k))*Rv_interp_uw_top(k))*(((V(i,j+1,k)-V(i,j-1,k))*inv_dz__2))+&
                             (V(i,j,k))*((V(i,j,k-1)*Rv_ddy_a(k)+V(i,j,k)*Rv_ddy_b(k)+V(i,j,k+1)*&
                             Rv_ddy_c(k))))+(((0.5*(((mu_arr(i,j,k-1)+mu_arr(i,j,k))))))*(((V(i+1,j,k)&
                             -2*V(i,j,k)+V(i-1,j,k))*inv_dx2)+((V(i,j+1,k)-2*V(i,j,k)+&
                             V(i,j-1,k))*inv_dz2)+(V(i,j,k-1)*Rv_d2dy2_a(k)+V(i,j,k)*Rv_d2dy2_b(k)&
                             +&
                             V(i,j,k+1)*Rv_d2dy2_c(k))+0.3333333333333333_dp*((div_now(i,j,k)-div_now(i,j,k-1))*inv_dy_u(k)))&
                             +(((((((V(i+1,j,k)))-((V(i-1,j,k))))*inv_dx__2)+(((0.5*(U(i,j,k)+&
                             U(i+1,j,k)))-(0.5*(U(i,j,k-1)+&
                             U(i+1,j,k-1))))*inv_dy_u(k))))*((((0.5*(mu_arr(i+1,j,k-1)+mu_arr(i+1,j,k)))-&
                             (0.5*(mu_arr(i-1,j,k-1)+mu_arr(i-1,j,k))))*inv_dx__2))+&
                             (2*((((V(i,j,k-1)))*Rv_ddy_a(k)+((V(i,j,k)))*Rv_ddy_b(k)+&
                             ((V(i,j,k+1)))*Rv_ddy_c(k)))-0.6666666666666666_dp*(0.5*(((div_now(i,j,k-1)&
                             +div_now(i,j,k))))))*(((((mu_arr(i,j,k)))-((mu_arr(i,j,k-1))))*inv_dy_u(k)))+&
                             ((((((V(i,j+1,k)))-((V(i,j-1,k))))*inv_dz__2)+(((0.5*(W(i,j,k)+&
                             W(i,j+1,k)))-(0.5*(W(i,j,k-1)+&
                             W(i,j+1,k-1))))*inv_dy_u(k))))*((((0.5*(mu_arr(i,j+1,k-1)+mu_arr(i,j+1,k)))-&
                             (0.5*(mu_arr(i,j-1,k-1)+mu_arr(i,j-1,k))))*inv_dz__2)))+&
                             (Gb_y*(0.5*(buoyancy_arr(i,j,k)+buoyancy_arr(i,j,k-1)))))/((0.5*(((rho_arr(i,j,k-1)+&
                             rho_arr(i,j,k)))))))
                V_hat(i,j,k)  =  V(i,j,k) + dtmax * (1.5*Rv_now_ijk - 0.5*Rv_old(i,j,k))
                Rv_old(i,j,k) =  Rv_now_ijk
              end do
            end do
          end do
        end if
      end block

      ! --------------------- build_W_hat ---------------------
      block
        real(dp) :: Rw_now_ijk
        if (use_incompressible) then
          !$acc parallel loop  collapse(3) default(present) private(Rw_now_ijk)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Rw_now_ijk = (-(0.25*(U(i,j,k)+U(i,j-1,k)+U(i+1,j,k)+&
                U(i+1,j-1,k)))*((W(i+1,j,k)-W(i-1,j,k))*inv_dx__2)-&
                (0.25*(V(i,j,k)+V(i,j-1,k)+V(i,j,k+1)+&
                V(i,j-1,k+1)))*(W(i,j,k-1)*Ruw_ddy_a(k)+W(i,j,k)*&
                Ruw_ddy_b(k)+W(i,j,k+1)*Ruw_ddy_c(k))-&
                (W(i,j,k))*((W(i,j+1,k)-W(i,j-1,k))*inv_dz__2)+&
                nu*((W(i+1,j,k)-2*W(i,j,k)+W(i-1,j,k))*inv_dx2+&
                W(i,j,k-1)*Ruw_d2dy2_a(k)+W(i,j,k)*&
                Ruw_d2dy2_b(k)+W(i,j,k+1)*Ruw_d2dy2_c(k)+&
                (W(i,j+1,k)-2*W(i,j,k)+W(i,j-1,k))*inv_dz2)+Sf_z)
                W_hat(i,j,k)  =  W(i,j,k) + dtmax * (1.5*Rw_now_ijk - 0.5*Rw_old(i,j,k))
                Rw_old(i,j,k) =  Rw_now_ijk
              end do
            end do
          end do
        else
          !$acc parallel loop  collapse(3) default(present) private(Rw_now_ijk)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Rw_now_ijk = (-((0.25*(U(i,j,k)+U(i,j-1,k)+U(i+1,j,k)+U(i+1,j-1,k)))*(((W(i+1,j,k)-&
                             W(i-1,j,k))*inv_dx__2))+(W(i,j,k))*(((W(i,j+1,k)-W(i,j-1,k))*inv_dz__2))&
                             +(0.25*(V(i,j,k)+V(i,j-1,k)+V(i,j,k+1)+V(i,j-1,k+1)))*((W(i,j,k-1)*&
                             Ruw_ddy_a(k)+W(i,j,k)*Ruw_ddy_b(k)+W(i,j,k+1)*Ruw_ddy_c(k))))+&
                             (((0.5*(((mu_arr(i,j-1,k)+mu_arr(i,j,k))))))*(((W(i+1,j,k)-2*W(i,j,k)+&
                             W(i-1,j,k))*inv_dx2)+((W(i,j+1,k)-2*W(i,j,k)+W(i,j-1,k))*inv_dz2)+(&
                             W(i,j,k-1)*Ruw_d2dy2_a(k)+W(i,j,k)*Ruw_d2dy2_b(k)+&
                             W(i,j,k+1)*Ruw_d2dy2_c(k))+0.3333333333333333_dp*((div_now(i,j,k)-div_now(i,j-1,k))*inv_dz))&
                             +(((((((W(i+1,j,k)))-((W(i-1,j,k))))*inv_dx__2)+(((0.5*(U(i,j,k)+&
                             U(i+1,j,k)))-(0.5*(U(i,j-1,k)+&
                             U(i+1,j-1,k))))*inv_dz)))*((((0.5*(mu_arr(i+1,j-1,k)+mu_arr(i+1,j,k)))-&
                             (0.5*(mu_arr(i-1,j-1,k)+mu_arr(i-1,j,k))))*inv_dx__2))+&
                             (((((W(i,j,k-1)))*Ruw_ddy_a(k)+((W(i,j,k)))*Ruw_ddy_b(k)+&
                             ((W(i,j,k+1)))*Ruw_ddy_c(k))+(((0.5*(V(i,j,k)+V(i,j,k+1)))-&
                             (0.5*(V(i,j-1,k)+V(i,j-1,k+1))))*inv_dz)))*(((0.5*(mu_arr(i,j-1,k-1)+&
                             mu_arr(i,j,k-1)))*Ruw_ddy_a(k)+(0.5*(mu_arr(i,j-1,k)+mu_arr(i,j,k)))*Ruw_ddy_b(k)+&
                             (0.5*(mu_arr(i,j-1,k+1)+mu_arr(i,j,k+1)))*Ruw_ddy_c(k)))+(2*(((((W(i,j+1,k)))-&
                             ((W(i,j-1,k))))*inv_dz__2))-0.6666666666666666_dp*(0.5*(((div_now(i,j-1,k)+&
                             div_now(i,j,k))))))*(((((mu_arr(i,j,k)))-((mu_arr(i,j-1,k))))*inv_dz)))+Sf_z+&
                             (Gb_z*(0.5*(buoyancy_arr(i,j,k)+buoyancy_arr(i,j-1,k)))))/((0.5*(((rho_arr(i,j-1,k)+&
                             rho_arr(i,j,k)))))))
                W_hat(i,j,k)  =  W(i,j,k) + dtmax * (1.5*Rw_now_ijk - 0.5*Rw_old(i,j,k))
                Rw_old(i,j,k) =  Rw_now_ijk
              end do
            end do
          end do
        end if
      end block
      if ((.not.use_incompressible).and.(use_utau_work)) then
        call build_utau_work(U, V, W, Ruw_pad_ddy_a, Ruw_pad_ddy_b, Ruw_pad_ddy_c, &
                             mu_arr, rho_arr, T_ext, utau_now, &
                             inv_dx, inv_dz, inv_dy_p, &
                             wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, &
                             wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, &
                             inv_dx__2, inv_dz__2, work_tr, buffer_ibm, &
                             nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, &
                             hl, obj_ibm_div, use_incompressible, use_utau_work, use_rough_surf)
      end if
      associate(Rt_now => temp_arr_n_xzy) ! Rt_now(i,j,k) is now replaced by P, only to save space
        ! --------------------- Advance_T ---------------------
        if (use_incompressible) then
          !$acc parallel loop  collapse(3) default(present)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Rt_now(i,j,k) = (-(0.5*(U(i,j,k)+U(i+1,j,k)))*((T(i+1,j,k)-&
                T(i-1,j,k))*inv_dx__2)-(0.5*(V(i,j,k)+&
                V(i,j,k+1)))*(T(i,j,k-1)*Ruw_ddy_a(k)+T(i,j,k)*&
                Ruw_ddy_b(k)+T(i,j,k+1)*Ruw_ddy_c(k))-&
                (0.5*(W(i,j,k)+W(i,j+1,k)))*((T(i,j+1,k)-&
                T(i,j-1,k))*inv_dz__2)+cond*((T(i+1,j,k)-2*T(i,j,k)&
                +T(i-1,j,k))*inv_dx2+T(i,j,k-1)*Ruw_d2dy2_a(k)&
                +T(i,j,k)*Ruw_d2dy2_b(k)+T(i,j,k+1)*&
                Ruw_d2dy2_c(k)+(T(i,j+1,k)-2*T(i,j,k)+T(i,j-1,k))&
                *inv_dz2)+Sq)
              end do
            end do
          end do
        else
          !$acc parallel loop  collapse(3) default(present)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                Rt_now(i,j,k) = (-(((0.5*(U(i,j,k)+U(i+1,j,k))))*(((T(i+1,j,k)-T(i-1,j,k))*inv_dx__2))+&
                                ((0.5*(W(i,j,k)+W(i,j+1,k))))*(((T(i,j+1,k)-T(i,j-1,k))*inv_dz__2))+&
                                ((0.5*(V(i,j,k)+V(i,j,k+1))))*(&
                                (T(i,j,k-1)*Ruw_ddy_a(k)+T(i,j,k)*Ruw_ddy_b(k)+T(i,j,k+1)*Ruw_ddy_c(k))&
                                ))+RHS_energy(i,j,k)/(rho_arr(i,j,k)*cp_arr(i,j,k)))
              end do
            end do
          end do
        end if
        !$acc parallel loop  collapse(3) default(present)
        do     k = 0, ny_loc -1
          do   j = 0, nz_loc -1
            do i = 0, nx_glob-1
              T(i,j,k)       =      T(i,j,k) + dtmax*( 1.5*Rt_now(i,j,k) - 0.5*Rt_old(i,j,k) )
              Rt_old(i,j,k)  =  Rt_now(i,j,k)
            end do
          end do
        end do
      end associate

        ! --------------------- update rough surface ---------------------
    if (use_rough_surf) then
      call  ibm_part_2(U_hat, V_hat, W_hat, T, &
                       rough_U_bot_inds, rough_U_top_inds, rough_W_bot_inds, rough_W_top_inds, &
                       rough_V_bot_inds, rough_V_top_inds, rough_T_bot_inds, rough_T_top_inds, &
                       wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                       wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                       nx_glob       , nz_loc       , ny_loc                     )
      if (use_rough_surf_quady) then 
        call cfd_periodicity(hl, U_hat, V_hat, W_hat, T, work_tr, mpi_pos_y, mpi_divs_y, &
                             wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                             wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
        call ibm_part_1(obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div, U_hat, V_hat, W_hat, T, div_now, buffer_ibm, &
                        wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                        wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                        mpi_pos_y    , mpi_divs_y   , .false.) ! disable div_now update
      end if
    else 
        call cfd_periodicity(hl, U_hat, V_hat, W_hat, T, work_tr, mpi_pos_y, mpi_divs_y, &
                             wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                             wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
    end if

        ! T periodicity
        block 
          integer :: ii_hl
          do ii_hl=0,2
            call diezdecomp_halos_execute_generic(hl(ii_hl), T, work_tr)
          end do 
        end block
        
        if (mpi_pos_y == 0             )  call set_val_slice(T, 0,  0, wall_BC_T_bot)
        if (mpi_pos_y == (mpi_divs_y-1))  call set_val_slice(T,-1, -1, wall_BC_T_top)

        if (.not.use_incompressible) then
          if (force_Sf_x_pressure) then
            !$acc parallel loop  collapse(3) default(present)
            do     k = -1, ny_loc
              do   j = -1, nz_loc
                do i = -1, nx_glob
                  rho_arr_old(i,j,k)  =  rho_arr(i,j,k)
                end do
              end do
            end do
          end if

          !$acc wait
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
          block
            integer :: i,j,k,mpi_ierr
            real(dp)  :: temp_real
            rho_min = 1e20
#if !defined(_BIN_REDUCTION)
            !$acc parallel loop collapse(3) default(present) reduction(min:rho_min)
            do      k=0, ny_loc -1
              do    j=0, nz_loc -1
                do  i=0, nx_glob-1
                    rho_min  =  min(rho_min, rho_arr(i,j,k))
                end do
              end do
            end do
#else
            call set_val_slice(temp_arr_n_xzy,0,-1, 9.9e20)
            !$acc parallel loop collapse(3) default(present)
            do      k=0, ny_loc -1
              do    j=0, nz_loc -1
                do  i=0, nx_glob-1
                    temp_arr_n_xzy(i,j,k)  =  ( rho_arr(i,j,k))
                end do
              end do
            end do
            call quick_binary_scalar_reduction(temp_arr_n_xzy, rho_min, nx_glob*nz_loc*ny_loc, -1, &
                                               utils_red_scalar)
#endif
            !$acc wait
            temp_real = rho_min
            call MPI_Allreduce(temp_real, rho_min, 1, MPI_DOUBLE_PRECISION, mpi_min, mpi_comm_world, mpi_ierr)
          end block
        end if

        ! --------------------- periodicity U_hat, W_hat, V_hat ---------------------

        call diezdecomp_halos_execute_generic(hl(0), U_hat, work_tr)
        call diezdecomp_halos_execute_generic(hl(1), W_hat, work_tr)
        call diezdecomp_halos_execute_generic(hl(2), V_hat, work_tr)

      if (mpi_pos_y == 0       ) then
        !$acc parallel loop  collapse(2) default(present)
        do   j = 0, nz_loc -1
          do i = 0, nx_glob-1
            V_hat(i, j, 0) = wall_BC_V_bot
          end do
        end do
      end if

      if (mpi_pos_y == (mpi_divs_y-1)) then
        !$acc parallel loop  collapse(2) default(present)
        do   j = 0, nz_loc -1
          do i = 0, nx_glob-1
            V_hat(i, j, ny_loc) = wall_BC_V_top
          end do
        end do
      end if

      ! --------------------- build_P_rhs ---------------------
      if (use_incompressible) then
        !$acc parallel loop  collapse(3) default(present)
        do     k = 0, ny_loc -1
          do   j = 0, nz_loc -1
            do i = 0, nx_glob-1
              P(i,j,k)  = inv_dtmax*((U_hat(i+1, j  , k  ) - U_hat(i,j,k))*inv_dx   + &
                                     (W_hat(i  , j+1, k  ) - W_hat(i,j,k))*inv_dz   + &
                                     (V_hat(i  , j  , k+1) - V_hat(i,j,k))*inv_dy_p(k))
            end do
          end do
        end do
      else
        block
          real(dp)  :: Ps_00,Ps_xp,Ps_xm,Ps_zp,Ps_zm,Ps_yp,Ps_ym, &
                       rho_xp,rho_xm,rho_zp,rho_zm,rho_yp,rho_ym
          !$acc parallel loop  collapse(3) default(present) &
          !$acc private(Ps_00,Ps_xp,Ps_xm,Ps_zp,Ps_zm,Ps_yp,Ps_ym) &
          !$acc private(rho_xp,rho_xm,rho_zp,rho_zm,rho_yp,rho_ym)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                  Ps_00   =  2*P_now(i  ,j  ,k  ) - P_old(i  ,j  ,k  )
                  Ps_xp   =  2*P_now(i+1,j  ,k  ) - P_old(i+1,j  ,k  )
                  Ps_xm   =  2*P_now(i-1,j  ,k  ) - P_old(i-1,j  ,k  )
                  Ps_zp   =  2*P_now(i  ,j+1,k  ) - P_old(i  ,j+1,k  )
                  Ps_zm   =  2*P_now(i  ,j-1,k  ) - P_old(i  ,j-1,k  )
                  Ps_yp   =  2*P_now(i  ,j  ,k+1) - P_old(i  ,j  ,k+1)
                  Ps_ym   =  2*P_now(i  ,j  ,k-1) - P_old(i  ,j  ,k-1)

                  rho_xp  =  1.0_dp - rho_min*(0.5_dp/rho_arr(i+1,j  ,k  ) + 0.5_dp/rho_arr(i,j,k))
                  rho_xm  =  1.0_dp - rho_min*(0.5_dp/rho_arr(i-1,j  ,k  ) + 0.5_dp/rho_arr(i,j,k))
                  rho_zp  =  1.0_dp - rho_min*(0.5_dp/rho_arr(i  ,j+1,k  ) + 0.5_dp/rho_arr(i,j,k))
                  rho_zm  =  1.0_dp - rho_min*(0.5_dp/rho_arr(i  ,j-1,k  ) + 0.5_dp/rho_arr(i,j,k))
                  rho_yp  =  1.0_dp - rho_min*(0.5_dp/rho_arr(i  ,j  ,k+1) + 0.5_dp/rho_arr(i,j,k))
                  rho_ym  =  1.0_dp - rho_min*(0.5_dp/rho_arr(i  ,j  ,k-1) + 0.5_dp/rho_arr(i,j,k))

                P(i,j,k)  = inv_dtmax*((U_hat(i+1, j  , k  ) - U_hat(i,j,k))*inv_dx   + &
                                       (W_hat(i  , j+1, k  ) - W_hat(i,j,k))*inv_dz   + &
                                       (V_hat(i  , j  , k+1) - V_hat(i,j,k))*inv_dy_p(k))*rho_min &
                          - inv_dtmax*rho_min*div_now(i,j,k) + &
                  ( (rho_xp*(Ps_xp - Ps_00)                   - rho_xm*(Ps_00 - Ps_xm)                    )*inv_dx2 + &
                    (rho_zp*(Ps_zp - Ps_00)                   - rho_zm*(Ps_00 - Ps_zm)                    )*inv_dz2 + &
                    (rho_yp*(Ps_yp - Ps_00)*inv_dy_p_pad(k+1) - rho_ym*(Ps_00 - Ps_ym) *inv_dy_p_pad(k))*inv_dy_p(k) )
              end do
            end do
          end do
        end block
      end if

#if defined(_EXTRA_OUTPUT)
      block
        character(len=36)  ::  my_file
        write(my_file,'(A23,I0.9,A4)') './ini/array_U_hat_ini0_',iters_full,'.dat'
        call cfd_write_arr(my_file, U_hat, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
      end block
      block
        character(len=36)  ::  my_file
        write(my_file,'(A23,I0.9,A4)') './ini/array_W_hat_ini0_',iters_full,'.dat'
        call cfd_write_arr(my_file, W_hat, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
      end block
      block
        character(len=36)  ::  my_file
        write(my_file,'(A23,I0.9,A4)') './ini/array_V_hat_ini0_',iters_full,'.dat'
        call cfd_write_arr(my_file, V_hat, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
      end block
      block
        character(len=36)  ::  my_file
        write(my_file,'(A23,I0.9,A4)') './ini/arr_rhs_hat_ini0_',iters_full,'.dat'
        call cfd_write_arr(my_file, P, 0, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
      end block
      block
        character(len=36)  ::  my_file
        write(my_file,'(A23,I0.9,A4)') './ini/array_T_hat_ini0_',iters_full,'.dat'
        call cfd_write_arr(my_file, T, -1, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
      end block
      block
        character(len=36)  ::  my_file
        write(my_file,'(A23,I0.9,A4)') './ini/arr_Rt_old_ini0__',iters_full,'.dat'
        call cfd_write_arr(my_file, Rt_old, 0, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                           irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
      end block
#endif

    end block
  end subroutine

  subroutine build_rhs_energy(T, cond_arr, div_now, div_jacobian_arr, RHS_energy, &
                              inv_dx2, inv_dz2, inv_dx__2, inv_dz__2, &
                              Ruw_d2dy2_a, Ruw_d2dy2_b, Ruw_d2dy2_c, Ruw_ddy_a, Ruw_ddy_b, Ruw_ddy_c, Sq, &
                              utau_now, utau_old, nx_glob, nz_loc, ny_loc, &
                              rough_T_bot_inds, rough_T_top_inds, &
                              hl, work_tr, buffer_ibm, mpi_pos_y, mpi_divs_y, wall_BC_T_bot, wall_BC_T_top, &
                              use_utau_work, use_incompressible, use_rough_surf, obj_ibm_T)
    use diezdecomp_api_generic, only: diezdecomp_props_halo, diezdecomp_halos_execute_generic
    use diezDecomp_api_ibm    , only: diezDecomp_ibm_type, diezdecomp_ibm_exec
    implicit none
    real(dp), contiguous ::  T(-1:,-1:,-1:), cond_arr(-1:,-1:,-1:), div_now(-1:,-1:,-1:), div_jacobian_arr(-1:,-1:,-1:), &
                             RHS_energy(0:,0:,0:), &
                             Ruw_d2dy2_a(0:), Ruw_d2dy2_b(0:), Ruw_d2dy2_c(0:), Ruw_ddy_a(0:), Ruw_ddy_b(0:), Ruw_ddy_c(0:), &
                             utau_now(0:,0:,0:), utau_old(0:,0:,0:), work_tr(0:), buffer_ibm(0:)
    real(dp)             ::  inv_dx2, inv_dz2, inv_dx__2, inv_dz__2, Sq, wall_BC_T_bot, wall_BC_T_top
    integer              ::  nx_glob, nz_loc, ny_loc, rough_T_bot_inds(0:,0:), rough_T_top_inds(0:,0:), mpi_pos_y, mpi_divs_y
    logical              ::  use_utau_work, use_incompressible, use_rough_surf
    type(diezdecomp_props_halo)  ::  hl(0:2)
    type(diezDecomp_ibm_type)    ::  obj_ibm_T

    if (use_incompressible) error stop 'ERROR: called build_rhs_energy.f90 with an incompressible flow'

    ! update T (rough surface)
    if (use_rough_surf) then
      call diezdecomp_ibm_exec(obj_ibm_T, T, buffer_ibm)
    end if
      if (mpi_pos_y ==  0            )  call set_val_slice(T, 0,  0, wall_BC_T_bot)
      if (mpi_pos_y == (mpi_divs_y-1))  call set_val_slice(T,-1, -1, wall_BC_T_top)

    ! get RHS_energy
      block
        integer   ::  i,j,k

        !$acc parallel loop  collapse(3) default(present)
        do     k = 0, ny_loc -1
          do   j = 0, nz_loc -1
            do i = 0, nx_glob-1
              RHS_energy(i,j,k)  =  ((cond_arr(i,j,k))*(((T(i+1,j,k)-2*T(i,j,k)+T(i-1,j,k))*inv_dx2)+&
                                    ((T(i,j+1,k)-2*T(i,j,k)+T(i,j-1,k))*inv_dz2)+(&
                                    T(i,j,k-1)*Ruw_d2dy2_a(k)+T(i,j,k)*Ruw_d2dy2_b(k)+&
                                    T(i,j,k+1)*Ruw_d2dy2_c(k)))+((((cond_arr(i+1,j,k)-&
                                    cond_arr(i-1,j,k))*inv_dx__2))*(((T(i+1,j,k)-T(i-1,j,k))*inv_dx__2))+&
                                    (((cond_arr(i,j+1,k)-cond_arr(i,j-1,k))*inv_dz__2))*(((T(i,j+1,k)-&
                                    T(i,j-1,k))*inv_dz__2))+((cond_arr(i,j,k-1)*Ruw_ddy_a(k)+cond_arr(i,j,k)*&
                                    Ruw_ddy_b(k)+cond_arr(i,j,k+1)*Ruw_ddy_c(k)))*((T(i,j,k-1)*Ruw_ddy_a(k)+&
                                    T(i,j,k)*Ruw_ddy_b(k)+T(i,j,k+1)*Ruw_ddy_c(k))))+Sq)
            end do
          end do
        end do

        if (use_utau_work) then
          !$acc parallel loop  collapse(3) default(present)
          do     k = 0, ny_loc -1
            do   j = 0, nz_loc -1
              do i = 0, nx_glob-1
                RHS_energy(i,j,k) =  RHS_energy(i,j,k) + 2*utau_now(i,j,k) - utau_old(i,j,k)
                utau_old(i,j,k)   =                        utau_now(i,j,k)
              end do
            end do
          end do
        end if
      end block


    ! get div_now
     block
       integer :: i,j,k
       !$acc parallel loop  collapse(3) default(present)
       do     k = 0, ny_loc -1
         do   j = 0, nz_loc -1
           do i = 0, nx_glob-1
             div_now(i,j,k) = div_jacobian_arr(i,j,k)*RHS_energy(i,j,k)
           end do
         end do
       end do
     end block

    if (use_rough_surf) then
      block
        integer :: i,j,k
        ! rough surf div_now
        !$acc parallel loop collapse(2) default(present)
        do     j  = 0, nz_loc -1
          do   i  = 0, nx_glob-1
            do k  = 0, rough_T_bot_inds(i,j), 1
              div_now(i,j,k) = 0.
            end do
          end do
        end do
        !$acc parallel loop collapse(2) default(present)
        do     j  = 0, nz_loc -1
          do   i  = 0, nx_glob-1
            do k  = rough_T_top_inds(i,j), ny_loc-1, 1
              div_now(i,j,k) = 0.
            end do
          end do
        end do
      end block
    end if

    block 
      integer :: ii_hl
      do ii_hl=0,2
        call diezdecomp_halos_execute_generic(hl(ii_hl), div_now, work_tr)
      end do 
    end block

    block
      integer :: i,j
     if (mpi_pos_y == 0) then
       !$acc parallel loop  collapse(2) default(present)
       do   j = -1, nz_loc
         do i = -1, nx_glob
           div_now(i,j,-1) = 0.
         end do
       end do
     end if

     if (mpi_pos_y == (mpi_divs_y -1)) then
       !$acc parallel loop  collapse(2) default(present)
       do   j = -1, nz_loc
         do i = -1, nx_glob
           div_now(i,j,ny_loc) = 0.
         end do
       end do
     end if
    end block
  end subroutine

  subroutine build_utau_work(U, V, W, Ruw_pad_ddy_a, Ruw_pad_ddy_b, Ruw_pad_ddy_c, &
                             mu_arr, rho_arr, T_ext, utau_now, &
                             inv_dx, inv_dz, inv_dy_p, &
                             wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, &
                             wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, &
                             inv_dx__2, inv_dz__2, work_tr, buffer_ibm, &
                             nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, &
                             hl, obj_ibm_div, use_incompressible, use_utau_work, use_rough_surf)
    use diezdecomp_api_generic, only: diezdecomp_props_halo, diezdecomp_halos_execute_generic
    use diezDecomp_api_ibm    , only: diezDecomp_ibm_type, diezdecomp_ibm_exec
    implicit none
    real(dp)                     ::  U(-1:,-1:,-1:), V(-1:,-1:,-1:), W(-1:,-1:,-1:), &
                                     Ruw_pad_ddy_a(0:), Ruw_pad_ddy_b(0:), Ruw_pad_ddy_c(0:), &
                                     mu_arr(-1:,-1:,-1:), rho_arr(-1:,-1:,-1:), T_ext(-1:,-1:,-1:), &
                                     utau_now(0:,0:,0:), inv_dx, inv_dz, inv_dy_p(0:), &
                                     wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, &
                                     wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, &
                                     inv_dx__2, inv_dz__2, work_tr(0:), buffer_ibm(0:)
    integer                      ::  nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y
    logical                      ::  use_incompressible, use_utau_work, use_rough_surf
    type(diezdecomp_props_halo)  ::  hl(0:2)
    type(diezDecomp_ibm_type)    ::  obj_ibm_div

    error stop 'ERROR: build_utau_work is a removed legacy subroutine.'

    if (use_incompressible) error stop 'ERROR: called build_utau_work with (use_incompressible)'
    if (.not.use_utau_work) error stop 'ERROR: called build_utau_work with (.NOT.use_utau_work)'
  end subroutine

  subroutine update_props(T_ext, T, work_tr, wall_BC_T_bot, wall_BC_T_top, &
                          mu_arr, cond_arr, rho_arr, buoyancy_arr, enth_arr, cp_arr , div_jacobian_arr, &
                          rough_T_bot_inds, rough_T_top_inds, nx_glob, nz_loc, ny_loc, &
                          lo_z_blocks_mpi, lo_y_blocks_mpi, mpi_pos_z, mpi_pos_y, &
                          use_incompressible, hl, use_rough_surf)
    use diezdecomp_api_generic, only: diezdecomp_props_halo, diezdecomp_halos_execute_generic
    implicit none
    real(dp)                     ::  wall_BC_T_bot, wall_BC_T_top
    real(dp), contiguous         ::  T_ext(-1:,-1:,-1:), T(-1:,-1:,-1:), work_tr(0:), &
                                     mu_arr          (-1:,-1:,-1:), cond_arr(-1:,-1:,-1:), rho_arr(-1:,-1:,-1:), &
                                     buoyancy_arr    (-1:,-1:,-1:), enth_arr(-1:,-1:,-1:), cp_arr (-1:,-1:,-1:), &
                                     div_jacobian_arr(-1:,-1:,-1:)
    integer                      ::  rough_T_bot_inds(0:,0:), rough_T_top_inds(0:,0:), nx_glob, nz_loc, ny_loc, &
                                     lo_z_blocks_mpi(0:), lo_y_blocks_mpi(0:), mpi_pos_z, mpi_pos_y
    logical                      ::  use_incompressible, use_rough_surf
    type(diezdecomp_props_halo)  ::  hl(0:2)
      if (use_incompressible) error stop 'ERROR: called update_props.f90 with an incompressible flow'

     block
      integer :: i,j,k
      !$acc parallel loop  collapse(3) default(present)
      do     k = -1, ny_loc
        do   j = -1, nz_loc
          do i = -1, nx_glob
            T_ext(i,j,k) = T(i,j,k)
          end do
        end do
      end do

      if (use_rough_surf) then
        ! ------------------------------ T ------------------------------
        !$acc parallel loop collapse(2) default(present)
        do      j  = 0, nz_loc-1
          do    i  = 0, nx_glob-1
            do k  = 0, rough_T_bot_inds(i,j), 1
              T_ext(i,j,k) = wall_BC_T_bot
            end do
          end do
        end do

        !$acc parallel loop collapse(2) default(present)
        do     j  = 0, nz_loc-1
          do   i  = 0, nx_glob-1
            do k  = rough_T_top_inds(i,j), ny_loc-1, 1
              T_ext(i,j,k) = wall_BC_T_top
            end do
          end do
        end do

        block 
          integer :: ii_hl
          do ii_hl=0,2
            call diezdecomp_halos_execute_generic(hl(ii_hl), T_ext, work_tr)
          end do 
        end block

      end if
     end block

    block
      integer  :: i,j,k,di_glob,dj_glob,dk_glob
      real(dp) :: T_ext_ijk, mu_ijk, cond_ijk, rho_ijk, cp_ijk, drho_dT_ijk, enth_ijk
      di_glob = 0
      dj_glob = lo_z_blocks_mpi(mpi_pos_z)
      dk_glob = lo_y_blocks_mpi(mpi_pos_y)
      !$acc parallel loop  collapse(3) private(T_ext_ijk,mu_ijk,cond_ijk,rho_ijk,cp_ijk,drho_dT_ijk,enth_ijk) default(present)
      do     k = -1, ny_loc
        do   j = -1, nz_loc
          do i = -1, nx_glob
           T_ext_ijk = T_ext(i,j,k)
           include "thermophysical_properties/expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90"
           ! enth_arr(i,j,k)         = ... T_ext(i,j,k) ...
           ! rho_arr(i,j,k)          = ... T_ext(i,j,k) ...
           ! mu_arr(i,j,k)           = ... T_ext(i,j,k) ...
           ! cond_arr(i,j,k)         = ... T_ext(i,j,k) ...
           ! buoyancy_arr(i,j,k)     = ... T_ext(i,j,k) ...
           ! cp_arr(i,j,k)           = ... T_ext(i,j,k) ...
           ! div_jacobian_arr(i,j,k) = ... T_ext(i,j,k) ...
          end do
        end do
      end do
    end block
  end subroutine

  subroutine cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                             wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                             wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
    use diezdecomp_api_generic, only: diezdecomp_props_halo, diezdecomp_halos_execute_generic
    implicit none
    type(diezdecomp_props_halo)  ::  hl(0:2)
    real(dp), contiguous  ::  work_tr(0:)   , &
                              U(-1:,-1:,-1:), &
                              V(-1:,-1:,-1:), &
                              W(-1:,-1:,-1:), &
                              T(-1:,-1:,-1:)
    real(dp)              ::  wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                              wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top
    integer               ::  mpi_pos_y, mpi_divs_y

    block 
      integer :: ii_hl
      do ii_hl=0,2
        call diezdecomp_halos_execute_generic(hl(ii_hl), U, work_tr)
        call diezdecomp_halos_execute_generic(hl(ii_hl), V, work_tr)
        call diezdecomp_halos_execute_generic(hl(ii_hl), W, work_tr)
        call diezdecomp_halos_execute_generic(hl(ii_hl), T, work_tr)
      end do 
    end block

   if (mpi_pos_y == 0) then
     call set_val_slice(U,0, 0, wall_BC_U_bot) ! set_val_slice is zero-indexed
     call set_val_slice(V,0, 1, wall_BC_V_bot)
     call set_val_slice(W,0, 0, wall_BC_W_bot)
     call set_val_slice(T,0, 0, wall_BC_T_bot)
   end if

   if (mpi_pos_y == (mpi_divs_y-1)) then
     call set_val_slice(U,-1, -1, wall_BC_U_top) ! set_val_slice is zero-indexed (negative indexes follow python-conventions)
     call set_val_slice(V,-1, -1, wall_BC_V_top)
     call set_val_slice(W,-1, -1, wall_BC_W_top)
     call set_val_slice(T,-1, -1, wall_BC_T_top)
   end if

  end subroutine

  subroutine filtering_UVWT(iters_full, nx_glob, nz_loc, ny_loc, kmin_Vhat, irank_mpi, temp_arr_n_xzy, U, V, W, T)
    integer  :: iters_full, nx_glob, nz_loc, ny_loc, kmin_Vhat, irank_mpi
    real(dp) :: temp_arr_n_xzy(0:,0:,0:), U(-1:,-1:,-1:), V(-1:,-1:,-1:), W(-1:,-1:,-1:), T(-1:,-1:,-1:)

    if (irank_mpi == 0) then
      write(6,'(A60,I10)') 'INFO: Filtering temperature and velocity fields, iteration: ', iters_full;flush(6)
    end if

    block
      integer :: i,j,k
      real(dp)  :: div_7
      div_7  =  1.d0/7.d0
        ! ------------------ FILTER U ------------------
        !$acc parallel loop  collapse(3) default(present)
        do         k = 0, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    temp_arr_n_xzy(i, j, k) = (U(i  , j  , k  ) + &
                                  U(i+1, j  , k  ) + &
                                  U(i-1, j  , k  ) + &
                                  U(i  , j+1, k  ) + &
                                  U(i  , j-1, k  ) + &
                                  U(i  , j  , k+1) + &
                                  U(i  , j  , k-1) )*div_7
                end do
            end do
        end do
        !$acc parallel loop  collapse(3) default(present)
        do         k = 0, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    U(i, j, k) = temp_arr_n_xzy(i, j, k)
                end do
            end do
        end do

        ! ------------------ FILTER V ------------------
        !$acc parallel loop  collapse(3) default(present)
        do         k = kmin_Vhat, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    temp_arr_n_xzy(i, j, k) = (V(i  , j  , k  ) + &
                                  V(i+1, j  , k  ) + &
                                  V(i-1, j  , k  ) + &
                                  V(i  , j+1, k  ) + &
                                  V(i  , j-1, k  ) + &
                                  V(i  , j  , k+1) + &
                                  V(i  , j  , k-1) )*div_7
                end do
            end do
        end do
        !$acc parallel loop  collapse(3) default(present)
        do         k = kmin_Vhat, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    V(i, j, k) = temp_arr_n_xzy(i, j, k)
                end do
            end do
        end do

        ! ------------------ FILTER W ------------------
        !$acc parallel loop  collapse(3) default(present)
        do         k = 0, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    temp_arr_n_xzy(i, j, k) = (W(i  , j  , k  ) + &
                                  W(i+1, j  , k  ) + &
                                  W(i-1, j  , k  ) + &
                                  W(i  , j+1, k  ) + &
                                  W(i  , j-1, k  ) + &
                                  W(i  , j  , k+1) + &
                                  W(i  , j  , k-1) )*div_7
                end do
            end do
        end do
        !$acc parallel loop  collapse(3) default(present)
        do         k = 0, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    W(i, j, k) = temp_arr_n_xzy(i, j, k)
                end do
            end do
        end do

        ! ------------------ FILTER T ------------------
        !$acc parallel loop  collapse(3) default(present)
        do         k = 0, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    temp_arr_n_xzy(i, j, k) = (T(i  , j  , k  ) + &
                                  T(i+1, j  , k  ) + &
                                  T(i-1, j  , k  ) + &
                                  T(i  , j+1, k  ) + &
                                  T(i  , j-1, k  ) + &
                                  T(i  , j  , k+1) + &
                                  T(i  , j  , k-1) )*div_7
                end do
            end do
        end do
        !$acc parallel loop  collapse(3) default(present)
        do         k = 0, ny_loc-1
            do     j = 0, nz_loc-1
                do i = 0, nx_glob-1
                    T(i, j, k) = temp_arr_n_xzy(i, j, k)
                end do
            end do
        end do
    end block
  end subroutine

  subroutine periodicity_pressure(hl, P_now, work_tr, nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, &
                                  use_incompressible)
    use diezdecomp_api_generic, only   :  diezdecomp_props_halo, diezdecomp_halos_execute_generic
    implicit none
    type(diezdecomp_props_halo)       ::  hl(0:2)
    real(dp), contiguous              ::  P_now(-1:,-1:,-1:), work_tr(0:)
    integer                           ::  nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y
    logical                           ::  use_incompressible
    if (use_incompressible) error stop 'ERROR: called periodicity_pressure.f90 with an incompressible flow'
    block
      integer :: i,j
      if (mpi_pos_y == 0) then
       !$acc parallel loop  collapse(2) default(present)
       do   j = -1, nz_loc
         do i = -1, nx_glob
           P_now(i,j,-1) = P_now(i,j,0)
         end do
       end do
      end if

      if (mpi_pos_y == (mpi_divs_y-1)) then
       !$acc parallel loop  collapse(2) default(present)
       do   j = -1, nz_loc
         do i = -1, nx_glob
           P_now(i,j,ny_loc) = P_now(i,j,ny_loc-1)
         end do
       end do
      end if
    end block
    block 
      integer :: ii_hl
      do ii_hl=0,2
        call diezdecomp_halos_execute_generic(hl(ii_hl), P_now, work_tr)
      end do 
    end block
  end subroutine

  subroutine update_vars(nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, kmin_Vhat, &
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
    use diezdecomp_api_generic, only: diezdecomp_props_halo, diezdecomp_halos_execute_generic
    implicit none
    integer                      ::  nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, kmin_Vhat, &
                                     bulk_U_bot_inds(0:,0:), bulk_U_top_inds(0:,0:)
    real(dp), contiguous         ::  P(0:,0:,0:), P_per(-1:,-1:,-1:), P_now(-1:,-1:,-1:), P_old(-1:,-1:,-1:), work_tr(0:), &
                                     U(-1:,-1:,-1:)    , V(-1:,-1:,-1:)    , W(-1:,-1:,-1:), &
                                     U_hat(-1:,-1:,-1:), V_hat(-1:,-1:,-1:), W_hat(-1:,-1:,-1:), &
                                     rho_arr(-1:,-1:,-1:), mu_arr(-1:,-1:,-1:), &
                                     dPdy_ddy_hi(0:), dPdy_ddy_lo(0:), inv_rho_dt_scale_Sf_x(0:,0:,0:), &
                                     dy_cell(0:), temp_arr_n_xzy(0:,0:,0:)
    real(dp), contiguous         ::  utils_red_scalar(0:)
    real(dp)                     ::  dtmax, inv_dx, inv_dz, &
                                     wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, &
                                     wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, &
                                     bulk_target_Ub, vol_U_bulk, &
                                     min_delta_U_bulk, max_delta_U_bulk, Ly, bulk_target_Reb, rho_min, &
                                     Sf_x_mid, Sf_x_mid_prev, min_Sf_x_mid, max_Sf_x_mid
    logical                      ::  use_incompressible, &
                                     bulk_use_rho_weighting, force_Sf_x_accel, force_Sf_x_pressure, bulk_use_target_Reb, &
                                     bulk_use_target_Ub, use_rough_surf, init_Sf_x_mid_prev
    type(diezdecomp_props_halo)  ::  hl(0:2)

    block
      integer :: i,j,k
      !$acc parallel loop  collapse(3) default(present)
      do     k = 0, ny_loc -1
        do   j = 0, nz_loc -1
          do i = 0, nx_glob-1
            P_per(i, j, k) = P(i,j,k)
          end do
        end do
      end do
      block 
        integer :: ii_hl
        do ii_hl=0,2
          call diezdecomp_halos_execute_generic(hl(ii_hl), P_per, work_tr)
        end do 
      end block

      block
        integer :: i,j
        if (mpi_pos_y == 0) then
         !$acc parallel loop  collapse(2) default(present)
         do   j = -1, nz_loc
           do i = -1, nx_glob
             P_per(i,j,-1) = P_per(i,j,0)
           end do
         end do
        end if

        if (mpi_pos_y == (mpi_divs_y-1)) then
         !$acc parallel loop  collapse(2) default(present)
         do   j = -1, nz_loc
           do i = -1, nx_glob
             P_per(i,j,ny_loc) = P_per(i,j,ny_loc-1)
           end do
         end do
        end if
      end block


      if (use_incompressible) then
          !$acc parallel loop  collapse(3) default(present)
          do     k = 0, ny_loc-1
            do   j = 0, nz_loc-1
              do i = 0, nx_glob-1
                U(i,j,k) = U_hat(i,j,k) - dtmax * (P_per(i,j,k) - P_per(i-1, j, k))*inv_dx
              end do
            end do
          end do

          !$acc parallel loop  collapse(3) default(present)
          do     k = 0, ny_loc-1
            do   j = 0, nz_loc-1
              do i = 0, nx_glob-1
                W(i,j,k) = W_hat(i,j,k) - dtmax * (P_per(i,j,k) - P_per(i, j-1, k))*inv_dz
              end do
            end do
          end do


          !$acc parallel loop  collapse(3) default(present)
          do     k = kmin_Vhat, ny_loc-1
            do   j = 0, nz_loc-1
              do i = 0, nx_glob-1
                V(i,j,k) = V_hat(i,j,k) - dtmax * (dPdy_ddy_hi(k) * P_per(i, j, k)   + &
                                                   dPdy_ddy_lo(k) * P_per(i, j, k-1) )
              end do
            end do
          end do

      else
        block
          real(dp)  :: P_fwd,P_bwd,inv_rho_min,inv_rho_now
          inv_rho_min = 1.d0/rho_min
          !$acc parallel loop  collapse(3) private(P_fwd,P_bwd,inv_rho_now) default(present)
          do     k = 0, ny_loc-1
            do   j = 0, nz_loc-1
              do i = 0, nx_glob-1
                inv_rho_now    =  0.5_dp/rho_arr(i,j,k) + 0.5_dp/rho_arr(i-1,j,k)
                P_fwd     =  P_per(i  ,j,k)*inv_rho_min + (inv_rho_now-inv_rho_min)*(2*P_now(i  ,j,k) - P_old(i  ,j,k))
                P_bwd     =  P_per(i-1,j,k)*inv_rho_min + (inv_rho_now-inv_rho_min)*(2*P_now(i-1,j,k) - P_old(i-1,j,k))
                U(i,j,k)  =  U_hat(i,j,k) - dtmax * (P_fwd - P_bwd)*inv_dx
              end do
            end do
          end do


          !$acc parallel loop  collapse(3) private(P_fwd,P_bwd,inv_rho_now) default(present)
          do     k = 0, ny_loc-1
            do   j = 0, nz_loc-1
              do i = 0, nx_glob-1
                inv_rho_now   =  0.5_dp/rho_arr(i,j,k) + 0.5_dp/rho_arr(i,j-1,k)
                P_fwd    =  P_per(i,j  ,k)*inv_rho_min + (inv_rho_now-inv_rho_min)*(2*P_now(i,j  ,k) - P_old(i,j  ,k))
                P_bwd    =  P_per(i,j-1,k)*inv_rho_min + (inv_rho_now-inv_rho_min)*(2*P_now(i,j-1,k) - P_old(i,j-1,k))
                W(i,j,k) =  W_hat(i,j,k) - dtmax * (P_fwd - P_bwd)*inv_dz
              end do
            end do
          end do


          !$acc parallel loop  collapse(3) private(P_fwd,P_bwd,inv_rho_now) default(present)
          do     k = kmin_Vhat, ny_loc-1
            do   j = 0, nz_loc-1
              do i = 0, nx_glob-1
                inv_rho_now    =  0.5_dp/rho_arr(i,j,k) + 0.5_dp/rho_arr(i,j,k-1)
                P_fwd     =  P_per(i,j,k  )*inv_rho_min + (inv_rho_now-inv_rho_min)*(2*P_now(i,j,k  ) - P_old(i,j,k  ))
                P_bwd     =  P_per(i,j,k-1)*inv_rho_min + (inv_rho_now-inv_rho_min)*(2*P_now(i,j,k-1) - P_old(i,j,k-1))
                V(i,j,k)  =  V_hat(i,j,k) - dtmax * (dPdy_ddy_hi(k) * P_fwd + dPdy_ddy_lo(k) * P_bwd )
              end do
            end do
          end do

        end block
      end if

      if (.not.use_incompressible) then
        !$acc parallel loop  collapse(3) default(present)
        do     k = -1, ny_loc
          do   j = -1, nz_loc
            do i = -1, nx_glob
              P_old(i,j,k) = P_now(i,j,k)
            end do
          end do
        end do

        !$acc parallel loop  collapse(3) default(present)
        do     k = -1, ny_loc
          do   j = -1, nz_loc
            do i = -1, nx_glob
              P_now(i,j,k) = P_per(i,j,k)
            end do
          end do
        end do
        call update_bulk_target(U, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, mu_arr, dy_cell, temp_arr_n_xzy, &
                          bulk_U_bot_inds, bulk_U_top_inds, min_delta_U_bulk, max_delta_U_bulk, Ly, bulk_target_Reb, &
                          bulk_use_rho_weighting, force_Sf_x_accel, force_Sf_x_pressure, bulk_use_target_Reb, &
                          bulk_use_target_Ub, use_rough_surf, init_Sf_x_mid_prev, utils_red_scalar, &
                          dtmax, inv_rho_dt_scale_Sf_x, Sf_x_mid, Sf_x_mid_prev, min_Sf_x_mid, max_Sf_x_mid, bulk_target_Ub)
      end if

      end block
      block 
        integer :: ii_hl
        do ii_hl=0,2
          call diezdecomp_halos_execute_generic(hl(ii_hl), U, work_tr)
          call diezdecomp_halos_execute_generic(hl(ii_hl), V, work_tr)
          call diezdecomp_halos_execute_generic(hl(ii_hl), W, work_tr)
        end do 
      end block

      if (mpi_pos_y == 0) then
        call set_val_slice(U,0, 0, wall_BC_U_bot) ! set_val_slice is zero-indexed
        call set_val_slice(V,0, 1, wall_BC_V_bot)
        call set_val_slice(W,0, 0, wall_BC_W_bot)
      end if

      if (mpi_pos_y == (mpi_divs_y-1)) then
        call set_val_slice(U,-1, -1, wall_BC_U_top) ! set_val_slice is zero-indexed (negative indexes follow python-conventions)
        call set_val_slice(V,-1, -1, wall_BC_V_top)
        call set_val_slice(W,-1, -1, wall_BC_W_top)
      end if
  end subroutine

  subroutine update_bulk_target(U, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, mu_arr, dy_cell, temp_arr_n_xzy, &
                          bulk_U_bot_inds, bulk_U_top_inds, min_delta_U_bulk, max_delta_U_bulk, Ly, bulk_target_Reb, &
                          bulk_use_rho_weighting, force_Sf_x_accel, force_Sf_x_pressure, bulk_use_target_Reb, &
                          bulk_use_target_Ub, use_rough_surf, init_Sf_x_mid_prev, utils_red_scalar, &
                          dtmax, inv_rho_dt_scale_Sf_x, Sf_x_mid, Sf_x_mid_prev, min_Sf_x_mid, max_Sf_x_mid, bulk_target_Ub)
    use mod_get_bulk, only: get_U_bulk, get_rho_mu_bulk
    implicit none
    integer  :: nx_glob, nz_loc, ny_loc, bulk_U_bot_inds(0:,0:), bulk_U_top_inds(0:,0:)
    real(dp) :: inv_rho_dt_scale_Sf_x(0:,0:,0:), U(-1:,-1:,-1:), vol_U_bulk, rho_arr(-1:,-1:,-1:), mu_arr(-1:,-1:,-1:), &
                dy_cell(0:), temp_arr_n_xzy(0:,0:,0:), min_delta_U_bulk, max_delta_U_bulk, Ly, bulk_target_Reb, dtmax, &
                Sf_x_mid, Sf_x_mid_prev, min_Sf_x_mid, max_Sf_x_mid, bulk_target_Ub
    real(dp), contiguous :: utils_red_scalar(0:)
    logical  :: bulk_use_rho_weighting, force_Sf_x_accel, force_Sf_x_pressure, bulk_use_target_Reb, bulk_use_target_Ub, &
                use_rough_surf, init_Sf_x_mid_prev
    block
      real(dp) :: delta_U_bulk
      if (bulk_use_target_Ub) then
        block
          real(dp)  ::  U_bulk
          call get_U_bulk(U, U_bulk, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, dy_cell, utils_red_scalar, &
                          temp_arr_n_xzy, use_rough_surf, bulk_U_bot_inds, bulk_U_top_inds, bulk_use_rho_weighting)
          delta_U_bulk      =  bulk_target_Ub - U_bulk
          min_delta_U_bulk  =  min(min_delta_U_bulk, delta_U_bulk)
          max_delta_U_bulk  =  max(max_delta_U_bulk, delta_U_bulk)
        end block
      end if

      if (bulk_use_target_Reb) then
        block
          real(dp)  ::  U_bulk, rho_bulk, mu_bulk
          call get_U_bulk(U, U_bulk, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, dy_cell, utils_red_scalar, &
                          temp_arr_n_xzy, use_rough_surf, bulk_U_bot_inds, bulk_U_top_inds, bulk_use_rho_weighting)
          call get_rho_mu_bulk(rho_bulk, mu_bulk, vol_U_bulk, nx_glob, nz_loc, ny_loc, rho_arr, mu_arr, dy_cell, &
                               utils_red_scalar, temp_arr_n_xzy, use_rough_surf, bulk_U_bot_inds, bulk_U_top_inds)
          delta_U_bulk      =  bulk_target_Reb/(rho_bulk*Ly/mu_bulk) - U_bulk  !  bulk_target_Ub - U_bulk
          min_delta_U_bulk  =  min(min_delta_U_bulk, delta_U_bulk)
          max_delta_U_bulk  =  max(max_delta_U_bulk, delta_U_bulk)
        end block
      end if

      if (force_Sf_x_accel) then
        block
          integer  :: i,j,k
          !$acc parallel loop  collapse(3) default(present)
          do     k = -1, ny_loc
            do   j = -1, nz_loc
              do i = -1, nx_glob
                U(i,j,k) = U(i,j,k) + delta_U_bulk
              end do
            end do
          end do
        end block
      end if

      if (force_Sf_x_pressure) then
        block
          real(dp) :: delta_Sf_x_mid, bulk_inv_rho_dt_scale_Sf_x
          !$acc wait
          call get_bulk_inv_rho_dt_scale_Sf_x(use_rough_surf, bulk_inv_rho_dt_scale_sf_x, inv_rho_dt_scale_Sf_x, dy_cell, &
                                      temp_arr_n_xzy, vol_U_bulk, nx_glob, nz_loc, ny_loc, bulk_U_bot_inds, bulk_U_top_inds,&
                                      utils_red_scalar)
          !$acc wait
          if (abs(dtmax)>0) then
            if (abs(bulk_inv_rho_dt_scale_Sf_x)<1e-20) then
              error stop 'abs(bulk_inv_rho_dt_scale_Sf_x)<1e-20'
            end if
            delta_Sf_x_mid = delta_U_bulk/bulk_inv_rho_dt_scale_Sf_x
            block
              integer  :: i,j,k
              !$acc parallel loop  collapse(3) default(present)
              do     k = 0, ny_loc-1
                do   j = 0, nz_loc-1
                  do i = 0, nx_glob-1
                    U(i,j,k) = U(i,j,k) + delta_Sf_x_mid*inv_rho_dt_scale_Sf_x(i,j,k)
                  end do
                end do
              end do
            end block
          end if
          block
            real(dp) :: Sf_x_mid_start
            Sf_x_mid        =  Sf_x_mid + delta_Sf_x_mid
            Sf_x_mid_start  =  Sf_x_mid
            if (.not.init_Sf_x_mid_prev) then
              init_Sf_x_mid_prev  =  .true.
              Sf_x_mid_prev       =  Sf_x_mid
            end if
            Sf_x_mid        =  2*Sf_x_mid - Sf_x_mid_prev
            Sf_x_mid_prev   =  Sf_x_mid_start
          end block
        end block
          min_Sf_x_mid  =  min(min_Sf_x_mid, Sf_x_mid)
          max_Sf_x_mid  =  max(max_Sf_x_mid, Sf_x_mid)
      end if
    end block
  end subroutine

  subroutine get_bulk_inv_rho_dt_scale_Sf_x(use_rough_surf, bulk_inv_rho_dt_scale_sf_x, inv_rho_dt_scale_Sf_x, dy_cell, &
                                      temp_arr_n_xzy, vol_U_bulk, nx_glob, nz_loc, ny_loc, bulk_U_bot_inds, bulk_U_top_inds, &
                                      utils_red_scalar)
    implicit none
    logical  ::  use_rough_surf
    real(dp) ::  bulk_inv_rho_dt_scale_sf_x, inv_rho_dt_scale_Sf_x(0:,0:,0:), dy_cell(0:), temp_arr_n_xzy(0:,0:,0:), &
                 vol_U_bulk
    real(dp), contiguous :: utils_red_scalar(0:)
    integer  ::  nx_glob, nz_loc, ny_loc, bulk_U_bot_inds(0:,0:), bulk_U_top_inds(0:,0:)
    block
      integer :: i, j, k, mpi_ierr
      ! ---------------- process bulk_inv_rho_dt_scale_Sf_x ----------------
      bulk_inv_rho_dt_scale_Sf_x = 0
      if (use_rough_surf) then
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(2) default(present) reduction(+:bulk_inv_rho_dt_scale_Sf_x)
        do      j = 0, nz_loc -1
          do    i = 0, nx_glob-1
            do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
              bulk_inv_rho_dt_scale_Sf_x = bulk_inv_rho_dt_scale_Sf_x + inv_rho_dt_scale_Sf_x(i,j,k)*dy_cell(k)
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
        !$acc parallel loop collapse(2) default(present)
        do      j = 0, nz_loc -1
          do    i = 0, nx_glob-1
            do  k = bulk_U_bot_inds(i,j),bulk_U_top_inds(i,j)
              temp_arr_n_xzy(i,j,k) = inv_rho_dt_scale_Sf_x(i,j,k)*dy_cell(k)
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, bulk_inv_rho_dt_scale_Sf_x, nx_glob*nz_loc*ny_loc, 0, &
                                           utils_red_scalar)
#endif
      else
#if !defined(_BIN_REDUCTION)
        !$acc parallel loop collapse(3) default(present) reduction(+:bulk_inv_rho_dt_scale_Sf_x)
        do      k  = 0, ny_loc -1
          do    j  = 0, nz_loc -1
            do  i  = 0, nx_glob-1
              bulk_inv_rho_dt_scale_Sf_x = bulk_inv_rho_dt_scale_Sf_x + inv_rho_dt_scale_Sf_x(i,j,k)*dy_cell(k)
            end do
          end do
        end do
#else
        call set_val_slice(temp_arr_n_xzy,0,-1, 0.d0)
        !$acc parallel loop collapse(3) default(present)
        do      k  = 0, ny_loc -1
          do    j  = 0, nz_loc -1
            do  i  = 0, nx_glob-1
              temp_arr_n_xzy(i,j,k) = inv_rho_dt_scale_Sf_x(i,j,k)*dy_cell(k)
            end do
          end do
        end do
        call quick_binary_scalar_reduction(temp_arr_n_xzy, bulk_inv_rho_dt_scale_Sf_x, nx_glob*nz_loc*ny_loc, 0, &
                                           utils_red_scalar)
#endif
      end if
      !$acc wait
      block
        real(dp)  ::  temp_real
        temp_real  =  bulk_inv_rho_dt_scale_Sf_x/vol_U_bulk
        call MPI_Allreduce(temp_real, bulk_inv_rho_dt_scale_Sf_x, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
      end block
      !$acc wait
    end block
  end subroutine
end module