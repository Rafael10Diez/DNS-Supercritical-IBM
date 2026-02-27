
! allocate 1-D arrays
allocate( Ruw_ddy_a          (0:ny_loc-1) ,&
          Ruw_ddy_b          (0:ny_loc-1) ,&
          Ruw_ddy_c          (0:ny_loc-1) ,&
          Ruw_d2dy2_a        (0:ny_loc-1) ,&
          Ruw_d2dy2_b        (0:ny_loc-1) ,&
          Ruw_d2dy2_c        (0:ny_loc-1) ,&
          Ruw_pad_d2dy2_a    (0:ny_loc-1) ,&
          Ruw_pad_d2dy2_b    (0:ny_loc-1) ,&
          Ruw_pad_d2dy2_c    (0:ny_loc-1) ,&
          Ruw_pad_ddy_a      (0:ny_loc-1) ,&
          Ruw_pad_ddy_b      (0:ny_loc-1) ,&
          Ruw_pad_ddy_c      (0:ny_loc-1) ,&
          inv_dy_p           (0:ny_loc-1) ,&
          inv_dy_u           (0:ny_loc  ) ,&
          inv_dy_p_pad       (0:ny_loc  ) ,&
          dy_cell            (0:ny_loc-1) ,&
          dPdy_ddy_lo        (0:ny_loc-1) ,&
          dPdy_ddy_hi        (0:ny_loc-1) ,&
          Rv_interp_uw_bottom(0:ny_loc  ) ,&
          Rv_interp_uw_top   (0:ny_loc  ) ,&
          Rv_ddy_a           (0:ny_loc  ) ,&
          Rv_ddy_b           (0:ny_loc  ) ,&
          Rv_ddy_c           (0:ny_loc  ) ,&
          Rv_d2dy2_a         (0:ny_loc  ) ,&
          Rv_d2dy2_b         (0:ny_loc  ) ,&
          Rv_d2dy2_c         (0:ny_loc  ) ,&
          buffer_1d          (0:max(ny_glob, nx_glob, nz_glob)))

! fill 1-D arrays
block
    integer :: ny_prev
    ny_prev = lo_y_blocks_mpi(mpi_pos_y)
    call read_column_r8 (Ruw_pad_ddy_a       , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_pad_ddy_abc.dat'  , 0, irank_mpi)
    call read_column_r8 (Ruw_pad_ddy_b       , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_pad_ddy_abc.dat'  , 1, irank_mpi)
    call read_column_r8 (Ruw_pad_ddy_c       , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_pad_ddy_abc.dat'  , 2, irank_mpi)

    call read_column_r8 (Ruw_ddy_a           , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_ddy_abc.dat'  , 0, irank_mpi)
    call read_column_r8 (Ruw_ddy_b           , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_ddy_abc.dat'  , 1, irank_mpi)
    call read_column_r8 (Ruw_ddy_c           , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_ddy_abc.dat'  , 2, irank_mpi)

    call read_column_r8 (Ruw_d2dy2_a         , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_d2dy2_abc.dat', 0, irank_mpi)
    call read_column_r8 (Ruw_d2dy2_b         , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_d2dy2_abc.dat', 1, irank_mpi)
    call read_column_r8 (Ruw_d2dy2_c         , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_d2dy2_abc.dat', 2, irank_mpi)

    call read_column_r8 (Ruw_pad_d2dy2_a     , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_pad_d2dy2_abc.dat', 0, irank_mpi)
    call read_column_r8 (Ruw_pad_d2dy2_b     , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_pad_d2dy2_abc.dat', 1, irank_mpi)
    call read_column_r8 (Ruw_pad_d2dy2_c     , ny_loc, ny_prev, buffer_1d, './geom_data/Ruw_pad_d2dy2_abc.dat', 2, irank_mpi)

    call read_column_r8 (dPdy_ddy_lo         , ny_loc, ny_prev, buffer_1d, './geom_data/dPdy_ddy_lh.dat'  , 0, irank_mpi)
    call read_column_r8 (dPdy_ddy_hi         , ny_loc, ny_prev, buffer_1d, './geom_data/dPdy_ddy_lh.dat'  , 1, irank_mpi)

    call read_column_r8 (P_band_a            , ny_loc, ny_prev, buffer_1d, './geom_data/P_band_abc.dat'   , 0, irank_mpi)
    call read_column_r8 (P_band_b            , ny_loc, ny_prev, buffer_1d, './geom_data/P_band_abc.dat'   , 1, irank_mpi)
    call read_column_r8 (P_band_c            , ny_loc, ny_prev, buffer_1d, './geom_data/P_band_abc.dat'   , 2, irank_mpi)
    call read_column_r8 (P_band_a_00         , ny_loc, ny_prev, buffer_1d, './geom_data/P_band_abc_00.dat', 0, irank_mpi)
    call read_column_r8 (P_band_b_00         , ny_loc, ny_prev, buffer_1d, './geom_data/P_band_abc_00.dat', 1, irank_mpi)
    call read_column_r8 (P_band_c_00         , ny_loc, ny_prev, buffer_1d, './geom_data/P_band_abc_00.dat', 2, irank_mpi)

    call read_column_r8 (inv_dy_p            , ny_loc  , ny_prev, buffer_1d, './geom_data/inv_dy_p.dat'    , 0, irank_mpi)
    call read_column_r8 (inv_dy_u            , ny_loc+1, ny_prev, buffer_1d, './geom_data/inv_dy_u.dat'    , 0, irank_mpi)
    call read_column_r8 (inv_dy_p_pad        , ny_loc+1, ny_prev, buffer_1d, './geom_data/inv_dy_p_pad.dat', 0, irank_mpi)

  call read_column_r8 (Rv_interp_uw_bottom,ny_loc+1,ny_prev,buffer_1d,'./geom_data/Rv_interp_uw_bottom_top.dat',0,irank_mpi)
  call read_column_r8 (Rv_interp_uw_top   ,ny_loc+1,ny_prev,buffer_1d,'./geom_data/Rv_interp_uw_bottom_top.dat',1,irank_mpi)

    call read_column_r8 (Rv_ddy_a            , ny_loc+1, ny_prev, buffer_1d, './geom_data/Rv_ddy_abc.dat'  , 0, irank_mpi)
    call read_column_r8 (Rv_ddy_b            , ny_loc+1, ny_prev, buffer_1d, './geom_data/Rv_ddy_abc.dat'  , 1, irank_mpi)
    call read_column_r8 (Rv_ddy_c            , ny_loc+1, ny_prev, buffer_1d, './geom_data/Rv_ddy_abc.dat'  , 2, irank_mpi)

    call read_column_r8 (Rv_d2dy2_a          , ny_loc+1, ny_prev, buffer_1d, './geom_data/Rv_d2dy2_abc.dat', 0, irank_mpi)
    call read_column_r8 (Rv_d2dy2_b          , ny_loc+1, ny_prev, buffer_1d, './geom_data/Rv_d2dy2_abc.dat', 1, irank_mpi)
    call read_column_r8 (Rv_d2dy2_c          , ny_loc+1, ny_prev, buffer_1d, './geom_data/Rv_d2dy2_abc.dat', 2, irank_mpi)

end block
dy_cell = 1.d0/inv_dy_p
! move to GPU 1-D arrays

!$acc enter data copyin(Ruw_ddy_a, Ruw_ddy_b, Ruw_ddy_c, Ruw_d2dy2_a, Ruw_d2dy2_b, Ruw_d2dy2_c)
!$acc enter data copyin(inv_dy_p, inv_dy_u, inv_dy_p_pad, dy_cell, dPdy_ddy_lo, dPdy_ddy_hi)
!$acc enter data copyin(Rv_interp_uw_bottom, Rv_interp_uw_top, Rv_ddy_a, Rv_ddy_b, Rv_ddy_c)
!$acc enter data copyin(Rv_d2dy2_a, Rv_d2dy2_b, Rv_d2dy2_c)

!$acc enter data copyin(Ruw_pad_d2dy2_a, Ruw_pad_d2dy2_b, Ruw_pad_d2dy2_c)
!$acc enter data copyin(Ruw_pad_ddy_a, Ruw_pad_ddy_b, Ruw_pad_ddy_c)

! 3-D arrays
allocate( Ru_old      ( 0:nx_loc-1,  0:nz_loc-1,  0:ny_loc-1) ,&
          Rw_old      ( 0:nx_loc-1,  0:nz_loc-1,  0:ny_loc-1) ,&
          Rv_old      ( 0:nx_loc-1,  0:nz_loc-1,  0:ny_loc-1) ,&
          Rt_old      ( 0:nx_loc-1,  0:nz_loc-1,  0:ny_loc-1) ,&
          T           (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          T_ext       (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          U           (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          W           (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          V           (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          U_hat       (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          W_hat       (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          V_hat       (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
          P_per       (-1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) )

Ru_old       = 0
Rw_old       = 0
Rv_old       = 0
Rt_old       = 0
U            = 0
W            = 0
V            = 0
T            = 0
T_ext        = 0
U_hat        = 0
W_hat        = 0
V_hat        = 0

!$acc enter data copyin(Ru_old, Rw_old, Rv_old, Rt_old, U, W, V, T, T_ext, U_hat, W_hat, V_hat, P_per)

if (.not.use_incompressible) then
  allocate(div_now          ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           P_now            ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           P_old            ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           RHS_energy       (  0:nx_loc-1,  0:nz_loc-1,  0:ny_loc-1) ,&
           mu_arr           ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           cond_arr         ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           rho_arr          ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           buoyancy_arr     ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           enth_arr         ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           cp_arr           ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) ,&
           div_jacobian_arr ( -1:nx_loc  , -1:nz_loc  , -1:ny_loc  ) )
  div_now          = 0
  P_now            = 0
  P_old            = 0
  RHS_energy       = 0
  mu_arr           = 0
  cond_arr         = 0
  rho_arr          = 0
  buoyancy_arr     = 0
  enth_arr         = 0
  cp_arr           = 0
  div_jacobian_arr = 0
  !$acc enter data copyin(div_now, P_now, P_old, RHS_energy, mu_arr, cond_arr, rho_arr, buoyancy_arr, enth_arr, cp_arr)
  !$acc enter data copyin(div_jacobian_arr)
  if (force_Sf_x_pressure) then
    allocate(rho_arr_old          ( -1:nx_loc   , -1:nz_loc   , -1:ny_loc   ) ,&
             inv_rho_dt_scale_Sf_x(  0:nx_loc-1 ,  0:nz_loc-1 ,  0:ny_loc-1 ) )
    !$acc enter data copyin(rho_arr_old, inv_rho_dt_scale_Sf_x)
    allocate(rho_profile_y(0:ny_loc-1))
    rho_profile_y = 1
    !$acc enter data copyin(rho_profile_y)
  end if
  if (use_utau_work) then
    allocate(utau_old         (  0:nx_loc-1,  0:nz_loc-1,  0:ny_loc-1) ,&
             utau_now         (  0:nx_loc-1,  0:nz_loc-1,  0:ny_loc-1) )
    utau_old         = 0
    utau_now         = 0
    !$acc enter data copyin(utau_old, utau_now)
  end if
end if

! --------------------- bulk parameters ---------------------
 vol_U_bulk  = -1 ! defined later
 vol_T_bulk  = -1
 vol_div_avg = -1

block
  integer :: mpi_ierr

! --------------------- avg parameters ---------------------
!__[PYTHON_HOLDER_AVG_SET]
! allocate({x}(0:nx-1, 0:nz-1, 0:ny-1))
! {x} = 0
! !$acc enter data copyin({x})
! !$acc wait
! if (avg_bin_nbins>0) then
!   allocate(binary_{x}(0:nx-1, 0:nz-1, 0:ny-1, 0:avg_bin_nbins-1))
!   binary_{x} = 0
! end if"

if (avg_bin_nbins>0) then
  allocate(binary_cycle(0:avg_bin_nbins-1)                        )
  binary_cycle = 0
end if

 avg_N_now     =  0
 avg_inv_nraw  =  1.d0/(0.d0 + avg_nraw)

 end block

min_delta_U_bulk     =  1e20
max_delta_U_bulk     = -1e20
min_Sf_x_mid         =  1e20
max_Sf_x_mid         = -1e20

init_Sf_x_mid_prev   =  .false.
Sf_x_mid             =  0.

allocate(temp_arr_n_xzy(0:nx_loc-1, 0:nz_loc-1, 0:ny_loc-1))
temp_arr_n_xzy = 0
!$acc enter data copyin(temp_arr_n_xzy)

!__[PYTHON_HOLDER_PROF1D_SET]
! allocate({x}(0:ny-1))
! {x} = 0
! !$acc enter data copyin({x})

block
  integer :: ny_loc_max
  ny_loc_max  =  maxval(ny_blocks_mpi(0:mpi_divs_y-1))
  if (ny_loc>ny_loc_max) error stop 'ny_loc>ny_loc_max'
  allocate(prof_1d_temp(0:ny_loc_max-1))
  prof_1d_temp   = 0
  !$acc enter data copyin(prof_1d_temp)
end block

prof_1d_avg_N = 0

block
  integer     ::  lo_ref(0:2)
  allocate(lo_00k_blocks_mpi(0:nproc_mpi-1), &
             n00k_blocks_mpi(0:nproc_mpi-1))
  call  get_lo_n_bounds_1d(nproc_mpi, ny_glob, lo_00k_blocks_mpi, n00k_blocks_mpi)
  mpi_pos_00k  =  mpi_pos_z + mpi_pos_y*mpi_divs_z
  if (mpi_pos_00k.ne.irank_mpi) error stop 'mpi_pos_00k.ne.irank_mpi'
  ny_loc_00k   =  n00k_blocks_mpi(mpi_pos_00k)
  lo_ref       =  (/ 0, 0, mpi_pos_00k /)
  call diezdecomp_track_mpi_decomp(lo_ref, obj_ranks_00k, irank_mpi, nproc_mpi)
end block

block
  integer  ::  lo_ref(0:2)
  lo_ref   =  (/ 0, mpi_pos_z, mpi_pos_y /)
  call diezdecomp_track_mpi_decomp(lo_ref, obj_ranks_0zy, irank_mpi, nproc_mpi)
end block


block
  integer(i8) ::  wsize_tr
  wsize_tr  =  1
  block
    logical     ::  allow_alltoallv, allow_autotune_reorder
    integer     ::  offset6_in(0:2,0:1), offset6_out(0:2,0:1), &
                    order_in(0:2), order_out(0:2), order_intermediate(0:2), ii, jj
    allow_alltoallv         =  .false.
    allow_autotune_reorder  =  .false.
    offset6_in              =  0
    offset6_out             =  0
    order_in                =  (/0, 1, 2/)
    order_out               =  (/0, 1, 2/)
    order_intermediate      =  (/0, 1, 2/)
    ii                      =  -1
    jj                      =  -1

    block
      integer     :: sp_in(0:2), sp_out(0:2)
      integer(i8) :: wsize
      sp_in                   =  (/ 1 ,          1, ny_loc     /) ! matches order_in
      sp_out                  =  (/ 1 , mpi_divs_z, ny_loc_00k /) ! matches order_out
      call diezdecomp_generic_fill_tr_obj(tr_red1d_fwd, obj_ranks_0zy, obj_ranks_00k, &
                                          sp_in, offset6_in,   &
                                          sp_out, offset6_out, &
                                          order_in, order_out, &
                                          order_intermediate, allow_alltoallv, ii, jj, &
                                          wsize, allow_autotune_reorder)
      wsize_tr  =  max(wsize_tr, wsize)
      call diezdecomp_generic_fill_tr_obj(tr_red1d_bwd, obj_ranks_00k, obj_ranks_0zy, &
                                          sp_out, offset6_out, &
                                          sp_in, offset6_in,   &
                                          order_out, order_in, &
                                          order_intermediate, allow_alltoallv, jj, ii, &
                                          wsize, allow_autotune_reorder)
      wsize_tr  =  max(wsize_tr, wsize)
      allocate(tr_red1d_B(0:0,0:mpi_divs_z-1,0:ny_loc_00k-1), tr_red1d_B_1d(0:ny_loc_00k-1))
      !$acc enter data create(tr_red1d_B, tr_red1d_B_1d)
    end block

    block
      integer     :: sp_in(0:2), sp_out(0:2)
      integer(i8) :: wsize
      sp_in                   =  (/ nx_glob, nz_loc,  ny_loc     /) ! matches order_in
      sp_out                  =  (/ nx_glob, nz_glob, ny_loc_00k /) ! matches order_out
      call diezdecomp_generic_fill_tr_obj(tr_io_fwd, obj_ranks_0zy, obj_ranks_00k, &
                                          sp_in, offset6_in,   &
                                          sp_out, offset6_out, &
                                          order_in, order_out, &
                                          order_intermediate, allow_alltoallv, ii, jj, &
                                          wsize, allow_autotune_reorder)
      wsize_tr  =  max(wsize_tr, wsize)
      call diezdecomp_generic_fill_tr_obj(tr_io_bwd, obj_ranks_00k, obj_ranks_0zy, &
                                          sp_out, offset6_out, &
                                          sp_in, offset6_in,   &
                                          order_out, order_in, &
                                          order_intermediate, allow_alltoallv, jj, ii, &
                                          wsize, allow_autotune_reorder)
      wsize_tr  =  max(wsize_tr, wsize)
      block
        integer :: ny_loc_00k_max
        ny_loc_00k_max  =  maxval(n00k_blocks_mpi(0:nproc_mpi-1))
        allocate(slab_00k(0:nx_glob-1, 0:nz_glob-1, 0:ny_loc_00k_max-1))
        !$acc enter data create(slab_00k)
      end block
    end block

    ! tr_io_fwd%send_autotuned       = 2 ; tr_io_fwd%recv_autotuned       = 2
    ! tr_io_fwd%send_mode_op_batched = 6 ; tr_io_fwd%recv_mode_op_batched = 6
    ! tr_io_fwd%send_mode_op_simul   = 6 ; tr_io_fwd%recv_mode_op_simul   = 6
    !
    ! tr_io_bwd%send_autotuned       = 2 ; tr_io_bwd%recv_autotuned       = 2
    ! tr_io_bwd%send_mode_op_batched = 6 ; tr_io_bwd%recv_mode_op_batched = 6
    ! tr_io_bwd%send_mode_op_simul   = 6 ; tr_io_bwd%recv_mode_op_simul   = 6
    !
    ! tr_red1d_fwd%send_autotuned       = 2 ; tr_red1d_fwd%recv_autotuned       = 2
    ! tr_red1d_fwd%send_mode_op_batched = 6 ; tr_red1d_fwd%recv_mode_op_batched = 6
    ! tr_red1d_fwd%send_mode_op_simul   = 6 ; tr_red1d_fwd%recv_mode_op_simul   = 6
    !
    ! tr_red1d_bwd%send_autotuned       = 2 ; tr_red1d_bwd%recv_autotuned       = 2
    ! tr_red1d_bwd%send_mode_op_batched = 6 ; tr_red1d_bwd%recv_mode_op_batched = 6
    ! tr_red1d_bwd%send_mode_op_simul   = 6 ; tr_red1d_bwd%recv_mode_op_simul   = 6
    !
    ! tr_poi_yz%send_autotuned       = 2 ; tr_poi_yz%recv_autotuned       = 2
    ! tr_poi_yz%send_mode_op_batched = 6 ; tr_poi_yz%recv_mode_op_batched = 6
    ! tr_poi_yz%send_mode_op_simul   = 6 ; tr_poi_yz%recv_mode_op_simul   = 6
    !
    ! tr_poi_zy%send_autotuned       = 2 ; tr_poi_zy%recv_autotuned       = 2
    ! tr_poi_zy%send_mode_op_batched = 6 ; tr_poi_zy%recv_mode_op_batched = 6
    ! tr_poi_zy%send_mode_op_simul   = 6 ; tr_poi_zy%recv_mode_op_simul   = 6
    !
    ! tr_poi_xy%send_autotuned       = 2 ; tr_poi_xy%recv_autotuned       = 2
    ! tr_poi_xy%send_mode_op_batched = 6 ; tr_poi_xy%recv_mode_op_batched = 6
    ! tr_poi_xy%send_mode_op_simul   = 6 ; tr_poi_xy%recv_mode_op_simul   = 6
    !
    ! tr_poi_yx%send_autotuned       = 2 ; tr_poi_yx%recv_autotuned       = 2
    ! tr_poi_yx%send_mode_op_batched = 6 ; tr_poi_yx%recv_mode_op_batched = 6
    ! tr_poi_yx%send_mode_op_simul   = 6 ; tr_poi_yx%recv_mode_op_simul   = 6

  end block

  block
    integer     ::  force_halo_sync, autotuned_pack
    integer     ::  A_shape(0:2), offset6(0:2,0:1), ii, nh_xyz(0:2), order_halo(0:2)
    logical     ::  periodic_xyz(0:2)
    integer(i8) ::  wsize
    if (nskip.ne.0) error stop 'skip.ne.0'
    if (nhalo.ne.1) error stop 'nhalo.ne.1'
    A_shape          =  (/ nx_glob, nz_loc, ny_loc /) + 2*nskip + 2*nhalo
    offset6          =  nskip
    nh_xyz           =  nhalo
    order_halo       =  (/0, 1, 2/)
    force_halo_sync  =  2
    autotuned_pack   =  2

    do ii=0,2
      periodic_xyz = .false.
      if (ii<2) periodic_xyz(ii) = .true.
      call diezdecomp_generic_fill_hl_obj(hl(ii), obj_ranks_0zy, A_shape, offset6, ii, nh_xyz, order_halo, periodic_xyz, wsize,&
                                              force_halo_sync, autotuned_pack)
      wsize_tr  =  max(wsize_tr, wsize)
    end do
  end block

  allocate(work_tr(0:wsize_tr-1))
  !$acc enter data create(work_tr)
end block

allocate(utils_red_scalar(0:0))
!$acc enter data create(utils_red_scalar)
!$acc wait