! ------------------ begin: Poisson solver declare ------------------

  integer                ::  nx_glob , nz_glob , ny_glob , n_fft_x  , n_fft_z, &
                             nx_loc  , nz_loc  , ny_loc  , nx_transp, nz_transp_py, &
                             mpi_divs_y, mpi_divs_z, mpi_pos_y, mpi_pos_z

  integer, allocatable   ::  nxt_blocks_mpi     (:), & ! mpi_divs_x == mpi_divs_z
                             nz_blocks_mpi      (:), &
                             ny_blocks_mpi      (:), &
                             nzt_py_blocks_mpi  (:), &
                             lo_xt_blocks_mpi   (:), &
                             lo_z_blocks_mpi    (:), &
                             lo_y_blocks_mpi    (:), &
                             lo_zt_py_blocks_mpi(:), &
                             poi_is_00x         (:), &
                             poi_is_00z         (:)

  real(dp), allocatable  ::  ts_cp         (:,:,:) ,&
                             ts_band_a     (:,:,:) ,&
                             ts_band_b     (:,:,:) ,&
                             P             (:,:,:) ,&
                             P_z           (:,:,:) ,&
                             a_x           (:)     ,&
                             a_z           (:)     ,&
                             poi_A         (:,:,:) ,&
                             poi_C         (:,:,:) ,&
                             inv_poi_B     (:,:,:) ,&
                             P_band_a      (:)     ,&
                             P_band_b      (:)     ,&
                             P_band_c      (:)     ,&
                             P_band_a_00   (:)     ,&
                             P_band_b_00   (:)     ,&
                             P_band_c_00   (:)     ,&
                             buffer_transp (:)     ,&
                             poi_slabs_tips(:,:,:) ,&
                             trispec(:,:,:)        ,&
                             spec_x (:,:,:)        ,&
                             spec_z (:,:,:)

#if defined(_USE_FFTW)
    type(C_PTR)           ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#else 
#if defined(_USE_HIPFFT)
    type(planinfo_hipfft) ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#else
    integer               ::  plan_fwd_x, plan_bwd_x, plan_fwd_z, plan_bwd_z
#endif
#endif
    integer(i8)  ::  wsize_transp_poi
! diezDecomp handles
  type(diezdecomp_props_transp)     ::  tr_poi_yz, tr_poi_zy, &
                                        tr_poi_xy, tr_poi_yx
  type(diezdecomp_parsed_mpi_ranks) ::  obj_poi_rank_txy_in, obj_poi_rank_txy_out, &
                                        obj_poi_rank_tyz_in, obj_poi_rank_tyz_out

! ------------------ end: Poisson solver declare ------------------