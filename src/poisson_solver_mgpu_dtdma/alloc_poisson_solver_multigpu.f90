
  ! ----------------------------- Begin: MPI setup finish -----------------------------
  n_fft_x = 2*(nx_glob / 2 + 1)
  n_fft_z = 2*(nz_glob / 2 + 1)
  if (nproc_mpi.ne.(mpi_divs_y*mpi_divs_z)) error stop 'nproc_mpi.ne.(mpi_divs_y*mpi_divs_z)'
  
  allocate(nxt_blocks_mpi     (0:mpi_divs_z-1), & ! mpi_divs_x == mpi_divs_z
           nz_blocks_mpi      (0:mpi_divs_z-1), &
           ny_blocks_mpi      (0:mpi_divs_y-1), &
           nzt_py_blocks_mpi  (0:mpi_divs_y-1), &
           lo_xt_blocks_mpi   (0:mpi_divs_z-1), &
           lo_z_blocks_mpi    (0:mpi_divs_z-1), &
           lo_y_blocks_mpi    (0:mpi_divs_y-1), &
           lo_zt_py_blocks_mpi(0:mpi_divs_y-1)  )

  mpi_pos_y = irank_mpi/mpi_divs_z
  mpi_pos_z = irank_mpi - mpi_pos_y*mpi_divs_z

  call  get_lo_n_bounds_1d(mpi_divs_z, n_fft_x, lo_xt_blocks_mpi   , nxt_blocks_mpi   ) 
  call  get_lo_n_bounds_1d(mpi_divs_z, nz_glob, lo_z_blocks_mpi    , nz_blocks_mpi    ) 
  call  get_lo_n_bounds_1d(mpi_divs_y, ny_glob, lo_y_blocks_mpi    , ny_blocks_mpi    ) 
  call  get_lo_n_bounds_1d(mpi_divs_y, n_fft_z, lo_zt_py_blocks_mpi, nzt_py_blocks_mpi) 

  ! ----------------------------- End: MPI setup finish -----------------------------


  ! ------------------ begin: Poisson solver allocate ------------------
  nx_loc        =  nx_glob
  nx_transp     =  nxt_blocks_mpi(mpi_pos_z)
  nz_loc        =  nz_blocks_mpi(mpi_pos_z)
  ny_loc        =  ny_blocks_mpi(mpi_pos_y)
  nz_transp_py  =  nzt_py_blocks_mpi(mpi_pos_y)

  if (ny_glob <2) error stop 'ny_glob <2' 
  if (ny_loc  <2) error stop 'ny_loc  <2' 

  if (nx_loc.ne.nx_glob)            error stop 'nx_loc.ne.nx_glob'
  if (n_fft_x.ne.(2*(nx_glob/2+1))) error stop 'n_fft_x.ne.(2*(nx_glob/2+1))'
  
  allocate(P             (0:nx_glob-1, 0:nz_loc   -1, 0:ny_loc-1), &
           spec_x        (0:n_fft_x-1, 0:nz_loc   -1, 0:ny_loc-1), &
           P_z           (0:nz_glob-1, 0:nx_transp-1, 0:ny_loc-1), &
           spec_z        (0:n_fft_z-1, 0:nx_transp-1, 0:ny_loc-1), &
           a_x           (0:nx_transp-1)                         , &
           a_z           (0:n_fft_z-1)                           , &
           poi_is_00x    (0:nx_transp-1)                         , &
           poi_is_00z    (0:n_fft_z-1)                           , &
           poi_C         (0:n_fft_z-1, 0:nx_transp-1, 0:ny_loc-1), &
           P_band_a      (0:ny_loc-1), &
           P_band_b      (0:ny_loc-1), &
           P_band_c      (0:ny_loc-1), &
           P_band_a_00   (0:ny_loc-1), &
           P_band_b_00   (0:ny_loc-1), &
           P_band_c_00   (0:ny_loc-1))
  P           =  0
  spec_x      =  0
  P_z         =  0
  spec_z      =  0
  a_x         =  0
  a_z         =  0
  poi_is_00x  =  0
  poi_is_00z  =  0
  poi_C       =  0
  !$acc enter data copyin(P, spec_x, P_z, spec_z, a_x, a_z, poi_is_00x, poi_is_00z, poi_C)
  !$acc enter data create(P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00)

  if (nproc_mpi>1) then 
    allocate(poi_slabs_tips(0:n_fft_z-1, 0:nx_transp-1, 0:1 ), &
             ts_cp         (0:nz_transp_py-1, 0:nx_transp-1,  0:(2*mpi_divs_y-1)), &
             ts_band_a     (0:nz_transp_py-1, 0:nx_transp-1,  0:(2*mpi_divs_y-1)), &
             ts_band_b     (0:nz_transp_py-1, 0:nx_transp-1,  0:(2*mpi_divs_y-1)), &
             trispec       (0:nz_transp_py-1, 0:nx_transp-1,  0:(2*mpi_divs_y-1)), &
             poi_A         (0:n_fft_z-1, 0:nx_transp-1, 0:ny_loc-1), &
             inv_poi_B     (0:n_fft_z-1, 0:nx_transp-1, 0:ny_loc-1))
    !$acc enter data create(poi_slabs_tips, ts_cp, ts_band_a, ts_band_b, trispec, poi_A, inv_poi_B)
  end if

  ! ------------------ end: Poisson solver allocate ------------------
  
  block
    logical     ::  allow_alltoallv, allow_autotune_reorder
    allow_alltoallv         =  .false.
    allow_autotune_reorder  =  .false.
    wsize_transp_poi        =  1

    ! --------------------------- tr_yz/zy ---------------------------
    block
      integer(i8) ::  wsize
      integer     ::  lo_ref_in(0:2), lo_ref_out(0:2), offset6_in(0:5), offset6_out(0:5), &
                      order_in(0:2), order_out(0:2), order_intermediate(0:2), ii, jj
  
      ! yz: (0:n_fft_z-1, 0:nx_transp-1, 0:1) -> (0:nz_transp_py-1, 0:nx_transp-1,  0:2*mpi_divs_y-1)
      lo_ref_in   = (/mpi_pos_z,  0        , mpi_pos_y /) ! xzy order
      lo_ref_out = (/ mpi_pos_z, mpi_pos_y , 0         /)
      call diezdecomp_track_mpi_decomp(lo_ref_in , obj_poi_rank_tyz_in , irank_mpi, nproc_mpi)
      call diezdecomp_track_mpi_decomp(lo_ref_out, obj_poi_rank_tyz_out, irank_mpi, nproc_mpi)
      offset6_in              =  0
      offset6_out             =  0
      order_in                =  (/1, 0, 2/)
      order_out               =  (/1, 0, 2/)
      order_intermediate      =  (/1, 0, 2/)
      ii                      =  -1
      jj                      =  -1
      block 
        integer  ::  sp_in(0:2), sp_out(0:2)
        sp_in                   =  (/ n_fft_z      , nx_transp, 2            /) ! matches order_in
        sp_out                  =  (/ nz_transp_py , nx_transp, 2*mpi_divs_y /) ! matches order_out
        call diezdecomp_generic_fill_tr_obj(tr_poi_yz, obj_poi_rank_tyz_in, obj_poi_rank_tyz_out, &
                                            sp_in, offset6_in,   &
                                            sp_out, offset6_out, &
                                            order_in, order_out, &
                                            order_intermediate, allow_alltoallv, ii, jj, &
                                            wsize, allow_autotune_reorder)
        wsize_transp_poi = max(wsize_transp_poi, wsize)
        call diezdecomp_generic_fill_tr_obj(tr_poi_zy, obj_poi_rank_tyz_out, obj_poi_rank_tyz_in, &
                                            sp_out, offset6_out, &
                                            sp_in, offset6_in,   &
                                            order_out, order_in, &
                                            order_intermediate, allow_alltoallv, jj, ii, &
                                            wsize, allow_autotune_reorder)
        wsize_transp_poi = max(wsize_transp_poi, wsize)
      end block
    end block 
  
    ! --------------------------- tr_xy/yx ---------------------------
    block
      integer(i8) ::  wsize
      integer     ::  lo_ref_in(0:2), lo_ref_out(0:2), offset6_in(0:5), offset6_out(0:5), &
                      order_in(0:2), order_out(0:2), order_intermediate(0:2), sp_in(0:2), sp_out(0:2), ii, jj
  
      ! xy: (0:n_fft_x-1, 0:nz_loc   -1, 0:ny_loc-1) => (0:nz_glob-1, 0:nx_transp-1, 0:ny_loc-1)
      lo_ref_in   =  (/ 0        , mpi_pos_z, mpi_pos_y /)
      lo_ref_out  =  (/ mpi_pos_z, 0        , mpi_pos_y /)
      call diezdecomp_track_mpi_decomp(lo_ref_in , obj_poi_rank_txy_in , irank_mpi, nproc_mpi)
      call diezdecomp_track_mpi_decomp(lo_ref_out, obj_poi_rank_txy_out, irank_mpi, nproc_mpi)
      offset6_in              =  0
      offset6_out             =  0
      order_in                =  (/0, 1, 2/)
      order_out               =  (/1, 0, 2/)
      order_intermediate      =  (/1, 0, 2/)
      ii                      =  -1
      jj                      =  -1
      sp_in                   =  (/ n_fft_x , nz_loc   , ny_loc /)
      sp_out                  =  (/ nz_glob , nx_transp, ny_loc /)
      call diezdecomp_generic_fill_tr_obj(tr_poi_xy, obj_poi_rank_txy_in, obj_poi_rank_txy_out, &
                                          sp_in, offset6_in,   &
                                          sp_out, offset6_out, &
                                          order_in, order_out, &
                                          order_intermediate, allow_alltoallv, ii, jj, &
                                          wsize, allow_autotune_reorder)
      wsize_transp_poi = max(wsize_transp_poi, wsize)
      call diezdecomp_generic_fill_tr_obj(tr_poi_yx, obj_poi_rank_txy_out, obj_poi_rank_txy_in, &
                                          sp_out, offset6_out, &
                                          sp_in, offset6_in,   &
                                          order_out, order_in, &
                                          order_intermediate, allow_alltoallv, jj, ii, &
                                          wsize, allow_autotune_reorder)
      wsize_transp_poi = max(wsize_transp_poi, wsize)
    end block 
    
      allocate(buffer_transp(0:wsize_transp_poi-1))
      buffer_transp = 0
      !$acc enter data copyin(buffer_transp)
  end block
