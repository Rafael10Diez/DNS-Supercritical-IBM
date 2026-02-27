  ! ------------------ begin: Poisson deallocate ------------------

  deallocate(nxt_blocks_mpi     , & ! mpi_divs_x == mpi_divs_z
             nz_blocks_mpi      , &
             ny_blocks_mpi      , &
             nzt_py_blocks_mpi  , &
             lo_xt_blocks_mpi   , &
             lo_z_blocks_mpi    , &
             lo_y_blocks_mpi    , &
             lo_zt_py_blocks_mpi, &
             P                                    , &
             spec_x                               , &
             P_z                                  , &
             spec_z                               , &
             a_x                                  , &
             a_z                                  , &
             poi_is_00x                           , &
             poi_is_00z                           , &
             poi_C                                , &
             P_band_a   , P_band_b   , P_band_c   , &
             P_band_a_00, P_band_b_00, P_band_c_00)
  !$acc exit data delete(P, spec_x, P_z, spec_z, a_x, a_z, poi_is_00x, poi_is_00z, poi_C)
  !$acc exit data delete(P_band_a, P_band_b, P_band_c, P_band_a_00, P_band_b_00, P_band_c_00)
  if (nproc_mpi>1) then 
    deallocate(poi_slabs_tips, &
             ts_cp           , &
             ts_band_a       , &
             ts_band_b       , &
             trispec         , &
             poi_A           , &
             inv_poi_B       )
    !$acc exit data delete(poi_slabs_tips, ts_cp, ts_band_a, ts_band_b, trispec, poi_A, inv_poi_B)
  end if
  
  if (wsize_transp_poi>=1) then 
    deallocate(buffer_transp)
    !$acc exit data delete(buffer_transp)
  end if
! ------------------ end: Poisson deallocate ------------------
