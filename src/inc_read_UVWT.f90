
block 
  character(len=24)  ::  my_file
    my_file  =  '../input/array_ini_U.dat'
    call cfd_read_arr(my_file, U, -1, tr_io_bwd, slab_00k, temp_arr_n_xzy, work_tr, &
                      irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
    my_file  =  '../input/array_ini_V.dat'
    call cfd_read_arr(my_file, V, -1, tr_io_bwd, slab_00k, temp_arr_n_xzy, work_tr, &
                      irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
    my_file  =  '../input/array_ini_W.dat'
    call cfd_read_arr(my_file, W, -1, tr_io_bwd, slab_00k, temp_arr_n_xzy, work_tr, &
                      irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
    my_file  =  '../input/array_ini_T.dat'
    call cfd_read_arr(my_file, T, -1, tr_io_bwd, slab_00k, temp_arr_n_xzy, work_tr, &
                      irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
end block

  call cfd_periodicity(hl, U, V, W, T, work_tr, mpi_pos_y, mpi_divs_y, &
                       wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                       wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top)
                       
block 
  character(len=28)  ::  my_file
  if (.not.use_incompressible) then 
    my_file  =  '../input/array_ini_P_old.dat'
    call cfd_read_arr(my_file, P_now, -1, tr_io_bwd, slab_00k, temp_arr_n_xzy, work_tr, &
                      irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
    call periodicity_pressure(hl, P_now, work_tr, &
                              nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, use_incompressible)
    P_old = P_now
    block 
      integer :: i,j,k
      !$acc parallel loop  collapse(3) default(present)
      do     k = -1, ny_loc
        do   j = -1, nz_loc
          do i = -1, nx_glob 
            P_old(i,j,k) = P_now(i,j,k)
          end do
        end do
      end do
    end block 
  
    my_file  =  '../input/array_ini_P_now.dat'
    call cfd_read_arr(my_file, P_now, -1, tr_io_bwd, slab_00k, temp_arr_n_xzy, work_tr, &
                      irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
    call periodicity_pressure(hl, P_now, work_tr, &
                              nx_glob, nz_loc, ny_loc, mpi_pos_y, mpi_divs_y, use_incompressible)
  end if 
end block
