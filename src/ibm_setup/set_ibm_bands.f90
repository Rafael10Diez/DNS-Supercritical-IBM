  block 
    integer  ::  iter_outer
    logical  ::  only_size
    do iter_outer=0,1
      only_size  =  (iter_outer==0)
      call ibm_read_coeffs('./geom_data/ibm_coeffs_U.dat', all_ijk_U, all_c_U, n_ghosts_U, n_fluids_U, &
                           nx_glob, nz_glob, ny_glob, only_size, irank_mpi, nproc_mpi)
      call ibm_read_coeffs('./geom_data/ibm_coeffs_V.dat', all_ijk_V, all_c_V, n_ghosts_V, n_fluids_V, &
                           nx_glob, nz_glob, ny_glob, only_size, irank_mpi, nproc_mpi)
      call ibm_read_coeffs('./geom_data/ibm_coeffs_W.dat', all_ijk_W, all_c_W, n_ghosts_W, n_fluids_W, &
                           nx_glob, nz_glob, ny_glob, only_size, irank_mpi, nproc_mpi)
      call ibm_read_coeffs('./geom_data/ibm_coeffs_T.dat', all_ijk_T, all_c_T, n_ghosts_T, n_fluids_T, &
                           nx_glob, nz_glob, ny_glob, only_size, irank_mpi, nproc_mpi)
      if (only_size) then 
        allocate(all_ijk_U(0:2,0:2,0:n_ghosts_U-1),&
                 all_ijk_V(0:2,0:2,0:n_ghosts_V-1),&
                 all_ijk_W(0:2,0:2,0:n_ghosts_W-1),&
                 all_ijk_T(0:2,0:2,0:n_ghosts_T-1),&
                 all_c_U  (    0:2,0:n_ghosts_U-1),&
                 all_c_V  (    0:2,0:n_ghosts_V-1),&
                 all_c_W  (    0:2,0:n_ghosts_W-1),&
                 all_c_T  (    0:2,0:n_ghosts_T-1))
        if (n_fluids_U.ne.2) error stop 'n_fluids_U.ne.2'
        if (n_fluids_V.ne.2) error stop 'n_fluids_V.ne.2'
        if (n_fluids_W.ne.2) error stop 'n_fluids_W.ne.2'
        if (n_fluids_T.ne.2) error stop 'n_fluids_T.ne.2'
      end if 
    end do 
  end block 

if (.not.use_incompressible) then 
  ! define all_ijk_div, all_c_div
  allocate(all_ijk_div, mold = all_ijk_T)
  allocate(all_c_div  , mold = all_c_T  )
  all_ijk_div  =  all_ijk_T
  all_c_div    =  0
  n_ghosts_div = n_ghosts_T
  n_fluids_div = n_fluids_T
  block 
    integer   :: i, mpi_ierr, temp_int
    real(dp)  :: dx, dz, d_ab, d_ta
    dx       =  1/inv_dx
    dz       =  1/inv_dz
    allocate(yp_full(0:ny_glob-1))
    if (irank_mpi == 0) then 
      open(19,file="./geom_data/yp.dat")
        read(19,*) temp_int
        if (temp_int /= ny_glob) error stop 'ERROR: Reading yp_full: temp_int /= ny_glob'
        do i=0,ny_glob-1
          read(19,*) yp_full(i)
        end do 
      close(19)
    end if 
    call MPI_BCAST(yp_full, ny_glob, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
    block 
      integer :: it,jt,kt,ia,ja,ka,ib,jb,kb
      do i=0, n_ghosts_div-1
        it =            all_ijk_div(0,0,i)
        jt =            all_ijk_div(1,0,i)
        kt =            all_ijk_div(2,0,i)
        ia = quick_pfix(all_ijk_div(0,1,i), it, nx_glob)
        ja = quick_pfix(all_ijk_div(1,1,i), jt, nz_glob)
        ka =            all_ijk_div(2,1,i)
        ib = quick_pfix(all_ijk_div(0,2,i), it, nx_glob)
        jb = quick_pfix(all_ijk_div(1,2,i), jt, nz_glob)
        kb =            all_ijk_div(2,2,i)
        d_ab = dsqrt( (dx*(ib -          ia) )**2 + &
                      (dz*(jb -          ja) )**2 + &
                  (yp_full(kb) - yp_full(ka) )**2 )
        d_ta = dsqrt( (dx*(it -          ia) )**2 + &
                      (dz*(jt -          ja) )**2 + &
                  (yp_full(kt) - yp_full(ka) )**2 )
        ! y = (yb - ya)*(x - xa)/(xb - xa) + ya
        ! y = (yb - ya)*-d_ta/d_ab + ya
        all_c_div(1,i) = 1 + d_ta/d_ab
        all_c_div(2,i) =    -d_ta/d_ab
      end do 
    end block 
    deallocate(yp_full)
  end block
end if 

block 
  integer :: loc_shape(0:2), aux_lo_00(0:0), aux_nx_00(0:0), pos_mpi_ref(0:2)
  loc_shape    =  (/ nx_loc+2*nhalo+2*nskip, nz_loc+2*nhalo+2*nskip, ny_loc+2*nhalo+2*nskip /)
  pos_mpi_ref  =  (/ 0, mpi_pos_z, mpi_pos_y /)
  aux_lo_00    =  0
  aux_nx_00    =  nx_glob
  call diezdecomp_ibm_init(obj_ibm_U, all_ijk_U, all_c_U, n_ghosts_U, n_fluids_U, loc_shape, &
                          aux_lo_00, lo_z_blocks_mpi, lo_y_blocks_mpi, &
                          aux_nx_00,   nz_blocks_mpi,   ny_blocks_mpi, pos_mpi_ref, &
                          nx_glob, nz_glob, ny_glob, nhalo, nskip, irank_mpi, nproc_mpi)
  call diezdecomp_ibm_init(obj_ibm_V, all_ijk_V, all_c_V, n_ghosts_V, n_fluids_V, loc_shape, &
                          aux_lo_00, lo_z_blocks_mpi, lo_y_blocks_mpi, &
                          aux_nx_00,   nz_blocks_mpi,   ny_blocks_mpi, pos_mpi_ref, &
                          nx_glob, nz_glob, ny_glob, nhalo, nskip, irank_mpi, nproc_mpi)
  call diezdecomp_ibm_init(obj_ibm_W, all_ijk_W, all_c_W, n_ghosts_W, n_fluids_W, loc_shape, &
                          aux_lo_00, lo_z_blocks_mpi, lo_y_blocks_mpi, &
                          aux_nx_00,   nz_blocks_mpi,   ny_blocks_mpi, pos_mpi_ref, &
                          nx_glob, nz_glob, ny_glob, nhalo, nskip, irank_mpi, nproc_mpi)
  call diezdecomp_ibm_init(obj_ibm_T, all_ijk_T, all_c_T, n_ghosts_T, n_fluids_T, loc_shape, &
                          aux_lo_00, lo_z_blocks_mpi, lo_y_blocks_mpi, &
                          aux_nx_00,   nz_blocks_mpi,   ny_blocks_mpi, pos_mpi_ref, &
                          nx_glob, nz_glob, ny_glob, nhalo, nskip, irank_mpi, nproc_mpi)
  if (.not.use_incompressible) then 
    call diezdecomp_ibm_init(obj_ibm_div, all_ijk_div, all_c_div, n_ghosts_div, n_fluids_div, loc_shape, &
                            aux_lo_00, lo_z_blocks_mpi, lo_y_blocks_mpi, &
                            aux_nx_00,   nz_blocks_mpi,   ny_blocks_mpi, pos_mpi_ref, &
                            nx_glob, nz_glob, ny_glob, nhalo, nskip, irank_mpi, nproc_mpi)
  end if
end block

allocate(aux_rough_slab(0:nx_glob-1, 0:nz_glob-1))
allocate(rough_U_bot_inds(0:nx_loc-1,0:nz_loc-1))
allocate(rough_V_bot_inds, mold=rough_U_bot_inds)
allocate(rough_W_bot_inds, mold=rough_U_bot_inds)
allocate(rough_T_bot_inds, mold=rough_U_bot_inds)
allocate(rough_U_top_inds, mold=rough_U_bot_inds)
allocate(rough_V_top_inds, mold=rough_U_bot_inds)
allocate(rough_W_top_inds, mold=rough_U_bot_inds)
allocate(rough_T_top_inds, mold=rough_U_bot_inds)

allocate(bulk_U_bot_inds, mold=rough_U_bot_inds)
allocate(bulk_U_top_inds, mold=rough_U_bot_inds)
allocate(bulk_T_bot_inds, mold=rough_U_bot_inds)
allocate(bulk_T_top_inds, mold=rough_U_bot_inds)

block 
  integer :: di_rank, dj_rank, dk_rank 
  di_rank = 0
  dj_rank = lo_z_blocks_mpi(mpi_pos_z)
  dk_rank = lo_y_blocks_mpi(mpi_pos_y)
  call ibm_fill_slabs_2d(rough_U_bot_inds, dk_rank, './geom_data/ibm_inds_U_bot.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
  call ibm_fill_slabs_2d(rough_V_bot_inds, dk_rank, './geom_data/ibm_inds_V_bot.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
  call ibm_fill_slabs_2d(rough_W_bot_inds, dk_rank, './geom_data/ibm_inds_W_bot.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
  call ibm_fill_slabs_2d(rough_T_bot_inds, dk_rank, './geom_data/ibm_inds_T_bot.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
  call ibm_fill_slabs_2d(rough_U_top_inds, dk_rank, './geom_data/ibm_inds_U_top.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
  call ibm_fill_slabs_2d(rough_V_top_inds, dk_rank, './geom_data/ibm_inds_V_top.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
  call ibm_fill_slabs_2d(rough_W_top_inds, dk_rank, './geom_data/ibm_inds_W_top.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
  call ibm_fill_slabs_2d(rough_T_top_inds, dk_rank, './geom_data/ibm_inds_T_top.dat', aux_rough_slab, nx_loc, nz_loc,&
                         di_rank, dj_rank, nx_glob, nz_glob, irank_mpi, nproc_mpi)
end block 

block
  integer :: i,j
  ! clip bounds bulk
  do   i=0,size(rough_U_bot_inds,1)-1
    do j=0,size(rough_U_bot_inds,2)-1
      bulk_U_bot_inds(i,j) = max(0,rough_U_bot_inds(i,j)+1)
    end do
  end do
  do   i=0,size(rough_T_bot_inds,1)-1
    do j=0,size(rough_T_bot_inds,2)-1
      bulk_T_bot_inds(i,j) = max(0,rough_T_bot_inds(i,j)+1)
    end do
  end do

  do   i=0,size(rough_U_top_inds,1)-1
    do j=0,size(rough_U_top_inds,2)-1
      bulk_U_top_inds(i,j) = min(ny_loc-1,rough_U_top_inds(i,j)-1)
    end do
  end do
  do   i=0,size(rough_T_top_inds,1)-1
    do j=0,size(rough_T_top_inds,2)-1
      bulk_T_top_inds(i,j) = min(ny_loc-1,rough_T_top_inds(i,j)-1)
    end do
  end do

  ! clip final bounds IBM
  do   i=0,size(rough_U_bot_inds,1)-1
    do j=0,size(rough_U_bot_inds,2)-1
      rough_U_bot_inds(i,j) = min(rough_U_bot_inds(i,j),ny_loc-1)
    end do
  end do
  do   i=0,size(rough_V_bot_inds,1)-1
    do j=0,size(rough_V_bot_inds,2)-1
      rough_V_bot_inds(i,j) = min(rough_V_bot_inds(i,j),ny_loc-1)
    end do
  end do
  do   i=0,size(rough_W_bot_inds,1)-1
    do j=0,size(rough_W_bot_inds,2)-1
      rough_W_bot_inds(i,j) = min(rough_W_bot_inds(i,j),ny_loc-1)
    end do
  end do
  do   i=0,size(rough_T_bot_inds,1)-1
    do j=0,size(rough_T_bot_inds,2)-1
      rough_T_bot_inds(i,j) = min(rough_T_bot_inds(i,j),ny_loc-1)
    end do
  end do

  do   i=0,size(rough_U_top_inds,1)-1
    do j=0,size(rough_U_top_inds,2)-1
      rough_U_top_inds(i,j) = max(rough_U_top_inds(i,j),0)
    end do
  end do
  do   i=0,size(rough_V_top_inds,1)-1
    do j=0,size(rough_V_top_inds,2)-1
      rough_V_top_inds(i,j) = max(rough_V_top_inds(i,j),0)
    end do
  end do
  do   i=0,size(rough_W_top_inds,1)-1
    do j=0,size(rough_W_top_inds,2)-1
      rough_W_top_inds(i,j) = max(rough_W_top_inds(i,j),0)
    end do
  end do
  do   i=0,size(rough_T_top_inds,1)-1
    do j=0,size(rough_T_top_inds,2)-1
      rough_T_top_inds(i,j) = max(rough_T_top_inds(i,j),0)
    end do
  end do
end block

allocate(buffer_ibm(0:max(obj_ibm_U  %wsize_buf_mpi,&
                          obj_ibm_V  %wsize_buf_mpi,&
                          obj_ibm_W  %wsize_buf_mpi,&
                          obj_ibm_T  %wsize_buf_mpi,&
                          obj_ibm_div%wsize_buf_mpi)-1))
!$acc enter data create(buffer_ibm)

!$acc enter data copyin(rough_T_bot_inds, rough_T_top_inds)
!$acc enter data copyin(rough_U_bot_inds, rough_U_top_inds)
!$acc enter data copyin(rough_V_bot_inds, rough_V_top_inds)
!$acc enter data copyin(rough_W_bot_inds, rough_W_top_inds)


!$acc enter data copyin(bulk_U_bot_inds, bulk_T_bot_inds)
!$acc enter data copyin(bulk_U_top_inds, bulk_T_top_inds)

deallocate(aux_rough_slab)